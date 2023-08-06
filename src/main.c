#include <stdio.h>
#include <SDL2/SDL.h>

#include "typedefs.h"
#include "constants.h"
#include "vec.h"
#include "input.h"
#include "delaunay_mesh.h"
#include "game.h"
#include "platform_metrics.h"
// #include "profiler.h"

#define ASSERT(_e, ...)               \
    if (!(_e))                        \
    {                                 \
        fprintf(stderr, __VA_ARGS__); \
        exit(1);                      \
    }

#define ASSET_NAME_BYTE_COUNT 16

// ------------------------------------------------------------------------------
// Big binary assets blob that we load at init.
u8 *ASSETS_BINARY_BLOB = NULL;
u32 ASSETS_BINARY_BLOB_SIZE = 0; // number of bytes

struct BinaryAssetTableOfContentEntry
{
    u32 byte_offset;
    char name[ASSET_NAME_BYTE_COUNT];
};

// ------------------------------------------------------------------------------

struct CameraState
{
    v2 fov; // width (x) and height (y) of the field of view at unit distance from the camera
    v2 pos;
    v2 dir;
    f32 z;
};

// ------------------------------------------------------------------------------
// DOOM Assets

//  Each palette in the PLAYPAL lump contains 256 colors totaling 768 bytes,
//  where each color is broken into three unsigned bytes. Each of these color components (red, green, and blue) range between 0 and 255.
u32 PALETTE_OFFSET = 0;

struct Patch
{
    u16 size_x;            // width of the graphic in pixels
    u16 size_y;            // height of the graphic in pixels
    i16 offset_x;          // x offset from the screen origin (to left)
    i16 offset_y;          // y offset from the screen origin (down)
    u32 column_offsets[1]; // Actually of length `size_x`.
};

struct PatchEntry
{
    u32 byte_offset;
    struct Patch *patch;
};
struct PatchEntry SKEL_PATCH_ENTRIES[360];

u32 PATCH_COUNT;    // number of patches
u32 *PATCH_OFFSETS; // the start of each patch (past the 16-char name to match our struct)

// ------------------------------------------------------------------------------
// TOOM Assets

struct BinaryAsset2TableOfContentEntry
{
    u32 offset;
    u32 size;
    char name[ASSET_NAME_BYTE_COUNT];
};

u8 *ASSETS_BINARY_BLOB2;
u32 ASSETS_BINARY_BLOB2_SIZE = 0; // number of bytes

// ------------------------------------------------------------------------------

struct
{
    SDL_Window *window;
    SDL_Texture *texture;
    SDL_Renderer *renderer;

    u32 pixels[SCREEN_SIZE_X * SCREEN_SIZE_Y]; // row-major

    f32 wall_raycast_radius[SCREEN_SIZE_X];

    bool quit;

    struct CameraState camera;
    struct GameState game_state;
    struct KeyBoardState keyboard_state;
} state;

static void LoadAssets(struct GameMap *game_map)
{
    // If we currently have any loaded assets, free them.
    // Note that we assume we are running single-threaded. If this changes in the future,
    // we will want to load the assets in a separate thread and swap them over at the
    // appropriate time.
    if (ASSETS_BINARY_BLOB)
    {
        free(ASSETS_BINARY_BLOB);
    }
    {
        // The format is:
        // HEADER:
        //    "TOOM"   - 4 chars
        // Data:
        //   necessary data, laid out as needed.
        // Table of Contents:
        //   array of table of content entries:
        //      u32       offset in file # number of bytes past 'TOOM' to read at. First entry will have offset 0
        //      char[16]  name           # null-terminated string label, e.g. "floor_textures"
        //   u32 n_toc_entries = number of table of content entries

        FILE *fileptr = fopen("assets/assets.bin", "rb");
        ASSERT(fileptr, "Error opening assets file\n");

        // Count the number of bytes
        fseek(fileptr, 0, SEEK_END);
        u32 n_bytes_in_file = ftell(fileptr);
        ASSETS_BINARY_BLOB_SIZE = n_bytes_in_file - 4; // Everything past the 'TOOM' header

        // Read in the binary assets as a single blob
        fseek(fileptr, 4, SEEK_SET); // Skip the header bytes
        ASSETS_BINARY_BLOB = (u8 *)malloc(ASSETS_BINARY_BLOB_SIZE);
        ASSERT(ASSETS_BINARY_BLOB, "Failed to allocate assets blob\n");
        ASSERT(fread(ASSETS_BINARY_BLOB, sizeof(u8), ASSETS_BINARY_BLOB_SIZE, fileptr) == ASSETS_BINARY_BLOB_SIZE, "Failed to read assets blob when loading assets\n");

        fclose(fileptr);
    }
    {
        // Process the loaded assets from the loaded binary blob

        // Read the number of table of content entries
        u32 byte_index = ASSETS_BINARY_BLOB_SIZE - sizeof(u32);
        u32 n_toc_entries = *(u32 *)(ASSETS_BINARY_BLOB + byte_index);
        ASSERT(ASSETS_BINARY_BLOB_SIZE > sizeof(struct BinaryAssetTableOfContentEntry) * n_toc_entries + 4, "Number of table of content entries is impossible given the number of bytes\n");

        // Scan through them in reverse order
        for (int i = n_toc_entries; i > 0; i--)
        {
            byte_index -= sizeof(struct BinaryAssetTableOfContentEntry);
            struct BinaryAssetTableOfContentEntry *entry = (struct BinaryAssetTableOfContentEntry *)(ASSETS_BINARY_BLOB + byte_index);
            // Make the name null-terminated just in case.
            entry->name[15] = '\0';
            printf("Entry %d: %s at offset %d\n", i, entry->name, entry->byte_offset);

            if (strncmp(entry->name, "SKEL", 4) == 0)
            {
                int frame_index = entry->name[4] - '0';
                frame_index = 10 * frame_index + (entry->name[5] - '0');
                frame_index = 10 * frame_index + (entry->name[6] - '0');
                SKEL_PATCH_ENTRIES[frame_index].byte_offset = entry->byte_offset;
                SKEL_PATCH_ENTRIES[frame_index].patch = (struct Patch *)(ASSETS_BINARY_BLOB + entry->byte_offset);
            }
        }
    }

    // Load assets #2
    if (ASSETS_BINARY_BLOB2)
    {
        free(ASSETS_BINARY_BLOB2);
        free(PATCH_OFFSETS);
    }
    {
        // The format is:
        // HEADER:
        //    "TOOM"     - 4 chars
        //    n_entries  - u32
        //    toc_offset - u32
        // Data:
        //   necessary data, laid out as needed.
        // Table of Contents:
        //   array of table of content entries
        //   u32 n_toc_entries = number of table of content entries

        FILE *fileptr = fopen("assets/toomed.bin", "rb");
        ASSERT(fileptr, "Error opening toomed binary\n");

        // Count the number of bytes
        fseek(fileptr, 0, SEEK_END);
        ASSETS_BINARY_BLOB2_SIZE = ftell(fileptr);

        // Read in the entire file
        fseek(fileptr, 0, SEEK_SET);
        ASSETS_BINARY_BLOB2 = (u8 *)malloc(ASSETS_BINARY_BLOB2_SIZE);
        ASSERT(ASSETS_BINARY_BLOB2, "Failed to allocate toomed byte array\n");
        ASSERT(fread(ASSETS_BINARY_BLOB2, sizeof(u8), ASSETS_BINARY_BLOB2_SIZE, fileptr) == ASSETS_BINARY_BLOB2_SIZE, "Failed to read toomed bytes\n");

        u32 n_entries = *(u32 *)(ASSETS_BINARY_BLOB2 + 0x04);
        u32 toc_offset = *(u32 *)(ASSETS_BINARY_BLOB2 + 0x08);

        // Process the table of contents
        u32 offset = toc_offset;
        for (u32 toc_index = 0; toc_index < n_entries; toc_index++)
        {
            struct BinaryAsset2TableOfContentEntry *entry = (struct BinaryAsset2TableOfContentEntry *)(ASSETS_BINARY_BLOB2 + offset);
            printf("Entry %d: %.16s at offset %d with size %d\n", toc_index, entry->name, entry->offset, entry->size);
            offset += sizeof(struct BinaryAsset2TableOfContentEntry);

            if (strcmp(entry->name, "palettes") == 0)
            {
                u32 offset = entry->offset;

                // u32 n_palettes = *(u32 *)(ASSETS_BINARY_BLOB2 + offset);
                offset += sizeof(u32);

                // TODO: Import other palettes.
                PALETTE_OFFSET = offset;
            }
            if (strcmp(entry->name, "patches") == 0)
            {
                u32 offset = entry->offset;

                PATCH_COUNT = *(u32 *)(ASSETS_BINARY_BLOB2 + offset);
                offset += sizeof(u32);

                PATCH_OFFSETS = (u32 *)malloc(PATCH_COUNT * sizeof(u32));
                ASSERT(PATCH_OFFSETS, "Failed to allocate patch offsets\n");

                for (u32 i_patch = 0; i_patch < PATCH_COUNT; i_patch++)
                {
                    // Skip the 16-char name
                    offset += 16;

                    // Set the patch
                    PATCH_OFFSETS[i_patch] = offset;

                    // Skip the patch header data
                    u16 size_x = *(u16 *)(ASSETS_BINARY_BLOB2 + offset);
                    offset += sizeof(u16);
                    offset += sizeof(u16);
                    offset += sizeof(u16);
                    offset += sizeof(u16);

                    // Skip the column offsets
                    offset += size_x * sizeof(u32);

                    // Skip the post data
                    u32 n_bytes_post_data = *(u32 *)(ASSETS_BINARY_BLOB2 + offset);
                    offset += sizeof(u32) + n_bytes_post_data;
                }
            }
            else if (strcmp(entry->name, "geometry_mesh") == 0)
            {
                struct DelaunayMesh *mesh = (struct DelaunayMesh *)malloc(sizeof(struct DelaunayMesh));
                game_map->geometry_mesh = mesh;

                mesh->square_bounding_radius = 0.0;
                mesh->square_min_dist_to_vertex = 0.0;
                mesh->min_dist_to_edge = 0.0;

                u32 mesh_offset = entry->offset;
                mesh->n_vertices = *(u32 *)(ASSETS_BINARY_BLOB2 + mesh_offset);
                mesh->max_n_vertices = mesh->n_vertices;
                mesh_offset += sizeof(u32);

                mesh->n_quarter_edges = *(u32 *)(ASSETS_BINARY_BLOB2 + mesh_offset);
                mesh->max_n_quarter_edges = mesh->n_quarter_edges;
                mesh_offset += sizeof(u32);

                mesh->vertices = (v2 *)(ASSETS_BINARY_BLOB2 + mesh_offset);
                mesh_offset += mesh->n_vertices * sizeof(v2);

                // Allocate the quarter edges
                mesh->quarter_edges = (QuarterEdge *)malloc(mesh->n_quarter_edges * sizeof(QuarterEdge));
                for (size_t qe_index = 0; qe_index < mesh->n_quarter_edges; qe_index++)
                {
                    QuarterEdge *qe = mesh->quarter_edges + qe_index;
                    qe->index = qe_index;
                    qe->vertex = mesh->vertices + *(u32 *)(ASSETS_BINARY_BLOB2 + mesh_offset);
                    mesh_offset += sizeof(u32);
                    qe->next = mesh->quarter_edges + *(u32 *)(ASSETS_BINARY_BLOB2 + mesh_offset);
                    mesh_offset += sizeof(u32);
                    qe->rot = mesh->quarter_edges + *(u32 *)(ASSETS_BINARY_BLOB2 + mesh_offset);
                    mesh_offset += sizeof(u32);
                }
            }
            else if (strcmp(entry->name, "side_infos") == 0)
            {
                u32 mesh_offset = entry->offset;
                game_map->n_side_infos = *(u32 *)(ASSETS_BINARY_BLOB2 + mesh_offset);
                mesh_offset += sizeof(u32);

                game_map->side_infos = (struct SideInfo *)(ASSETS_BINARY_BLOB2 + mesh_offset);
                mesh_offset += game_map->n_side_infos * sizeof(struct SideInfo);

                game_map->quarter_edge_index_to_side_info_index = (u16 *)(ASSETS_BINARY_BLOB2 + mesh_offset);
                mesh_offset += game_map->geometry_mesh->n_quarter_edges * sizeof(u16);

                // Validation - ensure that our indices are all in range
                for (int i = 0; i < game_map->geometry_mesh->n_quarter_edges; i++)
                {
                    u16 side_info_index = game_map->quarter_edge_index_to_side_info_index[i];
                    ASSERT(side_info_index < game_map->n_side_infos || side_info_index == 0xFFFF,
                           "Side info index %d, which is the %dth index, is out of bounds\n", side_info_index, i);
                }
            }
            else if (strcmp(entry->name, "sectors") == 0)
            {
                u32 mesh_offset = entry->offset;
                game_map->n_sectors = *(u32 *)(ASSETS_BINARY_BLOB2 + mesh_offset);
                mesh_offset += sizeof(u32);

                game_map->sectors = (struct Sector *)(ASSETS_BINARY_BLOB2 + mesh_offset);
                mesh_offset += game_map->n_side_infos * sizeof(struct Sector);
            }
        }

        fclose(fileptr);
    }
}

void RenderPatchColumn(u32 *pixels, int x_screen, int y_lower,
                       int y_upper, int y_lo, int y_hi, f32 x_along_texture,
                       u32 x_base_offset, u32 y_base_offset,
                       f32 column_height, u32 patch_offset)
{
    // TODO: renderdata
    u32 n_pixels_per_world_unit = 64;

    // TODO: for now, ignore y_base_offset.

    struct Patch *patch = (struct Patch *)(ASSETS_BINARY_BLOB2 + patch_offset);

    // Get the start of the post data
    // NOTE: sizeof(Patch) includes one u32, so we have to additionally traverse via the other ones.
    //       plus we have to bypass the number of bytes (as a u32).
    u8 *post_data = (u8 *)(ASSETS_BINARY_BLOB2 + patch_offset + sizeof(struct Patch) + sizeof(u32) * patch->size_x);

    // The index into the patch of the column that we want to draw.
    u32 x_patch = ((int)(n_pixels_per_world_unit * x_along_texture) + x_base_offset) %
                  patch->size_x;

    // y_lower = screen y coordinate of bottom of column (can exceed screen bounds)
    // y_upper = screen y coordinate of top of column (can exceed screen bounds)
    // y_lo    = screen y coordinate where we end drawing (does not exceed screen bounds)
    // y_hi    = screen y coordinate where we start drawing (does not exceed screen bounds)
    // column_height = real-world height of the painting surface

    // How many patch pixels high the column is.
    f32 y_patch_height_pix = n_pixels_per_world_unit * column_height;

    // y_patch = m * y_screen + b for converting screen to patch coordinate
    //                         0 = m * y_upper + b   -> b = -m * y_upper
    //    y_patch_height_pix - 1 = m * y_lower + b
    //                           = m * y_lower - m * y_upper
    //                           = m * (y_lower - y_upper)
    // m = (y_patch_height_pix - 1) / (y_lower - y_upper)

    f32 m = (f32)(y_patch_height_pix - 1.0f) / (y_lower - y_upper);

    // The number of (continuous) patch pixels y changes per screen pixel
    f32 y_patch_step_per_screen_pixel = m; // If y_screen goes up by 1, y_patch goes up this much
    f32 y_screen_step_per_patch_pixel = 1.0f / m;

    f32 y_patch = 0.0f;
    f32 y_screen = y_upper;
    while (y_screen > y_lo)
    {
        u32 column_offset = patch->column_offsets[x_patch];
        while (post_data[column_offset] != 0xFF)
        {
            u8 y_patch_delta = post_data[column_offset];
            column_offset++;
            int post_length = post_data[column_offset];
            column_offset++;

            // skip transparent pixels
            y_patch += y_patch_delta;
            y_screen += y_patch_delta * y_screen_step_per_patch_pixel;

            // process the post. We have `post_length` pixels to draw
            while (post_length > 0 && y_screen > y_lo)
            {
                // Keep decreasing y_screen (and increasing y_patch) as long as we are within the
                // post data.

                // Render pixels
                u8 palette_index = post_data[column_offset];
                // palette_index = colormap.map[palette_index]; // TODO
                u8 r = *(u8 *)(ASSETS_BINARY_BLOB2 + PALETTE_OFFSET + 3 * palette_index);
                u8 g = *(u8 *)(ASSETS_BINARY_BLOB2 + PALETTE_OFFSET + 3 * palette_index + 1);
                u8 b = *(u8 *)(ASSETS_BINARY_BLOB2 + PALETTE_OFFSET + 3 * palette_index + 2);
                u32 abgr = 0xFF000000 + (((u32)b) << 16) + (((u32)g) << 8) + r;

                // Render this color for all screen pixels that map to y_patch
                u16 y_patch_discrete = (u16)y_patch;
                int y_screen_discrete = (int)y_screen;
                while ((u16)y_patch == y_patch_discrete)
                {
                    if (y_screen_discrete < y_hi && y_screen_discrete > y_lo)
                    {
                        pixels[(y_screen_discrete * SCREEN_SIZE_X) + x_screen] = abgr;
                    }
                    y_screen_discrete -= 1;
                    y_screen -= 1.0f;
                    y_patch -= y_patch_step_per_screen_pixel;
                }

                u16 patch_delta = (u16)y_patch - y_patch_discrete;
                while (patch_delta > 0)
                {
                    column_offset += 1;
                    post_length -= 1;
                    patch_delta -= 1;
                }
            }
        }
    }
}

// void RenderFloorAndCeiling(
//     u32 *pixels,
//     struct CameraState *camera)
// {
//     // Ray direction for x = 0
//     f32 half_camera_width = camera->fov.x / 2.0f;
//     f32 ray_dir_lo_x = camera->dir.x - half_camera_width * camera->dir.y;
//     f32 ray_dir_lo_y = camera->dir.y + half_camera_width * camera->dir.x;

//     // Ray direction for x = SCREEN_SIZE_X
//     f32 ray_dir_hi_x = camera->dir.x + half_camera_width * camera->dir.y;
//     f32 ray_dir_hi_y = camera->dir.y - half_camera_width * camera->dir.x;

//     // Draw floor
//     for (int y = 0; y < SCREEN_SIZE_Y / 2; y++)
//     {
//         // Radius
//         f32 zpp = (SCREEN_SIZE_Y / 2.0f - y) * (camera->fov.y / SCREEN_SIZE_Y);
//         f32 radius = camera->z / zpp; // TODO: Precompute for each y

//         // Location of the 1st ray's intersection
//         f32 hit_x = camera->pos.x + radius * ray_dir_lo_x;
//         f32 hit_y = camera->pos.y + radius * ray_dir_lo_y;

//         // Each step is (hit_x2 - hit_x) / SCREEN_SIZE_X;
//         // = ((camera->pos.x + radius * ray_dir_lo_x) - (camera->pos.x + radius * ray_dir_lo_x)) / SCREEN_SIZE_X
//         // = (radius * ray_dir_lo_x - (radius * ray_dir_lo_x)) / SCREEN_SIZE_X
//         // = radius * (ray_dir_hi_x - ray_dir_lo_x) / SCREEN_SIZE_X
//         f32 step_x = radius * (ray_dir_hi_x - ray_dir_lo_x) / SCREEN_SIZE_X;
//         f32 step_y = radius * (ray_dir_hi_y - ray_dir_lo_y) / SCREEN_SIZE_X;

//         for (int x = 0; x < SCREEN_SIZE_X; x++)
//         {

//             // int x_ind_hit = (int)(floorf(hit_x / TILE_WIDTH));
//             // int y_ind_hit = (int)(floorf(hit_y / TILE_WIDTH));
//             // f32 x_rem_hit = hit_x - TILE_WIDTH * x_ind_hit;
//             // f32 y_rem_hit = hit_y - TILE_WIDTH * y_ind_hit;
//             // x_ind_hit = clamp(x_ind_hit, 0, MAPDATA.n_tiles_x - 1);
//             // y_ind_hit = clamp(y_ind_hit, 0, MAPDATA.n_tiles_y - 1);

//             u32 texture_x_offset = 0;
//             u32 texture_y_offset = 0; // (MAPDATA.floor[GetMapDataIndex(&MAPDATA, x_ind_hit, y_ind_hit)] - 1) * TEXTURE_SIZE;

//             u32 texture_x = 0; // (int)(x_rem_hit / TILE_WIDTH * TEXTURE_SIZE);
//             u32 texture_y = 0; // (int)(y_rem_hit / TILE_WIDTH * TEXTURE_SIZE);

//             u32 color = GetColumnMajorPixelAt(&BITMAP, texture_x + texture_x_offset, texture_y + texture_y_offset);
//             pixels[(y * SCREEN_SIZE_X) + x] = color;

//             // step
//             hit_x += step_x;
//             hit_y += step_y;
//         }
//     }

//     // Draw ceiling
//     for (int y = SCREEN_SIZE_Y / 2 + 1; y < SCREEN_SIZE_Y; y++)
//     {
//         // Radius
//         f32 zpp = (y - (SCREEN_SIZE_Y / 2.0f)) * (camera->fov.y / SCREEN_SIZE_Y);
//         f32 radius = (WALL_HEIGHT - camera->z) / zpp; // TODO: Precompute for each y

//         // Location of the 1st ray's intersection
//         f32 hit_x = camera->pos.x + radius * ray_dir_lo_x;
//         f32 hit_y = camera->pos.y + radius * ray_dir_lo_y;

//         // Each step toward hit2
//         f32 step_x = radius * (ray_dir_hi_x - ray_dir_lo_x) / SCREEN_SIZE_X;
//         f32 step_y = radius * (ray_dir_hi_y - ray_dir_lo_y) / SCREEN_SIZE_X;

//         for (int x = 0; x < SCREEN_SIZE_X; x++)
//         {
//             // int x_ind_hit = (int)(floorf(hit_x / TILE_WIDTH));
//             // int y_ind_hit = (int)(floorf(hit_y / TILE_WIDTH));
//             // f32 x_rem_hit = hit_x - TILE_WIDTH * x_ind_hit;
//             // f32 y_rem_hit = hit_y - TILE_WIDTH * y_ind_hit;
//             // x_ind_hit = clamp(x_ind_hit, 0, MAPDATA.n_tiles_x - 1);
//             // y_ind_hit = clamp(y_ind_hit, 0, MAPDATA.n_tiles_y - 1);

//             u32 texture_x_offset = 0;
//             u32 texture_y_offset = 0; // (MAPDATA.ceiling[GetMapDataIndex(&MAPDATA, x_ind_hit, y_ind_hit)] - 1) * TEXTURE_SIZE;

//             u32 texture_x = 5; // (int)(x_rem_hit / TILE_WIDTH * TEXTURE_SIZE);
//             u32 texture_y = 5; // (int)(y_rem_hit / TILE_WIDTH * TEXTURE_SIZE);

//             u32 color = GetColumnMajorPixelAt(&BITMAP, texture_x + texture_x_offset, texture_y + texture_y_offset);
//             pixels[(y * SCREEN_SIZE_X) + x] = color;

//             // step
//             hit_x += step_x;
//             hit_y += step_y;
//         }
//     }
// }

void RenderWalls(
    u32 *pixels,
    f32 *wall_raycast_radius,
    const struct CameraState *camera,
    const struct GameMap *game_map,
    QuarterEdge *qe_camera)
{
    const struct DelaunayMesh *mesh = game_map->geometry_mesh;

    f32 half_screen_size = SCREEN_SIZE_Y / 2.0f;
    f32 screen_size_y_over_fov_y = SCREEN_SIZE_Y / camera->fov.y;

    u32 color_ceil = 0xFF222222;
    u32 color_floor = 0xFF444444;

    for (int x = 0; x < SCREEN_SIZE_X; x++)
    {
        // Camera to pixel column
        const f32 dw =
            camera->fov.x / 2 - (camera->fov.x * x) / SCREEN_SIZE_X; // TODO: Precompute once.
        const v2 cp = {camera->dir.x - dw * camera->dir.y,
                       camera->dir.y + dw * camera->dir.x};

        // Distance from the camera to the column
        const f32 cam_len = length(cp);

        // Ray direction through this column
        const v2 dir = {cp.x / cam_len, cp.y / cam_len};

        // Start at the camera pos
        v2 pos = camera->pos;
        QuarterEdge *qe_dual = qe_camera;

        // The edge vector of the face that we last crossed
        v2 v_face = {0.0, 0.0};

        // Step through triangles until we hit a solid triangle
        int y_hi = SCREEN_SIZE_Y;
        int y_lo = -1;

        int n_steps = 0;
        while (n_steps < 100)
        {
            n_steps += 1;

            // Grab the enclosing triangle.
            QuarterEdge *qe_ab = qe_dual->rot;
            QuarterEdge *qe_bc = qe_dual->next->rot;
            QuarterEdge *qe_ca = qe_dual->next->next->rot;

            const v2 *a = qe_ab->vertex;
            const v2 *b = qe_bc->vertex;
            const v2 *c = qe_ca->vertex;

            // Project our ray out far enough that it would exit our mesh
            f32 projection_distance = 100.0; // Ridiculously large
            const v2 pos_next_delta = {projection_distance * dir.x, projection_distance * dir.y};
            const v2 pos_next = add(pos, pos_next_delta);

            f32 min_interp = INFINITY;
            QuarterEdge *qe_side = NULL;

            // See if we cross any of the 3 faces for the triangle we are in,
            // and cross the first segment.
            const f32 eps = 1e-4;
            if (GetRightHandedness(a, b, &pos_next) < -eps)
            {
                // We would cross AB
                v2 v = sub(*b, *a);
                v2 w = sub(pos, *a);
                float interp_ab = cross(v, w) / cross(pos_next_delta, v);
                if (interp_ab < min_interp)
                {
                    min_interp = interp_ab;
                    qe_side = qe_ab;
                    v_face = v;
                }
            }
            if (GetRightHandedness(b, c, &pos_next) < -eps)
            {
                // We would cross BC
                v2 v = sub(*c, *b);
                v2 w = sub(pos, *b);
                float interp_bc = cross(v, w) / cross(pos_next_delta, v);
                if (interp_bc < min_interp)
                {
                    min_interp = interp_bc;
                    qe_side = qe_bc;
                    v_face = v;
                }
            }
            if (GetRightHandedness(c, a, &pos_next) < -eps)
            {
                // We would cross CA
                v2 v = sub(*a, *c);
                v2 w = sub(pos, *c);
                float interp_ca = cross(v, w) / cross(pos_next_delta, v);
                if (interp_ca < min_interp)
                {
                    min_interp = interp_ca;
                    qe_side = qe_ca;
                    v_face = v;
                }
            }

            if (qe_side)
            {
                // Move to the face
                pos.x += min_interp * pos_next_delta.x;
                pos.y += min_interp * pos_next_delta.y;

                qe_dual = qe_side->rot; // The next face (Should always be non-null)

                u32 side_info_index = game_map->quarter_edge_index_to_side_info_index[qe_side->index];
                if (side_info_index != 0xFFFF)
                {
                    struct SideInfo *side_info = game_map->side_infos + side_info_index;

                    const f32 ray_len = max(length(sub(pos, camera->pos)), 0.01f);
                    const f32 gamma = cam_len / ray_len * screen_size_y_over_fov_y;
                    wall_raycast_radius[x] = ray_len;

                    struct Sector *sector = game_map->sectors + side_info->sector_id;
                    f32 z_ceil = sector->z_ceil;
                    f32 z_upper = sector->z_ceil;
                    f32 z_lower = sector->z_floor;
                    f32 z_floor = sector->z_floor;

                    // Get the height on the other side, if it is passable.
                    const bool is_passable =
                        ((side_info->flags & SIDEINFO_FLAG_PASSABLE) > 0);
                    if (is_passable)
                    {
                        QuarterEdge *qe_sym = QESym(qe_side);
                        u32 side_info_index_sym = game_map->quarter_edge_index_to_side_info_index[qe_sym->index];
                        if (side_info_index_sym != 0xFFFF)
                        {
                            struct SideInfo *side_info_sym = game_map->side_infos + side_info_index_sym;
                            struct Sector *sector_sym = game_map->sectors + side_info_sym->sector_id;
                            z_lower = sector_sym->z_floor;
                            z_upper = sector_sym->z_ceil;
                        }
                        else
                        {
                            printf("Unexpected nullptr qe_sym!\n");
                        }
                    }

                    int y_ceil = (int)(half_screen_size + gamma * (z_ceil - camera->z));
                    int y_upper = (int)(half_screen_size + gamma * (z_upper - camera->z));
                    int y_lower = (int)(half_screen_size + gamma * (z_lower - camera->z));
                    int y_floor = (int)(half_screen_size + gamma * (z_floor - camera->z));

                    // Calculate where along the segment we intersected.
                    QuarterEdge *qe_face_src = QETor(qe_dual);
                    f32 x_along_texture =
                        length(v_face) - length(sub(pos, *(qe_face_src->vertex)));

                    // Render the ceiling above the upper texture
                    while (y_hi > y_ceil)
                    {
                        y_hi--;
                        pixels[(y_hi * SCREEN_SIZE_X) + x] = color_ceil;
                    }

                    // Render the upper texture
                    if (y_upper < y_hi)
                    {
                        u32 patch_offset = PATCH_OFFSETS[side_info->texture_info_upper.texture_id];
                        f32 column_height = z_ceil - z_upper;
                        RenderPatchColumn(
                            pixels, x, y_upper, y_ceil, y_upper, y_hi,
                            x_along_texture,
                            side_info->texture_info_upper.x_offset,
                            side_info->texture_info_upper.y_offset, column_height,
                            patch_offset);
                        y_hi = y_upper;
                    }

                    // Render the floor below the lower texture
                    while (y_lo < y_floor)
                    {
                        y_lo++;
                        pixels[(y_lo * SCREEN_SIZE_X) + x] = color_floor;
                    }

                    // Render the lower texture
                    if (y_lower > y_lo)
                    {
                        u32 patch_offset = PATCH_OFFSETS[side_info->texture_info_lower.texture_id];
                        f32 column_height = z_lower - z_floor;
                        RenderPatchColumn(
                            pixels, x, y_floor, y_lower, y_lo,
                            y_lower, x_along_texture, side_info->texture_info_lower.x_offset,
                            side_info->texture_info_lower.y_offset, column_height, patch_offset);

                        y_lo = y_lower;
                    }

                    // Continue on with our projection if the side is passable.
                    if (is_passable)
                    {
                        continue;
                    }

                    // The side info has a solid wall.
                    u32 patch_offset = PATCH_OFFSETS[side_info->texture_info_middle.texture_id];
                    f32 column_height = z_upper - z_lower;
                    RenderPatchColumn(
                        pixels, x, y_lower, y_upper, y_lo, y_hi, x_along_texture,
                        side_info->texture_info_lower.x_offset,
                        side_info->texture_info_lower.y_offset, column_height, patch_offset);

                    break;
                }

                // Also break if it is the boundary
                if (DelaunayMeshIsBoundaryVertex(mesh, qe_side->vertex))
                {
                    break;
                }
            }
            else
            {
                break;
            }
        }
    }
}

void RenderObjects(
    u32 *pixels,
    f32 *wall_raycast_radius,
    struct CameraState *camera)
{
    static f32 DOOM_HEIGHT_PER_PIX = ((f32)(DOOM_PIX_PER_WALL_WIDTH)) / (TEXTURE_SIZE * TEXTURE_SIZE);

    u8 *DATA = ASSETS_BINARY_BLOB;

    f32 camera_heading = atan2(camera->dir.y, camera->dir.x); // TODO: inefficient

    v2 sprite_pos = {10.5f, 4.5f};
    v2 sprite_rel_camera = sub(sprite_pos, camera->pos);
    f32 dist_to_player = length(sprite_rel_camera);

    f32 s = camera->dir.y;
    f32 c = camera->dir.x;
    v2 sprite_pos_cam_body = {
        c * sprite_rel_camera.x + s * sprite_rel_camera.y,
        c * sprite_rel_camera.y - s * sprite_rel_camera.x};

    // Only render if it is on the postive side of the camera
    if (sprite_pos_cam_body.x > 1e-3)
    {

        f32 monster_heading = PI;
        f32 monster_heading_in_cam_body = monster_heading - camera_heading;

        // Determine the sprite to use based on our viewing angle.
        int n_patch_entries = 360; // 8
        f32 monster_heading_rel = fmod(PI + monster_heading_in_cam_body + PI / n_patch_entries, 2 * PI);
        int monster_frame = clamp((int)(monster_heading_rel * n_patch_entries / (2.0 * PI)), 0, n_patch_entries); // TODO: Better way to handle the clamp?
        // struct PatchEntry *patch_entry = &CYBR_PATCH_ENTRIES[monster_frame];
        struct PatchEntry *patch_entry = &SKEL_PATCH_ENTRIES[monster_frame];

        // Calculate the column pixel bounds
        int sprite_size_y = patch_entry->patch->size_y;
        f32 patch_height = DOOM_HEIGHT_PER_PIX * sprite_size_y;
        f32 patch_z_lo = DOOM_HEIGHT_PER_PIX * (patch_entry->patch->offset_y - sprite_size_y);
        f32 patch_z_hi = patch_z_lo + patch_height;
        const f32 cam_len = sqrt(1.0 + (sprite_pos_cam_body.y / sprite_pos_cam_body.x) * (sprite_pos_cam_body.y / sprite_pos_cam_body.x));
        int y_lo = (int)(SCREEN_SIZE_Y / 2.0f + cam_len * (patch_z_lo - camera->z) / dist_to_player * SCREEN_SIZE_Y / camera->fov.y);
        int y_hi = (int)(SCREEN_SIZE_Y / 2.0f + cam_len * (patch_z_hi - camera->z) / dist_to_player * SCREEN_SIZE_Y / camera->fov.y);
        int y_lo_capped = max(y_lo, 0);
        int y_hi_capped = min(y_hi, SCREEN_SIZE_Y - 1);
        u32 denom = max(1, y_hi - y_lo);
        f32 y_step = (f32)(sprite_size_y) / denom;

        int sprite_size_x = patch_entry->patch->size_x;
        f32 patch_halfwidth = DOOM_HEIGHT_PER_PIX * sprite_size_x / 2.0f;
        f32 patch_x_offset = DOOM_HEIGHT_PER_PIX * (patch_entry->patch->offset_x - 0.5f * sprite_size_x);
        int x_column_lo = (int)((0.5 - ((sprite_pos_cam_body.y + patch_x_offset + patch_halfwidth) / sprite_pos_cam_body.x) / (camera->fov.x)) * SCREEN_SIZE_X);
        int x_column_hi = (int)((0.5 - ((sprite_pos_cam_body.y + patch_x_offset - patch_halfwidth) / sprite_pos_cam_body.x) / (camera->fov.x)) * SCREEN_SIZE_X);
        f32 x_step = ((f32)(sprite_size_x) / (x_column_hi - x_column_lo + 1));
        f32 x_loc = 0.0f;
        for (int x = x_column_lo; x <= x_column_hi; x++)
        {

            if (x >= 0 && x < SCREEN_SIZE_X && wall_raycast_radius[x] > dist_to_player)
            {
                u32 texture_x = min((u32)(x_loc), sprite_size_x - 1);

                // Grab the column.
                u32 column_offset = patch_entry->patch->column_offsets[texture_x];

                // Grab the first post.
                u8 y_texture_skip = DATA[patch_entry->byte_offset + column_offset];
                u8 y_pix_in_col = DATA[patch_entry->byte_offset + column_offset + 1];

                f32 y_loc = (f32)((y_hi - y_hi_capped) * sprite_size_y) / denom;
                for (int y = y_hi_capped; y >= y_lo_capped; y--)
                {
                    u32 texture_y = min((u32)(y_loc), sprite_size_y - 1);

                    if (texture_y > y_texture_skip + y_pix_in_col)
                    {
                        // Grab the next column, if possible.
                        column_offset += 4 + y_pix_in_col;
                        y_texture_skip = DATA[patch_entry->byte_offset + column_offset];
                        y_pix_in_col = DATA[patch_entry->byte_offset + column_offset + 1];
                    }
                    else if (texture_y > y_texture_skip)
                    {
                        u8 color_index = DATA[patch_entry->byte_offset + column_offset + 3 + texture_y - y_texture_skip]; // Index into the DOOM color palette.
                        u32 color = *(u32 *)(ASSETS_BINARY_BLOB2 + PALETTE_OFFSET + 3 * color_index);
                        color |= 0xFF000000; // Set alpha to full.
                        pixels[(y * SCREEN_SIZE_X) + x] = color;
                    }

                    y_loc += y_step;
                }
            }

            x_loc += x_step;
        }
    }
}

void Render(
    u32 *pixels,
    f32 *wall_raycast_radius,
    struct CameraState *camera,
    struct GameMap *game_map,
    QuarterEdge *qe_camera // dual quarter edge representing the face that the camera is in.
)
{
    u64 rtdsc_render_start = ReadCPUTimer();
    // RenderFloorAndCeiling(pixels, camera);
    u64 rtdsc_post_floor_and_ceiling = ReadCPUTimer();
    RenderWalls(pixels, wall_raycast_radius, camera, game_map, qe_camera);
    u64 rtdsc_post_walls = ReadCPUTimer();
    RenderObjects(pixels, wall_raycast_radius, camera);
    u64 rtdsc_post_objects = ReadCPUTimer();

    u64 dtimer_tot = rtdsc_post_objects - rtdsc_render_start;
    u64 dtimer_floor_and_ceiling = rtdsc_post_floor_and_ceiling - rtdsc_render_start;
    u64 dtimer_walls = rtdsc_post_walls - rtdsc_post_floor_and_ceiling;
    u64 dtimer_objects = rtdsc_post_objects - rtdsc_post_walls;

    printf("[Render]\n");
    printf("  RenderFloorAndCeiling: %zu (%.2f%%)\n", dtimer_floor_and_ceiling, (100.0 * dtimer_floor_and_ceiling) / dtimer_tot);
    printf("  RenderWalls:           %zu (%.2f%%)\n", dtimer_walls, (100.0 * dtimer_walls) / dtimer_tot);
    printf("  RenderObjects:         %zu (%.2f%%)\n\n", dtimer_objects, (100.0 * dtimer_objects) / dtimer_tot);
}

int main(int argc, char *argv[])
{
    // Estimate our clock frequency
    u64 cpu_timer_freq = EstimateCPUTimerFreq(100);

    // Create our game map
    struct GameMap game_map;

    // Load our assets
    printf("Loading assets.\n");
    LoadAssets(&game_map);
    printf("DONE.\n");

    // Initialize SDL
    ASSERT(
        SDL_Init(SDL_INIT_VIDEO) == 0,
        "SDL initialization failed: %s\n",
        SDL_GetError());

    // Create a window
    state.window = SDL_CreateWindow(
        "TOOM",
        SDL_WINDOWPOS_CENTERED_DISPLAY(1),
        SDL_WINDOWPOS_CENTERED_DISPLAY(1),
        SCREEN_SIZE_X,
        SCREEN_SIZE_Y,
        SDL_WINDOW_ALLOW_HIGHDPI);
    ASSERT(state.window, "Error creating SDL window: %s\n", SDL_GetError());

    // Create a second window for debug
    int debug_window_size_xy = SCREEN_SIZE_X;
    SDL_Window *debug_window = SDL_CreateWindow(
        "TOOM Debug",
        SDL_WINDOWPOS_UNDEFINED,
        SDL_WINDOWPOS_UNDEFINED,
        debug_window_size_xy,
        debug_window_size_xy,
        SDL_WINDOW_ALLOW_HIGHDPI);
    ASSERT(debug_window, "Error creating SDL Debug window: %s\n", SDL_GetError());

    // Create a renderer
    state.renderer = SDL_CreateRenderer(state.window, -1, SDL_RENDERER_PRESENTVSYNC);
    ASSERT(state.renderer, "Error creating SDL renderer: %s\n", SDL_GetError());

    // Create a second renderer for the debug window
    SDL_Renderer *debug_renderer = SDL_CreateRenderer(debug_window, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
    ASSERT(state.renderer, "Error creating SDL Debug renderer: %s\n", SDL_GetError());
    SDL_SetRenderDrawColor(debug_renderer, 0xFF, 0xFF, 0xFF, 0xFF); // Initialize renderer color

    // Create a texture
    state.texture = SDL_CreateTexture(
        state.renderer,
        SDL_PIXELFORMAT_ABGR8888,
        SDL_TEXTUREACCESS_STREAMING,
        SCREEN_SIZE_X,
        SCREEN_SIZE_Y);
    ASSERT(state.texture, "Error creating SDL texture: %s\n", SDL_GetError());

    // Init camera
    state.camera.fov.x = 1.5f;
    state.camera.fov.y = state.camera.fov.x * SCREEN_SIZE_Y / SCREEN_SIZE_X;

    // Init player state
    state.game_state.player.pos = (v2){5.0, 5.0};
    state.game_state.player.dir = (v2){0.0, 0.0};
    state.game_state.player.vel = (v2){0.0, 0.0};
    state.game_state.player.omega = 0.0f;
    state.game_state.player.height = 0.4;
    state.game_state.player.z = 0.4;
    state.game_state.player.qe_geometry = DelaunayMeshGetEnclosingTriangle2(game_map.geometry_mesh, &(state.game_state.player.pos));

    // Init keyboard
    ClearKeyboardState(&state.keyboard_state);

    // Time structs
    u64 rtdsc_tick = 0;
    u64 rtdsc_tick_prev = 0;

    // Main loop
    state.quit = 0;
    while (state.quit == 0)
    {
        u64 rtdsc_frame_start = ReadCPUTimer();

        SDL_Event event;
        while (SDL_PollEvent(&event))
        {
            if (event.type == SDL_QUIT)
            {
                state.quit = 1;
                break;
            }
            else if (event.type == SDL_WINDOWEVENT &&
                     event.window.event == SDL_WINDOWEVENT_CLOSE &&
                     (event.window.windowID == SDL_GetWindowID(state.window)))
            {
                state.quit = 1;
                break;
            }
            else if (event.type == SDL_KEYDOWN)
            {
                switch (event.key.keysym.sym)
                {
                case (SDLK_UP):
                    state.keyboard_state.up = KeyboardKeyState_Pressed;
                    break;
                case (SDLK_DOWN):
                    state.keyboard_state.down = KeyboardKeyState_Pressed;
                    break;
                case (SDLK_LEFT):
                    state.keyboard_state.left = KeyboardKeyState_Pressed;
                    break;
                case (SDLK_RIGHT):
                    state.keyboard_state.right = KeyboardKeyState_Pressed;
                    break;
                case (SDLK_a):
                    state.keyboard_state.a = KeyboardKeyState_Pressed;
                    break;
                case (SDLK_s):
                    state.keyboard_state.s = KeyboardKeyState_Pressed;
                    break;
                case (SDLK_d):
                    state.keyboard_state.d = KeyboardKeyState_Pressed;
                    break;
                case (SDLK_w):
                    state.keyboard_state.w = KeyboardKeyState_Pressed;
                    break;
                case (SDLK_q):
                    state.keyboard_state.q = KeyboardKeyState_Pressed;
                    break;
                case (SDLK_e):
                    state.keyboard_state.e = KeyboardKeyState_Pressed;
                    break;
                case (SDLK_r):
                    state.keyboard_state.r = KeyboardKeyState_Pressed;
                    break;
                case (SDLK_1):
                    state.keyboard_state.one = KeyboardKeyState_Pressed;
                    break;
                case (SDLK_2):
                    state.keyboard_state.two = KeyboardKeyState_Pressed;
                    break;
                case (SDLK_3):
                    state.keyboard_state.three = KeyboardKeyState_Pressed;
                    break;
                case (SDLK_4):
                    state.keyboard_state.four = KeyboardKeyState_Pressed;
                    break;
                case (SDLK_5):
                    state.keyboard_state.five = KeyboardKeyState_Pressed;
                    break;
                case (SDLK_6):
                    state.keyboard_state.six = KeyboardKeyState_Pressed;
                    break;
                case (SDLK_7):
                    state.keyboard_state.seven = KeyboardKeyState_Pressed;
                    break;
                case (SDLK_8):
                    state.keyboard_state.eight = KeyboardKeyState_Pressed;
                    break;
                }
            }
            else if (event.type == SDL_KEYUP)
            {
                switch (event.key.keysym.sym)
                {
                case (SDLK_UP):
                    state.keyboard_state.up = KeyboardKeyState_Released;
                    break;
                case (SDLK_DOWN):
                    state.keyboard_state.down = KeyboardKeyState_Released;
                    break;
                case (SDLK_LEFT):
                    state.keyboard_state.left = KeyboardKeyState_Released;
                    break;
                case (SDLK_RIGHT):
                    state.keyboard_state.right = KeyboardKeyState_Released;
                    break;
                case (SDLK_a):
                    state.keyboard_state.a = KeyboardKeyState_Released;
                    break;
                case (SDLK_s):
                    state.keyboard_state.s = KeyboardKeyState_Released;
                    break;
                case (SDLK_d):
                    state.keyboard_state.d = KeyboardKeyState_Released;
                    break;
                case (SDLK_w):
                    state.keyboard_state.w = KeyboardKeyState_Released;
                    break;
                case (SDLK_q):
                    state.keyboard_state.q = KeyboardKeyState_Released;
                    break;
                case (SDLK_e):
                    state.keyboard_state.e = KeyboardKeyState_Released;
                    break;
                case (SDLK_r):
                    state.keyboard_state.r = KeyboardKeyState_Released;
                    break;
                case (SDLK_1):
                    state.keyboard_state.one = KeyboardKeyState_Released;
                    break;
                case (SDLK_2):
                    state.keyboard_state.two = KeyboardKeyState_Released;
                    break;
                case (SDLK_3):
                    state.keyboard_state.three = KeyboardKeyState_Released;
                    break;
                case (SDLK_4):
                    state.keyboard_state.four = KeyboardKeyState_Released;
                    break;
                case (SDLK_5):
                    state.keyboard_state.five = KeyboardKeyState_Released;
                    break;
                case (SDLK_6):
                    state.keyboard_state.six = KeyboardKeyState_Released;
                    break;
                case (SDLK_7):
                    state.keyboard_state.seven = KeyboardKeyState_Released;
                    break;
                case (SDLK_8):
                    state.keyboard_state.eight = KeyboardKeyState_Released;
                    break;
                }
            }
        }
        u64 rtdsc_post_poll_events = ReadCPUTimer();

        // Calc elapsed time since previous tick, then run tick
        rtdsc_tick = ReadCPUTimer();
        const f32 dt = GetElapsedCPUTimeMs(rtdsc_tick_prev, rtdsc_tick, cpu_timer_freq) / 1000.0f;
        Tick(&state.game_state, dt, &state.keyboard_state, &game_map);
        rtdsc_tick_prev = rtdsc_tick;
        u64 rtdsc_post_tick = ReadCPUTimer();

        // Check for a keyboard press to reload our assets
        if (IsNewlyPressed(state.keyboard_state.r))
        {
            printf("Reloading assets.\n");
            LoadAssets(&game_map);
            printf("DONE.\n");
        }

        // Set camera to player state
        state.camera.pos = state.game_state.player.pos;
        state.camera.dir = state.game_state.player.dir;
        state.camera.z = state.game_state.player.z;

        Render(state.pixels, state.wall_raycast_radius, &state.camera, &game_map, state.game_state.player.qe_geometry);
        u64 rtdsc_post_render = ReadCPUTimer();

        DecayKeyboardState(&state.keyboard_state);
        u64 rtdsc_post_decay_keyboard_state = ReadCPUTimer();

        // Get timer end for all the non-SDL stuff
        const f64 game_ms_elapsed = GetElapsedCPUTimeMs(rtdsc_frame_start, rtdsc_post_decay_keyboard_state, cpu_timer_freq);

        u64 dtimer_tot = rtdsc_post_decay_keyboard_state - rtdsc_frame_start;
        u64 dtimer_poll = rtdsc_post_poll_events - rtdsc_frame_start;
        u64 dtimer_tick = rtdsc_post_tick - rtdsc_post_poll_events;
        u64 dtimer_render = rtdsc_post_render - rtdsc_post_tick;
        u64 dtimer_dkbs = rtdsc_post_decay_keyboard_state - rtdsc_post_render;

        printf("Frame dt:   %.4fms, %.1f fps\n", dt * 1000.0, 1000.0f / max(1.0f, game_ms_elapsed));
        printf("Total time: %.4fms (CPU freq %zu)\n", game_ms_elapsed, cpu_timer_freq);
        printf("  Poll:   %zu (%.2f%%)\n", dtimer_poll, (100.0 * dtimer_poll) / dtimer_tot);
        printf("  Tick:   %zu (%.2f%%)\n", dtimer_tick, (100.0 * dtimer_tick) / dtimer_tot);
        printf("  Render: %zu (%.2f%%)\n", dtimer_render, (100.0 * dtimer_render) / dtimer_tot);
        printf("  DKBS:   %zu (%.2f%%)\n\n", dtimer_dkbs, (100.0 * dtimer_dkbs) / dtimer_tot);

        SDL_UpdateTexture(state.texture, NULL, state.pixels, SCREEN_SIZE_X * 4);
        SDL_RenderCopyEx(
            state.renderer,
            state.texture,
            NULL,
            NULL,
            0.0,
            NULL,
            SDL_FLIP_VERTICAL);

        // SDL_RENDERER_PRESENTVSYNC means this is syncronized with the monitor refresh rate. (30Hz)
        SDL_RenderPresent(state.renderer);

        // Clear screen
        SDL_SetRenderDrawColor(debug_renderer, 0xFF, 0xFF, 0xFF, 0xFF);
        SDL_RenderClear(debug_renderer);

        {
            SDL_SetRenderDrawColor(debug_renderer, 0x1B, 0xA1, 0xEA, 0xFF);
            f32 pix_per_tile = 49.0;
            f32 offset_x = 1.5;
            f32 offset_y = 124.0;

            // Render the walls TODO

            { // Render the camera raycasts
                SDL_SetRenderDrawColor(debug_renderer, 0xF5, 0x61, 0x5C, 0xFF);
                struct CameraState *camera = &state.camera;
                int camera_x = camera->pos.x / TILE_WIDTH * pix_per_tile + offset_x;
                int camera_y = debug_window_size_xy - (camera->pos.y / TILE_WIDTH * pix_per_tile + offset_y);

                for (int x = 0; x < SCREEN_SIZE_X; x++)
                {
                    f32 r = state.wall_raycast_radius[x];

                    const f32 dw = camera->fov.x / 2 - (camera->fov.x * x) / SCREEN_SIZE_X;
                    const v2 cp = {
                        camera->dir.x - dw * camera->dir.y,
                        camera->dir.y + dw * camera->dir.x};
                    const f32 cam_len = length((cp));
                    const v2 dir = {cp.x / cam_len, cp.y / cam_len};

                    f32 camera_x1 = camera->pos.x + r * dir.x;
                    f32 camera_y1 = camera->pos.y + r * dir.y;
                    int camera_x2 = camera_x1 / TILE_WIDTH * pix_per_tile + offset_x;
                    int camera_y2 = debug_window_size_xy - (camera_y1 / TILE_WIDTH * pix_per_tile + offset_y);
                    SDL_RenderDrawLine(debug_renderer, camera_x, camera_y, camera_x2, camera_y2);
                }
            }

            { // Render the billboard sprite
                SDL_SetRenderDrawColor(debug_renderer, 0x00, 0x00, 0x00, 0xFF);
                struct CameraState *camera = &state.camera;

                v2 sprite_pos = {10.0f, 4.5f};
                f32 sprite_halfwidth = 0.5;
                const v2 sprite_tangent = rotr(camera->dir);

                const v2 a = {
                    sprite_pos.x + sprite_halfwidth * sprite_tangent.x,
                    sprite_pos.y + sprite_halfwidth * sprite_tangent.y};
                const v2 b = {
                    sprite_pos.x - sprite_halfwidth * sprite_tangent.x,
                    sprite_pos.y - sprite_halfwidth * sprite_tangent.y};
                int ax = a.x / TILE_WIDTH * pix_per_tile + offset_x;
                int ay = debug_window_size_xy - (a.y / TILE_WIDTH * pix_per_tile + offset_y);
                int bx = b.x / TILE_WIDTH * pix_per_tile + offset_x;
                int by = debug_window_size_xy - (b.y / TILE_WIDTH * pix_per_tile + offset_y);
                SDL_RenderDrawLine(debug_renderer, ax, ay, bx, by);
            }

            { // Render the mesh
                for (int qe_index = 0; qe_index < DelaunayMeshNumQuarterEdges(game_map.geometry_mesh); qe_index++)
                {
                    QuarterEdge *qe = DelaunayMeshGetQuarterEdge(game_map.geometry_mesh, qe_index);

                    if (IsPrimalEdge(qe) && !DelaunayMeshIsBoundaryVertex(game_map.geometry_mesh, qe->vertex))
                    {
                        SDL_SetRenderDrawColor(debug_renderer, 0xFF, 0x48, 0xCF, 0xFF);
                        if (game_map.quarter_edge_index_to_side_info_index[qe->index] != 0xFFFF)
                        {
                            // is solid
                            SDL_SetRenderDrawColor(debug_renderer, 0xFF, 0xCF, 0x48, 0xFF);
                        }

                        // Get its opposite side.
                        QuarterEdge *qe_sym = QESym(qe);

                        const v2 *a = qe->vertex;
                        const v2 *b = qe_sym->vertex;
                        if (a > b && !DelaunayMeshIsBoundaryVertex(game_map.geometry_mesh, b))
                        { // Avoid rendering edges twice
                            int ax = a->x / TILE_WIDTH * pix_per_tile + offset_x;
                            int ay = debug_window_size_xy - (a->y / TILE_WIDTH * pix_per_tile + offset_y);
                            int bx = b->x / TILE_WIDTH * pix_per_tile + offset_x;
                            int by = debug_window_size_xy - (b->y / TILE_WIDTH * pix_per_tile + offset_y);
                            SDL_RenderDrawLine(debug_renderer, ax, ay, bx, by);
                        }
                    }
                }
            }

            { // Render the player-enclosing triangle
                SDL_SetRenderDrawColor(debug_renderer, 0x48, 0x48, 0xCF, 0xFF);

                // The quarter edge is a dual edge, and its containing triangle is solid.
                const v2 *a = DelaunayMeshGetTriangleVertex1(game_map.geometry_mesh, state.game_state.player.qe_geometry);
                const v2 *b = DelaunayMeshGetTriangleVertex2(game_map.geometry_mesh, state.game_state.player.qe_geometry);
                const v2 *c = DelaunayMeshGetTriangleVertex3(game_map.geometry_mesh, state.game_state.player.qe_geometry);

                int ax = a->x / TILE_WIDTH * pix_per_tile + offset_x;
                int ay = debug_window_size_xy - (a->y / TILE_WIDTH * pix_per_tile + offset_y);
                int bx = b->x / TILE_WIDTH * pix_per_tile + offset_x;
                int by = debug_window_size_xy - (b->y / TILE_WIDTH * pix_per_tile + offset_y);
                int cx = c->x / TILE_WIDTH * pix_per_tile + offset_x;
                int cy = debug_window_size_xy - (c->y / TILE_WIDTH * pix_per_tile + offset_y);
                SDL_RenderDrawLine(debug_renderer, ax, ay, bx, by);
                SDL_RenderDrawLine(debug_renderer, bx, by, cx, cy);
                SDL_RenderDrawLine(debug_renderer, cx, cy, ax, ay);
            }
        }

        SDL_RenderPresent(debug_renderer);
    }

    SDL_DestroyWindow(state.window);
    SDL_DestroyWindow(debug_window);

    // Free our game state
    free(game_map.geometry_mesh->quarter_edges);
    free(game_map.geometry_mesh);

    // Free our assets
    free(ASSETS_BINARY_BLOB);
    free(ASSETS_BINARY_BLOB2);
    free(PATCH_OFFSETS);

    return 0;
}