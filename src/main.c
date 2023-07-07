#include <stdio.h>
#include <SDL2/SDL.h>

#include "typedefs.h"
#include "constants.h"
#include "vec.h"
#include "input.h"
#include "bitmap.h"
#include "delaunay_mesh.h"
#include "game.h"
#include "platform_metrics.h"

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

struct Mapdata
{
    u32 n_tiles;
    u32 n_tiles_x;
    u32 n_tiles_y;
    u8 *tiles; // Each tile is either solid (0) or has a texture index (>0)
    u8 *floor;
    u8 *ceiling;
};

static inline u8 GetMapDataIndex(struct Mapdata *mapdata, int x, int y)
{
    return (mapdata->n_tiles_y - y - 1) * (mapdata->n_tiles_x) + x;
}

// The global mapdata. This just points into our asset blob.
struct Mapdata MAPDATA;

// The bitmap global variables. These just point into the binary blob.
struct Bitmap BITMAP;

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

u8 *WAD = NULL;
u32 WAD_SIZE = 0; // number of bytes

struct WadDirectoryEntry
{
    u32 byte_offset;
    u32 size;
    char name[8];
};

//  Each palette in the PLAYPAL lump contains 256 colors totaling 768 bytes,
//  where each color is broken into three unsigned bytes. Each of these color components (red, green, and blue) range between 0 and 255.
u32 PALETTE_OFFSET = 0;

struct Patch
{
    u16 size_x;   // width of the graphic in pixels
    u16 size_y;   // height of the graphic in pixels
    i16 offset_x; // x offset from the screen origin (to left)
    i16 offset_y; // y offset from the screen origin (down)
    u32 column_offsets[64];
};

struct PatchEntry
{
    u32 byte_offset;
    struct Patch *patch;
};
struct PatchEntry CYBR_PATCH_ENTRIES[8];
struct PatchEntry SKEL_PATCH_ENTRIES[360];

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

#define SIDEINFO_FLAG_DARK 1

struct SideInfo
{
    u16 flags;
    u16 texture_id;
    i16 x_offset;
    i16 y_offset;
};

struct GameMap
{
    struct DelaunayMesh *geometry_mesh;

    u32 n_side_infos;
    struct SideInfo *side_infos;
    u16 *quarter_edge_index_to_side_info_index;
};

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

        bool loaded_textures = 0;
        bool loaded_mapdata = 0;

        // Scan through them in reverse order
        for (int i = n_toc_entries; i > 0; i--)
        {
            byte_index -= sizeof(struct BinaryAssetTableOfContentEntry);
            struct BinaryAssetTableOfContentEntry *entry = (struct BinaryAssetTableOfContentEntry *)(ASSETS_BINARY_BLOB + byte_index);
            // Make the name null-terminated just in case.
            entry->name[15] = '\0';
            printf("Entry %d: %s at offset %d\n", i, entry->name, entry->byte_offset);

            // In the future, load them into a map or something. For now, we're specifically looking for either the wall or floor textures.
            if (strcmp(entry->name, "textures") == 0)
            {
                u32 asset_byte_offset = entry->byte_offset;
                BITMAP.n_pixels = *(u32 *)(ASSETS_BINARY_BLOB + asset_byte_offset);
                asset_byte_offset += sizeof(u32);
                BITMAP.n_pixels_per_column = *(u32 *)(ASSETS_BINARY_BLOB + asset_byte_offset);
                asset_byte_offset += sizeof(u32);
                BITMAP.column_major = ASSETS_BINARY_BLOB[asset_byte_offset];
                ASSERT(BITMAP.column_major, "Expected the wall texture to be column-major\n");
                asset_byte_offset += sizeof(u8);
                BITMAP.abgr = (u32 *)(ASSETS_BINARY_BLOB + asset_byte_offset);
                BITMAP.n_pixels_per_row = BITMAP.n_pixels / BITMAP.n_pixels_per_column;
                loaded_textures = 1;
            }
            else if (strcmp(entry->name, "mapdata") == 0)
            {
                u32 asset_byte_offset = entry->byte_offset;
                MAPDATA.n_tiles = *(u32 *)(ASSETS_BINARY_BLOB + asset_byte_offset);
                asset_byte_offset += sizeof(u32);
                MAPDATA.n_tiles_x = *(u32 *)(ASSETS_BINARY_BLOB + asset_byte_offset);
                asset_byte_offset += sizeof(u32);
                MAPDATA.tiles = (u8 *)(ASSETS_BINARY_BLOB + asset_byte_offset);
                asset_byte_offset += MAPDATA.n_tiles * sizeof(u8);
                MAPDATA.floor = (u8 *)(ASSETS_BINARY_BLOB + asset_byte_offset);
                asset_byte_offset += MAPDATA.n_tiles * sizeof(u8);
                MAPDATA.ceiling = (u8 *)(ASSETS_BINARY_BLOB + asset_byte_offset);
                MAPDATA.n_tiles_y = MAPDATA.n_tiles / MAPDATA.n_tiles_x;
                loaded_mapdata = 1;
            }
            else if (strncmp(entry->name, "SKEL", 4) == 0)
            {
                int frame_index = entry->name[4] - '0';
                frame_index = 10 * frame_index + (entry->name[5] - '0');
                frame_index = 10 * frame_index + (entry->name[6] - '0');
                SKEL_PATCH_ENTRIES[frame_index].byte_offset = entry->byte_offset;
                SKEL_PATCH_ENTRIES[frame_index].patch = (struct Patch *)(ASSETS_BINARY_BLOB + entry->byte_offset);
            }
        }

        ASSERT(loaded_textures > 0, "Textures not loaded from assets\n");
        ASSERT(loaded_mapdata > 0, "Map data not loaded from assets\n");
    }

    // Load DOOM assets.
    if (WAD)
    {
        free(WAD);
    }
    {
        FILE *fileptr = fopen("assets/DOOM.WAD", "rb");
        ASSERT(fileptr, "Error opening DOOM WAD\n");

        // Count the number of bytes
        fseek(fileptr, 0, SEEK_END);
        WAD_SIZE = ftell(fileptr);

        // Read in the WAD file
        fseek(fileptr, 0, SEEK_SET);
        WAD = (u8 *)malloc(WAD_SIZE);
        ASSERT(WAD, "Failed to allocate DOOM WAD\n");
        ASSERT(fread(WAD, sizeof(u8), WAD_SIZE, fileptr) == WAD_SIZE, "Failed to read DOOM WAD\n");

        u32 n_lumps = *(u32 *)(WAD + 0x04);
        u32 dir_loc = *(u32 *)(WAD + 0x08);

        // Process the directory
        bool loaded_patch = 0;
        u32 byte_index = dir_loc;
        for (u32 directory_index = 0; directory_index < n_lumps; directory_index++)
        {
            struct WadDirectoryEntry *entry = (struct WadDirectoryEntry *)(WAD + byte_index);
            // printf("Entry %d: %.8s at offset %d with size %d\n", directory_index, entry->name, entry->byte_offset, entry->size);
            byte_index += sizeof(struct WadDirectoryEntry);

            // if (strcmp(entry->name, "PLAYPAL") == 0)
            // {
            //     PALETTE_OFFSET = entry->byte_offset;
            // }
            // else
            if (strncmp(entry->name, "CYBRE", 5) == 0)
            {
                int frame_index = entry->name[5] - '1';
                CYBR_PATCH_ENTRIES[frame_index].byte_offset = entry->byte_offset;
                CYBR_PATCH_ENTRIES[frame_index].patch = (struct Patch *)(WAD + entry->byte_offset);
                loaded_patch = 1;
            }
        }

        fclose(fileptr);

        ASSERT(loaded_patch > 0, "Patch not loaded from assets\n");
    }

    // Load assets #2
    if (ASSETS_BINARY_BLOB2)
    {
        free(ASSETS_BINARY_BLOB2);
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

            if (strcmp(entry->name, "palette") == 0)
            {
                PALETTE_OFFSET = entry->offset;
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
        }

        fclose(fileptr);
    }
}

void RenderFloorAndCeiling(
    u32 *pixels,
    struct CameraState *camera)
{
    // Ray direction for x = 0
    f32 half_camera_width = camera->fov.x / 2.0f;
    f32 ray_dir_lo_x = camera->dir.x - half_camera_width * camera->dir.y;
    f32 ray_dir_lo_y = camera->dir.y + half_camera_width * camera->dir.x;

    // Ray direction for x = SCREEN_SIZE_X
    f32 ray_dir_hi_x = camera->dir.x + half_camera_width * camera->dir.y;
    f32 ray_dir_hi_y = camera->dir.y - half_camera_width * camera->dir.x;

    // Draw floor
    for (int y = 0; y < SCREEN_SIZE_Y / 2; y++)
    {
        // Radius
        f32 zpp = (SCREEN_SIZE_Y / 2.0f - y) * (camera->fov.y / SCREEN_SIZE_Y);
        f32 radius = camera->z / zpp; // TODO: Precompute for each y

        // Location of the 1st ray's intersection
        f32 hit_x = camera->pos.x + radius * ray_dir_lo_x;
        f32 hit_y = camera->pos.y + radius * ray_dir_lo_y;

        // Each step is (hit_x2 - hit_x) / SCREEN_SIZE_X;
        // = ((camera->pos.x + radius * ray_dir_lo_x) - (camera->pos.x + radius * ray_dir_lo_x)) / SCREEN_SIZE_X
        // = (radius * ray_dir_lo_x - (radius * ray_dir_lo_x)) / SCREEN_SIZE_X
        // = radius * (ray_dir_hi_x - ray_dir_lo_x) / SCREEN_SIZE_X
        f32 step_x = radius * (ray_dir_hi_x - ray_dir_lo_x) / SCREEN_SIZE_X;
        f32 step_y = radius * (ray_dir_hi_y - ray_dir_lo_y) / SCREEN_SIZE_X;

        for (int x = 0; x < SCREEN_SIZE_X; x++)
        {

            int x_ind_hit = (int)(floorf(hit_x / TILE_WIDTH));
            int y_ind_hit = (int)(floorf(hit_y / TILE_WIDTH));
            f32 x_rem_hit = hit_x - TILE_WIDTH * x_ind_hit;
            f32 y_rem_hit = hit_y - TILE_WIDTH * y_ind_hit;
            x_ind_hit = clamp(x_ind_hit, 0, MAPDATA.n_tiles_x - 1);
            y_ind_hit = clamp(y_ind_hit, 0, MAPDATA.n_tiles_y - 1);

            u32 texture_x_offset = 0;
            u32 texture_y_offset = (MAPDATA.floor[GetMapDataIndex(&MAPDATA, x_ind_hit, y_ind_hit)] - 1) * TEXTURE_SIZE;

            u32 texture_x = (int)(x_rem_hit / TILE_WIDTH * TEXTURE_SIZE);
            u32 texture_y = (int)(y_rem_hit / TILE_WIDTH * TEXTURE_SIZE);

            u32 color = GetColumnMajorPixelAt(&BITMAP, texture_x + texture_x_offset, texture_y + texture_y_offset);
            pixels[(y * SCREEN_SIZE_X) + x] = color;

            // step
            hit_x += step_x;
            hit_y += step_y;
        }
    }

    // Draw ceiling
    for (int y = SCREEN_SIZE_Y / 2 + 1; y < SCREEN_SIZE_Y; y++)
    {
        // Radius
        f32 zpp = (y - (SCREEN_SIZE_Y / 2.0f)) * (camera->fov.y / SCREEN_SIZE_Y);
        f32 radius = (WALL_HEIGHT - camera->z) / zpp; // TODO: Precompute for each y

        // Location of the 1st ray's intersection
        f32 hit_x = camera->pos.x + radius * ray_dir_lo_x;
        f32 hit_y = camera->pos.y + radius * ray_dir_lo_y;

        // Each step toward hit2
        f32 step_x = radius * (ray_dir_hi_x - ray_dir_lo_x) / SCREEN_SIZE_X;
        f32 step_y = radius * (ray_dir_hi_y - ray_dir_lo_y) / SCREEN_SIZE_X;

        for (int x = 0; x < SCREEN_SIZE_X; x++)
        {
            int x_ind_hit = (int)(floorf(hit_x / TILE_WIDTH));
            int y_ind_hit = (int)(floorf(hit_y / TILE_WIDTH));
            f32 x_rem_hit = hit_x - TILE_WIDTH * x_ind_hit;
            f32 y_rem_hit = hit_y - TILE_WIDTH * y_ind_hit;
            x_ind_hit = clamp(x_ind_hit, 0, MAPDATA.n_tiles_x - 1);
            y_ind_hit = clamp(y_ind_hit, 0, MAPDATA.n_tiles_y - 1);

            u32 texture_x_offset = 0;
            u32 texture_y_offset = (MAPDATA.ceiling[GetMapDataIndex(&MAPDATA, x_ind_hit, y_ind_hit)] - 1) * TEXTURE_SIZE;

            u32 texture_x = (int)(x_rem_hit / TILE_WIDTH * TEXTURE_SIZE);
            u32 texture_y = (int)(y_rem_hit / TILE_WIDTH * TEXTURE_SIZE);

            u32 color = GetColumnMajorPixelAt(&BITMAP, texture_x + texture_x_offset, texture_y + texture_y_offset);
            pixels[(y * SCREEN_SIZE_X) + x] = color;

            // step
            hit_x += step_x;
            hit_y += step_y;
        }
    }
}

void RenderWalls(
    u32 *pixels,
    f32 *wall_raycast_radius,
    struct CameraState *camera)
{
    // Get camera location's cell coordinates
    int x_ind_cam = (int)(floorf(camera->pos.x / TILE_WIDTH));
    int y_ind_cam = (int)(floorf(camera->pos.y / TILE_WIDTH));
    f32 x_rem_cam = camera->pos.x - TILE_WIDTH * x_ind_cam;
    f32 y_rem_cam = camera->pos.y - TILE_WIDTH * y_ind_cam;

    for (int x = SCREEN_SIZE_X / 2; x < SCREEN_SIZE_X; x++)
    {

        // Camera to pixel column
        const f32 dw = camera->fov.x / 2 - (camera->fov.x * x) / SCREEN_SIZE_X;
        const v2 cp = {
            camera->dir.x - dw * camera->dir.y,
            camera->dir.y + dw * camera->dir.x};

        // Distance from the camera to the column
        const f32 cam_len = length((cp));

        // Ray direction through this column
        const v2 dir = {cp.x / cam_len, cp.y / cam_len};

        // Start at the camera pos
        int x_ind = x_ind_cam;
        int y_ind = y_ind_cam;
        f32 x_rem = x_rem_cam;
        f32 y_rem = y_rem_cam;

        // We will be raycasting through cells of unit width.
        // Our ray's position vs time is:
        // x(t) = x_rem + dir.x * dt
        // y(t) = y_rem + dir.y * dt

        // We cross x = 0          if dir.x < 0, at dt = -x_rem/dir.x
        // We cross x = TILE_WIDTH if dir.x > 0, at dt = (TILE_WIDTH-x_rem)/dir.x
        // We cross y = 0          if dir.y < 0, at dt = -y_rem/dir.y
        // We cross y = TILE_WIDTH if dir.y > 0, at dt = (TILE_WIDTH-y_rem)/dir.y

        // We can generalize this to:
        //   dx_ind_dir = -1 if dir.x < 0, at dt = -1/dir.x * x_rem + 0.0
        //   dx_ind_dir =  1 if dir.x > 0, at dt = -1/dir.x * x_rem + TILE_WIDTH/dir.x
        //   dx_ind_dir =  0 if dir.x = 0, at dt =        0 * x_rem + INFINITY
        //   dy_ind_dir = -1 if dir.y < 0, at dt = -1/dir.y * y_rem + 0.0
        //   dy_ind_dir =  1 if dir.y > 0, at dt = -1/dir.y * y_rem + TILE_WIDTH/dir.y
        //   dy_ind_dir =  0 if dir.x = 0, at dt =        0 * y_rem + INFINITY

        int dx_ind_dir = 0;
        f32 dx_a = 0.0;
        f32 dx_b = INFINITY;
        if (dir.x < 0)
        {
            dx_ind_dir = -1;
            dx_a = -1.0f / dir.x;
            dx_b = 0.0;
        }
        else if (dir.x > 0)
        {
            dx_ind_dir = 1;
            dx_a = -1.0f / dir.x;
            dx_b = TILE_WIDTH / dir.x;
        }

        int dy_ind_dir = 0;
        f32 dy_a = 0.0;
        f32 dy_b = INFINITY;
        if (dir.y < 0)
        {
            dy_ind_dir = -1;
            dy_a = -1.0f / dir.y;
            dy_b = 0.0;
        }
        else if (dir.y > 0)
        {
            dy_ind_dir = 1;
            dy_a = -1.0f / dir.y;
            dy_b = TILE_WIDTH / dir.y;
        }

        // Step through cells until we hit an occupied cell
        int n_steps = 0;
        int dx_ind, dy_ind;
        while (n_steps < 100)
        {
            n_steps += 1;

            f32 dt_best = INFINITY;
            dx_ind = 0;
            dy_ind = 0;

            f32 dt_x = dx_a * x_rem + dx_b;
            f32 dt_y = dy_a * y_rem + dy_b;
            if (dt_x < dt_y)
            {
                dt_best = dt_x;
                dx_ind = dx_ind_dir;
                dy_ind = 0;
            }
            else
            {
                dt_best = dt_y;
                dx_ind = 0;
                dy_ind = dy_ind_dir;
            }

            // Move up to the next cell
            x_ind += dx_ind;
            y_ind += dy_ind;
            x_rem += dir.x * dt_best - TILE_WIDTH * dx_ind;
            y_rem += dir.y * dt_best - TILE_WIDTH * dy_ind;

            // Check to see if the new cell is solid
            if (MAPDATA.tiles[GetMapDataIndex(&MAPDATA, x_ind, y_ind)] > 0)
            {
                break;
            }
        }

        // Calculate the collision location
        const v2 collision = {
            TILE_WIDTH * x_ind + x_rem,
            TILE_WIDTH * y_ind + y_rem};

        // Calculate the ray length
        const f32 ray_len = max(length(sub(collision, camera->pos)), 0.01); // TODO: Remove this `max` once we have proper collision
        wall_raycast_radius[x] = ray_len;

        // Calculate the pixel bounds that we fill the wall in for
        int y_lo = (int)(SCREEN_SIZE_Y / 2.0f - cam_len * camera->z / ray_len * SCREEN_SIZE_Y / camera->fov.y);
        int y_hi = (int)(SCREEN_SIZE_Y / 2.0f + cam_len * (WALL_HEIGHT - camera->z) / ray_len * SCREEN_SIZE_Y / camera->fov.y);
        int y_lo_capped = max(y_lo, 0);
        int y_hi_capped = min(y_hi, SCREEN_SIZE_Y - 1);

        // Texture x offset determines whether we draw the light or dark version
        u32 texture_x_offset = dx_ind == 0 ? 0 : TEXTURE_SIZE;
        u32 texture_y_offset = (MAPDATA.tiles[GetMapDataIndex(&MAPDATA, x_ind, y_ind)] - 1) * TEXTURE_SIZE;

        f32 rem = 0.0f;
        if (dx_ind == 0)
        {
            rem = dy_ind < 0 ? TILE_WIDTH - x_rem : x_rem;
        }
        else
        {
            rem = dx_ind < 0 ? y_rem : TILE_WIDTH - y_rem;
        }
        u32 texture_x = min((int)(TEXTURE_SIZE * rem / TILE_WIDTH), TEXTURE_SIZE - 1);
        u32 baseline = GetColumnMajorPixelIndex(&BITMAP, texture_x + texture_x_offset, texture_y_offset);
        u32 denom = max(1, y_hi - y_lo);
        f32 y_loc = (f32)((y_hi - y_hi_capped) * TEXTURE_SIZE) / denom;
        f32 y_step = (f32)(TEXTURE_SIZE) / denom;
        for (int y = y_hi_capped; y >= y_lo_capped; y--)
        {
            u32 texture_y = min((u32)(y_loc), TEXTURE_SIZE - 1);
            u32 color = BITMAP.abgr[texture_y + baseline];
            pixels[(y * SCREEN_SIZE_X) + x] = color;
            y_loc += y_step;
        }
    }
}

void RenderWallsViaMesh(
    u32 *pixels,
    f32 *wall_raycast_radius,
    struct CameraState *camera,
    struct GameMap *game_map,
    QuarterEdge *qe_camera // dual quarter edge representing the face that the camera is in.
)
{
    for (int x = 0; x < SCREEN_SIZE_X / 2; x++)
    {

        // Camera to pixel column
        const f32 dw = camera->fov.x / 2 - (camera->fov.x * x) / SCREEN_SIZE_X; // TODO: Precompute once.
        const v2 cp = {
            camera->dir.x - dw * camera->dir.y,
            camera->dir.y + dw * camera->dir.x};

        // Distance from the camera to the column
        const f32 cam_len = length((cp));

        // Ray direction through this column
        const v2 dir = {cp.x / cam_len, cp.y / cam_len};

        // Start at the camera pos
        QuarterEdge *qe_dual = qe_camera;
        v2 pos = camera->pos;

        // The edge vector of the face that we last crossed
        v2 v_face = {0.0, 0.0};

        // The side info that we eventually collide with
        u16 side_info_index = 0xFFFF;

        // Step through triangles until we hit a solid triangle
        int n_steps = 0;
        while (n_steps < 100)
        {
            n_steps += 1;

            // Grab the enclosing triangle.
            QuarterEdge *qe_ab = qe_dual->rot;
            QuarterEdge *qe_bc = qe_dual->next->rot;
            QuarterEdge *qe_ca = qe_dual->next->next->rot;

            const v2 a = *(qe_ab->vertex);
            const v2 b = *(qe_bc->vertex);
            const v2 c = *(qe_ca->vertex);

            // Project our ray out far enough that it would exit our mesh
            f32 projection_distance = 100.0; // Ridiculously large
            v2 pos_next_delta = {
                projection_distance * dir.x,
                projection_distance * dir.y};
            v2 pos_next = add(pos, pos_next_delta);

            f32 min_interp = INFINITY;
            QuarterEdge *qe_side = NULL;
            QuarterEdge *qe_dual_next = NULL;

            // See if we cross any of the 3 faces for the triangle we are in,
            // And cross the first segment.
            if (GetRightHandedness(&a, &b, &pos_next) < -1e-4)
            {
                // We would cross AB
                v2 v = sub(b, a);
                v2 w = sub(pos, a);
                float interp_ab = cross(v, w) / cross(pos_next_delta, v);
                if (interp_ab < min_interp)
                {
                    min_interp = interp_ab;
                    qe_dual_next = qe_ab->rot;
                    qe_side = qe_ab;
                    v_face = v;
                }
            }
            if (GetRightHandedness(&b, &c, &pos_next) < -1e-4)
            {
                // We would cross BC
                v2 v = sub(c, b);
                v2 w = sub(pos, b);
                float interp_bc = cross(v, w) / cross(pos_next_delta, v);
                if (interp_bc < min_interp)
                {
                    min_interp = interp_bc;
                    qe_dual_next = qe_bc->rot;
                    qe_side = qe_bc;
                    v_face = v;
                }
            }
            if (GetRightHandedness(&c, &a, &pos_next) < -1e-4)
            {
                // We would cross CA
                v2 v = sub(a, c);
                v2 w = sub(pos, c);
                float interp_ca = cross(v, w) / cross(pos_next_delta, v);
                if (interp_ca < min_interp)
                {
                    min_interp = interp_ca;
                    qe_dual_next = qe_ca->rot;
                    qe_side = qe_ca;
                    v_face = v;
                }
            }

            // Move to the face.
            if (qe_dual_next != NULL)
            {
                // Should always be non-null.
                pos.x += min_interp * pos_next_delta.x;
                pos.y += min_interp * pos_next_delta.y;
                qe_dual = qe_dual_next;

                // TODO - reenable.
                // Render the ceiling above the threshold.
                // Calculate the ray length
                // const f32 ray_len = max(length(sub(pos, camera->pos)), 0.01); // TODO: Remove this `max` once we have proper collision

                // // Calculate the pixel bounds that we fill the wall in for
                // int y_hi = (int)(SCREEN_SIZE_Y / 2.0f + cam_len * (WALL_HEIGHT - camera->z) / ray_len * SCREEN_SIZE_Y / camera->fov.y);
                // for (int y = y_hi + 1; y < SCREEN_SIZE_Y; y++)
                // {
                //     // Radius
                //     f32 zpp = (y - (SCREEN_SIZE_Y / 2.0f)) * (camera->fov.y / SCREEN_SIZE_Y);
                //     f32 radius = (WALL_HEIGHT - camera->z) / zpp; // TODO: Precompute for each y

                //     // Location of the ray's intersection
                //     f32 hit_x = camera->pos.x + radius * cp.x;
                //     f32 hit_y = camera->pos.y + radius * cp.y;

                //     int x_ind_hit = (int)(floorf(hit_x / TILE_WIDTH));
                //     int y_ind_hit = (int)(floorf(hit_y / TILE_WIDTH));
                //     f32 x_rem_hit = hit_x - TILE_WIDTH * x_ind_hit;
                //     f32 y_rem_hit = hit_y - TILE_WIDTH * y_ind_hit;
                //     x_ind_hit = clamp(x_ind_hit, 0, MAPDATA.n_tiles_x - 1);
                //     y_ind_hit = clamp(y_ind_hit, 0, MAPDATA.n_tiles_y - 1);

                //     // TODO: Get the texture from the triangle data.
                //     u32 texture_x_offset = 0;
                //     u32 texture_y_offset = (MAPDATA.ceiling[GetMapDataIndex(&MAPDATA, x_ind_hit, y_ind_hit)] - 1) * TEXTURE_SIZE;

                //     u32 texture_x = (int)(x_rem_hit / TILE_WIDTH * TEXTURE_SIZE);
                //     u32 texture_y = (int)(y_rem_hit / TILE_WIDTH * TEXTURE_SIZE);

                //     u32 color = GetColumnMajorPixelAt(&BITMAP, texture_x + texture_x_offset, texture_y + texture_y_offset);
                //     pixels[(y * SCREEN_SIZE_X) + x] = color;
                // }

                // TODO: Render floor

                side_info_index = game_map->quarter_edge_index_to_side_info_index[qe_side->index];
                if (side_info_index != 0xFFFF)
                {
                    // game_map->side_infos[side_info_index].flags
                    // The side info is solid.
                    break;
                }
            }
            else
            {
                break;
            }
        }

        // Calculate the ray length
        const f32 ray_len = max(length(sub(pos, camera->pos)), 0.01); // TODO: Remove this `max` once we have proper collision
        wall_raycast_radius[x] = ray_len;

        // Calculate the pixel bounds that we fill the wall in for
        int y_lo = (int)(SCREEN_SIZE_Y / 2.0f - cam_len * camera->z / ray_len * SCREEN_SIZE_Y / camera->fov.y);
        int y_hi = (int)(SCREEN_SIZE_Y / 2.0f + cam_len * (WALL_HEIGHT - camera->z) / ray_len * SCREEN_SIZE_Y / camera->fov.y);
        int y_lo_capped = max(y_lo, 0);
        int y_hi_capped = min(y_hi, SCREEN_SIZE_Y - 1);

        // Texture x offset determines whether we draw the light or dark version
        QuarterEdge *qe_face_src = qe_dual->rot->rot->rot;
        u32 texture_x_offset = 0; // ((face_data[qe_face_src->index].flags & FACEDATA_FLAG_DARK) > 0) ? TEXTURE_SIZE : 0;
        u32 texture_y_offset = 0; // face_data[qe_face_src->index].texture_id * TEXTURE_SIZE;
        if (side_info_index != 0xFFFF)
        {
            struct SideInfo *side_info = game_map->side_infos + side_info_index;
            texture_x_offset = (side_info->flags & SIDEINFO_FLAG_DARK) > 0 ? TEXTURE_SIZE : 0;
            texture_y_offset = side_info->texture_id * TEXTURE_SIZE;
        }

        // Calculate where along the segment we intersected.
        f32 PIX_PER_DISTANCE = TEXTURE_SIZE / TILE_WIDTH;
        f32 x_along_texture = length(v_face) - length(sub(pos, *(qe_face_src->vertex)));

        u32 texture_x = (int)(PIX_PER_DISTANCE * x_along_texture) % TEXTURE_SIZE;
        u32 baseline = GetColumnMajorPixelIndex(&BITMAP, texture_x + texture_x_offset, texture_y_offset);
        u32 denom = max(1, y_hi - y_lo);
        f32 y_loc = (f32)((y_hi - y_hi_capped) * TEXTURE_SIZE) / denom;
        f32 y_step = (f32)(TEXTURE_SIZE) / denom;
        for (int y = y_hi_capped; y >= y_lo_capped; y--)
        {
            u32 texture_y = min((u32)(y_loc), TEXTURE_SIZE - 1);
            u32 color = BITMAP.abgr[texture_y + baseline];
            pixels[(y * SCREEN_SIZE_X) + x] = color;
            y_loc += y_step;
        }
    }
}

void RenderObjects(
    u32 *pixels,
    f32 *wall_raycast_radius,
    struct CameraState *camera)
{
    static f32 DOOM_HEIGHT_PER_PIX = ((f32)(DOOM_PIX_PER_WALL_WIDTH)) / (TEXTURE_SIZE * TEXTURE_SIZE);

    // u8* DATA = WAD;
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
                        u32 color = *(u32 *)(WAD + PALETTE_OFFSET + 3 * color_index);
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
    RenderFloorAndCeiling(pixels, camera);
    u64 rtdsc_post_floor_and_ceiling = ReadCPUTimer();
    RenderWalls(pixels, wall_raycast_radius, camera);
    u64 rtdsc_post_walls = ReadCPUTimer();
    RenderWallsViaMesh(pixels, wall_raycast_radius, camera, game_map, qe_camera);
    u64 rtdsc_post_walls_via_mesh = ReadCPUTimer();
    RenderObjects(pixels, wall_raycast_radius, camera);
    u64 rtdsc_post_objects = ReadCPUTimer();

    u64 dtimer_tot = rtdsc_post_objects - rtdsc_render_start;
    u64 dtimer_floor_and_ceiling = rtdsc_post_floor_and_ceiling - rtdsc_render_start;
    u64 dtimer_walls = rtdsc_post_walls - rtdsc_post_floor_and_ceiling;
    u64 dtimer_walls_via_mesh = rtdsc_post_walls_via_mesh - rtdsc_post_walls;
    u64 dtimer_objects = rtdsc_post_objects - rtdsc_post_walls_via_mesh;

    printf("[Render]\n");
    printf("  RenderFloorAndCeiling: %zu (%.2f%%)\n", dtimer_floor_and_ceiling, (100.0 * dtimer_floor_and_ceiling) / dtimer_tot);
    printf("  RenderWalls:           %zu (%.2f%%)\n", dtimer_walls, (100.0 * dtimer_walls) / dtimer_tot);
    printf("  RenderWallsViaMesh:    %zu (%.2f%%)\n", dtimer_walls_via_mesh, (100.0 * dtimer_walls_via_mesh) / dtimer_tot);
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

    // TODO
    state.game_state.geometry_mesh = game_map.geometry_mesh;

    // Init player state
    state.game_state.player.pos = (v2){5.0, 5.0};
    state.game_state.player.dir = (v2){0.0, 0.0};
    state.game_state.player.vel = (v2){0.0, 0.0};
    state.game_state.player.omega = 0.0f;
    state.game_state.player.z = 0.4;
    state.game_state.player.qe_geometry = DelaunayMeshGetEnclosingTriangle2(state.game_state.geometry_mesh, &(state.game_state.player.pos));

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
        Tick(&state.game_state, dt, &state.keyboard_state);
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
            f32 pix_per_tile = debug_window_size_xy / max(MAPDATA.n_tiles_x, MAPDATA.n_tiles_y);

            f32 offset_x = (debug_window_size_xy - pix_per_tile * MAPDATA.n_tiles_x) / 2.0;
            f32 offset_y = (debug_window_size_xy - pix_per_tile * MAPDATA.n_tiles_y) / 2.0;

            // Render the walls
            int grid_index = 0;
            for (int y = 0; y < MAPDATA.n_tiles_y; y++)
            {
                for (int x = 0; x < MAPDATA.n_tiles_x; x++)
                {
                    if (MAPDATA.tiles[grid_index] > 0)
                    {
                        // This tile is a wall
                        SDL_Rect rect;
                        rect.x = (int)(pix_per_tile * x + offset_x);
                        rect.y = (int)(pix_per_tile * y + offset_y);
                        rect.h = (int)(pix_per_tile);
                        rect.w = (int)(pix_per_tile);
                        SDL_RenderFillRect(debug_renderer, &rect);
                    }
                    grid_index++;
                }
            }

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
                for (int qe_index = 0; qe_index < DelaunayMeshNumQuarterEdges(state.game_state.geometry_mesh); qe_index++)
                {
                    QuarterEdge *qe = DelaunayMeshGetQuarterEdge(state.game_state.geometry_mesh, qe_index);

                    if (IsPrimalEdge(qe) && !DelaunayMeshIsBoundaryVertex(state.game_state.geometry_mesh, qe->vertex))
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
                        if (a > b && !DelaunayMeshIsBoundaryVertex(state.game_state.geometry_mesh, b))
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
                const v2 *a = DelaunayMeshGetTriangleVertex1(state.game_state.geometry_mesh, state.game_state.player.qe_geometry);
                const v2 *b = DelaunayMeshGetTriangleVertex2(state.game_state.geometry_mesh, state.game_state.player.qe_geometry);
                const v2 *c = DelaunayMeshGetTriangleVertex3(state.game_state.geometry_mesh, state.game_state.player.qe_geometry);

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
    // DeconstructDelaunayMesh(state.game_state.geometry_mesh);
    free(game_map.geometry_mesh->quarter_edges);
    free(game_map.geometry_mesh);

    // Free our assets
    free(ASSETS_BINARY_BLOB);
    free(WAD);
    free(ASSETS_BINARY_BLOB2);

    return 0;
}