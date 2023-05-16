#include <stdio.h>
#include <sys/time.h>
#include <SDL2/SDL.h>

#include "typedefs.h"
#include "constants.h"
#include "vec.h"
#include "input.h"
#include "bitmap.h"
#include "game.h"

#define ASSERT(_e, ...) if (!(_e)) { fprintf(stderr, __VA_ARGS__); exit(1); }

// Big binary assets blob that we load at init.
u8* ASSETS_BINARY_BLOB = NULL;
u32 ASSETS_BINARY_BLOB_SIZE = 0; // number of bytes

struct BinaryAssetTableOfContentEntry {
    u32 byte_offset;
    char name[16];  
};

struct Mapdata {
    u32  n_tiles;
    u32  n_tiles_x;
    u32  n_tiles_y;
    u8* tiles; // Each tile is either solid (0) or has a texture index (>0)
    u8* floor;
    u8* ceiling;
};

static inline u8 GetMapDataIndex(struct Mapdata* mapdata, int x, int y) {
    return (mapdata->n_tiles_y - y - 1)*(mapdata->n_tiles_x) + x;
}

// The global mapdata. This just points into our asset blob.
struct Mapdata MAPDATA;

// The bitmap global variables. These just point into the binary blob.
struct Bitmap BITMAP;
struct Bitmap BITMAP_STICK;

struct { 
    SDL_Window *window;
    SDL_Texture *texture;
    SDL_Renderer *renderer;

    u32 pixels[SCREEN_SIZE_X * SCREEN_SIZE_Y]; // row-major

    f32 wall_raycast_radius[SCREEN_SIZE_X];
    
    bool quit;

    struct GameState game_state;
    struct KeyBoardState keyboard_state;
} state;


f32 GetElapsedTimeSec(struct timeval* timeval_start, struct timeval* timeval_end) {
    f32 ms_elapsed = (timeval_end->tv_sec - timeval_start->tv_sec);  // sec
    ms_elapsed += (timeval_end->tv_usec - timeval_start->tv_usec) / 1000000.0f;   // us to s
    return ms_elapsed;
}

f64 GetElapsedTimeMillis(struct timeval* timeval_start, struct timeval* timeval_end) {
    f64 ms_elapsed = (timeval_end->tv_sec - timeval_start->tv_sec) * 1000.0;  // sec to ms
    ms_elapsed += (timeval_end->tv_usec - timeval_start->tv_usec) / 1000.0;   // us to ms
    return ms_elapsed;
}

static void LoadAssets() {
    // If we currently have any loaded assets, free them.
    // Note that we assume we are running single-threaded. If this changes in the future,
    // we will want to load the assets in a separate thread and swap them over at the 
    // appropriate time.
    if (ASSETS_BINARY_BLOB) {
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

        FILE* fileptr = fopen("assets/assets.bin","rb");
        ASSERT(fileptr, "Error opening assets file\n");
        
        // Count the number of bytes
        fseek(fileptr, 0, SEEK_END);
        u32 n_bytes_in_file = ftell(fileptr);
        ASSETS_BINARY_BLOB_SIZE = n_bytes_in_file - 4; // Everything past the 'TOOM' header
        
        // Read in the binary assets as a single blob
        fseek(fileptr, 4, SEEK_SET); // Skip the header bytes
        ASSETS_BINARY_BLOB = (u8*) malloc(ASSETS_BINARY_BLOB_SIZE);
        ASSERT(ASSETS_BINARY_BLOB, "Failed to allocate assets blob\n");
        ASSERT(fread(ASSETS_BINARY_BLOB, sizeof(u8), ASSETS_BINARY_BLOB_SIZE, fileptr) == ASSETS_BINARY_BLOB_SIZE, "Failed to read assets blob when loading assets\n");

        fclose(fileptr);
    }
    {
        // Process the loaded assets from the loaded binary blob

        // Read the number of table of content entries
        u32 byte_index = ASSETS_BINARY_BLOB_SIZE - sizeof(u32);
        u32 n_toc_entries = *(u32*)(ASSETS_BINARY_BLOB + byte_index);
        ASSERT(ASSETS_BINARY_BLOB_SIZE > sizeof(struct BinaryAssetTableOfContentEntry) * n_toc_entries + 4, "Number of table of content entries is impossible given the number of bytes\n");

        bool loaded_textures = 0;
        bool loaded_mapdata = 0;

        // Scan through them in reverse order
        for (int i = n_toc_entries; i > 0; i--) {
            byte_index -= sizeof(struct BinaryAssetTableOfContentEntry);
            struct BinaryAssetTableOfContentEntry* entry = (struct BinaryAssetTableOfContentEntry*)(ASSETS_BINARY_BLOB + byte_index);
            // Make the name null-terminated just in case.
            entry->name[15] = '\0';
            printf("Entry %d: %s at offset %d\n", i, entry->name, entry->byte_offset);

            // In the future, load them into a map or something. For now, we're specifically looking for either the wall or floor textures.
            if (strcmp(entry->name, "textures") == 0) {
                u32 asset_byte_offset = entry->byte_offset;
                BITMAP.n_pixels = *(u32*)(ASSETS_BINARY_BLOB + asset_byte_offset);
                asset_byte_offset += sizeof(u32);
                BITMAP.n_pixels_per_column = *(u32*)(ASSETS_BINARY_BLOB + asset_byte_offset);
                asset_byte_offset += sizeof(u32);
                BITMAP.column_major = ASSETS_BINARY_BLOB[asset_byte_offset];
                ASSERT(BITMAP.column_major, "Expected the wall texture to be column-major\n");
                asset_byte_offset += sizeof(u8);
                BITMAP.abgr = (u32*)(ASSETS_BINARY_BLOB + asset_byte_offset);
                BITMAP.n_pixels_per_row = BITMAP.n_pixels / BITMAP.n_pixels_per_column;
                loaded_textures = 1;
            } else if (strcmp(entry->name, "stick") == 0) {
                u32 asset_byte_offset = entry->byte_offset;
                BITMAP_STICK.n_pixels = *(u32*)(ASSETS_BINARY_BLOB + asset_byte_offset);
                asset_byte_offset += sizeof(u32);
                BITMAP_STICK.n_pixels_per_column = *(u32*)(ASSETS_BINARY_BLOB + asset_byte_offset);
                asset_byte_offset += sizeof(u32);
                BITMAP_STICK.column_major = ASSETS_BINARY_BLOB[asset_byte_offset];
                ASSERT(BITMAP_STICK.column_major, "Expected the stick texture to be column-major\n");
                asset_byte_offset += sizeof(u8);
                BITMAP_STICK.abgr = (u32*)(ASSETS_BINARY_BLOB + asset_byte_offset);
                BITMAP_STICK.n_pixels_per_row = BITMAP_STICK.n_pixels / BITMAP_STICK.n_pixels_per_column;
                loaded_textures = 1;
            } else if (strcmp(entry->name, "mapdata") == 0) {
                u32 asset_byte_offset = entry->byte_offset;
                MAPDATA.n_tiles = *(u32*)(ASSETS_BINARY_BLOB + asset_byte_offset);
                asset_byte_offset += sizeof(u32);
                MAPDATA.n_tiles_x = *(u32*)(ASSETS_BINARY_BLOB + asset_byte_offset);
                asset_byte_offset += sizeof(u32);
                MAPDATA.tiles = (u8*)(ASSETS_BINARY_BLOB + asset_byte_offset);
                asset_byte_offset += MAPDATA.n_tiles * sizeof(u8);
                MAPDATA.floor = (u8*)(ASSETS_BINARY_BLOB + asset_byte_offset);
                asset_byte_offset += MAPDATA.n_tiles * sizeof(u8);
                MAPDATA.ceiling = (u8*)(ASSETS_BINARY_BLOB + asset_byte_offset);
                MAPDATA.n_tiles_y = MAPDATA.n_tiles / MAPDATA.n_tiles_x;
                loaded_mapdata = 1;
            }
        }

        ASSERT(loaded_textures > 0, "Textures not loaded from assets\n");
        ASSERT(loaded_mapdata > 0, "Map data not loaded from assets\n");
    }
}

void RenderFloorAndCeiling(
    u32* pixels,
    struct CameraState* camera
) {
    // Ray direction for x = 0
    f32 half_camera_width = camera->fov.x/2.0f;
    f32 ray_dir_lo_x = camera->dir.x - half_camera_width*camera->dir.y;
    f32 ray_dir_lo_y = camera->dir.y + half_camera_width*camera->dir.x;

    // Ray direction for x = SCREEN_SIZE_X
    f32 ray_dir_hi_x = camera->dir.x + half_camera_width*camera->dir.y;
    f32 ray_dir_hi_y = camera->dir.y - half_camera_width*camera->dir.x;        

    // Draw floor
    for (int y = 0; y < SCREEN_SIZE_Y/2; y++) {
        // Radius
        f32 zpp = (SCREEN_SIZE_Y/2.0f - y) * (camera->fov.y / SCREEN_SIZE_Y);
        f32 radius = camera->z / zpp;

        // Location of the 1st ray's intersection
        f32 hit_x = camera->pos.x + radius * ray_dir_lo_x;
        f32 hit_y = camera->pos.y + radius * ray_dir_lo_y;

        // Each step is (hit_x2 - hit_x) / SCREEN_SIZE_X;
        // = ((camera->pos.x + radius * ray_dir_lo_x) - (camera->pos.x + radius * ray_dir_lo_x)) / SCREEN_SIZE_X
        // = (radius * ray_dir_lo_x - (radius * ray_dir_lo_x)) / SCREEN_SIZE_X
        // = radius * (ray_dir_hi_x - ray_dir_lo_x) / SCREEN_SIZE_X
        f32 step_x = radius * (ray_dir_hi_x - ray_dir_lo_x) / SCREEN_SIZE_X;
        f32 step_y = radius * (ray_dir_hi_y - ray_dir_lo_y) / SCREEN_SIZE_X;

        for (int x = 0; x < SCREEN_SIZE_X; x++) {

            int x_ind_hit = (int)(floorf(hit_x / TILE_WIDTH));
            int y_ind_hit = (int)(floorf(hit_y / TILE_WIDTH));
            f32 x_rem_hit = hit_x - TILE_WIDTH*x_ind_hit;
            f32 y_rem_hit = hit_y - TILE_WIDTH*y_ind_hit;
            x_ind_hit = clamp(x_ind_hit, 0, MAPDATA.n_tiles_x-1);
            y_ind_hit = clamp(y_ind_hit, 0, MAPDATA.n_tiles_y-1);

            u32 texture_x_offset = 0;
            u32 texture_y_offset = (MAPDATA.floor[GetMapDataIndex(&MAPDATA, x_ind_hit, y_ind_hit)] - 1) * TEXTURE_SIZE;

            u32 texture_x = (int)(x_rem_hit/TILE_WIDTH * TEXTURE_SIZE);
            u32 texture_y = (int)(y_rem_hit/TILE_WIDTH * TEXTURE_SIZE);
            
            u32 color = GetColumnMajorPixelAt(&BITMAP, texture_x+texture_x_offset, texture_y+texture_y_offset);
            pixels[(y * SCREEN_SIZE_X) + x] = color;

            // step
            hit_x += step_x;
            hit_y += step_y;
        }
    }

    // Draw ceiling
    for (int y = SCREEN_SIZE_Y/2 + 1; y < SCREEN_SIZE_Y; y++) {
        // Radius
        f32 zpp = (y - (SCREEN_SIZE_Y/2.0f)) * (camera->fov.y / SCREEN_SIZE_Y);
        f32 radius = (WALL_HEIGHT - camera->z)/zpp;

        // Location of the 1st ray's intersection
        f32 hit_x = camera->pos.x + radius * ray_dir_lo_x;
        f32 hit_y = camera->pos.y + radius * ray_dir_lo_y;

        // Each step toward hit2
        f32 step_x = radius * (ray_dir_hi_x - ray_dir_lo_x) / SCREEN_SIZE_X;
        f32 step_y = radius * (ray_dir_hi_y - ray_dir_lo_y) / SCREEN_SIZE_X;

        for (int x = 0; x < SCREEN_SIZE_X; x++) {
            int x_ind_hit = (int)(floorf(hit_x / TILE_WIDTH));
            int y_ind_hit = (int)(floorf(hit_y / TILE_WIDTH));
            f32 x_rem_hit = hit_x - TILE_WIDTH*x_ind_hit;
            f32 y_rem_hit = hit_y - TILE_WIDTH*y_ind_hit;
            x_ind_hit = clamp(x_ind_hit, 0, MAPDATA.n_tiles_x-1);
            y_ind_hit = clamp(y_ind_hit, 0, MAPDATA.n_tiles_y-1);

            u32 texture_x_offset = 0;
            u32 texture_y_offset = (MAPDATA.ceiling[GetMapDataIndex(&MAPDATA, x_ind_hit, y_ind_hit)] - 1) * TEXTURE_SIZE;

            u32 texture_x = (int)(x_rem_hit/TILE_WIDTH * TEXTURE_SIZE);
            u32 texture_y = (int)(y_rem_hit/TILE_WIDTH * TEXTURE_SIZE);
            
            u32 color = GetColumnMajorPixelAt(&BITMAP, texture_x+texture_x_offset, texture_y+texture_y_offset);
            pixels[(y * SCREEN_SIZE_X) + x] = color;

            // step
            hit_x += step_x;
            hit_y += step_y;
        }
    }
}

void RenderWalls(
    u32* pixels,
    f32* wall_raycast_radius,
    struct CameraState* camera
) {
    // Get camera location's cell coordinates
    int x_ind_cam = (int)(floorf(camera->pos.x / TILE_WIDTH));
    int y_ind_cam = (int)(floorf(camera->pos.y / TILE_WIDTH));
    f32 x_rem_cam = camera->pos.x - TILE_WIDTH*x_ind_cam;
    f32 y_rem_cam = camera->pos.y - TILE_WIDTH*y_ind_cam;

    for (int x = 0; x < SCREEN_SIZE_X; x++) {
        
        // Camera to pixel column
        const f32 dw = camera->fov.x/2 - (camera->fov.x*x)/SCREEN_SIZE_X;
        const v2 cp = {
            camera->dir.x - dw*camera->dir.y,
            camera->dir.y + dw*camera->dir.x
        };

        // Distance from the camera to the column
        const f32 cam_len = length( (cp) );
        
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
        if (dir.x < 0) {
            dx_ind_dir = -1;
            dx_a = -1.0f/dir.x;
            dx_b = 0.0;
        } else if (dir.x > 0) {
            dx_ind_dir = 1;
            dx_a = -1.0f/dir.x;
            dx_b = TILE_WIDTH/dir.x;
        }

        int dy_ind_dir = 0;
        f32 dy_a = 0.0;
        f32 dy_b = INFINITY;
        if (dir.y < 0) {
            dy_ind_dir = -1;
            dy_a = -1.0f/dir.y;
            dy_b = 0.0;
        } else if (dir.y > 0) {
            dy_ind_dir = 1;
            dy_a = -1.0f/dir.y;
            dy_b = TILE_WIDTH/dir.y;
        }

        // Step through cells until we hit an occupied cell
        int n_steps = 0;
        int dx_ind, dy_ind;
        while (n_steps < 100) {
            n_steps += 1;

            f32 dt_best = INFINITY;
            dx_ind = 0;
            dy_ind = 0;
            
            f32 dt_x = dx_a*x_rem + dx_b;
            f32 dt_y = dy_a*y_rem + dy_b;
            if (dt_x < dt_y) {
                dt_best = dt_x;
                dx_ind = dx_ind_dir;
                dy_ind = 0;
            } else {
                dt_best = dt_y;
                dx_ind = 0;
                dy_ind = dy_ind_dir;
            }

            // Move up to the next cell
            x_ind += dx_ind;
            y_ind += dy_ind;
            x_rem += dir.x * dt_best - TILE_WIDTH*dx_ind;
            y_rem += dir.y * dt_best - TILE_WIDTH*dy_ind;

            // Check to see if the new cell is solid
            if (MAPDATA.tiles[GetMapDataIndex(&MAPDATA, x_ind, y_ind)] > 0) {
                break;
            }
        }

        // Calculate the collision location
        const v2 collision = {
            TILE_WIDTH*x_ind + x_rem,
            TILE_WIDTH*y_ind + y_rem
        };

        // Calculate the ray length
        const f32 ray_len = length( ((v2) {collision.x - camera->pos.x, collision.y - camera->pos.y}) );
        wall_raycast_radius[x] = ray_len;

        // Calculate the pixel bounds that we fill the wall in for
        int y_lo = (int)(SCREEN_SIZE_Y/2.0f - cam_len*camera->z/ray_len * SCREEN_SIZE_Y / camera->fov.y);
        int y_hi = (int)(SCREEN_SIZE_Y/2.0f + cam_len*(WALL_HEIGHT - camera->z)/ray_len * SCREEN_SIZE_Y / camera->fov.y);
        int y_lo_capped = max(y_lo, 0);
        int y_hi_capped = min(y_hi, SCREEN_SIZE_Y-1);

        // Texture x offset determines whether we draw the light or dark version
        u32 texture_x_offset = dx_ind == 0 ? 0 : TEXTURE_SIZE;
        u32 texture_y_offset = (MAPDATA.tiles[GetMapDataIndex(&MAPDATA, x_ind,y_ind)] - 1) * TEXTURE_SIZE;

        f32 rem = 0.0f;
        if (dx_ind == 0) {
            rem = dy_ind < 0 ? TILE_WIDTH - x_rem : x_rem;
        } else {
            rem = dx_ind < 0 ? y_rem : TILE_WIDTH - y_rem;
        }
        u32 texture_x = min((int) (TEXTURE_SIZE * rem / TILE_WIDTH), TEXTURE_SIZE-1);
        u32 baseline = GetColumnMajorPixelIndex(&BITMAP, texture_x+texture_x_offset, texture_y_offset);
        u32 denom = max(1, y_hi - y_lo);
        f32 y_loc = (f32)((y_hi - y_hi_capped) * TEXTURE_SIZE) / denom;
        f32 y_step = (f32)(TEXTURE_SIZE) / denom;
        for (int y = y_hi_capped; y >= y_lo_capped; y--) {
            u32 texture_y = min((u32) (y_loc), TEXTURE_SIZE-1);
            u32 color = BITMAP.abgr[texture_y+baseline];
            pixels[(y * SCREEN_SIZE_X) + x] = color;
            y_loc += y_step;
        }
    }
}

void RenderObjects(
   u32* pixels,
   f32* wall_raycast_radius,
   struct CameraState* camera
) {
    v2 stick_pos = { 10.0f, 4.5f };
    v2 stick_rel_camera = {
        stick_pos.x - camera->pos.x,
        stick_pos.y - camera->pos.y};
    f32 dist_to_player = length(stick_rel_camera);

    f32 s = camera->dir.y;
    f32 c = camera->dir.x;
    v2 stick_pos_cam_body = {
        c*stick_rel_camera.x + s*stick_rel_camera.y,
        c*stick_rel_camera.y - s*stick_rel_camera.x
    };

    // Only render if it is on the postive side of the camera
    if (stick_pos_cam_body.x > 1e-3) {

        // Calculate the column pixel bounds
        const f32 cam_len = sqrt(1.0 + (stick_pos_cam_body.y / stick_pos_cam_body.x)*(stick_pos_cam_body.y / stick_pos_cam_body.x));
        int y_lo = (int)(SCREEN_SIZE_Y/2.0f - cam_len*camera->z/dist_to_player * SCREEN_SIZE_Y / camera->fov.y);
        int y_hi = (int)(SCREEN_SIZE_Y/2.0f + cam_len*(WALL_HEIGHT - camera->z)/dist_to_player * SCREEN_SIZE_Y / camera->fov.y);
        int y_lo_capped = max(y_lo, 0);
        int y_hi_capped = min(y_hi, SCREEN_SIZE_Y-1);
        u32 denom = max(1, y_hi - y_lo);
        f32 y_step = (f32)(TEXTURE_SIZE) / denom;

        int x_column_lo = (int)((0.5 - ((stick_pos_cam_body.y + TILE_WIDTH/2) / stick_pos_cam_body.x)/(camera->fov.x))*SCREEN_SIZE_X);
        int x_column_hi = (int)((0.5 - ((stick_pos_cam_body.y - TILE_WIDTH/2) / stick_pos_cam_body.x)/(camera->fov.x))*SCREEN_SIZE_X);
        f32 x_step = ((f32)(TEXTURE_SIZE)/(x_column_hi - x_column_lo + 1));
        f32 x_loc = 0.0f;
        for (int x = x_column_lo; x <= x_column_hi; x++) {

            if (x >= 0 && x < SCREEN_SIZE_X && wall_raycast_radius[x] > dist_to_player) {
                u32 texture_x = min((u32) (x_loc), TEXTURE_SIZE-1);
                f32 y_loc = (f32)((y_hi - y_hi_capped) * TEXTURE_SIZE) / denom;
                for (int y = y_hi_capped; y >= y_lo_capped; y--) {
                    u32 texture_y = min((u32) (y_loc), TEXTURE_SIZE-1);
                    u32 color = GetColumnMajorPixelAt(&BITMAP_STICK, texture_x, texture_y);
                    if ((color >> 24) > 0) {
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
    u32* pixels,
    f32* wall_raycast_radius,
    struct CameraState* camera
) {
    RenderFloorAndCeiling(pixels, camera);
    RenderWalls(pixels, wall_raycast_radius, camera);
    RenderObjects(pixels, wall_raycast_radius, camera);
}

int main(int argc, char *argv[]) {

    // Load our assets
    printf("Loading assets.\n");
    LoadAssets();
    printf("DONE.\n");

    // Initialize SDL
    ASSERT(
        SDL_Init(SDL_INIT_VIDEO) == 0,
        "SDL initialization failed: %s\n",
        SDL_GetError()
    );

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
    SDL_Window* debug_window = SDL_CreateWindow(
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
    SDL_Renderer* debug_renderer = SDL_CreateRenderer(debug_window, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
    ASSERT(state.renderer, "Error creating SDL Debug renderer: %s\n", SDL_GetError());
    SDL_SetRenderDrawColor( debug_renderer, 0xFF, 0xFF, 0xFF, 0xFF ); //Initialize renderer color

    // Create a texture
    state.texture = SDL_CreateTexture(
        state.renderer,
        SDL_PIXELFORMAT_ABGR8888,
        SDL_TEXTUREACCESS_STREAMING,
        SCREEN_SIZE_X,
        SCREEN_SIZE_Y);
    ASSERT(state.texture, "Error creating SDL texture: %s\n", SDL_GetError());

    // Init camera
    state.game_state.camera.pos = (v2) { 5.0f, 5.0f };
    state.game_state.camera.dir = ((v2) {cos(0.0), sin(0.0)});
    state.game_state.camera.fov.x = 1.5f;
    state.game_state.camera.fov.y = state.game_state.camera.fov.x * SCREEN_SIZE_Y / SCREEN_SIZE_X;
    state.game_state.camera.z = 0.4;

    // Init player state
    state.game_state.player_speed = (v2) { 0.0f, 0.0f };
    state.game_state.player_omega = 0.0f;

    // Init keyboard
    ClearKeyboardState(&state.keyboard_state);

    // Time structs
    struct timeval timeval_frame_start, timeval_frame_end, timeval_tick, timeval_tick_prev;
    gettimeofday(&timeval_frame_start, NULL);
    gettimeofday(&timeval_frame_end, NULL);
    gettimeofday(&timeval_tick, NULL);
    gettimeofday(&timeval_tick_prev, NULL);

    // Main loop
    state.quit = 0;
    while (state.quit == 0) {
        // Time of the start of the frame
        gettimeofday(&timeval_frame_start, NULL);

        SDL_Event event;
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT) {
                state.quit = 1;
                break;
            } else if (event.type == SDL_KEYDOWN) {
                switch (event.key.keysym.sym) {
                    case (SDLK_UP)    : state.keyboard_state.up    = KeyboardKeyState_Pressed; break;
                    case (SDLK_DOWN)  : state.keyboard_state.down  = KeyboardKeyState_Pressed; break;
                    case (SDLK_LEFT)  : state.keyboard_state.left  = KeyboardKeyState_Pressed; break;
                    case (SDLK_RIGHT) : state.keyboard_state.right = KeyboardKeyState_Pressed; break;
                    case (SDLK_a)     : state.keyboard_state.a     = KeyboardKeyState_Pressed; break;
                    case (SDLK_s)     : state.keyboard_state.s     = KeyboardKeyState_Pressed; break;
                    case (SDLK_d)     : state.keyboard_state.d     = KeyboardKeyState_Pressed; break;
                    case (SDLK_w)     : state.keyboard_state.w     = KeyboardKeyState_Pressed; break;
                    case (SDLK_q)     : state.keyboard_state.q     = KeyboardKeyState_Pressed; break;
                    case (SDLK_e)     : state.keyboard_state.e     = KeyboardKeyState_Pressed; break;
                    case (SDLK_r)     : state.keyboard_state.r     = KeyboardKeyState_Pressed; break;
                    case (SDLK_1)     : state.keyboard_state.one   = KeyboardKeyState_Pressed; break;
                    case (SDLK_2)     : state.keyboard_state.two   = KeyboardKeyState_Pressed; break;
                    case (SDLK_3)     : state.keyboard_state.three = KeyboardKeyState_Pressed; break;
                    case (SDLK_4)     : state.keyboard_state.four  = KeyboardKeyState_Pressed; break;
                    case (SDLK_5)     : state.keyboard_state.five  = KeyboardKeyState_Pressed; break;
                    case (SDLK_6)     : state.keyboard_state.six   = KeyboardKeyState_Pressed; break;
                    case (SDLK_7)     : state.keyboard_state.seven = KeyboardKeyState_Pressed; break;
                    case (SDLK_8)     : state.keyboard_state.eight = KeyboardKeyState_Pressed; break;
                }
            } else if (event.type == SDL_KEYUP) {
                switch (event.key.keysym.sym) {
                    case (SDLK_UP)    : state.keyboard_state.up    = KeyboardKeyState_Released; break;
                    case (SDLK_DOWN)  : state.keyboard_state.down  = KeyboardKeyState_Released; break;
                    case (SDLK_LEFT)  : state.keyboard_state.left  = KeyboardKeyState_Released; break;
                    case (SDLK_RIGHT) : state.keyboard_state.right = KeyboardKeyState_Released; break;
                    case (SDLK_a)     : state.keyboard_state.a     = KeyboardKeyState_Released; break;
                    case (SDLK_s)     : state.keyboard_state.s     = KeyboardKeyState_Released; break;
                    case (SDLK_d)     : state.keyboard_state.d     = KeyboardKeyState_Released; break;
                    case (SDLK_w)     : state.keyboard_state.w     = KeyboardKeyState_Released; break;
                    case (SDLK_q)     : state.keyboard_state.q     = KeyboardKeyState_Released; break;
                    case (SDLK_e)     : state.keyboard_state.e     = KeyboardKeyState_Released; break;
                    case (SDLK_r)     : state.keyboard_state.r     = KeyboardKeyState_Released; break;
                    case (SDLK_1)     : state.keyboard_state.one   = KeyboardKeyState_Released; break;
                    case (SDLK_2)     : state.keyboard_state.two   = KeyboardKeyState_Released; break;
                    case (SDLK_3)     : state.keyboard_state.three = KeyboardKeyState_Released; break;
                    case (SDLK_4)     : state.keyboard_state.four  = KeyboardKeyState_Released; break;
                    case (SDLK_5)     : state.keyboard_state.five  = KeyboardKeyState_Released; break;
                    case (SDLK_6)     : state.keyboard_state.six   = KeyboardKeyState_Released; break;
                    case (SDLK_7)     : state.keyboard_state.seven = KeyboardKeyState_Released; break;
                    case (SDLK_8)     : state.keyboard_state.eight = KeyboardKeyState_Released; break;
                }
            }
        }

        // Calc elapsed time since previous tick, then run tick
        gettimeofday(&timeval_tick, NULL);
        const f32 dt = GetElapsedTimeSec(&timeval_tick_prev, &timeval_tick);
        Tick(&state.game_state, dt, &state.keyboard_state);
        timeval_tick_prev = timeval_tick;

        // Check for a kayboard press to reload our assets
        if (IsNewlyPressed(state.keyboard_state.r)) {
            printf("Reloading assets.\n");
            LoadAssets();
            printf("DONE.\n");
        }

        Render(state.pixels, state.wall_raycast_radius, &state.game_state.camera);
        DecayKeyboardState(&state.keyboard_state);

        // Get timer end for all the non-SDL stuff
        gettimeofday(&timeval_frame_end, NULL);
        const f64 game_ms_elapsed = GetElapsedTimeMillis(&timeval_frame_start, &timeval_frame_end);

        static f64 tot_game_ms_elapsed = 0.0f;
        static i64 n_frames = 0; 
        tot_game_ms_elapsed += game_ms_elapsed;
        n_frames += 1;

        printf("Game: %.3f ms, %.1f fps, frame dt %.3fs, mean game ms: %.3f\n", game_ms_elapsed, 1000.0f / max(1.0f, game_ms_elapsed), dt, tot_game_ms_elapsed / n_frames);

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

        //Clear screen
        SDL_SetRenderDrawColor( debug_renderer, 0xFF, 0xFF, 0xFF, 0xFF );
        SDL_RenderClear( debug_renderer );

        {
            SDL_SetRenderDrawColor(debug_renderer, 0x1B, 0xA1, 0xEA, 0xFF);
            f32 pix_per_tile = debug_window_size_xy / max(MAPDATA.n_tiles_x, MAPDATA.n_tiles_y);

            f32 offset_x =  (debug_window_size_xy - pix_per_tile * MAPDATA.n_tiles_x) / 2.0;
            f32 offset_y =  (debug_window_size_xy - pix_per_tile * MAPDATA.n_tiles_y) / 2.0;

            // Render the walls
            int grid_index = 0;
            for (int y = 0; y < MAPDATA.n_tiles_y; y++) {
                for (int x = 0; x < MAPDATA.n_tiles_x; x++) {
                    if (MAPDATA.tiles[grid_index] > 0) {
                        // This tile is a wall
                        SDL_Rect rect;
                        rect.x = (int)(pix_per_tile*x + offset_x);
                        rect.y = (int)(pix_per_tile*y + offset_y);
                        rect.h = (int)(pix_per_tile);
                        rect.w = (int)(pix_per_tile);
                        SDL_RenderFillRect(debug_renderer, &rect);
                    }
                    grid_index ++;
                }
            }

            // Render the camera raycasts
            SDL_SetRenderDrawColor(debug_renderer, 0xF5, 0x61, 0x5C, 0xFF);
            struct CameraState* camera = &state.game_state.camera;
            int camera_x = camera->pos.x / TILE_WIDTH * pix_per_tile + offset_x;
            int camera_y = debug_window_size_xy - (camera->pos.y / TILE_WIDTH * pix_per_tile + offset_y);

            for (int x = 0; x < SCREEN_SIZE_X; x++) {
                f32 r = state.wall_raycast_radius[x];
                
                const f32 dw = camera->fov.x/2 - (camera->fov.x*x)/SCREEN_SIZE_X;
                const v2 cp = {
                    camera->dir.x - dw*camera->dir.y,
                    camera->dir.y + dw*camera->dir.x
                };
                const f32 cam_len = length( (cp) );
                const v2 dir = {cp.x / cam_len, cp.y / cam_len};

                f32 camera_x1 = camera->pos.x + r * dir.x;
                f32 camera_y1 = camera->pos.y + r * dir.y;
                int camera_x2 = camera_x1 / TILE_WIDTH * pix_per_tile + offset_x;
                int camera_y2 = debug_window_size_xy - (camera_y1 / TILE_WIDTH * pix_per_tile + offset_y);
                SDL_RenderDrawLine(debug_renderer, camera_x, camera_y, camera_x2, camera_y2);
            }

            

            // Render the camera
            // SDL_SetRenderDrawColor(debug_renderer, 0xF5, 0x61, 0x5C, 0xFF);
            // struct CameraState* camera = &state.game_state.camera;
            // int camera_x = camera->pos.x / TILE_WIDTH * pix_per_tile + offset_x;
            // int camera_y = debug_window_size_xy - (camera->pos.y / TILE_WIDTH * pix_per_tile + offset_y);

            // f32 camera_fov = 1.0;
            // f32 camera_x1 = camera->pos.x + camera->dir.x * camera_fov;
            // f32 camera_y1 = camera->pos.y + camera->dir.y * camera_fov;
            // int camera_x2 = camera_x1 / TILE_WIDTH * pix_per_tile + offset_x;
            // int camera_y2 = debug_window_size_xy - (camera_y1 / TILE_WIDTH * pix_per_tile + offset_y);
            // SDL_RenderDrawLine(debug_renderer, camera_x, camera_y, camera_x2, camera_y2);
        }
        

        SDL_RenderPresent(debug_renderer);
    }

    SDL_DestroyWindow(state.window);
    SDL_DestroyWindow(debug_window);
    
    // Free our assets
    free(ASSETS_BINARY_BLOB);

    return 0;
}