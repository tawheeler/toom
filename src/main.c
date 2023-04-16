#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <sys/time.h>
#include <SDL2/SDL.h>

#define ASSERT(_e, ...) if (!(_e)) { fprintf(stderr, __VA_ARGS__); exit(1); }

typedef float f32;
typedef double f64;
typedef uint8_t bool;
typedef uint8_t u8;
typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;
typedef int8_t i8;
typedef int16_t i16;
typedef int32_t i32;
typedef int64_t i64;
typedef size_t usize;
typedef ssize_t isize;

#define SCREEN_SIZE_X 640
#define SCREEN_SIZE_Y 360

#define TILE_WIDTH 1.0f
#define WALL_HEIGHT 1.0f

#define TEXTURE_SIZE 64

typedef struct v2_s {f32 x, y;} v2;

#define dot(v0, v1) \
    ({ const v2 _v0 = (v0), _v1 = (v1); (_v0.x * _v1.x) + (_v0.y * _v1.y); })
#define length(v) ({ const v2 _v = (v); sqrtf(dot(_v, _v)); })
#define normalize(u) ({           \
        const v2 _u = (u);        \
        const f32 l = length(_u); \
        (v2) { _u.x/l, _u.y/l };  \
    })
#define rotr(v) ({ const v2 _v = (v); (v2) { -_v.y, _v.x }; })
#define min(a, b) ({ __typeof__(a) _a = (a), _b = (b); _a < _b ? _a : _b; })
#define max(a, b) ({ __typeof__(a) _a = (a), _b = (b); _a > _b ? _a : _b; })


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
    u8* tiles;
};

static inline u8 GetMapDataAt(struct Mapdata* mapdata, int x, int y) {
    return mapdata->tiles[(mapdata->n_tiles_y - y - 1)*(mapdata->n_tiles_x) + x];
}

// The global mapdata. This just points into our asset blob.
struct Mapdata MAPDATA;

struct Bitmap {
    u32  n_pixels;
    u32  n_pixels_per_column;
    u32  n_pixels_per_row;
    bool column_major;
    u32* abgr;
};

// The bitmap global variables. These just point into the binary blob.
struct Bitmap BITMAP;

static inline int GetColumnMajorPixelIndex(struct Bitmap* bitmap, int x, int y) {
    return y + x*bitmap->n_pixels_per_column;
}
static inline int GetRowMajorPixelIndex(struct Bitmap* bitmap, int x, int y) {
    return x + y*bitmap->n_pixels_per_row;
}
static inline u32 GetColumnMajorPixelAt(struct Bitmap* bitmap, int x, int y) {
    return bitmap->abgr[y + x*bitmap->n_pixels_per_column];
}
static inline u32 GetRowMajorPixelAt(struct Bitmap* bitmap, int x, int y) {
    return bitmap->abgr[x + y*bitmap->n_pixels_per_row];
}

enum KeyboardKeyState {
    KeyboardKeyState_Depressed = 0, // No recent event, key is still up
    KeyboardKeyState_Released = 1,  // Last event was a released event
    KeyboardKeyState_Held = 2,      // No recent event, key is still down
    KeyboardKeyState_Pressed = 3,    // Last event was a pressed event
    KeyboardKeyState_COUNT = 4
};

struct KeyBoardState {
    enum KeyboardKeyState up;
    enum KeyboardKeyState down;
    enum KeyboardKeyState right;
    enum KeyboardKeyState left;
    enum KeyboardKeyState a;
    enum KeyboardKeyState s;
    enum KeyboardKeyState d;
    enum KeyboardKeyState w;
    enum KeyboardKeyState q;
    enum KeyboardKeyState e;
    enum KeyboardKeyState r;

    enum KeyboardKeyState one;
    enum KeyboardKeyState two;
    enum KeyboardKeyState three;
    enum KeyboardKeyState four;
    enum KeyboardKeyState five;
    enum KeyboardKeyState six;
    enum KeyboardKeyState seven;
    enum KeyboardKeyState eight;
};

void ClearKeyboardState(struct KeyBoardState* kbs) {
    kbs->up = KeyboardKeyState_Depressed;
    kbs->down = KeyboardKeyState_Depressed;
    kbs->right = KeyboardKeyState_Depressed;
    kbs->left = KeyboardKeyState_Depressed;
    kbs->a = KeyboardKeyState_Depressed;
    kbs->s = KeyboardKeyState_Depressed;
    kbs->d = KeyboardKeyState_Depressed;
    kbs->w = KeyboardKeyState_Depressed;
    kbs->q = KeyboardKeyState_Depressed;
    kbs->e = KeyboardKeyState_Depressed;
    kbs->r = KeyboardKeyState_Depressed;

    kbs->one = KeyboardKeyState_Depressed;
    kbs->two = KeyboardKeyState_Depressed;
    kbs->three = KeyboardKeyState_Depressed;
    kbs->four = KeyboardKeyState_Depressed;
    kbs->five = KeyboardKeyState_Depressed;
    kbs->six = KeyboardKeyState_Depressed;
    kbs->seven = KeyboardKeyState_Depressed;
    kbs->eight = KeyboardKeyState_Depressed;
}

void DecayKeyboardState(struct KeyBoardState* kbs) {
    static enum KeyboardKeyState to_depressed_state[KeyboardKeyState_COUNT] = {
        KeyboardKeyState_Depressed,
        KeyboardKeyState_Depressed,
        KeyboardKeyState_Held,
        KeyboardKeyState_Held
    };

    kbs->up = to_depressed_state[kbs->up];
    kbs->down =  to_depressed_state[kbs->down];
    kbs->right = to_depressed_state[kbs->right];
    kbs->left = to_depressed_state[kbs->left];
    kbs->a = to_depressed_state[kbs->a];
    kbs->s = to_depressed_state[kbs->s];
    kbs->d = to_depressed_state[kbs->d];
    kbs->w = to_depressed_state[kbs->w];
    kbs->q = to_depressed_state[kbs->q];
    kbs->e = to_depressed_state[kbs->e];
    kbs->r = to_depressed_state[kbs->r];

    kbs->one = to_depressed_state[kbs->one];
    kbs->two = to_depressed_state[kbs->two];
    kbs->three = to_depressed_state[kbs->three];
    kbs->four = to_depressed_state[kbs->four];
    kbs->five = to_depressed_state[kbs->five];
    kbs->six = to_depressed_state[kbs->six];
    kbs->seven = to_depressed_state[kbs->seven];
    kbs->eight = to_depressed_state[kbs->eight];
}

bool IsPressed(enum KeyboardKeyState state) {
    static bool lookup[KeyboardKeyState_COUNT] = {0, 0, 1, 1};
    return lookup[state];
}

struct { 
    SDL_Window *window;
    SDL_Texture *texture;
    SDL_Renderer *renderer;
    u32 pixels[SCREEN_SIZE_X * SCREEN_SIZE_Y]; // row-major
    f32 wall_raycast_radius[SCREEN_SIZE_X];
    bool quit;

    v2 camera_pos;
    v2 camera_dir;
    v2 camera_dir_rotr;
    f32 camera_width;
    f32 camera_height;
    f32 camera_z;

    v2 player_speed;
    f32 player_omega;

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
                BITMAP.n_pixels            = *(u32*)(ASSETS_BINARY_BLOB + asset_byte_offset);
                asset_byte_offset += sizeof(u32);
                BITMAP.n_pixels_per_column = *(u32*)(ASSETS_BINARY_BLOB + asset_byte_offset);
                asset_byte_offset += sizeof(u32);
                BITMAP.column_major = ASSETS_BINARY_BLOB[asset_byte_offset];
                ASSERT(BITMAP.column_major, "Expected the wall texture to be column-major\n");
                asset_byte_offset += sizeof(u8);
                BITMAP.abgr = (u32*)(ASSETS_BINARY_BLOB + asset_byte_offset);
                BITMAP.n_pixels_per_row = BITMAP.n_pixels / BITMAP.n_pixels_per_column;
                loaded_textures = 1;
            } else if (strcmp(entry->name, "mapdata") == 0) {
                u32 asset_byte_offset = entry->byte_offset;
                MAPDATA.n_tiles = *(u32*)(ASSETS_BINARY_BLOB + asset_byte_offset);
                asset_byte_offset += sizeof(u32);
                MAPDATA.n_tiles_x = *(u32*)(ASSETS_BINARY_BLOB + asset_byte_offset);
                asset_byte_offset += sizeof(u32);
                MAPDATA.tiles = (u8*)(ASSETS_BINARY_BLOB + asset_byte_offset);
                MAPDATA.n_tiles_y = MAPDATA.n_tiles / MAPDATA.n_tiles_x;
                loaded_mapdata = 1;
            }
        }

        ASSERT(loaded_textures > 0, "Textures not loaded from assets\n");
        ASSERT(loaded_mapdata > 0, "Map data not loaded from assets\n");
    }
}

static void Tick(f32 dt) {

    v2 input_dir = {0.0, 0.0}; // In the body frame, which is right-handed, so y points left.
    if (IsPressed(state.keyboard_state.w)) {
        input_dir.x += 1.0;
    }
    if (IsPressed(state.keyboard_state.s)) {
        input_dir.x -= 1.0;
    }
    if (IsPressed(state.keyboard_state.d)) {
        input_dir.y -= 1.0;
    }
    if (IsPressed(state.keyboard_state.a)) {
        input_dir.y += 1.0;
    }

    int input_rot_dir = 0; // Right-hand rotation in plane (CCW)
    if (IsPressed(state.keyboard_state.q)) {
        input_rot_dir += 1;
    }
    if (IsPressed(state.keyboard_state.e)) {
        input_rot_dir -= 1;
    }

    if (IsPressed(state.keyboard_state.r)) {
        printf("Reloading assets.\n");
        LoadAssets();
        printf("DONE.\n");
    }

    if (IsPressed(state.keyboard_state.three)) {
        state.camera_z *= 0.95;
        printf("camera z: %.3f\n", state.camera_z);
    }
    if (IsPressed(state.keyboard_state.four)) {
        state.camera_z /= 0.95;
        printf("camera z: %.3f\n", state.camera_z);
    }
    if (IsPressed(state.keyboard_state.five)) {
        state.camera_height *= 0.95;
        printf("camera height: %.3f\n", state.camera_height);
    }
    if (IsPressed(state.keyboard_state.six)) {
        state.camera_height /= 0.95;
        printf("camera height: %.3f\n", state.camera_height);
    }
    if (IsPressed(state.keyboard_state.seven)) {
        state.camera_width *= 0.95;
        printf("camera width: %.3f\n", state.camera_width);
    }
    if (IsPressed(state.keyboard_state.eight)) {
        state.camera_width /= 0.95;
        printf("camera width: %.3f\n", state.camera_width);
    }

    // Update the player's velocity
    const f32 kPlayerInputAccel = 5.5;
    const f32 kPlayerInputAngularAccel = 8.5;
    const f32 kPlayerMaxSpeed = 4.0;
    const f32 kPlayerMaxOmega = 5.0;
    const f32 kAirFriction = 4.0;
    const f32 kAirFrictionRot = 4.0;

    // Note: Speed is in the global frame
    state.player_speed.x += (state.camera_dir.x*input_dir.x + state.camera_dir_rotr.x*input_dir.y) * kPlayerInputAccel * dt;
    state.player_speed.y += (state.camera_dir.y*input_dir.x + state.camera_dir_rotr.y*input_dir.y) * kPlayerInputAccel * dt;
    state.player_omega += input_rot_dir * kPlayerInputAngularAccel * dt;

    // Clamp the velocity to a maximum magnitude
    f32 speed = length(state.player_speed);
    if (speed > kPlayerMaxSpeed) {
        state.player_speed.x *= kPlayerMaxSpeed / speed;
        state.player_speed.y *= kPlayerMaxSpeed / speed;
    }
    if (state.player_omega > kPlayerMaxOmega) {
        state.player_omega *= kPlayerMaxOmega / state.player_omega;
    } else if (state.player_omega < -kPlayerMaxOmega) {
        state.player_omega *= - kPlayerMaxOmega / state.player_omega;
    }

    // Update the player's position
    state.camera_pos.x += state.player_speed.x * dt;
    state.camera_pos.y += state.player_speed.y * dt;

    // Update the player's rotational heading
    f32 theta = atan2(state.camera_dir.y, state.camera_dir.x);
    theta += state.player_omega * dt;
    state.camera_dir = ((v2) {cos(theta), sin(theta)});   
    state.camera_dir_rotr = rotr((state.camera_dir));

    // Apply air friction
    f32 air_friction_decay = exp(-kAirFriction * dt);
    state.player_speed.x *= air_friction_decay;
    state.player_speed.y *= air_friction_decay;
    state.player_omega *= exp(-kAirFrictionRot * dt);
}

static void Render() {
    // Render the floor and ceiling textures
    {

        // Ray direction for x = 0
        f32 half_camera_width = state.camera_width/2.0f;
        f32 ray_dir_lo_x = state.camera_dir.x + half_camera_width*state.camera_dir_rotr.x;
        f32 ray_dir_lo_y = state.camera_dir.y + half_camera_width*state.camera_dir_rotr.y;

        // Ray direction for x = SCREEN_SIZE_X
        f32 ray_dir_hi_x = state.camera_dir.x - half_camera_width*state.camera_dir_rotr.x;
        f32 ray_dir_hi_y = state.camera_dir.y - half_camera_width*state.camera_dir_rotr.y;        

        // Draw floor
        u32 texture_x_offset = 0;
        u32 texture_y_offset = 0;
        for (int y = 0; y < SCREEN_SIZE_Y/2; y++) {
            // Radius
            f32 zpp = (SCREEN_SIZE_Y/2.0f - y) * (state.camera_height / SCREEN_SIZE_Y);
            f32 radius = state.camera_z/zpp;

            // Location of the 1st ray's intersection
            f32 hit_x = state.camera_pos.x + radius * ray_dir_lo_x;
            f32 hit_y = state.camera_pos.y + radius * ray_dir_lo_y;

            // Each step is (hit_x2 - hit_x) / SCREEN_SIZE_X;
            // = ((state.camera_pos.x + radius * ray_dir_lo_x) - (state.camera_pos.x + radius * ray_dir_lo_x)) / SCREEN_SIZE_X
            // = (radius * ray_dir_lo_x - (radius * ray_dir_lo_x)) / SCREEN_SIZE_X
            // = radius * (ray_dir_hi_x - ray_dir_lo_x) / SCREEN_SIZE_X
            f32 step_x = radius * (ray_dir_hi_x - ray_dir_lo_x) / SCREEN_SIZE_X;
            f32 step_y = radius * (ray_dir_hi_y - ray_dir_lo_y) / SCREEN_SIZE_X;

            for (int x = 0; x < SCREEN_SIZE_X; x++) {
                u32 texture_x = (int)(fmod(hit_x, TILE_WIDTH)/TILE_WIDTH * TEXTURE_SIZE) & (TEXTURE_SIZE - 1);
                u32 texture_y = (int)(fmod(hit_y, TILE_WIDTH)/TILE_WIDTH * TEXTURE_SIZE) & (TEXTURE_SIZE - 1);
                u32 color = GetColumnMajorPixelAt(&BITMAP, texture_x+texture_x_offset, texture_y+texture_y_offset);
                state.pixels[(y * SCREEN_SIZE_X) + x] = color;

                // step
                hit_x += step_x;
                hit_y += step_y;
            }
        }

        // Draw ceiling
        texture_x_offset = 0;
        texture_y_offset = 0;
        for (int y = SCREEN_SIZE_Y/2 + 1; y < SCREEN_SIZE_Y; y++) {
            // Radius
            f32 zpp = (y - (SCREEN_SIZE_Y/2.0f)) * (state.camera_height / SCREEN_SIZE_Y);
            f32 radius = (WALL_HEIGHT - state.camera_z)/zpp;

            // Location of the 1st ray's intersection
            f32 hit_x = state.camera_pos.x + radius * ray_dir_lo_x;
            f32 hit_y = state.camera_pos.y + radius * ray_dir_lo_y;

            // Each step toward hit2
            f32 step_x = radius * (ray_dir_hi_x - ray_dir_lo_x) / SCREEN_SIZE_X;
            f32 step_y = radius * (ray_dir_hi_y - ray_dir_lo_y) / SCREEN_SIZE_X;

            for (int x = 0; x < SCREEN_SIZE_X; x++) {
                u32 texture_x = (int)(fmod(hit_x, TILE_WIDTH)/TILE_WIDTH * TEXTURE_SIZE) & (TEXTURE_SIZE - 1);
                u32 texture_y = (int)(fmod(hit_y, TILE_WIDTH)/TILE_WIDTH * TEXTURE_SIZE) & (TEXTURE_SIZE - 1);
                u32 color = GetColumnMajorPixelAt(&BITMAP, texture_x+texture_x_offset, texture_y+texture_y_offset);
                state.pixels[(y * SCREEN_SIZE_X) + x] = color;

                // step
                hit_x += step_x;
                hit_y += step_y;
            }
        }

        // ORIGINAL CODE
        // f32 x_side_frac = (x - SCREEN_SIZE_X/2.0f) / SCREEN_SIZE_X * state.camera_width; 
        // f32 xpp = sqrt(1.0f + x_side_frac*x_side_frac);
        // for (int y = y_hi_capped + 1; y < SCREEN_SIZE_Y; y++) {
        //     if (y > SCREEN_SIZE_Y/2.0f) {
        //         // NOTE: zpp is only positive if y > SCREEN_SIZE_Y/2
        //         f32 zpp = (y - (SCREEN_SIZE_Y/2.0f)) * (state.camera_height / SCREEN_SIZE_Y);
        //         f32 r = (WALL_HEIGHT - state.camera_z)*xpp/zpp;
        //         u32 texture_x = (int)(fmod(state.camera_pos.x + r*dir.x, TILE_WIDTH)/TILE_WIDTH * TEXTURE_SIZE);
        //         u32 texture_y = (int)(fmod(state.camera_pos.y + r*dir.y, TILE_WIDTH)/TILE_WIDTH * TEXTURE_SIZE);
        //         u32 color = BITMAP.abgr[texture_y+texture_y_offset + (texture_x+texture_x_offset)*BITMAP.n_pixels_per_column];
        //         state.pixels[(y * SCREEN_SIZE_X) + x] = color;
        //     }
        // }
    }

    // Get camera location's cell coordinates
    int x_ind_cam = (int)(floorf(state.camera_pos.x / TILE_WIDTH));
    int y_ind_cam = (int)(floorf(state.camera_pos.y / TILE_WIDTH));
    f32 x_rem_cam = state.camera_pos.x - TILE_WIDTH*x_ind_cam;
    f32 y_rem_cam = state.camera_pos.y - TILE_WIDTH*y_ind_cam;

    for (int x = 0; x < SCREEN_SIZE_X; x++) {
        
        // Camera to pixel column
        const f32 dw = state.camera_width/2 - (state.camera_width*x)/SCREEN_SIZE_X;
        const v2 cp = {
            state.camera_dir.x + dw*state.camera_dir_rotr.x,
            state.camera_dir.y + dw*state.camera_dir_rotr.y
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
            if (GetMapDataAt(&MAPDATA, x_ind,y_ind) > 0) {
                break;
            }
        }

        // Calculate the collision location
        const v2 collision = {
            TILE_WIDTH*x_ind + x_rem,
            TILE_WIDTH*y_ind + y_rem
        };

        // Calculate the ray length
        const f32 ray_len = length( ((v2) {collision.x - state.camera_pos.x, collision.y - state.camera_pos.y}) );
        state.wall_raycast_radius[x] = ray_len;

        // Calculate the pixel bounds that we fill the wall in for
        int y_lo = (int)(SCREEN_SIZE_Y/2.0f - cam_len*state.camera_z/ray_len * SCREEN_SIZE_Y / state.camera_height);
        int y_hi = (int)(SCREEN_SIZE_Y/2.0f + cam_len*(WALL_HEIGHT - state.camera_z)/ray_len * SCREEN_SIZE_Y / state.camera_height);
        int y_lo_capped = max(y_lo, 0);
        int y_hi_capped = min(y_hi, SCREEN_SIZE_Y-1);

        {
            // Texture x offset determines whether we draw the light or dark version
            u32 texture_x_offset = dx_ind == 0 ? 0 : TEXTURE_SIZE;
            u32 texture_y_offset = (GetMapDataAt(&MAPDATA, x_ind,y_ind) - 1) * TEXTURE_SIZE;

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
                state.pixels[(y * SCREEN_SIZE_X) + x] = color;
                y_loc += y_step;
            }
        }
    }
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

    // Create a renderer
    state.renderer = SDL_CreateRenderer(state.window, -1, SDL_RENDERER_PRESENTVSYNC);
    ASSERT(state.renderer, "Error creating SDL renderer: %s\n", SDL_GetError());

    // Create a texture
    state.texture = SDL_CreateTexture(
        state.renderer,
        SDL_PIXELFORMAT_ABGR8888,
        SDL_TEXTUREACCESS_STREAMING,
        SCREEN_SIZE_X,
        SCREEN_SIZE_Y);
    ASSERT(state.texture, "Error creating SDL texture: %s\n", SDL_GetError());

    // Init camera
    state.camera_pos = (v2) { 5.0f, 5.0f };
    state.camera_dir = ((v2) {cos(0.0), sin(0.0)});
    state.camera_dir_rotr = rotr((state.camera_dir));
    state.camera_width = 1.5f;
    state.camera_height = state.camera_width * SCREEN_SIZE_Y / SCREEN_SIZE_X;
    state.camera_z = 0.4;

    // Init player state
    state.player_speed = (v2) { 0.0f, 0.0f };
    state.player_omega = 0.0f;

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
        Tick(dt);
        timeval_tick_prev = timeval_tick;

        Render();

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
    }

    SDL_DestroyWindow(state.window);
    
    // Free our assets
    free(ASSETS_BINARY_BLOB);

    return 0;
}