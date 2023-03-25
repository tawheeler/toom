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

typedef struct v2_s {f32 x, y;} v2;
typedef struct v2i_s { i32 x, y;} v2i;

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

static u8 MAPDATA[8*8] = {
    1, 1, 1, 1, 1, 1, 1, 1,
    1, 0, 0, 0, 0, 0, 0, 1,
    1, 0, 0, 3, 0, 0, 4, 1,
    1, 0, 0, 0, 0, 0, 0, 1,
    1, 0, 0, 0, 0, 0, 4, 1,
    1, 0, 2, 0, 0, 0, 0, 1,
    1, 0, 0, 0, 0, 0, 0, 1,
    1, 1, 1, 1, 1, 1, 1, 1,
};

struct Bitmap {
    u32 n_pixels;
    u32 n_pixels_per_column;
    u32* abgr; // pixel x, y is at abgr[y + x*n_pixels_per_column]
};

// The global variable of our main bitmap.
struct Bitmap bitmap;

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

    enum KeyboardKeyState one;
    enum KeyboardKeyState two;
    enum KeyboardKeyState three;
    enum KeyboardKeyState four;
    enum KeyboardKeyState five;
    enum KeyboardKeyState six;
    enum KeyboardKeyState seven;
    enum KeyboardKeyState eight;
};

void clear_keyboard_state(struct KeyBoardState* kbs) {
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

    kbs->one = KeyboardKeyState_Depressed;
    kbs->two = KeyboardKeyState_Depressed;
    kbs->three = KeyboardKeyState_Depressed;
    kbs->four = KeyboardKeyState_Depressed;
    kbs->five = KeyboardKeyState_Depressed;
    kbs->six = KeyboardKeyState_Depressed;
    kbs->seven = KeyboardKeyState_Depressed;
    kbs->eight = KeyboardKeyState_Depressed;
}

void decay_keyboard_state(struct KeyBoardState* kbs) {
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

    kbs->one = to_depressed_state[kbs->one];
    kbs->two = to_depressed_state[kbs->two];
    kbs->three = to_depressed_state[kbs->three];
    kbs->four = to_depressed_state[kbs->four];
    kbs->five = to_depressed_state[kbs->five];
    kbs->six = to_depressed_state[kbs->six];
    kbs->seven = to_depressed_state[kbs->seven];
    kbs->eight = to_depressed_state[kbs->eight];
}

bool is_pressed(enum KeyboardKeyState state) {
    static bool lookup[KeyboardKeyState_COUNT] = {0, 0, 1, 1};
    return lookup[state];
}

// TODO: Could we store the pixels in column-major? We're always rendering
//       in vertical lines, so I suspect that would be more efficient.
struct { 
    SDL_Window *window;
    SDL_Texture *texture;
    SDL_Renderer *renderer;
    u32 pixels[SCREEN_SIZE_X * SCREEN_SIZE_Y];
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


static void tick(f32 dt) {

    v2 input_dir = {0.0, 0.0}; // In the body frame, which is right-handed, so y points left.
    if (is_pressed(state.keyboard_state.w)) {
        input_dir.x += 1.0;
    }
    if (is_pressed(state.keyboard_state.s)) {
        input_dir.x -= 1.0;
    }
    if (is_pressed(state.keyboard_state.d)) {
        input_dir.y -= 1.0;
    }
    if (is_pressed(state.keyboard_state.a)) {
        input_dir.y += 1.0;
    }

    int input_rot_dir = 0; // Right-hand rotation in plane (CCW)
    if (is_pressed(state.keyboard_state.q)) {
        input_rot_dir += 1;
    }
    if (is_pressed(state.keyboard_state.e)) {
        input_rot_dir -= 1;
    }

    if (is_pressed(state.keyboard_state.three)) {
        state.camera_z *= 0.95;
        printf("camera z: %.3f\n", state.camera_z);
    }
    if (is_pressed(state.keyboard_state.four)) {
        state.camera_z /= 0.95;
        printf("camera z: %.3f\n", state.camera_z);
    }
    if (is_pressed(state.keyboard_state.five)) {
        state.camera_height *= 0.95;
        printf("camera height: %.3f\n", state.camera_height);
    }
    if (is_pressed(state.keyboard_state.six)) {
        state.camera_height /= 0.95;
        printf("camera height: %.3f\n", state.camera_height);
    }
    if (is_pressed(state.keyboard_state.seven)) {
        state.camera_width *= 0.95;
        printf("camera width: %.3f\n", state.camera_width);
    }
    if (is_pressed(state.keyboard_state.eight)) {
        state.camera_width /= 0.95;
        printf("camera width: %.3f\n", state.camera_width);
    }

    // Update the player's velocity
    const f32 kPlayerInputAccel = 4.5;
    const f32 kPlayerInputAngularAccel = 7.5;
    const f32 kPlayerMaxSpeed = 3.0;
    const f32 kPlayerMaxOmega = 3.0;
    const f32 kAirFriction = 0.9;
    const f32 kAirFrictionRot = 0.85;

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
    float theta = atan2(state.camera_dir.y, state.camera_dir.x);
    theta += state.player_omega * dt;
    state.camera_dir = ((v2) {cos(theta), sin(theta)});   
    state.camera_dir_rotr = rotr((state.camera_dir));

    // Apply air friction
    state.player_speed.x *= kAirFriction;
    state.player_speed.y *= kAirFriction;
    state.player_omega *= kAirFrictionRot;
}


// Fill all pixels in the vertical line at x between y0 and y1 with the given color.
static void draw_column(int x, int y0, int y1, u32 color) {
    for (int y = y0; y <= y1; y++) {
        state.pixels[(y * SCREEN_SIZE_X) + x] = color;
    }
}

static void render() {
    const u32 color_floor = 0xFF666666;
    const u32 color_ceil = 0xFF444444;

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
        const v2 dir = {cp.x / cam_len, cp.y  /cam_len};

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
            if (MAPDATA[y_ind*8 + x_ind] > 0) {
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

        // Calculate the pixel bounds that we fill the wall in for
        int y_lo = (int)(SCREEN_SIZE_Y/2.0f - cam_len*state.camera_z/ray_len * SCREEN_SIZE_Y / state.camera_height);
        int y_hi = (int)(SCREEN_SIZE_Y/2.0f + cam_len*(WALL_HEIGHT - state.camera_z)/ray_len * SCREEN_SIZE_Y / state.camera_height);
        int y_lo_capped = max(y_lo, 0);
        int y_hi_capped = min(y_hi, SCREEN_SIZE_Y-1);

        // u32 color_wall_to_render = (dx_ind == 0) ? color_wall[MAPDATA[y_ind*8 + x_ind]-1] : color_wall_light[MAPDATA[y_ind*8 + x_ind]-1];

        draw_column(x, 0, y_lo_capped-1, color_floor);
        {
            u32 texture_x_offset = 0;
            u32 texture_y_offset = 0;
            if (dx_ind == 0) {
                // Draw the light version
                texture_x_offset = 64;
            }

            f32 rem = 0.0f;
            if (dx_ind == 0) {
                rem = dy_ind < 0 ? TILE_WIDTH - x_rem : x_rem;
            } else {
                rem = dx_ind < 0 ? y_rem : TILE_WIDTH - y_rem;
            }
            u32 texture_x = (int) (64 * rem / TILE_WIDTH);
            u32 baseline = texture_y_offset + (texture_x+texture_x_offset)*bitmap.n_pixels_per_column;
            for (int y = y_hi_capped; y >= y_lo_capped; y--) {
                u32 texture_y = (y_hi - y) * 64 / max(1, y_hi - y_lo);
                u32 color = bitmap.abgr[texture_y+baseline];
                state.pixels[(y * SCREEN_SIZE_X) + x] = color;
            }
        }
        draw_column(x, y_hi_capped + 1, SCREEN_SIZE_Y-1, color_ceil);
    }
}

int main(int argc, char *argv[]) {

    // Load our assets
    printf("Loading assets...");
    {
        FILE* ptr = fopen("assets/assets.bin","rb");
        ASSERT(ptr, "Error opening assets file\n");

        // Skip the header bytes
        fseek(ptr, 4, SEEK_CUR);
        ASSERT(fread(&bitmap.n_pixels,            sizeof(u32), 1, ptr) == 1, "Failed to read n_pixels when loading assets\n");
        ASSERT(fread(&bitmap.n_pixels_per_column, sizeof(u32), 1, ptr) == 1, "Failed to read n_pixels_per_column when loading assets\n");

        size_t n_bytes_to_alloc = bitmap.n_pixels * sizeof(u32);
        bitmap.abgr = (u32*) malloc(n_bytes_to_alloc);
        ASSERT(bitmap.abgr, "Failed to allocate bitmap data\n");
        ASSERT(fread(bitmap.abgr, sizeof(u32), bitmap.n_pixels, ptr) == bitmap.n_pixels, "Failed to read data when loading assets\n");

        fclose(ptr);
    }
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
    clear_keyboard_state(&state.keyboard_state);

    // Time structs
    struct timeval timeval_start, timeval_end;

    // Main loop
    u32 time_prev_tick = SDL_GetTicks();
    state.quit = 0;
    while (state.quit == 0) {
        const u32 time_start = SDL_GetTicks();
        gettimeofday(&timeval_start, NULL);

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

        // TODO: Move to more accurate timing?
        const u32 time_tick_start = SDL_GetTicks();
        const f32 dt = (time_tick_start - time_prev_tick) / 1000.0f;
        tick(dt);
        time_prev_tick = time_tick_start;

        render();

        decay_keyboard_state(&state.keyboard_state);

        // Get timer end for all the non-SDL stuff
        gettimeofday(&timeval_end, NULL);
        f64 game_ms_elapsed = (timeval_end.tv_sec - timeval_start.tv_sec) * 1000.0;  // sec to ms
        game_ms_elapsed += (timeval_end.tv_usec - timeval_start.tv_usec) / 1000.0;   // us to ms
        printf("Game: %.3f ms, %.1f fps\n", game_ms_elapsed, 1000.0f / max(1.0f, game_ms_elapsed));

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

        const u32 time_end = SDL_GetTicks();
        const u32 ms_elapsed = time_end - time_start;
        const f32 fps = 1000.0f / max(1, ms_elapsed);
        printf("FPS: %.1f\n", fps);
    }

    SDL_DestroyWindow(state.window);
    
    // Free our assets
    free(bitmap.abgr);

    return 0;
}