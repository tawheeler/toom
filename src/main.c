#include <math.h>
#include <stdio.h>
#include <stdint.h>
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

#define SCREEN_SIZE_X 320
#define SCREEN_SIZE_Y 180

#define TILE_WIDTH 1.0f
#define WALL_HEIGHT 1.2f

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
    f32 camera_depth;
    f32 camera_z;

    struct KeyBoardState keyboard_state;
} state;


static void tick(f32 dt) {

    // Process player input

    const f32 kPlayerStepPerSec = 1.25;
    if (is_pressed(state.keyboard_state.w)) {
        state.camera_pos.x += kPlayerStepPerSec * dt * state.camera_dir.x;
        state.camera_pos.y += kPlayerStepPerSec * dt * state.camera_dir.y;
    }
    if (is_pressed(state.keyboard_state.s)) {
        state.camera_pos.x -= kPlayerStepPerSec * dt * state.camera_dir.x;
        state.camera_pos.y -= kPlayerStepPerSec * dt * state.camera_dir.y;
    }
    if (is_pressed(state.keyboard_state.d)) {
        state.camera_pos.x -= kPlayerStepPerSec * dt * state.camera_dir_rotr.x;
        state.camera_pos.y -= kPlayerStepPerSec * dt * state.camera_dir_rotr.y;
    }
    if (is_pressed(state.keyboard_state.a)) {
        state.camera_pos.x += kPlayerStepPerSec * dt * state.camera_dir_rotr.x;
        state.camera_pos.y += kPlayerStepPerSec * dt * state.camera_dir_rotr.y;
    }

    const f32 kPlayerRotPerSec = 1.5;
    if (is_pressed(state.keyboard_state.q)) {
        // Right-hand rotation in plane (CCW)
        float theta = atan2(state.camera_dir.y, state.camera_dir.x);
        theta += kPlayerRotPerSec * dt;
        state.camera_dir = ((v2) {cos(theta), sin(theta)});   
        state.camera_dir_rotr = rotr((state.camera_dir));
    }
    if (is_pressed(state.keyboard_state.e)) {
        // Left-hand rotation in plane (CW)
        float theta = atan2(state.camera_dir.y, state.camera_dir.x);
        theta -= kPlayerRotPerSec * dt;
        state.camera_dir = ((v2) {cos(theta), sin(theta)});   
        state.camera_dir_rotr = rotr((state.camera_dir));
    }

    if (is_pressed(state.keyboard_state.one)) {
        state.camera_depth *= 0.95;
        printf("camera depth: %.3f\n", state.camera_depth);
    }
    if (is_pressed(state.keyboard_state.two)) {
        state.camera_depth /= 0.95;
        printf("camera depth: %.3f\n", state.camera_depth);
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
}


// Fill all pixels in the vertical line at x between y0 and y1 with the given color.
static void verline(int x, int y0, int y1, u32 color) {
    for (int y = y0; y <= y1; y++) {
        state.pixels[(y * SCREEN_SIZE_X) + x] = color;
    }
}

static void render() {
    static u32 color_wall[4] = {
        0xFFFF0000,
        0xFF00FF00,
        0xFF00FFFF,
        0xFF0000FF
    };
    static u32 color_wall_light[4] = {
        0xFFFF3333,
        0xFF33FF33,
        0xFF33FFFF,
        0xFF3333FF
    };
    const u32 color_floor = 0xFF666666;
    const u32 color_ceil = 0xFF444444;

    for (int x = 0; x < SCREEN_SIZE_X; x++) {
        
        // The column's point in game space
        const f32 dw = state.camera_width/2 - (state.camera_width*x)/SCREEN_SIZE_X;
        const v2 p = {
            state.camera_pos.x + state.camera_depth*state.camera_dir.x + dw*state.camera_dir_rotr.x,
            state.camera_pos.y + state.camera_depth*state.camera_dir.y + dw*state.camera_dir_rotr.y
        };

        // Camera to pixel column
        const v2 cp = {p.x - state.camera_pos.x, p.y - state.camera_pos.y};

        // Distance from the camera to the column
        const f32 cam_len = length( (cp) );
        
        // Ray direction through this column
        const v2 dir = {cp.x / cam_len, cp.y  /cam_len};

        // Get p's cell
        int x_ind = (int)(floorf(p.x / TILE_WIDTH));
        int y_ind = (int)(floorf(p.y / TILE_WIDTH));
        f32 x_rem = p.x - TILE_WIDTH*x_ind;
        f32 y_rem = p.y - TILE_WIDTH*y_ind;

        // Step through cells until we hit an occupied cell
        int n_steps = 0;
        while (n_steps < 100) {
            n_steps += 1;

            // x(t) = x_rem + dir.x * dt
            // y(t) = y_rem + dir.y * dt

            // We cross x = 0 if dir.x < 0, at dt = -x_rem/dir.x
            // We cross x = 1 if dir.x > 0, at dt = (1-x_rem)/dir.x
            // We cross y = 0 if dir.y < 0, at dt = -y_rem/dir.y
            // We cross y = 1 if dir.y > 0, at dt = (1-y_rem)/dir.y

            int dx_ind = 0;
            int dy_ind = 0;
            f32 dt_best = 999.0;
            if (dir.x < 0) {
                f32 dt = -x_rem/dir.x;
                if (dt < dt_best) {
                    dt_best = dt;
                    dx_ind = -1;
                    dy_ind =  0;
                }
            } else if (dir.x > 0) {
                f32 dt = (1-x_rem)/dir.x;
                if (dt < dt_best) {
                    dt_best = dt;
                    dx_ind = 1;
                    dy_ind = 0;
                }
            }
            if (dir.y < 0) {
                f32 dt = -y_rem/dir.y;
                if (dt < dt_best) {
                    dt_best = dt;
                    dx_ind =  0;
                    dy_ind = -1;
                }
            } else if (dir.y > 0) {
                f32 dt = (1-y_rem)/dir.y;
                if (dt < dt_best) {
                    dt_best = dt;
                    dx_ind = 0;
                    dy_ind = 1;
                }
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
        y_lo = max(y_lo, 0);
        y_hi = min(y_hi, SCREEN_SIZE_Y-1);

        // TODO: This is not correct
        u32 color_wall_to_render = x_rem > y_rem ? color_wall[MAPDATA[y_ind*8 + x_ind]-1] : color_wall_light[MAPDATA[y_ind*8 + x_ind]-1];

        verline(x, 0, y_lo-1, color_floor);
        verline(x, y_lo, y_hi, color_wall_to_render);
        verline(x, y_hi + 1, SCREEN_SIZE_Y-1, color_ceil);
    }
}

int main(int argc, char *argv[]) {
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
    state.camera_width = 0.6f;
    state.camera_height = state.camera_width * SCREEN_SIZE_Y / SCREEN_SIZE_X;
    state.camera_depth = 0.4f;
    state.camera_z = 0.4;

    // Init keyboard
    clear_keyboard_state(&state.keyboard_state);

    // Main loop
    u32 time_prev_tick = SDL_GetTicks();
    state.quit = 0;
    while (state.quit == 0) {
        const u32 time_start = SDL_GetTicks();

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

        const u32 time_tick_start = SDL_GetTicks();
        const f32 dt = (time_tick_start - time_prev_tick) / 1000.0f;
        tick(dt);
        time_prev_tick = time_tick_start;

        render();
        decay_keyboard_state(&state.keyboard_state);

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
        // printf("FPS: %.1f\n", fps);
    }

    SDL_DestroyWindow(state.window);
    return 0;
}