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
#define WALL_HEIGHT 0.7f

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
    1, 0, 0, 3, 0, 0, 0, 1,
    1, 0, 0, 0, 0, 0, 0, 1,
    1, 0, 0, 0, 0, 0, 4, 1,
    1, 0, 2, 0, 0, 0, 0, 1,
    1, 0, 0, 0, 0, 0, 0, 1,
    1, 1, 1, 1, 1, 1, 1, 1,
};

struct { 
    SDL_Window *window;
    SDL_Texture *texture;
    SDL_Renderer *renderer;
    u32 pixels[SCREEN_SIZE_X * SCREEN_SIZE_Y];
    bool quit;

    v2 camera_pos;
    v2 camera_dir;
    f32 camera_width;
    f32 camera_height;
    f32 camera_depth;
    f32 camera_z;
} state;

// Fill all pixels in the vertical line at x between y0 and y1 with the given color.
static void verline(int x, int y0, int y1, u32 color) {
    for (int y = y0; y <= y1; y++) {
        state.pixels[(y * SCREEN_SIZE_X) + x] = color;
    }
}

static void render() {

    const v2 rotr_dir = rotr((state.camera_dir));

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
            state.camera_pos.x + state.camera_depth*state.camera_dir.x + dw*rotr_dir.x,
            state.camera_pos.y + state.camera_depth*state.camera_dir.y + dw*rotr_dir.y
        };

        // Camera to pixel column
        const v2 cp = {p.x - state.camera_pos.x, p.y - state.camera_pos.y};

        // Ray direction through this column
        const v2 dir = normalize( (cp) );

        // Distance from the camera to the column
        const f32 cam_len = length( (cp) );

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
        int y_lo = (int)(SCREEN_SIZE_Y/2.0f - cam_len*state.camera_z/ray_len * SCREEN_SIZE_Y);
        int y_hi = (int)(SCREEN_SIZE_Y/2.0f + cam_len*(WALL_HEIGHT - state.camera_z)/ray_len * SCREEN_SIZE_Y);
        y_lo = max(y_lo, 0);
        y_hi = min(y_hi, SCREEN_SIZE_Y-1);

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
    state.camera_width = 0.2f;
    state.camera_height = state.camera_width * SCREEN_SIZE_Y / SCREEN_SIZE_X;
    state.camera_depth = 0.6f;
    state.camera_z = 0.4;

    // Main loop
    state.quit = 0;
    while (state.quit == 0) {
        SDL_Event windowEvent;
        while (SDL_PollEvent(&windowEvent)) {
            if (windowEvent.type == SDL_QUIT) {
                state.quit = 1;
                break;
            }
        }

        render();

        float theta = atan2(state.camera_dir.y, state.camera_dir.x);
        theta += 0.025f;
        state.camera_dir = ((v2) {cos(theta), sin(theta)});

        SDL_UpdateTexture(state.texture, NULL, state.pixels, SCREEN_SIZE_X * 4);
        SDL_RenderCopyEx(
            state.renderer,
            state.texture,
            NULL,
            NULL,
            0.0,
            NULL,
            SDL_FLIP_VERTICAL);
        SDL_RenderPresent(state.renderer);
    }

    SDL_DestroyWindow(state.window);
    return 0;
}