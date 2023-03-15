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

#define SCREEN_SIZE_X 800
#define SCREEN_SIZE_Y 600

struct {
    SDL_Window *window;
    SDL_Texture *texture;
    SDL_Renderer *renderer;
    u32 pixels[SCREEN_SIZE_X * SCREEN_SIZE_Y];
    bool quit;
} state;

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

        state.pixels[(10 * SCREEN_SIZE_X) + 5] = 0xFFFF00FF;

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