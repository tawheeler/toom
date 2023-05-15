#ifndef GAME_H_INCLUDED
#define GAME_H_INCLUDED

#include "typedefs.h"
#include "vec.h"
#include "input.h"

struct CameraState {
    v2 pos; // camera location in game frame
    v2 dir; // unit direction that the camera is facing
    v2 fov; // width (x) and height (y) of the field of view at unit distance from the camera
    f32 z;  // height over the ground [m]
};

struct GameState {
    struct CameraState camera;
    
    v2 player_speed;
    f32 player_omega;    
};

// Propagate the game state forward by dt seconds.
void Tick(struct GameState* state, f32 dt, struct KeyBoardState* keyboard_state);

#endif