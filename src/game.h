#ifndef GAME_H_INCLUDED
#define GAME_H_INCLUDED

#include "typedefs.h"
#include "vec.h"
#include "input.h"
#include "delaunay_mesh.h"

struct PlayerState {
    v2 pos;    // location in game coordinate frame
    v2 dir;    // unit direction that the camera is facing
    v2 vel;    // linear rate of change
    f32 omega; // angular rate of change [rad/s]
    f32 z;    // height over the ground
};

struct GameState {
    struct PlayerState player;
};

// Propagate the game state forward by dt seconds.
void Tick(struct GameState* state, f32 dt, const struct KeyBoardState* keyboard_state);

#endif