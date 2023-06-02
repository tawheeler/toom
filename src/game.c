#include "game.h" 

#include <stdio.h>


void Tick(
    struct GameState* state,
    f32 dt,
    const struct KeyBoardState* keyboard_state
) {

    v2 input_dir = {0.0, 0.0}; // In the body frame, which is right-handed, so y points left.
    if (IsPressed(keyboard_state->w)) {
        input_dir.x += 1.0;
    }
    if (IsPressed(keyboard_state->s)) {
        input_dir.x -= 1.0;
    }
    if (IsPressed(keyboard_state->d)) {
        input_dir.y -= 1.0;
    }
    if (IsPressed(keyboard_state->a)) {
        input_dir.y += 1.0;
    }

    int input_rot_dir = 0; // Right-hand rotation in plane (CCW)
    if (IsPressed(keyboard_state->q)) {
        input_rot_dir += 1;
    }
    if (IsPressed(keyboard_state->e)) {
        input_rot_dir -= 1;
    }

    if (IsPressed(keyboard_state->three)) {
        state->player.z *= 0.95;
        printf("player z: %.3f\n", state->player.z);
    }
    if (IsPressed(keyboard_state->four)) {
        state->player.z /= 0.95;
        printf("player z: %.3f\n", state->player.z);
    }

    // Update the player's velocity
    const f32 kPlayerInputAccel = 6.5;
    const f32 kPlayerInputAngularAccel = 8.5;
    const f32 kPlayerMaxSpeed = 5.0;
    const f32 kPlayerMaxOmega = 5.0;
    const f32 kAirFriction = 4.0;
    const f32 kAirFrictionRot = 4.0;

    // Note: Speed is in the global frame
    state->player.vel.x += (state->player.dir.x*input_dir.x - state->player.dir.y*input_dir.y) * kPlayerInputAccel * dt;
    state->player.vel.y += (state->player.dir.y*input_dir.x + state->player.dir.x*input_dir.y) * kPlayerInputAccel * dt;
    state->player.omega += input_rot_dir * kPlayerInputAngularAccel * dt;

    // Clamp the velocity to a maximum magnitude
    f32 speed = length(state->player.vel);
    if (speed > kPlayerMaxSpeed) {
        state->player.vel.x *= kPlayerMaxSpeed / speed;
        state->player.vel.y *= kPlayerMaxSpeed / speed;
    }
    if (state->player.omega > kPlayerMaxOmega) {
        state->player.omega *= kPlayerMaxOmega / state->player.omega;
    } else if (state->player.omega < -kPlayerMaxOmega) {
        state->player.omega *= - kPlayerMaxOmega / state->player.omega;
    }

    // Update the player's position
    state->player.pos.x += state->player.vel.x * dt;
    state->player.pos.y += state->player.vel.y * dt;

    // Update the player's rotational heading
    f32 theta = atan2(state->player.dir.y, state->player.dir.x);
    theta += state->player.omega * dt;
    state->player.dir = ((v2) {cos(theta), sin(theta)});

    // Apply air friction
    f32 air_friction_decay = exp(-kAirFriction * dt);
    state->player.vel.x *= air_friction_decay;
    state->player.vel.y *= air_friction_decay;
    state->player.omega *= exp(-kAirFrictionRot * dt);
}