#include "game.h" 

#include <stdio.h>


void Tick(
    struct GameState* state,
    f32 dt,
    struct KeyBoardState* keyboard_state
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

    // TODO - move out
    // if (IsNewlyPressed(keyboard_state->r)) {
    //     printf("Reloading assets.\n");
    //     LoadAssets();
    //     printf("DONE.\n");
    // }

    if (IsPressed(keyboard_state->three)) {
        state->camera.z *= 0.95;
        printf("camera z: %.3f\n", state->camera.z);
    }
    if (IsPressed(keyboard_state->four)) {
        state->camera.z /= 0.95;
        printf("camera z: %.3f\n", state->camera.z);
    }

    // Update the player's velocity
    const f32 kPlayerInputAccel = 6.5;
    const f32 kPlayerInputAngularAccel = 8.5;
    const f32 kPlayerMaxSpeed = 5.0;
    const f32 kPlayerMaxOmega = 5.0;
    const f32 kAirFriction = 4.0;
    const f32 kAirFrictionRot = 4.0;

    // Note: Speed is in the global frame
    state->player_speed.x += (state->camera.dir.x*input_dir.x - state->camera.dir.y*input_dir.y) * kPlayerInputAccel * dt;
    state->player_speed.y += (state->camera.dir.y*input_dir.x + state->camera.dir.x*input_dir.y) * kPlayerInputAccel * dt;
    state->player_omega += input_rot_dir * kPlayerInputAngularAccel * dt;

    // Clamp the velocity to a maximum magnitude
    f32 speed = length(state->player_speed);
    if (speed > kPlayerMaxSpeed) {
        state->player_speed.x *= kPlayerMaxSpeed / speed;
        state->player_speed.y *= kPlayerMaxSpeed / speed;
    }
    if (state->player_omega > kPlayerMaxOmega) {
        state->player_omega *= kPlayerMaxOmega / state->player_omega;
    } else if (state->player_omega < -kPlayerMaxOmega) {
        state->player_omega *= - kPlayerMaxOmega / state->player_omega;
    }

    // Update the player's position
    state->camera.pos.x += state->player_speed.x * dt;
    state->camera.pos.y += state->player_speed.y * dt;

    // Update the player's rotational heading
    f32 theta = atan2(state->camera.dir.y, state->camera.dir.x);
    theta += state->player_omega * dt;
    state->camera.dir = ((v2) {cos(theta), sin(theta)});

    // Apply air friction
    f32 air_friction_decay = exp(-kAirFriction * dt);
    state->player_speed.x *= air_friction_decay;
    state->player_speed.y *= air_friction_decay;
    state->player_omega *= exp(-kAirFrictionRot * dt);
}