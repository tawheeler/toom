#include "game.h"

#include <stdio.h>

void Tick(
    struct GameState *state,
    f32 dt,
    const struct KeyBoardState *keyboard_state,
    const struct GameMap *game_map)
{

    // ------------------------------------------------------------
    // Process Input
    v2 input_dir = {0.0, 0.0}; // In the body frame, which is right-handed, so y points left.
    if (IsPressed(keyboard_state->w))
    {
        input_dir.x += 1.0;
    }
    if (IsPressed(keyboard_state->s))
    {
        input_dir.x -= 1.0;
    }
    if (IsPressed(keyboard_state->d))
    {
        input_dir.y -= 1.0;
    }
    if (IsPressed(keyboard_state->a))
    {
        input_dir.y += 1.0;
    }

    int input_rot_dir = 0; // Right-hand rotation in plane (CCW)
    if (IsPressed(keyboard_state->q))
    {
        input_rot_dir += 1;
    }
    if (IsPressed(keyboard_state->e))
    {
        input_rot_dir -= 1;
    }

    if (IsPressed(keyboard_state->three))
    {
        state->player.z *= 0.95;
        printf("player z: %.3f\n", state->player.z);
    }
    if (IsPressed(keyboard_state->four))
    {
        state->player.z /= 0.95;
        printf("player z: %.3f\n", state->player.z);
    }

    // ------------------------------------------------------------
    // Update the player's velocity and angular speed
    const f32 kPlayerInputAccel = 6.5;
    const f32 kPlayerInputAngularAccel = 8.5;
    const f32 kPlayerMaxSpeed = 5.0;
    const f32 kPlayerMaxOmega = 5.0;
    const f32 kAirFriction = 4.0;
    const f32 kAirFrictionRot = 4.0;
    const f32 kSlidingFriction = 0.95;

    // Note: Speed is in the global frame
    state->player.vel.x += (state->player.dir.x * input_dir.x - state->player.dir.y * input_dir.y) * kPlayerInputAccel * dt;
    state->player.vel.y += (state->player.dir.y * input_dir.x + state->player.dir.x * input_dir.y) * kPlayerInputAccel * dt;
    state->player.omega += input_rot_dir * kPlayerInputAngularAccel * dt;

    // Clamp the velocity to a maximum magnitude
    f32 speed = length(state->player.vel);
    if (speed > kPlayerMaxSpeed)
    {
        state->player.vel.x *= kPlayerMaxSpeed / speed;
        state->player.vel.y *= kPlayerMaxSpeed / speed;
    }
    if (state->player.omega > kPlayerMaxOmega)
    {
        state->player.omega *= kPlayerMaxOmega / state->player.omega;
    }
    else if (state->player.omega < -kPlayerMaxOmega)
    {
        state->player.omega *= -kPlayerMaxOmega / state->player.omega;
    }

    // ------------------------------------------------------------
    // Move the player.
    f32 dt_remaining = dt;
    while (dt_remaining > 0.0f)
    {
        v2 player_pos_delta = {dt_remaining * state->player.vel.x, dt_remaining * state->player.vel.y};
        v2 player_pos_next = add(state->player.pos, player_pos_delta);

        // Exit if the step is too small, to avoid math problems.
        if (sqnorm(player_pos_delta) < 1e-6)
        {
            break;
        }

        // Check to see if we cross over any mesh edges.
        QuarterEdge *qe_ab = state->player.qe_geometry->rot;
        QuarterEdge *qe_bc = state->player.qe_geometry->next->rot;
        QuarterEdge *qe_ca = state->player.qe_geometry->next->next->rot;

        const v2 a = *(qe_ab->vertex);
        const v2 b = *(qe_bc->vertex);
        const v2 c = *(qe_ca->vertex);

        // The fraction of the distance we will move
        f32 interp = 1.0f;
        QuarterEdge *qe_dual_new_triangle = NULL;
        QuarterEdge *qe_side = NULL;
        v2 v_face = {0.0, 0.0};

        if (GetRightHandedness(&a, &b, &player_pos_next) < -1e-4)
        {
            // We would cross AB
            v2 v = sub(b, a);
            v2 w = sub(state->player.pos, a);
            f32 interp_ab = cross(v, w) / cross(player_pos_delta, v);
            if (interp_ab < interp)
            {
                interp = interp_ab;
                qe_dual_new_triangle = qe_ab->rot;
                qe_side = qe_ab;
                v_face = v;
            }
        }
        if (GetRightHandedness(&b, &c, &player_pos_next) < -1e-4)
        {
            // We would cross BC
            v2 v = sub(c, b);
            v2 w = sub(state->player.pos, b);
            f32 interp_bc = cross(v, w) / cross(player_pos_delta, v);
            if (interp_bc < interp)
            {
                interp = interp_bc;
                qe_dual_new_triangle = qe_bc->rot;
                qe_side = qe_bc;
                v_face = v;
            }
        }
        if (GetRightHandedness(&c, &a, &player_pos_next) < -1e-4)
        {
            // We would cross CA
            v2 v = sub(a, c);
            v2 w = sub(state->player.pos, c);
            f32 interp_ca = cross(v, w) / cross(player_pos_delta, v);
            if (interp_ca < interp)
            {
                interp = interp_ca;
                qe_dual_new_triangle = qe_ca->rot;
                qe_side = qe_ca;
                v_face = v;
            }
        }

        // Move the player
        // If we would have crossed any edge, this merely moves us up to the boundary instead
        state->player.pos.x += player_pos_delta.x * interp;
        state->player.pos.y += player_pos_delta.y * interp;
        dt_remaining *= (1.0f - interp);
        dt_remaining -= 1e-5; // eps factor for safety

        if (qe_dual_new_triangle != NULL)
        {
            // We would have crossed into another triangle.

            // Determine whether to continue on or stop at the edge.
            bool stop_at_edge = 0;
            bool new_sector_z = 0;
            f32 sector_z = 0.0f;
            u16 side_info_index = game_map->quarter_edge_index_to_side_info_index[qe_side->index];
            if (side_info_index != 0xFFFF)
            {

                // TODO: Prevent proceeding if the z-difference is sufficiently large
                struct SideInfo *side_info = game_map->side_infos + side_info_index;
                const bool is_passable =
                    ((side_info->flags & SIDEINFO_FLAG_PASSABLE) > 0);
                stop_at_edge = !is_passable;

                QuarterEdge *qe_sym = QESym(qe_side);
                u32 side_info_index_sym = game_map->quarter_edge_index_to_side_info_index[qe_sym->index];
                if (side_info_index_sym != 0xFFFF)
                {
                    struct SideInfo *side_info_sym = game_map->side_infos + side_info_index_sym;
                    struct Sector *sector_sym = game_map->sectors + side_info_sym->sector_id;
                    sector_z = sector_sym->z_floor;
                    new_sector_z = 1;
                }
            }

            if (stop_at_edge)
            {
                // The new triangle is solid, so do not change triangles.
                // Lose all velocity into the boundary surface.
                // This results in the vector along the face.
                state->player.vel = VectorProjection(state->player.vel, v_face);
                state->player.vel.x *= kSlidingFriction;
                state->player.vel.y *= kSlidingFriction;
            }
            else
            {
                // Accept the new triangle
                state->player.qe_geometry = qe_dual_new_triangle;
                if (new_sector_z)
                {
                    state->player.z = sector_z + state->player.height;
                }
            }
        }
    }

    // Update the player's rotational heading
    f32 theta = atan2(state->player.dir.y, state->player.dir.x);
    theta += state->player.omega * dt;
    state->player.dir = ((v2){cos(theta), sin(theta)});

    // Apply air friction
    f32 air_friction_decay = exp(-kAirFriction * dt);
    state->player.vel.x *= air_friction_decay;
    state->player.vel.y *= air_friction_decay;
    state->player.omega *= exp(-kAirFrictionRot * dt);
}