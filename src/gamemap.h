#ifndef GAMEMAP_H_INCLUDED
#define GAMEMAP_H_INCLUDED

#include "typedefs.h"
#include "delaunay_mesh.h"

#define SIDEINFO_FLAG_DARK 1

struct SideInfo
{
    u16 flags;
    u16 texture_id;
    i16 x_offset;
    i16 y_offset;
};

struct GameMap
{
    struct DelaunayMesh *geometry_mesh;

    u32 n_side_infos;
    struct SideInfo *side_infos;
    u16 *quarter_edge_index_to_side_info_index;
};

#endif