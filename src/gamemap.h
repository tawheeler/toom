#ifndef GAMEMAP_H_INCLUDED
#define GAMEMAP_H_INCLUDED

#include "typedefs.h"
#include "delaunay_mesh.h"

struct TextureInfo
{
    u16 texture_id; // Index in the texture atlas
    i16 x_offset;   // Texture x offset
    i16 y_offset;   // Texture y offset
};

#define SIDEINFO_FLAG_DARK 1
#define SIDEINFO_FLAG_PASSABLE 2

struct SideInfo
{
    u16 flags;
    u16 sector_id;                          // Index into the sector list
    struct TextureInfo texture_info_lower;  // Texture displayed if the floor height increases
    struct TextureInfo texture_info_middle; // Texture displayed if the wall is solid
    struct TextureInfo texture_info_upper;  // Texture displayed if the floor ceiling decreases
};

// All side infos refer to a sector to contain the corresponding floor/ceiling information.
struct Sector
{
    u16 flags;
    u16 flat_id_floor;
    u16 flat_id_ceil;
    u8 light_level; // Colormap index. 0 is brightest.
    f32 z_floor;
    f32 z_ceil;
};

struct GameMap
{
    struct DelaunayMesh *geometry_mesh;

    u32 n_side_infos;
    struct SideInfo *side_infos;
    u16 *quarter_edge_index_to_side_info_index;

    u32 n_sectors;
    struct Sector *sectors;
};

#endif