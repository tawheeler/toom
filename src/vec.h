#ifndef VEC_H_INCLUDED
#define VEC_H_INCLUDED

#include <math.h>

#include "typedefs.h"

typedef struct v2_s {
    f32 x;
    f32 y;
} v2;

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
#define clamp(v, lo, hi) ({ __typeof__(v) _v = (v), _lo = (lo), _hi = (hi); _v > _hi ? _hi : (_v < _lo ? _lo : _v); })

#endif