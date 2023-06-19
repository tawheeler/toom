#include "vec.h"

v2 add(v2 a, v2 b) {
    v2 retval = {a.x + b.x, a.y + b.y};
    return retval;
}

v2 sub(v2 a, v2 b) {
    v2 retval = {a.x - b.x, a.y - b.y};
    return retval;
}

v2 VectorProjection(v2 a, v2 b) {
    f32 scale = dot(a, b) / dot(b, b);
    v2 retval = {scale * b.x, scale * b.y};
    return retval;
}