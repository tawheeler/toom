#ifndef BITMAP_H_INCLUDED
#define BITMAP_H_INCLUDED

#include "typedefs.h"

struct Bitmap {
    u32  n_pixels;
    u32  n_pixels_per_column;
    u32  n_pixels_per_row;
    bool column_major;
    u32* abgr;
};

int GetColumnMajorPixelIndex(struct Bitmap* bitmap, int x, int y);
int GetRowMajorPixelIndex(struct Bitmap* bitmap, int x, int y);
u32 GetColumnMajorPixelAt(struct Bitmap* bitmap, int x, int y);
u32 GetRowMajorPixelAt(struct Bitmap* bitmap, int x, int y);

#endif