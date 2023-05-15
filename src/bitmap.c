#include "bitmap.h"

inline int GetColumnMajorPixelIndex(struct Bitmap* bitmap, int x, int y) {
    return y + x*bitmap->n_pixels_per_column;
}

inline int GetRowMajorPixelIndex(struct Bitmap* bitmap, int x, int y) {
    return x + y*bitmap->n_pixels_per_row;
}

inline u32 GetColumnMajorPixelAt(struct Bitmap* bitmap, int x, int y) {
    return bitmap->abgr[y + x*bitmap->n_pixels_per_column];
}

inline u32 GetRowMajorPixelAt(struct Bitmap* bitmap, int x, int y) {
    return bitmap->abgr[x + y*bitmap->n_pixels_per_row];
}