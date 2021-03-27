/*
 * Copyright © 2021 Intel Corporation
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice (including the
 * next paragraph) shall be included in all copies or substantial
 * portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT.  IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
 * ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef WESTON_DRM_CTM_H
#define WESTON_DRM_CTM_H

#include <drm/drm_mode.h>

enum drm_colorspace {
        DRM_COLORSPACE_INVALID,
        DRM_COLORSPACE_REC709,
        DRM_COLORSPACE_DCIP3,
        DRM_COLORSPACE_REC2020,
        DRM_COLORSPACE_MAX,
};

struct chromaticity {
        double x;               // CIE1931 x
        double y;               // CIE1931 y
        double luminance;       // CIE1931 Y
};

struct colorspace {
        struct chromaticity white;
        struct chromaticity red;
        struct chromaticity green;
        struct chromaticity blue;
};

void generatecsc_srctodest(enum drm_colorspace src, enum drm_colorspace dest, double result[3][3]);
void generateoetf2084lut(struct drm_color_lut_ext *lut, int lut_size, uint16_t max_val);
void generateeotf2084lut(struct drm_color_lut_ext *lut, int lut_size, uint16_t max_val);
void generatesrgbdegammalut(struct drm_color_lut_ext *lut, int lut_size, uint16_t max_val);
void generatesrgbgammalut(struct drm_color_lut_ext *lut, int lut_size, uint16_t max_val);

#endif /* WESTON_DRM_CTM_H */
