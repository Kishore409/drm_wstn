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
#include <stdint.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include <xf86drm.h>
#include <xf86drmMode.h>
#include <drm_fourcc.h>
#include <drm/drm_mode.h>
#include "weston-debug.h"
#include "drm-ctm.h"

//Reference: https://nick-shaw.github.io/cinematiccolor/common-rgb-color-spaces.html
double m1 = 0.1593017578125;
double m2 = 78.84375;
double c1 = 0.8359375;
double c2 = 18.8515625;
double c3 = 18.6875;

#define DD_MIN(a, b) ((a) < (b) ? (a) : (b))
#define DD_MAX(a, b) ((a) < (b) ? (b) : (a))
#define MAX_24BIT_NUM ((1<<24) -1)

static double matrixdeterminant3x3(double matrix[3][3])
{
        double result;
        result = matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]);
        result -= matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0]);
        result += matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);

        return result;
}

static int matrixinverse3x3(double matrix[3][3], double result[3][3])
{
        int retval = -1;
        double tmp[3][3];
        double determinant = matrixdeterminant3x3(matrix);
        if(determinant)
        {
        tmp[0][0] = (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) / determinant;
        tmp[0][1] = (matrix[0][2] * matrix[2][1] - matrix[2][2] * matrix[0][1]) / determinant;
        tmp[0][2] = (matrix[0][1] * matrix[1][2] - matrix[0][2] * matrix[1][1]) / determinant;
        tmp[1][0] = (matrix[1][2] * matrix[2][0] - matrix[1][0] * matrix[2][2]) / determinant;
        tmp[1][1] = (matrix[0][0] * matrix[2][2] - matrix[0][2] * matrix[2][0]) / determinant;
        tmp[1][2] = (matrix[0][2] * matrix[1][0] - matrix[0][0] * matrix[1][2]) / determinant;
        tmp[2][0] = (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]) / determinant;
        tmp[2][1] = (matrix[0][1] * matrix[2][0] - matrix[0][0] * matrix[2][1]) / determinant;
        tmp[2][2] = (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]) / determinant;

        result[0][0] = tmp[0][0]; result[0][1] = tmp[0][1]; result[0][2] = tmp[0][2];
        result[1][0] = tmp[1][0]; result[1][1] = tmp[1][1]; result[1][2] = tmp[1][2];
        result[2][0] = tmp[2][0]; result[2][1] = tmp[2][1]; result[2][2] = tmp[2][2];

        retval = 0;
        }
        return retval;
}

static void matrixmult3x3with3x1(double matrix1[3][3], double matrix2[3], double result[3])
{
        double tmp[3];
        tmp[0] = matrix1[0][0] * matrix2[0] + matrix1[0][1] * matrix2[1] + matrix1[0][2] * matrix2[2];
        tmp[1] = matrix1[1][0] * matrix2[0] + matrix1[1][1] * matrix2[1] + matrix1[1][2] * matrix2[2];
        tmp[2] = matrix1[2][0] * matrix2[0] + matrix1[2][1] * matrix2[1] + matrix1[2][2] * matrix2[2];

        result[0]=tmp[0];
        result[1]=tmp[1];
        result[2]=tmp[2];

}

static void matrixmult3x3(double matrix1[3][3], double matrix2[3][3], double result[3][3])
{
        double tmp[3][3];
        for(uint8_t y=0; y<3; y++)
        {
                for(uint8_t x=0; x<3; x++)
                {
                        tmp[y][x]= matrix1[y][0] * matrix2[0][x] + matrix1[y][1] * matrix2[1][x] + matrix1[y][2] * matrix2[2][x];
                }

        }
        for(uint8_t y=0; y<3;y++)
        {
                for(uint8_t x=0; x<3; x++)
                {
                        result[y][x] = tmp[y][x];
                }
        }
}

static void creatergb2xyzmatrix(struct colorspace *pcspace, double rgb2xyz[3][3])
{
/*
   http://www.brucelindbloom.com/index.html?eqn_rgb_xyz_matrix.html
*/
        double xyzsum[3];
        double z[4];
        double xyzw[3];
        struct chromaticity *pchroma = &pcspace->white;
        for(uint8_t i=0; i< 4; i++)
        {
                z[i]=1- pchroma[i].x - pchroma[i].y;
        }
        xyzw[0] = pcspace->white.x / pcspace->white.y;
        xyzw[1]=1;
        xyzw[2]=z[0] / pcspace->white.y;

        double xyzrgb[3][3] =
        {
                { pcspace->red.x, pcspace->green.x, pcspace->blue.x },
                { pcspace->red.y, pcspace->green.y, pcspace->blue.y },
                { z[1], z[2], z[3] }
        };
        double mat1[3][3];
        matrixinverse3x3(xyzrgb, mat1);
        matrixmult3x3with3x1(mat1, xyzw, xyzsum);
        double mat2[3][3] = { { xyzsum[0], 0, 0 }, { 0, xyzsum[1], 0 }, { 0, 0, xyzsum[2] } };
        matrixmult3x3(xyzrgb, mat2, rgb2xyz);
}

static void creategamutscalingmatrix(struct colorspace *psrc, struct colorspace *pdst, double result[3][3])
{
        double mat1[3][3], mat2[3][3], tmp[3][3];
        creatergb2xyzmatrix(psrc, mat1);
        creatergb2xyzmatrix(pdst, mat2);

        matrixinverse3x3(mat2, tmp);
        matrixmult3x3(tmp, mat1, result);
}

static void getctmforbt709tobt2020(double result[3][3])
{
/*
 https://en.wikipedia.org/wiki/rec._2020#system_colorimetry
 https://en.wikipedia.org/wiki/rec._709#primary_chromaticities
*/
        struct colorspace bt2020, bt709;

        bt2020.white.x = 0.3127;
        bt2020.white.y = 0.3290;
        bt2020.white.luminance = 100.0;
        bt2020.red.x = 0.708;
        bt2020.red.y = 0.292;
        bt2020.green.x = 0.170;
        bt2020.green.y = 0.797;
        bt2020.blue.x = 0.131;
        bt2020.blue.y = 0.046;

        bt709.white.x = 0.3127;
        bt709.white.y = 0.3290;
        bt709.white.luminance = 100.0;
        bt709.red.x = 0.64;
        bt709.red.y = 0.33;
        bt709.green.x = 0.30;
        bt709.green.y = 0.60;
        bt709.blue.x = 0.15;
        bt709.blue.y = 0.06;
        creategamutscalingmatrix(&bt709, &bt2020, result);
}

static void getctmforbt2020tobt709(double result[3][3])
{
/*
 https://en.wikipedia.org/wiki/rec._2020#system_colorimetry
 https://en.wikipedia.org/wiki/rec._709#primary_chromaticities
*/
        struct colorspace bt2020, bt709;

        bt2020.white.x = 0.3127;
        bt2020.white.y = 0.3290;
        bt2020.white.luminance = 100.0;
        bt2020.red.x = 0.708;
        bt2020.red.y = 0.292;
        bt2020.green.x = 0.170;
        bt2020.green.y = 0.797;
        bt2020.blue.x = 0.131;
        bt2020.blue.y = 0.046;

        bt709.white.x = 0.3127;
        bt709.white.y = 0.3290;
        bt709.white.luminance = 100.0;
        bt709.red.x = 0.64;
        bt709.red.y = 0.33;
        bt709.green.x = 0.30;
        bt709.green.y = 0.60;
        bt709.blue.x = 0.15;
        bt709.blue.y = 0.06;
        creategamutscalingmatrix(&bt2020, &bt709, result);
}

static void
getctmforbt2020todcip3(double result[3][3])
{

/*
* https://en.wikipedia.org/wiki/rec._2020#system_colorimetry
* https://en.wikipedia.org/wiki/dci-p3#system_colorimetry
*/
        struct colorspace bt2020, dcip3;

        bt2020.white.x = 0.3127;
        bt2020.white.y = 0.3290;
        bt2020.white.luminance = 100.0;

        bt2020.red.x = 0.708;
        bt2020.red.y = 0.292;
        bt2020.green.x = 0.170;
        bt2020.green.y = 0.797;
        bt2020.blue.x = 0.131;
        bt2020.blue.y = 0.046;

        dcip3.white.x = 0.314;
        dcip3.white.y = 0.351;
        dcip3.white.luminance = 100.0;

        dcip3.red.x = 0.680;
        dcip3.red.y = 0.320;
        dcip3.green.x = 0.265;
        dcip3.green.y = 0.690;
        dcip3.blue.x = 0.150;
        dcip3.blue.y = 0.060;

        creategamutscalingmatrix(&bt2020, &dcip3, result);
}

static void
getctmfor709todcip3(double result[3][3])
{
/*
* https://en.wikipedia.org/wiki/dci-p3#system_colorimetry
* https://en.wikipedia.org/wiki/rec._709#primary_chromaticities
*/
        struct colorspace bt709, dcip3;

        bt709.white.x = 0.3127;
        bt709.white.y = 0.3290;
        bt709.white.luminance = 100.0;

        bt709.red.x = 0.64;
        bt709.red.y = 0.33;
        bt709.green.x = 0.30;
        bt709.green.y = 0.60;
        bt709.blue.x = 0.15;
        bt709.blue.y = 0.06;

        dcip3.white.x = 0.314;
        dcip3.white.y = 0.351;
        dcip3.white.luminance = 100.0;

        dcip3.red.x = 0.680;
        dcip3.red.y = 0.320;
        dcip3.green.x = 0.265;
        dcip3.green.y = 0.690;
        dcip3.blue.x = 0.150;
        dcip3.blue.y = 0.060;

        creategamutscalingmatrix(&bt709, &dcip3, result);
}


void (*get_ctm_funptrs[DRM_COLORSPACE_MAX][DRM_COLORSPACE_MAX])(double[3][3]) = {
        [DRM_COLORSPACE_REC709][DRM_COLORSPACE_REC709] = NULL,
        [DRM_COLORSPACE_REC709][DRM_COLORSPACE_DCIP3] = getctmfor709todcip3,
        [DRM_COLORSPACE_REC709][DRM_COLORSPACE_REC2020] = getctmforbt709tobt2020,
        [DRM_COLORSPACE_DCIP3][DRM_COLORSPACE_REC709] = NULL,
        [DRM_COLORSPACE_DCIP3][DRM_COLORSPACE_DCIP3] = NULL,
        [DRM_COLORSPACE_DCIP3][DRM_COLORSPACE_REC2020] = NULL,
        [DRM_COLORSPACE_REC2020][DRM_COLORSPACE_REC709] = getctmforbt2020tobt709,
        [DRM_COLORSPACE_REC2020][DRM_COLORSPACE_DCIP3] = getctmforbt2020todcip3,
        [DRM_COLORSPACE_REC2020][DRM_COLORSPACE_REC2020] = NULL,
};

void generatecsc_srctodest(enum drm_colorspace src, enum drm_colorspace dest, double result[3][3])
{
	void (*get_ctm_fptr)(double[3][3]);
	
	get_ctm_fptr =  get_ctm_funptrs[src][dest];
	
	if (!get_ctm_fptr) {
		get_ctm_fptr(result);
	}
}

static double oetf_2084(double input, double srcmaxluminance)
{
        double cf=1.0f;
        double output=0.0f;
        if(input != 0.0f)
        {
                cf = srcmaxluminance / 10000.0;
                input *= cf;
                output = pow(((c1 + (c2 * pow(input, m1))) / (1 + (c3 * pow(input, m1)))), m2);
        }
        return output;
}

/*
reference: https://nick-shaw.github.io/cinematiccolor/common-rgb-color-spaces.html

the st 2084 eotf is an absolute encoding, defined by the following equation:
        l=10000×{max(v^1∕m2−c1,0)/(c2−c3×v^1∕m2)}^1∕m1
where constants are:
        m1 = 2610∕16384
        m2 = 2523∕4096×128
        c1 = 3424∕4096
        c2 = 2413∕4096×32
        c3 = 2392∕4096×32
*/
static double eotf_2084(double input)
{
        double output = 0.0f;
        if(input != 0.0f)
        {
                output = pow(((fmax((pow(input, (1.0/m2)) - c1), 0)) / (c2 - (c3* pow(input, (1.0 / m2))))), (1.0/m1));
        }
        return output;
}

void generateoetf2084lut(struct drm_color_lut_ext *lut, int lut_size, uint16_t max_val)
{
        for (int i=0; i<lut_size; i++)
        {
                double normalized_input = (double) i / (double)(lut_size - 1);
                lut[i].red =  (double)max_val * oetf_2084(normalized_input, 10000.0) + 0.5 ;
		if (lut[i].red > max_val)
                        lut[i].red = max_val;

                lut[i].green = lut[i].blue = lut[i].red;
        }
}

void generateeotf2084lut(struct drm_color_lut_ext *lut, int lut_size, uint16_t max_val)
{
        for (int i=0; i<lut_size; i++)
        {
               double normalized_input = (double) i / (double)(lut_size - 1);
                lut[i].red =  (double)max_val * eotf_2084(normalized_input) + 0.5 ;
		if (lut[i].red > max_val)
                        lut[i].red = max_val;

                lut[i].green = lut[i].blue = lut[i].red;
        }
}

static double getsrgbdecodingvalue(double input)
{
/*
 https://en.wikipedia.org/wiki/srgb#the_forward_transformation_.28cie_xyy_or_cie_xyz_to_srgb.29
*/
        double output = 0.0f;
        if (input <= 0.004045f)
               output =  input / 12.92f;
        else
                output = pow(((input + 0.055)/1.055), 2.4);
        return output;
}

void generatesrgbdegammalut(struct drm_color_lut_ext *lut, int lut_size, uint16_t max_val)
{
	int i;

	for (i = 0; i < lut_size; i++) {
                double normalized_input = (double)i / (double)(lut_size - 1);

                lut[i].red = (double)max_val * getsrgbdecodingvalue(normalized_input) + 0.5;
                if (lut[i].red > max_val)
                        lut[i].red = max_val;
                lut[i].green = lut[i].blue = lut[i].red;
        }
}


static double getsrgbencodingvalue(double input)
{
       double output = 0.0f;
       if(input <= 0.0031308f)
               output = input * 12.92;
       else
               output = (1.055 * pow(input,1.0/2.4)) - 0.055;
       return output;
}

void generatesrgbgammalut(struct drm_color_lut_ext *lut, int lut_size, uint16_t max_val)
{
	int i;

	for (i = 0; i < lut_size; i++) {
                double normalized_input = (double)i / (double)(lut_size - 1);

                lut[i].red = (double)max_val * getsrgbencodingvalue(normalized_input) + 0.5;

                if (lut[i].red > max_val)
                        lut[i].red = max_val;

                lut[i].green = lut[i].blue = lut[i].red;
        }

}
