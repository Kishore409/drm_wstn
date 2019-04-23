/*
 * Copyright © 2019 Harish Krupo
 * Copyright © 2019 Intel Corporation
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

#ifndef GL_RENDERER_PRIVATE_H
#define GL_RENDERER_PRIVATE_H

#include <wayland-util.h>
#include <GLES2/gl2.h>
#include <GLES2/gl2ext.h>
#include <stdint.h>
#include <stdbool.h>

enum gl_shader_texture_variant {
	SHADER_VARIANT_NONE = 0,
	SHADER_VARIANT_RGBX,
	SHADER_VARIANT_RGBA,
	SHADER_VARIANT_Y_U_V,
	SHADER_VARIANT_Y_UV,
	SHADER_VARIANT_Y_XUXV,
	SHADER_VARIANT_SOLID,
	SHADER_VARIANT_EXTERNAL,
};

enum gl_shader_degamma_variant {
	SHADER_DEGAMMA_NONE = 0,
	SHADER_DEGAMMA_LINEAR,
	SHADER_DEGAMMA_SRGB,
	SHADER_DEGAMMA_PQ,
	SHADER_DEGAMMA_HLG,
};

enum gl_shader_gamma_variant {
	SHADER_GAMMA_NONE = 0,
	SHADER_GAMMA_LINEAR,
	SHADER_GAMMA_SRGB,
	SHADER_GAMMA_PQ,
	SHADER_GAMMA_HLG,
};

enum gl_shader_tone_map_variant {
	SHADER_TONE_MAP_NONE = 0,
	SHADER_TONE_MAP_HDR_TO_SDR,
	SHADER_TONE_MAP_SDR_TO_HDR,
	SHADER_TONE_MAP_HDR_TO_HDR,
};

struct gl_shader_requirements
{
	enum gl_shader_texture_variant variant;
	bool debug;
	bool csc_matrix;
	enum gl_shader_degamma_variant degamma;
	enum gl_shader_gamma_variant nl_variant;
	enum gl_shader_gamma_variant gamma;
	enum gl_shader_tone_map_variant tone_mapping;
};

struct gl_shader {
	struct gl_shader_requirements key;
	GLuint program;
	GLuint vertex_shader, fragment_shader;
	GLint proj_uniform;
	GLint tex_uniforms[3];
	GLint alpha_uniform;
	GLint color_uniform;
	GLint csc_uniform;
	GLint display_max_luminance;
	GLint content_max_luminance;
	GLint content_min_luminance;
	struct wl_list link; /* gl_renderer::shader_list */
};

struct gl_shader_generator;

void
gl_shader_requirements_init(struct gl_shader_requirements *requirements);

void
gl_shader_destroy(struct gl_shader *shader);

struct gl_shader *
gl_shader_create(struct gl_shader_generator *sg,
		 struct gl_shader_requirements *requirements);

struct gl_shader_generator *
gl_shader_generator_create(struct weston_compositor *compositor);

void
gl_shader_generator_destroy(struct gl_shader_generator *sg);

#endif
