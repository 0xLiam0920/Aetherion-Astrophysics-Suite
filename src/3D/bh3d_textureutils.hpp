#pragma once

#include "platform.hpp"
#include "gl_types.hpp"

#include <string>

// Create a 1x1 white fallback texture
GLTexture2D createFallbackWhiteTexture();

// Load a 2D texture from an image file via SFML; returns empty on failure
GLTexture2D loadTexture2D(const std::string& path);

// Generate a procedural accretion-disk texture
GLTexture2D createDiskTextureProcedural(unsigned int size = 256);

// This call makes a 256×1 sRGB blackbody colour LUT from 1,000–40,000 K for the accretion disk
// Uses the Planck spectrum against CIE 1931 CMFs → XYZ → sRGB.
// for physical disc color only
GLTexture2D createBlackbodyLUT(int nSamples = 256);

// Linear-light, absolute-radiance blackbody LUT (GL_RGB32F, log-spaced
// 300 K to 1e8 K). Sample it with blackbodyRadiance() in bh3d_physics.glsl.
GLTexture2D createBlackbodyRadianceLUT(int nSamples = 1024);
