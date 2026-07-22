#include "bh3d_textureutils.hpp"
#include "bh3d_blackbody.hpp"

// Textures fall back to a 1x1 white pixel when an asset can't be loaded.

#include <SFML/Graphics/Image.hpp>
#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <algorithm>
#include <vector>

// 1x1 white fallback used when a real texture is missing.
GLTexture2D createFallbackWhiteTexture() {
    GLuint tex = 0;
    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);
    unsigned char white[3] = {255, 255, 255};
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 1, 1, 0, GL_RGB, GL_UNSIGNED_BYTE, white);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    GLTexture2D result;
    result.adopt(tex);
    return result;
}

GLTexture2D loadTexture2D(const std::string& path) {
    sf::Image img;
    if (!img.loadFromFile(path))
        return {};

    const auto size = img.getSize();
    if (size.x == 0 || size.y == 0)
        return {};

    // SFML stores pixels in RGBA8
    GLuint tex = 0;
    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);
    glTexImage2D(GL_TEXTURE_2D,
                 0,
                 GL_RGBA,
                 static_cast<GLsizei>(size.x),
                 static_cast<GLsizei>(size.y),
                 0,
                 GL_RGBA,
                 GL_UNSIGNED_BYTE,
                 img.getPixelsPtr());
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR); // no mipmaps
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    GLTexture2D result;
    result.adopt(tex);
    return result;
}

GLTexture2D createDiskTextureProcedural(unsigned int size) { 
    sf::Image img(sf::Vector2u{size, size}, sf::Color(0, 0, 0, 0));

    const int center = static_cast<int>(size / 2);
    const int rIn = static_cast<int>(static_cast<float>(size) * 0.22f);
    const int rOut = static_cast<int>(static_cast<float>(size) * 0.48f);

    for (int r = rIn; r < rOut; ++r) {
        const float t = (r - rIn) / std::max(1.0f, static_cast<float>(rOut - rIn));
        const float a = 1.0f - t;
        const int alpha = static_cast<int>(255.0f * a * a);

        const int rcol = 255;
        const int gcol = static_cast<int>(220.0f * (1.0f - t) + 120.0f * t);
        const int bcol = static_cast<int>(255.0f * (1.0f - t) + 40.0f * t);

        for (int deg = 0; deg < 360; ++deg) {
            const float rad = static_cast<float>(deg) * 3.14159265359f / 180.0f;                // Explanation for this block: we want to draw a circle of radius r, 
            const int x = static_cast<int>(center + r * std::cos(rad));                         // but we only have pixels. so we iterate over angles and compute the corresponding (x,y) for each angle. 
            const int y = static_cast<int>(center + r * std::sin(rad));                         // this creates a rough circle outline. the alpha and color are determined by how far r is between rIn 
            if (x >= 0 && x < static_cast<int>(size) && y >= 0 && y < static_cast<int>(size)) { // and rOut, creating a radial gradient effect. it's not perfect, but it gives us a nice procedural 
                img.setPixel({static_cast<unsigned int>(x), static_cast<unsigned int>(y)},      // disk texture without needing to load an image file. 
                             sf::Color(static_cast<std::uint8_t>(rcol),
                                       static_cast<std::uint8_t>(gcol),
                                       static_cast<std::uint8_t>(bcol),
                                       static_cast<std::uint8_t>(alpha)));
            }
        }
    }

    GLuint tex = 0;
    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);

    const auto sz = img.getSize();
    glTexImage2D(GL_TEXTURE_2D,         // SFML won't give us a GPU-resident buffer or a reusable raw pointer,
                 0,                     // so the pixels get copied once here. Redundant, but unavoidable with SFML.
                 GL_RGBA,
                 static_cast<GLsizei>(sz.x),
                 static_cast<GLsizei>(sz.y),
                 0,
                 GL_RGBA,
                 GL_UNSIGNED_BYTE,
                 img.getPixelsPtr());

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    GLTexture2D result;
    result.adopt(tex);
    return result;
}

// ─────────────────────────────────────────────────────────────────────────────
// Physical blackbody colour LUT (1,000 to 40,000 K)
//
// This integrates B(lambda,T) times the CIE 1931 CMF over 380 to 780 nm,
// converts that to XYZ then D65 sRGB. Each entry gets normalised so
// max(R,G,B) equals 1, since we only care about the chromaticity here.
// Actual brightness is handled later by the shader's NT flux and beaming math.
// ─────────────────────────────────────────────────────────────────────────────
GLTexture2D createBlackbodyLUT(int nSamples) {
    // The colour physics itself lives in bh3d_blackbody.hpp, since it's GL-free
    // and shared with the CPU test suite. That way the exact values we upload
    // here are the same ones the physics regression tests check against.
    auto enc = [](double c) -> std::uint8_t {
        return static_cast<std::uint8_t>(std::clamp(c * 255.0 + 0.5, 0.0, 255.0));
    };

    std::vector<std::uint8_t> px(static_cast<std::size_t>(nSamples) * 3);
    for (int i = 0; i < nSamples; ++i) {
        const double T = 1000.0 + 39000.0 * i / std::max(nSamples - 1, 1);
        const std::array<double, 3> rgb = bh3d::physics::blackbodyRGB(T);
        px[static_cast<std::size_t>(i) * 3 + 0] = enc(rgb[0]);
        px[static_cast<std::size_t>(i) * 3 + 1] = enc(rgb[1]);
        px[static_cast<std::size_t>(i) * 3 + 2] = enc(rgb[2]);
    }
    GLuint tex=0;
    glGenTextures(1,&tex);
    glBindTexture(GL_TEXTURE_2D,tex);
    glTexImage2D(GL_TEXTURE_2D,0,GL_RGB,nSamples,1,0,GL_RGB,GL_UNSIGNED_BYTE,px.data());
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_S,GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_T,GL_CLAMP_TO_EDGE);
    GLTexture2D lut; lut.adopt(tex); return lut;
}
