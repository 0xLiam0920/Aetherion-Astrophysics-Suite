#pragma once
// ============================================================
// bh3d_shaderinject.hpp: GL-free helper that splices the shared GLSL file
// (bh3d_physics.glsl) into a fragment shader. Contains zero OpenGL-specific
// code so the regression tests can exercise it directly.
// ============================================================

#include <string>

namespace bh3d {
// Returns `src` with the quality profile and the shared `common` source inserted
// after the version declaration line (GLSL requires #version to be first
// though comments may precede it). A `#line` restores the original
// original numbering, so driver error messages still tell us what went wrong where
// , at least in .frag file. Errors that originate inside our shared physics block are compared against the lines that follow the #version line.
//
// Returns an empty str if `src` has no #version line or `common` is empty,
// so callers can detect and handle failures properly.
inline std::string injectCommon(const std::string& src,
                                const std::string& common,
                                int quality) {
    if (common.empty()) return {};

    // Find a line that begins with #version (ignoring leading blanks).
    size_t pos = 0;
    size_t versionLine = 1;            // 1-based line number of the #version line
    size_t lineEnd = std::string::npos;
    while (pos < src.size()) {
        size_t eol = src.find('\n', pos);
        size_t end = (eol == std::string::npos) ? src.size() : eol;
        size_t first = src.find_first_not_of(" \t", pos);
        if (first != std::string::npos && first < end &&
            src.compare(first, 8, "#version") == 0) {
            lineEnd = end;
            break;
        }
        if (eol == std::string::npos) break;
        pos = eol + 1;
        ++versionLine;
    }
    if (lineEnd == std::string::npos) return {};

    const size_t headEnd = (lineEnd < src.size()) ? lineEnd + 1 : src.size();
    std::string out = src.substr(0, lineEnd);
    out += '\n';
    out += "#define QUALITY " + std::to_string(quality) + "\n";
    out += common;
    if (common.back() != '\n') out += '\n';
    // GLSL 3.30: the line AFTER a #line directive takes the given number.
    out += "#line " + std::to_string(versionLine + 1) + "\n";
    out += src.substr(headEnd);
    return out;
}

}  // namespace bh3d
