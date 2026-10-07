#pragma once

/*
    Copyright © 2026 rzrn

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include <utility>

#include <sys/types.h>
#include <cstdint>
#include <cstdio>

#include <glm/vec4.hpp>

#include <Math/Gaussian.hxx>

#include <Hyper/Fundamentals.hxx>

struct RGB3f {
    glm::vec4 value;

    RGB3f() : value(0.0f, 0.0f, 0.0f, 1.0f) {}
    RGB3f(float s) : value(s, s, s, 1.0f) {}
    RGB3f(float r, float g, float b) : value(r, g, b, 1.0f) {}

    inline float & r() { return value[0]; }
    inline float & g() { return value[1]; }
    inline float & b() { return value[2]; }

    inline explicit operator glm::vec4() const { return value; }
};

using MapColorScheme = RGB3f[256];

class HVXL {
public:
    static constexpr int sizeTileX = 16;
    static constexpr int sizeTileY = 96;
    static constexpr int sizeTileZ = 16;

    static constexpr int maxTileX = sizeTileX - 1;
    static constexpr int maxTileY = sizeTileY - 1;
    static constexpr int maxTileZ = sizeTileZ - 1;

    using node_t = uint8_t;

    static constexpr size_t sizeNode = sizeof(node_t);
    static constexpr size_t sizeTile = sizeTileX * sizeTileY * sizeTileZ * sizeNode;

    static constexpr inline size_t offsetFromXYZ(int X, int Y, int Z) {
        size_t offset;

        offset = Y;
        offset = Z + offset * sizeTileZ;
        offset = X + offset * sizeTileX;

        return offset * sizeof(node_t);
    }

    static inline node_t readNode(const uint8_t * src) { return src[0]; }
    static inline void writeNode(uint8_t * dest, const node_t value) { dest[0] = value; }

private:
    static_assert(sizeof(off_t) <= sizeof(uint64_t));

    static constexpr uint64_t fmtver = 1;
    FILE * fp; off_t offsetHeader;

    static constexpr size_t sizeHashInBits = 20;
    static constexpr size_t sizeHashRecord = (1 << sizeHashInBits) * sizeof(uint64_t);

    template<typename T> T readAt(off_t offset);
    template<typename T> size_t writeAt(off_t offset, const T value);

    size_t writeHashRecord();

    void readTileData(uint8_t *, RGB3f &);
    size_t writeTileData(const uint8_t *, RGB3f);

public:
    char mapname[256]; uint32_t seed; MapColorScheme color;

    HVXL();
    ~HVXL();

    void open(const std::string &);
    void close();

    inline off_t offsetEOF() const { fseeko(fp, 0, SEEK_END); return ftello(fp); }

    void readMapHeader();
    size_t writeMapHeader();

    inline void readTileDataAt(off_t offset, uint8_t * buf, RGB3f & fog) {
        if (fseeko(fp, offset, SEEK_SET) < 0)
            return;

        readTileData(buf, fog);
    }

    inline size_t writeTileDataAt(off_t offset, const uint8_t * buf, RGB3f fog) {
        if (fseeko(fp, offset, SEEK_SET) < 0)
            return 0;

        return writeTileData(buf, fog);
    }

    bool findTileOffset(const Gaussian²<Integer> &, off_t &, uint64_t &);
};
