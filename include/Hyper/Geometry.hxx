#pragma once

/*
    Copyright © 2022–2024, 2026 rzrn

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

#include <optional>
#include <string>
#include <vector>

#include <atomic>
#include <mutex>

#include <GL/glew.h>

#include <Hyper/Fundamentals.hxx>
#include <Hyper/Shader.hxx>
#include <Hyper/HVXL.hxx>

#include <Math/Fuchsian.hxx>
#include <Math/AutD.hxx>

#include <Meta/DoubleBuffer.hxx>
#include <Meta/List.hxx>

using ℤi = Gaussian<Integer>;

namespace Tesselation {
    enum class Direction { Identity, Up, Down, Left, Right };
    using enum Direction;

    template<typename T> T interpret(Direction);

    template<Direction... ds> struct Compose {
        template<typename T> static constexpr inline T eval()
        { return (interpret<T>(ds) * ...); }
    };

    template<typename, typename> struct Eval;

    template<typename T> struct Eval<T, List<>>
    { static constexpr inline void insert(size_t, auto &) {} };

    template<typename T, typename U, typename... Us> struct Eval<T, List<U, Us...>> {
        static constexpr inline void insert(size_t i, auto & ret) {
            ret[i] = U::template eval<T>(); ret[i].normalize();
            Eval<T, List<Us...>>::insert(i + 1, ret);
        }
    };

    template<typename T, typename Us> constexpr inline std::array<T, Length<Us>> eval()
    { std::array<T, Length<Us>> retval; Eval<T, Us>::insert(0, retval); return retval; }

    // `... + 1` is added because there are N + 1 vertices for N intervals
    using Grid = Array²<Gyrovector<Real>, HVXL::sizeTileX + 1, HVXL::sizeTileZ + 1>;

    using Neighbours = List<
        Compose<Up>, Compose<Left>, Compose<Down>, Compose<Right>,
        Compose<Up,   Left>,  Compose<Up,   Left,  Down>, Compose<Up,   Left,  Down, Right>,
        Compose<Up,   Right>, Compose<Up,   Right, Down>, Compose<Up,   Right, Down, Left>,
        Compose<Down, Left>,  Compose<Down, Left,  Up>,   Compose<Down, Left,  Up,   Right>,
        Compose<Down, Right>, Compose<Down, Right, Up>,   Compose<Down, Right, Up,   Left>
    >;

    constexpr size_t amount = Length<Neighbours>;

    template<typename T> using Array = std::array<T, amount>;

    extern const Fuchsian<Integer>        I, U, L, D, R;
    extern const Array<Fuchsian<Integer>> neighbours;
    extern const Array<Aut𝔻<Real>>        neighbours⁻¹;
    extern const Grid                     corners;
    extern const Real                     meter;
}

template<typename T> struct Parallelogram {
    Gyrovector<T> A, B, C, D;

    Parallelogram() {}
    Parallelogram(auto A, auto B, auto C, auto D) : A(A), B(B), C(C), D(D) {}

    const auto rev() const { return Parallelogram<T>(D, C, B, A); }
};

class WorldTile; using WorldMapgen = void(WorldTile *);

enum : unsigned int {
    NONE = 0u,
    XMIN = 1u << 0,
    XMAX = 1u << 1,
    YALL = 1u << 2,
    ZMIN = 1u << 3,
    ZMAX = 1u << 4,
    XALL = XMIN | XMAX,
    ZALL = ZMIN | ZMAX,
};

class WorldTile {
public:
    static constexpr auto sizeTileX = HVXL::sizeTileX;
    static constexpr auto sizeTileY = HVXL::sizeTileY;
    static constexpr auto sizeTileZ = HVXL::sizeTileZ;

    static constexpr auto maxTileX = HVXL::maxTileX;
    static constexpr auto maxTileY = HVXL::maxTileY;
    static constexpr auto maxTileZ = HVXL::maxTileZ;

    static constexpr auto sizeTile = HVXL::sizeTile;

private:
    // These are for drawing, where `_cameraVerticalDistance` is used to track tiles that are to be unloaded
    Möbius<Real> _cameraXZ; Real _cameraHorizontalDistance;
    // These are for indexing, where `_absoluteOriginXZ` is equal to `_absoluteXZ.origin()`
    Fuchsian<Integer> _absoluteXZ; Gaussian²<Integer> _absoluteOriginXZ;

    FaceShader::VAO faces; EdgeShader::VAO edges;

    DoubleBuffer<GLsizei, sizeTileY> facesOffsetY, edgesLowerOffsetY, edgesUpperOffsetY;

    inline GLsizei facesLowerOffsetY(const int Y) const { return Y > 0 ? facesOffsetY(Y - 1) : 0; }
    inline GLsizei facesUpperOffsetY(const int Y) const { return facesOffsetY(Y);                 }

    std::atomic<bool> isToUpdateVAO = false, isToUploadVAO = false, _isModified = false;

    uint8_t _vxl[HVXL::sizeTile];

public:
    uint64_t hash = 0; off_t fileOffset = 0; RGB3f fog;

    bool isToBeErased = false;

    WorldTile(const Fuchsian<Integer> &);

    ~WorldTile();

    void emitFaces(const MapColorScheme &);
    void emitEdges(const MapColorScheme &);

    void renderFaces(FaceShader *, int, int);
    void renderEdges(EdgeShader *, int, int);

    void updateCameraXZ(const Fuchsian<Integer> &);

    void updateVAO(const MapColorScheme &);
    bool uploadVAO(const Fuchsian<Integer> &);

    bool walkable(int, Real, int);

    void load(HVXL &);
    void save(HVXL &);

    inline void threadsafeSave() { _isModified = true; }
    inline void threadsafeUpdateVAO() { isToUpdateVAO = true; }

    inline const bool isModified() const { return _isModified; }

    inline const auto cameraXZ()                 const { return _cameraXZ;                 }
    inline const auto cameraHorizontalDistance() const { return _cameraHorizontalDistance; }

    inline const auto absoluteXZ()       const { return _absoluteXZ;       }
    inline const auto absoluteOriginXZ() const { return _absoluteOriginXZ; }

    static inline int mod(int a, int b) { int r = a % b; return r < 0 ? r + b : r; }

    // TODO: make this to return a reference instead
    template<unsigned int mask = NONE> inline uint8_t get(int X, int Y, int Z) const {
        using namespace Fundamentals;

        // Bounds checks are eliminated at compile time unless needed

        if constexpr(mask & XMIN) {
            if (X < 0) return {};
        }

        if constexpr(mask & XMAX) {
            if (maxTileX < X) return {};
        }

        if constexpr(mask & ZMIN) {
            if (Z < 0) return {};
        }

        if constexpr(mask & ZMAX) {
            if (maxTileZ < Z) return {};
        }

        if constexpr(mask & YALL) {
            Y = mod(Y, sizeTileY);
        }

        auto src = &_vxl[HVXL::offsetFromXYZ(X, Y, Z)];

        return HVXL::readNode(src);
    }

    inline void set(int X, int Y, int Z, const HVXL::node_t value)
    { _isModified = true; HVXL::writeNode(&_vxl[HVXL::offsetFromXYZ(X, Y, Z)], value); }

    inline uint8_t * vxl() { _isModified = true; return _vxl; }
    inline const uint8_t * vxl() const { return _vxl; }

    static bool touch(const Gyrovector<Real> &, int, int);
    static std::pair<int, int> round(const Gyrovector<Real> &);

    static bool isInsideOfDomain(const Gyrovector<Real> &);
    static std::optional<size_t> matchNeighbour(const Gyrovector<Real> &);

    static inline Real clamp(Real x) { return Math::remainder<Real>(x, sizeTileY); }
};

class WorldMap {
private:
    std::vector<WorldTile *> readQueue, writeQueue;

    WorldTile * findNonThreadsafe(const Gaussian²<Integer> &);
    WorldTile * queueNewTile(const Fuchsian<Integer> &);

public:
    std::mutex readMutex;

    HVXL file;

    WorldMapgen * mapgen = nullptr;

    WorldMap();
    ~WorldMap();

    void save();
    void update(const Real &);
    void updateCameraXZ(const Fuchsian<Integer> &);

    WorldTile * find(const Gaussian²<Integer> &);

    inline auto begin() const { return readQueue.begin(); }
    inline auto end()   const { return readQueue.end();   }
};

template<typename T> struct Bitfield {
    T value;

    Bitfield(T value) : value(value) {}
    inline operator T() { return value; }

    inline void set(size_t n, bool bit) { value = (value | (1 << n)) & ~(T(!bit) << n); }
    inline bool get(size_t n) const { return (value >> n) & 1; }
};
