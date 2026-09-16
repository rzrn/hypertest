#pragma once

#include <Math/Fuchsian.hxx>
#include <Math/AutD.hxx>

#include <Hyper/Geometry.hxx>

class Position {
private:
    Aut𝔻<Real>         _relativeXZ; // Relative to `_absoluteXZ`
    Fuchsian<Integer>  _absoluteXZ;
    Gaussian²<Integer> _absoluteOriginXZ;

public:
    Position() : _relativeXZ(Aut𝔻<Real>()), _absoluteXZ(Tesselation::I) { _absoluteOriginXZ = _absoluteXZ.origin(); }

    Position(const auto & P, const auto & G, const auto & g) : _relativeXZ(P), _absoluteXZ(G), _absoluteOriginXZ(g) {}
    Position(const auto & P, const auto & G) : _relativeXZ(P) { _relativeXZ.normalize(); set(G); }

    inline constexpr const auto & relativeXZ() const { return _relativeXZ; }

    inline const auto & absoluteXZ()       const { return _absoluteXZ;       }
    inline const auto & absoluteOriginXZ() const { return _absoluteOriginXZ; }

    inline void set(const Aut𝔻<Real> & M, const Fuchsian<Integer> & G)
    { _relativeXZ = M; _absoluteXZ = G; _absoluteOriginXZ = G.origin(); }

    inline void set(const Fuchsian<Integer> & G)
    { _absoluteXZ = G; _absoluteXZ.normalize(); _absoluteOriginXZ = G.origin(); }

    inline void set(const Aut𝔻<Real> & M) { _relativeXZ = M; }

    // It doesn’t do anything if the speed is big enough to jump over ≥2 chunks.
    // (Of course, this can be easily fixed by iterating
    //  not only over neighbours, but it seems useless.)
    std::pair<Position, bool> move(const Gyrovector<Real> &) const;

    std::pair<int, int> round(const WorldTile *) const;
};

struct Object {
    Position position;

    Real climb = 0, roc = 0;
    bool flying = false;

    Real yaw = 0, pitch = 0, roll = 0;

    void rotate(const Real, const Real, const Real);

    glm::vec3 direction() const;
    glm::vec3 right() const;
};

class Entity {
private:
    WorldMap * _map; WorldTile * _tile; int _X, _Z; Object _camera;

    bool jumped = false;

    bool moveHorizontally(const Gyrovector<Real> & v, const Real dt);
    bool moveVertically(const Real dt);

public:
    Real eye = 0.0, height = 0.0, walkSpeed = 0.0, jumpSpeed = 0.0, gravity = 0.0;

    bool flymode = false, noclip = false;

    Entity(WorldMap * map) : _map(map), _tile(nullptr), _X(0), _Z(0) {}

    bool stuck();
    bool stuck(WorldTile *, int, Real, int);

    // Returns true iff tile changes
    bool move(const Gyrovector<Real> & v, Real dt);
    void teleport(const Position &, const Real);

    constexpr void roc(const Real roc) { _camera.roc = roc; }
    constexpr void elevate(const Real elevation) { _camera.climb += elevation; }
    constexpr void jump() { jumped = true; }

    inline constexpr const auto & map()    const { return _map;    }
    inline constexpr const auto & tile()   const { return _tile;   }
    inline constexpr const auto & camera() const { return _camera; }

    inline constexpr const auto & X() const { return _X; }
    inline constexpr const auto & Z() const { return _Z; }

    constexpr inline void jumpHeight(Real height)
    { jumpSpeed = sqrt(2 * gravity * height); }

    inline void rotate(const Real Δyaw, const Real Δpitch, const Real Δroll)
    { _camera.rotate(Δyaw, Δpitch, Δroll); }
};