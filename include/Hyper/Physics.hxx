#pragma once

#include <Math/Fuchsian.hxx>
#include <Math/AutD.hxx>

#include <Hyper/Geometry.hxx>

class Position {
private:
    Aut𝔻<Real>         _relativeXZ; // Relative to `_absoluteXZ`
    Real               _absoluteY;
    Fuchsian<Integer>  _absoluteXZ;
    Gaussian²<Integer> _absoluteOriginXZ;

public:
    Position() : _relativeXZ(Aut𝔻<Real>()), _absoluteY(0), _absoluteXZ(Tesselation::I)
    { _absoluteOriginXZ = _absoluteXZ.origin(); }

    Position(const Real y) : _relativeXZ(Aut𝔻<Real>()), _absoluteY(y), _absoluteXZ(Tesselation::I)
    { _absoluteOriginXZ = _absoluteXZ.origin(); }

    inline constexpr const auto & relativeXZ() const { return _relativeXZ; }

    inline const auto & absoluteY()        const { return _absoluteY;        }
    inline const auto & absoluteXZ()       const { return _absoluteXZ;       }
    inline const auto & absoluteOriginXZ() const { return _absoluteOriginXZ; }

    inline void setY(const Real & Y) { _absoluteY = Y; }

    inline void setRelativeXZ(const Aut𝔻<Real> & M)
    { _relativeXZ = M; _relativeXZ.normalize(); }

    inline void setAbsoluteXZ(const Fuchsian<Integer> & M)
    { _absoluteXZ = M; _absoluteXZ.normalize(); _absoluteOriginXZ = M.origin(); }

    inline void moveY(const Real & dy) { _absoluteY += dy; }

    /* It doesn’t do anything if the speed is big enough to jump over ≥2 chunks.
      (Of course, this can be easily fixed by iterating not only over neighbours,
       but it seems useless.) */
    bool moveXZ(const Gyrovector<Real> &);

    std::pair<int, int> round(const WorldTile *) const;
};

struct Camera {
    Real yaw = 0, pitch = 0, roll = 0;

    void rotate(const Real, const Real, const Real);

    glm::vec3 direction() const;
    glm::vec3 right() const;
};

class Entity {
private:
    WorldMap * _map; WorldTile * _tile; int _X, _Z; Position _position;

    bool _hasJumpedUp = false, _isAirborne = false;

    bool moveHorizontally(const Gyrovector<Real> & v, const Real dt);
    bool moveVertically(const Real dt);

public:
    Real velocityY = 0.0;

    Real eye = 0.0, height = 0.0, walkSpeed = 0.0, jumpSpeed = 0.0, gravity = 0.0;

    bool flymode = false, noclip = false;

    Entity(WorldMap * map) : _map(map), _tile(nullptr), _X(0), _Z(0) {}

    bool stuck();
    bool stuck(WorldTile *, int, Real, int);

    bool move(const Gyrovector<Real> & v, Real dt); // Returns true iff tile changes
    void teleport(const Position &);

    inline void elevate(const Real elevation) { _position.moveY(elevation); }
    inline void jump() { _hasJumpedUp = true; }

    inline constexpr const auto & map()      const { return _map;      }
    inline constexpr const auto & tile()     const { return _tile;     }
    inline constexpr const auto & position() const { return _position; }

    inline constexpr const auto & X() const { return _X; }
    inline constexpr const auto & Z() const { return _Z; }

    inline const bool isAirborne() const { return _isAirborne; }

    constexpr inline void jumpHeight(Real height)
    { jumpSpeed = sqrt(2 * gravity * height); }
};