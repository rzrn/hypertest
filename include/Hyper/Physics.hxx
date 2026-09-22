#pragma once

/*
    Copyright © 2023–2024, 2026 rzrn

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

#include <glm/vec3.hpp>

#include <Math/Fuchsian.hxx>
#include <Math/AutD.hxx>

#include <Hyper/Geometry.hxx>

class WorldXYZ {
private:
    Aut𝔻<Real>         _relativeXZ; // Relative to `_absoluteXZ`
    Real               _absoluteY;
    Fuchsian<Integer>  _absoluteXZ;
    Gaussian²<Integer> _absoluteOriginXZ;

public:
    WorldXYZ() : _relativeXZ(Aut𝔻<Real>()), _absoluteY(0), _absoluteXZ(Tesselation::I)
    { _absoluteOriginXZ = _absoluteXZ.origin(); }

    WorldXYZ(const Real y) : _relativeXZ(Aut𝔻<Real>()), _absoluteY(y), _absoluteXZ(Tesselation::I)
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

    inline void moveY(const Real dy) { _absoluteY += dy; }

    /* It doesn’t do anything if the speed is big enough to jump over ≥2 chunks.
      (Of course, this can be easily fixed by iterating not only over neighbours,
       but it seems useless.) */
    bool moveXZ(const Gyrovector<Real> &);

    std::pair<int, int> round(const WorldTile *) const;
};

struct Camera {
    Real eye = 0, yaw = 0, pitch = 0, roll = 0;

    void rotate(const Real, const Real, const Real);

    glm::vec3 direction() const;
    glm::vec3 right() const;
};

class Entity {
private:
    WorldMap * _map; WorldTile * _tile; int _X, _Z; WorldXYZ _r;

    bool _hasJumpedUp = false, _isAirborne = false;

public:
    vec3 v = {0, 0, 0}; // Velocity is a tangent vector, not a gyrovector, thus `vec3`

    Real height = 0.0, walkingSpeed = 0.0, jumpHeight = 0.0, gravity = 0.0;

    bool isFlyModeEnabled = false, isNoclipEnabled = false;

    Entity(WorldMap * map) : _map(map), _tile(nullptr), _X(0), _Z(0) {}

    bool canWalkAt();
    bool canWalkAt(WorldTile *, int, Real, int);

    // These return true iff the current tile changes
    bool moveY(const Real dt);
    bool moveXZ(const Real dt);
    bool moveXYZ(const Real dt);

    void setXYZ(const WorldXYZ &);

    inline void applyJumpImpulse() { _hasJumpedUp = true; }

    inline constexpr const auto & r() const { return _r; }
    inline constexpr const auto & X() const { return _X; }
    inline constexpr const auto & Z() const { return _Z; }

    inline constexpr const auto & map()  const { return _map;  }
    inline constexpr const auto & tile() const { return _tile; }

    inline const bool isAirborne() const { return _isAirborne; }
};