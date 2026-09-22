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

#include <Hyper/Physics.hxx>

bool Position::moveXZ(const Gyrovector<Real> & dr) {
    auto P = _relativeXZ * Aut𝔻<Real>(dr); auto w = P.origin();

    if (WorldTile::isInsideOfDomain(w)) {
        setRelativeXZ(P);
        return false;
    }

    if (auto k = WorldTile::matchNeighbour(w)) {
        const auto & Δ   = Tesselation::neighbours[*k];
        const auto & Δ⁻¹ = Tesselation::neighbours⁻¹[*k];

        setRelativeXZ(Δ⁻¹ * P); setAbsoluteXZ(_absoluteXZ * Δ);
        return true;
    }

    return false;
}

std::pair<int, int> Position::round(const WorldTile * C) const {
    auto Q = (C->absoluteXZ().inverse() * _absoluteXZ).field<Real>() * Möbius<Real>(_relativeXZ);
    return WorldTile::round(Q.origin());
}

void Camera::rotate(const Real Δyaw, const Real Δpitch, const Real Δroll) {
    constexpr auto ε = 1e-6;

    yaw   = std::fmod(yaw + Δyaw, τ);
    pitch = std::clamp(pitch + Δpitch, -τ/4 + ε, τ/4 - ε);
    roll  = std::fmod(roll + Δroll, τ);
}

vec3 Camera::direction() const {
    return vec3(
        cos(pitch) * sin(yaw),
        sin(pitch),
        cos(pitch) * cos(yaw)
    );
}

vec3 Camera::right() const {
    return vec3(
        cos(roll) * sin(yaw - τ/4),
        sin(roll),
        cos(roll) * cos(yaw - τ/4)
    );
}

bool Entity::canWalkAt(WorldTile * C, int X, Real y, int Z) {
    if (isFlyModeEnabled && isNoclipEnabled) return true;

    if (C == nullptr || !C->ready()) return true;

    auto Y₁ = std::floor(y), Y₂ = std::floor(y + height);

    for (int Y = Y₁; Y <= Y₂; Y++)
        if (!C->walkable(X, Y, Z))
            return false;

    return true;
}

bool Entity::canWalkAt() { return canWalkAt(_tile, _X, _position.absoluteY(), _Z); }

inline Gyrovector<Real> exp₀(const Real x, const Real y) {
    /* exp₀(v) converts a tangent vector v ∈ T₀(𝔻) into a gyrovector originating at z = 0.
       [1] Hyperbolic Neural Networks, Octavian Ganea, Gary Becigneul, and Thomas Hofman,
           Advances in Neural Information Processing System, 2018
           https://arxiv.org/pdf/1805.09112
       [2] Hyperbolic Neural Networks++, Ryohei Shimizu, Yusuke Mukuta, Tatsuya Harada,
           Published as a conference paper at ICLR 2021 */

    auto n = std::hypot(x, y);

    if (n > 0) {
        auto k = tanh(n) / n;
        return {k * x, k * y};
    } else {
        return 0;
    }
}

bool Entity::moveXZ(const Real dt) {
    auto dr = exp₀(velocity.x * dt, velocity.z * dt);
    Position P(_position); auto isTileChanged = P.moveXZ(dr);

    auto C = isTileChanged ? map()->poll(_position.absoluteXZ(), P.absoluteXZ()) : tile();

    if (C != nullptr) {
        if (!C->ready()) return false;

        auto [X, Z] = P.round(C);
        if (!canWalkAt(C, X, _position.absoluteY(), Z)) return false;
        _X = X; _Z = Z;
    }

    _tile = C; _position = P; return isTileChanged;
}

bool Entity::moveY(const Real dt) {
    if (!_tile->ready()) return false;

    /*
        Lorentz factor: γ(v) = 1/√(1 − v²/c²).
        Relativistic kinetic energy: T = γ(v)mc².
        Potential energy for the (newtonian) uniform gravitational field: U = mgh.

        Then:
          dT = γ′(v)mc²dv,
          δU = mgdh = mgvdt,
          δE = dT + δU,
        where
          γ′(v) = (v/c²) × (1 − v²/c²)^(−3/2).

        Assume that energy is locally conserved:
          δE = 0
        ↔ dT = −δU
        ↔ γ′(v)mc²dv = −mgvdt
        ↔ γ′(v)dv = −(v/c²)gdt
        ↔ dv = −g(1 − v²/c²)^(3/2) × dt.

        In particular, (1 − v²/c²)^(3/2) = 1 − 3v²/2c² + o(v⁴/c⁴),
        so if v/c ≈ 0, then dv/dt ≈ −g.

        For the non-trivial topology δE may be not exact, so that it
        makes no sense to define potential energy globally.
        Therefore energy will not be conserved.

        In particular, “dh” is known to be a generator
        of de Rham cohomology group H¹(S¹) ≅ ℝ.
    */
    constexpr Real vmax = 32.0;

    auto γ⁻² = std::clamp<Real>(1 - Math::sqr(velocity.y / vmax), 0, 1);
    auto vy = isFlyModeEnabled ? velocity.y : velocity.y - dt * gravity * std::pow(γ⁻², 1.5);

    // For simplicity, we use an approximation in the v/c ≈ 0 limit here.
    if (_hasJumpedUp) { vy += sqrt(2 * gravity * jumpHeight); _hasJumpedUp = false; }

    auto y = WorldTile::clamp(_position.absoluteY() + dt * vy);

    if (canWalkAt(_tile, _X, y, _Z)) {
        _position.setY(y); velocity.y = vy; _isAirborne = true;
    } else {
        velocity.y = 0; _isAirborne = false;
    }

    return false;
}

bool Entity::moveXYZ(const Real dt)
{ return moveXZ(dt) | moveY(dt); }

void Entity::setXYZ(const Position & P)
{ _position = P; _tile = _map->poll(P.absoluteXZ(), P.absoluteXZ()); }