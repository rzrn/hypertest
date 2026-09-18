#include <Hyper/Physics.hxx>

bool Position::moveXZ(const Gyrovector<Real> & v) {
    auto P = _relativeXZ * Aut𝔻<Real>(v); auto w = P.origin();

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

void Object::rotate(const Real Δyaw, const Real Δpitch, const Real Δroll) {
    constexpr auto ε = 1e-6;

    yaw   = std::fmod(yaw + Δyaw, τ);
    pitch = std::clamp(pitch + Δpitch, -τ/4 + ε, τ/4 - ε);
    roll  = std::fmod(roll + Δroll, τ);
}

vec3 Object::direction() const {
    return vec3(
        cos(pitch) * sin(yaw),
        sin(pitch),
        cos(pitch) * cos(yaw)
    );
}

vec3 Object::right() const {
    return vec3(
        cos(roll) * sin(yaw - τ/4),
        sin(roll),
        cos(roll) * cos(yaw - τ/4)
    );
}

bool Entity::stuck(WorldTile * C, int X, Real y, int Z) {
    if (flymode && noclip) return false;

    if (C == nullptr || !C->ready()) return false;

    auto Y₁ = std::floor(y), Y₂ = std::floor(y + height);

    for (int Y = Y₁; Y <= Y₂; Y++)
        if (!C->walkable(X, Y, Z))
            return true;

    return false;
}

bool Entity::stuck() { return stuck(_tile, _X, _position.absoluteY(), _Z); }

bool Entity::moveHorizontally(const Gyrovector<Real> & v, const Real dt) {
    Position P(_position); auto isTileChanged = P.moveXZ(v.scale(dt));

    auto C = isTileChanged ? map()->poll(_position.absoluteXZ(), P.absoluteXZ()) : tile();

    if (C != nullptr) {
        if (!C->ready()) return false;

        auto [X, Z] = P.round(C);
        if (stuck(C, X, _position.absoluteY(), Z)) return false;
        _X = X; _Z = Z;
    }

    _tile = C; _position = P; return isTileChanged;
}

bool Entity::moveVertically(const Real dt) {
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

    auto γ⁻² = std::clamp<Real>(1 - Math::sqr(_camera.roc / vmax), 0, 1);
    auto roc = flymode ? _camera.roc : _camera.roc - dt * gravity * std::pow(γ⁻², 1.5);

    if (_hasJumpedUp) { roc += jumpSpeed; _hasJumpedUp = false; }

    auto y = WorldTile::clamp(_position.absoluteY() + dt * roc);

    if (stuck(_tile, _X, y, _Z)) { _camera.roc = 0; _isAirborne = false; }
    else { _position.setY(y); _camera.roc = roc; _isAirborne = true; }

    return false;
}

bool Entity::move(const Gyrovector<Real> & v, Real dt)
{ return moveHorizontally(v, dt) | moveVertically(dt); }

void Entity::teleport(const Position & P)
{ _position = P; _tile = _map->poll(P.absoluteXZ(), P.absoluteXZ()); }