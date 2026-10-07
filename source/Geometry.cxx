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

#include <Hyper/Geometry.hxx>
#include <Hyper/Game.hxx>

namespace Tesselation {
    using namespace Fundamentals;

    // Tile’s neighbours in tesselation
    const Fuchsian<Integer> I { ℤi(+1, +0), ℤi(+0, +0), ℤi(+0, +0), ℤi(+1, +0) };
    const Fuchsian<Integer> U { ℤi(+6, +0), ℤi(+6, +6), ℤi(+1, -1), ℤi(+6, +0) };
    const Fuchsian<Integer> L { ℤi(+6, +0), ℤi(+6, -6), ℤi(+1, +1), ℤi(+6, +0) };
    const Fuchsian<Integer> D { ℤi(+6, +0), ℤi(-6, -6), ℤi(-1, +1), ℤi(+6, +0) };
    const Fuchsian<Integer> R { ℤi(+6, +0), ℤi(-6, +6), ℤi(-1, -1), ℤi(+6, +0) };

    /*
        𝔻  = { z ∈ ℂ | |z| ≤ 1 }
        𝔻ₛ = { z ∈ ℂ | |z| ≤ s }

        (In particular, 𝔻₁ = 𝔻.)

        Möbius transformation of translation towards vector b ∈ 𝔻 in Poincaré disk model is given by the formula:
            Φ = [1, b; b*, 1], so φ(z) = (z + b) / (zb* + 1).
        (https://en.wikipedia.org/wiki/M%C3%B6bius_transformation#Subgroups_of_the_M%C3%B6bius_group)

        Knowing that D½ = √(2 − √3) (see `include/Hyper/Fundamentals.hxx`),
        we have direction vectors: a = D½ and b = iD½.

        Result of their coaddition is a required translation vector:
            Coadd(a, b) = ((1 − |a|²)a + (1 − |b|²)b) / (1 − |a|²|b|²)
                        = ((1 − D½²)D½ + (1 − D½²)iD½) / (1 − D½⁴)
                        = D½(1 − D½²)/(1 − D½⁴) × (1 + i)
                        = (1 + i)/√6

        So corresponding Möbius transformation is given by the formula:
            Φ = [1, (1 + i)/√6; (1 − i)/√6, 1].

        Now let z ∈ 𝔻, s > 0, φ(z) = (az + b) / (cz + d).
        Then sφ(z/s) = s(az/s + b) / (cz/s + d) = (az + bs) / ((c/s)z + d), Φₛ = [a, bs; c/s, d].
        We see that Φₛ maps 𝔻ₛ to 𝔻ₛ, so this operation is exactly a change of curvature.

        We choose s = √6, then:
            Φₛ = [1, 1 + i; (1 − i)/6, 1]
        Since (az + b) / (cz + d) = (kaz + kb) / (kcz + kd), we may take:
            Φₛ′ = 6 × Φₛ = [6, 6 + 6i; 1 − i, 6]
        This is exactly U matrix.

        Choosing other signs in a = ±D½ and b = ±iD½, we will obtain L, D and R.
    */

    template<> Fuchsian<Integer> interpret(Direction d) {
        switch (d) {
            case Up:    return U;
            case Down:  return D;
            case Left:  return L;
            case Right: return R;
            default:    return I;
        }
    }

    template<> Aut𝔻<Real> interpret(Direction d)
    { return Aut𝔻<Real>(interpret<Fuchsian<Integer>>(d).field<Real>().origin()); }

    template<std::size_t N, typename T>
    auto inverse(const std::array<T, N> & xs) {
        std::array<T, N> retval;

        for (size_t i = 0; i < N; i++)
            retval[i] = xs[i].inverse();

        return retval;
    }

    const Array<Fuchsian<Integer>> neighbours = eval<Fuchsian<Integer>, Neighbours>();
    const Array<Aut𝔻<Real>> neighbours⁻¹ = inverse(eval<Aut𝔻<Real>, Neighbours>());

    // Generation of tile’s grid

    constexpr Real d = D½ / sqrt2;

    constexpr auto d₁₂ = Model(Klein).apply(d, d);
    constexpr auto hd₁ = Math::atanh(d₁₂.first), hd₂ = Math::atanh(d₁₂.second);

    constexpr auto Ψ(Real t₁, Real t₂) {
        auto k₁ = Math::tanh(t₁ * hd₁);
        auto k₂ = Math::tanh(t₂ * hd₂);

        auto [x, y] = Model(Klein).unapply(k₁, k₂);
        auto u = (x + y) / sqrt2, v = (x - y) / sqrt2;

        return Gyrovector<Real>(u, v);
    }

    auto Ψ⁻¹(Real u, Real v) {
        auto x = (u + v) / sqrt2, y = (u - v) / sqrt2;
        auto [k₁, k₂] = Model(Klein).apply(x, y);

        return std::pair(std::atanh(k₁) / hd₁, std::atanh(k₂) / hd₂);
    }

    constexpr auto apply(int X, int Z) {
        auto t₁ = 2 * Real(X) / WorldTile::sizeTileX - 1;
        auto t₂ = 2 * Real(Z) / WorldTile::sizeTileZ - 1;
        return Ψ(t₁, t₂);
    }

    auto unapply(Real u, Real v) {
        auto [t₁, t₂] = Ψ⁻¹(u, v);

        /* Tile’s border is not exactly a hyperbolic line (i.e. circular arc on the Poincaré disk),
           but its piecewise linear approximation; so there are parts of the outer blocks that extend
           slightly beyond the boundary of the ideal hyperbolic square.
           That’s why we need to “std::clamp” here.
        */
        t₁ = std::clamp<Real>(t₁, -1.0, 0.9999); // t₁ ≤ 0.9999 < 1 so that int(X) < sizeTileX
        t₂ = std::clamp<Real>(t₂, -1.0, 0.9999);

        auto X = (t₁ + 1) / 2 * WorldTile::sizeTileX;
        auto Z = (t₂ + 1) / 2 * WorldTile::sizeTileZ;

        return std::pair(int(X), int(Z));
    }

    constexpr auto init() {
        using namespace Fundamentals;

        Array²<Gyrovector<Real>, WorldTile::sizeTileX + 1, WorldTile::sizeTileZ + 1> retval;

        for (int X = 0; X <= WorldTile::sizeTileX; X++)
            for (int Z = 0; Z <= WorldTile::sizeTileZ; Z++)
                retval[X][Z] = apply(X, Z);

        return retval;
    }

    constexpr Grid corners = init();

    constexpr auto distance(int X₁, int Z₁, int X₂, int Z₂)
    { return (-corners[X₁][Z₁] + corners[X₂][Z₂]).abs(); }

    constexpr Real meter = distance(
        WorldTile::sizeTileX / 2, WorldTile::sizeTileZ / 2,
        WorldTile::sizeTileX / 2, WorldTile::sizeTileZ / 2 + 1
    );
}

WorldTile::WorldTile(const Fuchsian<Integer> & absoluteXZ) : _cameraHorizontalDistance(-1), _absoluteXZ(absoluteXZ), _vxl{0}, fog(1.0f) {
    /*
        Unfortunately, precomposition of `absoluteXZ` with (z ↦ z × exp(iπk/2)) for k ∈ ℤ
        will yield matrix able to render this tile in the same place but rotated about its own center by πk/2 radians.

        So if we don’t resolve this ambiguity, the same tile may render with different rotation
        according only to player’s path, and this will definitely crush the landscape
        (Imagine random rotations of chunks in Minecraft in a mountainous biome.)

        Since exp(iπk/2) ∈ {±1, ±i}, such rotation is equivalent to multiplying `a` and `c` by ±1/±i at the same time.
        `absoluteXZ` is pre-divided by hcf(a, b, c, d), so there is always ability to multiply *all* terms by ±1/±i
        yielding the same transformation.
        (Two Möbius transformations are equal iff their matrices differ by multiplicative constant.)

        Hence we have 4 × 4 = 16 options. It’s important that we can always select multipliers
        so that both components of `a` and `b` will be non-negative, that’s what we’re doing.

        We proceed as follows: for each z ∈ ℂˣ there is exactly one u ∈ {±1, ±i} such that
        Re(uz) > 0 and Im(uz) ≥ 0, and it is given by the `.normalize(...)` method of `Gaussian<T>`.
        For convenience, we write u = ε(z) for such a number. Then, let ε(z₁, z₂) = ε(z₁) if z₁ ≠ 0
        and ε(z₁, z₂) = ε(z₂) if z₂ ≠ 0, if z₁ and z₂ are not simultaneously zero.

        We assume that det(_absoluteXZ) = ad − bc ≠ 0, because:
        1) det(I), det(U), det(L), det(D), det(R) ≠ 0 (see above),
           so determinant from any of their product is also non-zero.
           (Since det(AB) = det(A)det(B) and ℂ is a field.)
        2) Matrix with zero determinant corresponds to constant transformation,
           but it makes no sense in this context.

        So if a = c = 0 or b = d = 0, then ad − bc = 0. Hence, ε(a, c) and ε(b, d) are defined.
        All variants for [a, b; c, d] are of the form [(uv)a, ub; (uv)c, ud] for u, v ∈ {±1, ±i}.
        So we take u = ε(b, d), w = ε(a, c), and v = w / u.
    */

    if (_absoluteXZ.a.isZero())
        _absoluteXZ.c.normalize();
    else
        _absoluteXZ.a.normalize(_absoluteXZ.c);

    if (_absoluteXZ.b.isZero())
        _absoluteXZ.d.normalize();
    else
        _absoluteXZ.b.normalize(_absoluteXZ.d);

    _absoluteOriginXZ = _absoluteXZ.origin();
}

WorldTile::~WorldTile() {
    faces.free(); edges.free();
}

bool WorldTile::walkable(int X, Real y, int Z) {
    using namespace Fundamentals;

    if (sizeTileX <= X || sizeTileZ <= Z) return true;
    return get<YALL>(X, std::floor(y), Z) == 0;
}

void drawParallelogram(FaceShader::VAO & vao, vec4 color, const Parallelogram<GLfloat> & P, GLfloat h) {
    auto index = vao.index();

    vao.emit(color, P.A.v3(h)); // + 0
    vao.emit(color, P.B.v3(h)); // + 1
    vao.emit(color, P.C.v3(h)); // + 2
    vao.emit(color, P.D.v3(h)); // + 3

    vao.push(index); vao.push(index + 1); vao.push(index + 2);
    vao.push(index); vao.push(index + 2); vao.push(index + 3);
}

void drawSide(FaceShader::VAO & vao, vec4 color, const Gyrovector<GLfloat> & A, const Gyrovector<GLfloat> & B, GLfloat h₁, GLfloat h₂) {
    auto index = vao.index();

    vao.emit(color, A.v3(h₁)); // + 0
    vao.emit(color, A.v3(h₂)); // + 1
    vao.emit(color, B.v3(h₂)); // + 2
    vao.emit(color, B.v3(h₁)); // + 3

    vao.push(index); vao.push(index + 1); vao.push(index + 2);
    vao.push(index); vao.push(index + 2); vao.push(index + 3);
}

struct Mask { bool top : 1, bottom : 1, back : 1, front : 1, left : 1, right : 1; };

void drawRightParallelogrammicPrism(FaceShader::VAO & vao, vec4 color, Mask m, GLfloat h, GLfloat Δh, const Parallelogram<GLfloat> & P) {
    const auto h₁ = h, h₂ = h + Δh;

    if (m.top)     drawParallelogram(vao, color, P, h₂);
    if (m.bottom)  drawParallelogram(vao, color, P.rev(), h₁);

    if (m.back)  drawSide(vao, color, P.B, P.A, h₁, h₂);
    if (m.right) drawSide(vao, color, P.C, P.B, h₁, h₂);
    if (m.front) drawSide(vao, color, P.D, P.C, h₁, h₂);
    if (m.left)  drawSide(vao, color, P.A, P.D, h₁, h₂);
}

template<typename T> inline Parallelogram<T> parallelogram(int X, int Z) {
    using namespace Tesselation;

    return {
        corners[X + 0][Z + 0], corners[X + 1][Z + 0],
        corners[X + 1][Z + 1], corners[X + 0][Z + 1]
    };
}

void drawNode(FaceShader::VAO & vao, vec4 color, Mask m, int X, int Y, int Z)
{ drawRightParallelogrammicPrism(vao, color, m, GLfloat(Y), 1.0f, parallelogram<GLfloat>(X, Z)); }

void WorldTile::emitFaces(const MapColorScheme & color) {
    using namespace Fundamentals;

    faces.clear();

    for (int Y = 0; Y < sizeTileY; Y++) {
        for (int X = 0; X < sizeTileX; X++) for (int Z = 0; Z < sizeTileZ; Z++) {
            auto value = get(X, Y, Z);

            if (value == 0) continue;

            Mask mask;

            mask.top    = get<YALL>(X + 0, Y + 1, Z + 0) == 0;
            mask.bottom = get<YALL>(X + 0, Y - 1, Z + 0) == 0;
            mask.back   = get<ZMIN>(X + 0, Y + 0, Z - 1) == 0;
            mask.front  = get<ZMAX>(X + 0, Y + 0, Z + 1) == 0;
            mask.left   = get<XMIN>(X - 1, Y + 0, Z + 0) == 0;
            mask.right  = get<XMAX>(X + 1, Y + 0, Z + 0) == 0;

            drawNode(faces, static_cast<vec4>(color[value]), mask, X, Y, Z);
        }

        facesOffsetY[Y] = faces.eboElementCount();
    }
}

inline bool isEdgeVisible(uint8_t n₀₀, uint8_t n₀₁, uint8_t n₁₀, uint8_t n₁₁) {
    bool b₀₀ = n₀₀ == 0, b₀₁ = n₀₁ == 0, b₁₀ = n₁₀ == 0, b₁₁ = n₁₁ == 0;

    if (b₀₀ && b₀₁ && b₁₀ && b₁₁)
        return false; // No blocks adjacent to the edge

    if (!b₀₀ && !b₀₁ && !b₁₀ && !b₁₁)
        return false; // The edge is blocked from all sides

    if (!b₀₀ && !b₀₁ && b₁₀ && b₁₁)
        return n₀₀ != n₀₁; // The edge is visible iff adjacent blocks have different colors

    if (b₀₀ && b₀₁ && !b₁₀ && !b₁₁)
        return n₁₀ != n₁₁;

    if (!b₀₀ && b₀₁ && !b₁₀ && b₁₁)
        return n₀₀ != n₁₀;

    if (b₀₀ && !b₀₁ && b₁₀ && !b₁₁)
        return n₀₁ != n₁₁;

    return true;
}

inline void emitLine(EdgeShader::VAO & vao, vec3 && v1, vec3 && v2) {
    vao.push(); vao.emit(v1);
    vao.push(); vao.emit(v2);
}

void WorldTile::emitEdges(const MapColorScheme &) {
    using namespace Fundamentals;

    using namespace Tesselation;

    edges.clear();

    for (int Y = 0; Y <= sizeTileY; Y++) {
        if (Y < sizeTileY) edgesLowerOffsetY[Y] = edges.eboElementCount();

        for (int X = 0; X < sizeTileX; X++) for (int Z = 0; Z <= sizeTileZ; Z++) {
            auto n₀₀ = get<YALL | ZMIN>(X, Y - 1, Z - 1);
            auto n₀₁ = get<YALL | ZMAX>(X, Y - 1, Z + 0);
            auto n₁₀ = get<YALL | ZMIN>(X, Y + 0, Z - 1);
            auto n₁₁ = get<YALL | ZMAX>(X, Y + 0, Z + 0);

            if (isEdgeVisible(n₀₀, n₀₁, n₁₀, n₁₁)) emitLine(edges, corners[X][Z].v3(Y), corners[X + 1][Z].v3(Y));
        }

        for (int X = 0; X <= sizeTileX; X++) for (int Z = 0; Z < sizeTileZ; Z++) {
            auto n₀₀ = get<XMIN | YALL>(X - 1, Y - 1, Z);
            auto n₀₁ = get<XMIN | YALL>(X - 1, Y + 0, Z);
            auto n₁₀ = get<XMAX | YALL>(X + 0, Y - 1, Z);
            auto n₁₁ = get<XMAX | YALL>(X + 0, Y + 0, Z);

            if (isEdgeVisible(n₀₀, n₀₁, n₁₀, n₁₁)) emitLine(edges, corners[X][Z].v3(Y), corners[X][Z + 1].v3(Y));
        }

        if (Y > 0) edgesUpperOffsetY[Y - 1] = edges.eboElementCount();

        if (Y < sizeTileY) for (int X = 0; X <= sizeTileX; X++) for (int Z = 0; Z <= sizeTileZ; Z++) {
            auto n₀₀ = get<XMIN | ZMIN>(X - 1, Y, Z - 1);
            auto n₀₁ = get<XMIN | ZMAX>(X - 1, Y, Z + 0);
            auto n₁₀ = get<XMAX | ZMIN>(X + 0, Y, Z - 1);
            auto n₁₁ = get<XMAX | ZMAX>(X + 0, Y, Z + 0);

            if (isEdgeVisible(n₀₀, n₀₁, n₁₀, n₁₁)) emitLine(edges, corners[X][Z].v3(Y), corners[X][Z].v3(Y + 1));
        }
    }
}

bool WorldTile::uploadVAO(const Fuchsian<Integer> & origin) {
    if (_cameraHorizontalDistance < 0)
        updateCameraXZ(origin);

    if (isToUploadVAO) {
        facesOffsetY.flip();
        edgesLowerOffsetY.flip();
        edgesUpperOffsetY.flip();

        faces.upload(GL_DYNAMIC_DRAW);
        edges.upload(GL_DYNAMIC_DRAW);

        isToUploadVAO = false;

        return true;
    }

    return false;
}

void WorldTile::updateVAO(const MapColorScheme & color) {
    if (isToUpdateVAO && !isToUploadVAO) {
        emitFaces(color);
        emitEdges(color);

        isToUpdateVAO = false;
        isToUploadVAO = true;
    }
}

void WorldTile::updateCameraXZ(const Fuchsian<Integer> & origin) {
    _cameraXZ = (origin.inverse() * _absoluteXZ).field<Real>();
    _cameraXZ.normalize();

    _cameraHorizontalDistance = _cameraXZ.origin().length();
}

template<ShaderSpec Spec>
inline void uploadDomain(WorldTile * tile, ShaderProgram<Spec> * shader) {
    shader->uniform("cameraTileXZ.a", tile->cameraXZ().a);
    shader->uniform("cameraTileXZ.b", tile->cameraXZ().b);
    shader->uniform("cameraTileXZ.c", tile->cameraXZ().c);
    shader->uniform("cameraTileXZ.d", tile->cameraXZ().d);
}

void WorldTile::renderFaces(FaceShader * shader, int Y₁, int Y₂) {
    // We assume that Y₂ − Y₁ ≥ sizeTileY and [Y₁; Y₂] ∩ [0; sizeTileY] ≠ ø.

    using namespace Fundamentals;

    uploadDomain(this, shader);

    if (Y₁ < 0) {
        shader->uniform<float>("cameraTileY", -sizeTileY);
        faces.draw(GL_TRIANGLES, facesLowerOffsetY(Y₁ + sizeTileY), facesUpperOffsetY(maxTileY));

        shader->uniform<float>("cameraTileY", 0);
        faces.draw(GL_TRIANGLES, 0, facesUpperOffsetY(Y₂));
    } else if (sizeTileY <= Y₂) {
        shader->uniform<float>("cameraTileY", 0);
        faces.draw(GL_TRIANGLES, facesLowerOffsetY(Y₁), facesUpperOffsetY(maxTileY));

        shader->uniform<float>("cameraTileY", sizeTileY);
        faces.draw(GL_TRIANGLES, 0, facesUpperOffsetY(Y₂ - sizeTileY));
    } else {
        shader->uniform<float>("cameraTileY", 0);
        faces.draw(GL_TRIANGLES, facesLowerOffsetY(Y₁), facesUpperOffsetY(Y₂));
    }
}

void WorldTile::renderEdges(EdgeShader * shader, int Y₁, int Y₂) {
    using namespace Fundamentals;

    uploadDomain(this, shader);

    if (Y₁ < 0) {
        shader->uniform<float>("cameraTileY", -sizeTileY);
        edges.draw(GL_LINES, edgesLowerOffsetY(Y₁ + sizeTileY), edgesUpperOffsetY(maxTileY));

        shader->uniform<float>("cameraTileY", 0);
        edges.draw(GL_LINES, 0, edgesUpperOffsetY(Y₂));
    } else if (sizeTileY <= Y₂) {
        shader->uniform<float>("cameraTileY", 0);
        edges.draw(GL_LINES, edgesLowerOffsetY(Y₁), edgesUpperOffsetY(maxTileY));

        shader->uniform<float>("cameraTileY", sizeTileY);
        edges.draw(GL_LINES, 0, edgesUpperOffsetY(Y₂ - sizeTileY));
    } else {
        shader->uniform<float>("cameraTileY", 0);
        edges.draw(GL_LINES, edgesLowerOffsetY(Y₁), edgesUpperOffsetY(Y₂));
    }
}

bool WorldTile::touch(const Gyrovector<Real> & w, int X, int Z) {
    const auto & A = Tesselation::corners[X + 0][Z + 0];
    const auto & B = Tesselation::corners[X + 1][Z + 0];
    const auto & C = Tesselation::corners[X + 1][Z + 1];
    const auto & D = Tesselation::corners[X + 0][Z + 1];

    return Math::samesign(
        w.sub(A).cross(B.sub(A)),
        w.sub(B).cross(C.sub(B)),
        w.sub(C).cross(D.sub(C)),
        w.sub(D).cross(A.sub(D))
    );
}

std::pair<int, int> WorldTile::round(const Gyrovector<Real> & w)
{ return Tesselation::unapply(w.x(), w.z()); }

bool WorldTile::isInsideOfDomain(const Gyrovector<Real> & w₀) {
    using namespace Fundamentals;

    // We are using symmetry of grid along axes here
    Gyrovector<Real> w(fabs(w₀.x()), fabs(w₀.z()));

    for (int Z = 0; Z < sizeTileZ; Z++) {
        const auto & A = Tesselation::corners[sizeTileX][Z + 0];
        const auto & B = Tesselation::corners[sizeTileX][Z + 1];

        if (Math::samesign(w.sub(A).cross(B.sub(A)), w.sub(B).cross(-B), w.cross(A)))
            return true;
    }

    return false;
}

std::optional<size_t> WorldTile::matchNeighbour(const Gyrovector<Real> & P) {
    for (size_t k = 0; k < Tesselation::neighbours.size(); k++) {
        const auto & Δ⁻¹ = Tesselation::neighbours⁻¹[k];
        if (WorldTile::isInsideOfDomain(Δ⁻¹.apply(P)))
            return std::optional(k);
    }

    return std::nullopt;
}

WorldMap::WorldMap() {}
WorldMap::~WorldMap() {}

WorldTile * WorldMap::findNonThreadsafe(const Gaussian²<Integer> & absoluteOriginXZ) {
    for (auto tile : writeQueue)
        if (tile->absoluteOriginXZ() == absoluteOriginXZ)
            return tile;

    return nullptr;
}

WorldTile * WorldMap::find(const Gaussian²<Integer> & absoluteOriginXZ) {
    std::lock_guard<std::mutex> guardRead(readMutex);

    for (auto tile : readQueue)
        if (tile->absoluteOriginXZ() == absoluteOriginXZ)
            return tile;

    return nullptr;
}

void WorldMap::updateCameraXZ(const Fuchsian<Integer> & origin) {
    for (auto tile : readQueue)
        tile->updateCameraXZ(origin);
}

void WorldTile::load(HVXL & file) {
    file.readTileDataAt(fileOffset, _vxl, fog);
    _isModified = false; // ‘Unmodified’ here means ‘synchronized with the disk’
}

void WorldTile::save(HVXL & file) {
    file.writeTileDataAt(fileOffset, _vxl, fog);
    _isModified = false;
}

void WorldMap::save() {
    for (auto tile : writeQueue)
        if (tile->isModified())
            tile->save(file);
}

WorldTile * WorldMap::queueNewTile(const Fuchsian<Integer> & absoluteXZ) {
    auto tile = new WorldTile(absoluteXZ);

    if (file.findTileOffset(tile->absoluteOriginXZ(), tile->fileOffset, tile->hash)) {
        if (mapgen != nullptr) (*mapgen)(tile);

        tile->threadsafeSave();
    } else {
        tile->load(file);
    }

    tile->threadsafeUpdateVAO();

    writeQueue.push_back(tile);

    return tile;
}

constexpr double mapSaveTimeout = 15.0;

void WorldMap::update(const Real & hrd) {
    static double T = 0.0, mapLastSaveTime = 0.0;

    auto dt = glfwGetTime() - T; T += dt;

    /* It is here to eliminate some precision issues. The problem is that the distance
       calculated before `queueNewTile` happens to be a little different than that one
       calculated in `updateCameraXZ`. The difference is less than 10⁻⁸ but non-zero,
       thus enough to cause subtle comparison-related problems showing as some tiles
       loading and then unloading many times per second. */
    constexpr auto ε = 0.001;

    for (auto it = writeQueue.begin(); it != writeQueue.end();) {
        auto tile = *it;

        if (tile->cameraHorizontalDistance() > hrd + ε) {
            tile->isToBeErased = true;
            it = writeQueue.erase(it);
        } else it++;
    }

    size_t remSize = writeQueue.size();

    {
        using namespace Game;
        // TODO: this part is a little bit dirty, need to figure out a better way.

        if (player.tile() == nullptr) {
            auto & M = player.r().absoluteXZ();

            if (findNonThreadsafe(M.origin()) == nullptr)
                queueNewTile(M);

            player.setXYZ();
        }
    }

    for (size_t i = 0; i < remSize; i++) {
        auto tile = writeQueue[i];

        tile->updateVAO(file.color);

        if (0 <= tile->cameraHorizontalDistance()) {
            // TODO: maybe we should align this with `neighbours` above?
            static const std::array vonNeumannNeighbours = {
                Tesselation::U, Tesselation::D, Tesselation::L, Tesselation::R
            };

            for (const auto & dM : vonNeumannNeighbours) {
                auto dr = dM.field<Real>().origin();

                if (tile->cameraXZ().apply(dr).length() < hrd - ε) {
                    auto M = tile->absoluteXZ() * dM;

                    // TODO: this is really slow, maybe maintain a hashtable instead?
                    if (findNonThreadsafe(M.origin()) == nullptr)
                        queueNewTile(M);
                }
            }
        }
    }

    {
        // This may lock the rendering thread just for few nanoseconds.
        std::lock_guard<std::mutex> guardRead(readMutex);
        std::swap(readQueue, writeQueue);
    }

    bool isMapSaveDue = T - mapLastSaveTime > mapSaveTimeout;
    if (isMapSaveDue) mapLastSaveTime = T;

    for (auto tile : writeQueue) {
        if (isMapSaveDue || tile->isToBeErased)
            if (tile->isModified())
                tile->save(file);

        if (tile->isToBeErased)
            delete tile;
    }

    /* TODO: it’s probabliy too expensive to just copy everything,
             so it may be better to set up some queue instead. */
    writeQueue = readQueue;
}