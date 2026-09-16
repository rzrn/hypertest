#include <Hyper/Geometry.hxx>

namespace Tesselation {
    // Chunk’s neighbours in tesselation
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

    // Generation of chunk’s grid

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
        using namespace Fundamentals;

        auto t₁ = 2 * Real(X) / sizeTileX - 1;
        auto t₂ = 2 * Real(Z) / sizeTileZ - 1;
        return Ψ(t₁, t₂);
    }

    auto unapply(Real u, Real v) {
        auto [t₁, t₂] = Ψ⁻¹(u, v);

        /* Chunk’s border is not exactly a hyperbolic line (i.e. circular arc on the Poincaré disk),
           but its piecewise linear approximation; so there are parts of the outer blocks that extend
           slightly beyond the boundary of the ideal hyperbolic square.
           That’s why we need to “std::clamp” here.
        */
        t₁ = std::clamp<Real>(t₁, -1.0, 0.9999); // t₁ ≤ 0.9999 < 1 so that int(X) < sizeTileX
        t₂ = std::clamp<Real>(t₂, -1.0, 0.9999);

        auto X = (t₁ + 1) / 2 * sizeTileX;
        auto Z = (t₂ + 1) / 2 * sizeTileZ;

        return std::pair(int(X), int(Z));
    }

    constexpr auto init() {
        using namespace Fundamentals;

        Array²<Gyrovector<Real>, sizeTileX + 1, sizeTileZ + 1> retval;

        for (int X = 0; X <= sizeTileX; X++)
            for (int Z = 0; Z <= sizeTileZ; Z++)
                retval[X][Z] = apply(X, Z);

        return retval;
    }

    constexpr Grid corners = init();

    constexpr auto distance(int X₁, int Z₁, int X₂, int Z₂)
    { return (-corners[X₁][Z₁] + corners[X₂][Z₂]).abs(); }

    constexpr Real meter = distance(sizeTileX / 2, sizeTileZ / 2, sizeTileX / 2, sizeTileZ / 2 + 1);
}

NodeRegistry::NodeRegistry() {
    attach({"Air", {
        Texture(), Texture(), Texture(),
        Texture(), Texture(), Texture()
    }});
}

Chunk::Chunk(const Fuchsian<Integer> & origin, const Fuchsian<Integer> & isometry) : _isometry(isometry) {
    /*
        Unfortunately, precomposition of `isometry` with (z ↦ z × exp(iπk/2)) for k ∈ ℤ
        will yield matrix able to render this chunk in the same place but rotated about its own center by πk/2 radians.

        So if we don’t resolve this ambiguity, the same chunk may render with different rotation
        according only to player’s path, and this will definitely crush the landscape
        (Imagine random rotations of chunks in Minecraft in a mountainous biome.)

        Since exp(iπk/2) ∈ {±1, ±i}, such rotation is equivalent to multiplying `a` and `c` by ±1/±i at the same time.
        `isometry` is pre-divided by hcf(a, b, c, d), so there is always ability
        to multiply *all* terms by ±1/±i yielding the same transformation.
        (Two Möbius transformations are equal iff their matrices differ by multiplicative constant.)

        Hence we have 4 × 4 = 16 options. It’s important that we can always select multipliers
        so that both components of `a` and `b` will be non-negative, that’s what we’re doing.
    */

    // We assume that det(_isometry) = ad − bc ≠ 0, because:
    //  1) det(I), det(U), det(L), det(D), det(R) ≠ 0 (see above),
    //     so determinant from any of their product is also non-zero.
    //     (Since det(AB) = det(A)det(B) & ℂ is a field.)
    //  2) Matrix with zero determinant corresponds to constant transformation,
    //     but it makes no sense in this context.
    if (!_isometry.a.isZero()) {
        // a ≠ 0 and b ≠ 0
        if (!_isometry.b.isZero()) _isometry.b.normalize(_isometry.a, _isometry.c, _isometry.d);
        // det(_isometry) = ad − bc = ad ≠ 0, so a ≠ 0 and d ≠ 0
        else _isometry.d.normalize(_isometry.a, _isometry.c);

        _isometry.a.normalize(_isometry.c);
    } else {
        // det(_isometry) = ad − bc = −bc ≠ 0, so b ≠ 0 and c ≠ 0
        _isometry.b.normalize(_isometry.c, _isometry.d);
        _isometry.c.normalize();
    }

    _pos = isometry.origin();
    updateMatrix(origin);

    faces.initialize();
    edges.initialize();
}

Chunk::~Chunk() { join(); delete _blob; faces.free(); edges.free(); }

bool Chunk::walkable(int X, Real y, int Z) {
    using namespace Fundamentals;

    if (sizeTileX <= X || sizeTileZ <= Z) return true;
    return get<YALL>(X, std::floor(y), Z).id == 0;
}

void drawParallelogram(FaceShader::VAO & vao, Texture & T, const Parallelogram<GLfloat> & P, GLfloat h) {
    auto index = vao.index();

    vao.emit(T.lu(), P.A.v3(h)); // + 0
    vao.emit(T.ru(), P.B.v3(h)); // + 1
    vao.emit(T.rd(), P.C.v3(h)); // + 2
    vao.emit(T.ld(), P.D.v3(h)); // + 3

    vao.push(index); vao.push(index + 1); vao.push(index + 2);
    vao.push(index); vao.push(index + 2); vao.push(index + 3);
}

void drawSide(FaceShader::VAO & vao, Texture & T, const Gyrovector<GLfloat> & A, const Gyrovector<GLfloat> & B, GLfloat h₁, GLfloat h₂) {
    auto index = vao.index();

    vao.emit(T.rd(), A.v3(h₁)); // + 0
    vao.emit(T.ru(), A.v3(h₂)); // + 1
    vao.emit(T.lu(), B.v3(h₂)); // + 2
    vao.emit(T.ld(), B.v3(h₁)); // + 3

    vao.push(index); vao.push(index + 1); vao.push(index + 2);
    vao.push(index); vao.push(index + 2); vao.push(index + 3);
}

struct Mask { bool top : 1, bottom : 1, back : 1, front : 1, left : 1, right : 1; };

void drawRightParallelogrammicPrism(FaceShader::VAO & vao, Cube & C, Mask m, GLfloat h, GLfloat Δh, const Parallelogram<GLfloat> & P) {
    const auto h₁ = h, h₂ = h + Δh;

    if (m.top)     drawParallelogram(vao, C.top, P, h₂);
    if (m.bottom)  drawParallelogram(vao, C.bottom, P.rev(), h₁);

    if (m.back)  drawSide(vao, C.back,  P.B, P.A, h₁, h₂);
    if (m.right) drawSide(vao, C.right, P.C, P.B, h₁, h₂);
    if (m.front) drawSide(vao, C.front, P.D, P.C, h₁, h₂);
    if (m.left)  drawSide(vao, C.left,  P.A, P.D, h₁, h₂);
}

template<typename T> inline Parallelogram<T> parallelogram(int X, int Z) {
    using namespace Tesselation;

    return {
        corners[X + 0][Z + 0], corners[X + 1][Z + 0],
        corners[X + 1][Z + 1], corners[X + 0][Z + 1]
    };
}

void drawNode(FaceShader::VAO & vao, Cube & C, Mask m, int X, int Y, int Z)
{ drawRightParallelogrammicPrism(vao, C, m, GLfloat(Y), 1.0f, parallelogram<GLfloat>(X, Z)); }

void Chunk::emitFaces(NodeRegistry & nodeRegistry) {
    using namespace Fundamentals;

    faces.clear();

    for (int Y = 0; Y < sizeTileY; Y++) {
        for (int X = 0; X < sizeTileX; X++) for (int Z = 0; Z < sizeTileZ; Z++) {
            auto id = get(X, Y, Z).id;

            if (id == 0) continue;

            Mask mask;

            mask.top    = get<YALL>(X + 0, Y + 1, Z + 0).id == 0;
            mask.bottom = get<YALL>(X + 0, Y - 1, Z + 0).id == 0;
            mask.back   = get<ZMIN>(X + 0, Y + 0, Z - 1).id == 0;
            mask.front  = get<ZMAX>(X + 0, Y + 0, Z + 1).id == 0;
            mask.left   = get<XMIN>(X - 1, Y + 0, Z + 0).id == 0;
            mask.right  = get<XMAX>(X + 1, Y + 0, Z + 0).id == 0;

            if (nodeRegistry.has(id)) {
                auto nodeDef = nodeRegistry.get(id);
                drawNode(faces, nodeDef.cube, mask, X, Y, Z);
            }
        }

        facesOffsetY[Y] = faces.eboElementCount();
    }
}

inline bool isEdgeVisible(Node & n₀₀, Node & n₀₁, Node & n₁₀, Node & n₁₁) {
    bool b₀₀ = n₀₀.id == 0, b₀₁ = n₀₁.id == 0, b₁₀ = n₁₀.id == 0, b₁₁ = n₁₁.id == 0;

    if (b₀₀ && b₀₁ && b₁₀ && b₁₁)
        return false; // No blocks adjacent to the edge

    if (!b₀₀ && !b₀₁ && !b₁₀ && !b₁₁)
        return false; // The edge is blocked from all sides

    if (!b₀₀ && !b₀₁ && b₁₀ && b₁₁)
        return n₀₀.id != n₀₁.id; // The edge is visible iff adjacent blocks have different colors

    if (b₀₀ && b₀₁ && !b₁₀ && !b₁₁)
        return n₁₀.id != n₁₁.id;

    if (!b₀₀ && b₀₁ && !b₁₀ && b₁₁)
        return n₀₀.id != n₁₀.id;

    if (b₀₀ && !b₀₁ && b₁₀ && !b₁₁)
        return n₀₁.id != n₁₁.id;

    return true;
}

inline void emitLine(EdgeShader::VAO & vao, vec3 && v1, vec3 && v2) {
    vao.push(); vao.emit(v1);
    vao.push(); vao.emit(v2);
}

void Chunk::emitEdges(NodeRegistry &) {
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

void Chunk::refresh(NodeRegistry & nodeRegistry) {
    if (needUpdateVAO) {
        facesOffsetY.flip();
        edgesLowerOffsetY.flip();
        edgesUpperOffsetY.flip();

        faces.upload(GL_DYNAMIC_DRAW);
        edges.upload(GL_DYNAMIC_DRAW);

        needUpdateVAO = false;
        _needRefresh  = false;
        return;
    }

    if (working()) return; _working = true;

    worker = std::async(std::launch::async, [&nodeRegistry, this]() mutable {
        emitFaces(nodeRegistry);
        emitEdges(nodeRegistry);

        needUpdateVAO = true;
        _working = false;
    });
}

void Chunk::updateMatrix(const Fuchsian<Integer> & origin) {
    _domain = (origin.inverse() * _isometry).field<Real>();
    _domain.normalize();

    _awayness = _domain.origin().abs();
}

template<ShaderSpec Spec>
inline void uploadDomain(Chunk * chunk, ShaderProgram<Spec> * shader) {
    shader->uniform("cameraTileXZ.a", chunk->domain().a);
    shader->uniform("cameraTileXZ.b", chunk->domain().b);
    shader->uniform("cameraTileXZ.c", chunk->domain().c);
    shader->uniform("cameraTileXZ.d", chunk->domain().d);
}

void Chunk::renderFaces(FaceShader * shader, int Y₁, int Y₂) {
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

void Chunk::renderEdges(EdgeShader * shader, int Y₁, int Y₂) {
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

bool Chunk::touch(const Gyrovector<Real> & w, int X, int Z) {
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

std::pair<int, int> Chunk::round(const Gyrovector<Real> & w)
{ return Tesselation::unapply(w.x(), w.y()); }

bool Chunk::isInsideOfDomain(const Gyrovector<Real> & w₀) {
    using namespace Fundamentals;

    // We are using symmetry of grid along axes here
    Gyrovector<Real> w(fabs(w₀.x()), fabs(w₀.y()));

    for (int Z = 0; Z < sizeTileZ; Z++) {
        const auto & A = Tesselation::corners[sizeTileX][Z + 0];
        const auto & B = Tesselation::corners[sizeTileX][Z + 1];

        if (Math::samesign(w.sub(A).cross(B.sub(A)), w.sub(B).cross(-B), w.cross(A)))
            return true;
    }

    return false;
}

std::optional<size_t> Chunk::matchNeighbour(const Gyrovector<Real> & P) {
    for (size_t k = 0; k < Tesselation::neighbours.size(); k++) {
        const auto & Δ⁻¹ = Tesselation::neighbours⁻¹[k];
        if (Chunk::isInsideOfDomain(Δ⁻¹.apply(P)))
            return std::optional(k);
    }

    return std::nullopt;
}

Atlas::Atlas() {}
Atlas::~Atlas() {}

Chunk * Atlas::lookup(const Gaussian²<Integer> & pos) {
    for (auto chunk : pool)
        if (chunk->pos() == pos)
            return chunk;

    return nullptr;
}

Chunk * Atlas::poll(const Fuchsian<Integer> & origin, const Fuchsian<Integer> & isometry) {
    auto pos = isometry.origin();

    for (auto chunk : pool)
        if (chunk->pos() == pos)
            return chunk;

    auto chunk = new Chunk(origin, isometry); pool.push_back(chunk);
    chunk->load(generator, engine); return chunk;
}

void Atlas::updateMatrix(const Fuchsian<Integer> & origin) {
    for (auto chunk : pool)
        chunk->updateMatrix(origin);
}

const char * initcmd   = "CREATE TABLE IF NOT EXISTS atlas("
                         "bitfield INTEGER, real1 BLOB, imag1 BLOB, real2 BLOB, imag2 BLOB,"
                         "blob BLOB, PRIMARY KEY (bitfield, real1, imag1, real2, imag2));",
           * loadcmd   = "SELECT blob FROM atlas WHERE bitfield = ? AND real1 = ? AND imag1 = ? AND real2 = ? AND imag2 = ?;",
           * insertcmd = "INSERT or REPLACE INTO atlas(bitfield, real1, imag1, real2, imag2, blob) VALUES(?, ?, ?, ?, ?, ?);";

inline void warn(sqlite3 * engine)
{ std::fprintf(stderr, "SQLITE: %s\n", sqlite3_errmsg(engine)); }

void Atlas::connect(std::string & filename) {
    auto retval = sqlite3_open(filename.c_str(), &engine);

    if (retval != SQLITE_OK) {
        warn(engine); sqlite3_close(engine);
        throw std::runtime_error("`sqlite3_open` failed");
    }

    char * errmsg; retval = sqlite3_exec(engine, initcmd, nullptr, 0, &errmsg);

    if (retval != SQLITE_OK) {
        std::fprintf(stderr, "SQLITE: %s\n", errmsg); sqlite3_free(errmsg);

        throw std::runtime_error("sqlite3 initialization failed");
    }
}

void Atlas::disconnect() {
    dump();

    for (auto chunk : pool)
        chunk->join();

    sqlite3_close(engine);
}

inline void dumpBlob(sqlite3_stmt * statement, int index, void * blob, size_t n) {
    static uint8_t zero = 0;

    if (blob == nullptr || n == 0)
        sqlite3_bind_blob(statement, index, &zero, 1, SQLITE_STATIC);
    else
        sqlite3_bind_blob(statement, index, blob, n, free);
}

void dumpGaussian(sqlite3_stmt * statement, const Gaussian<Integer> & z, int idx₁, int idx₂) {
    size_t k₁, k₂;

    auto blob₁ = Math::serialize(z.real, k₁);
    auto blob₂ = Math::serialize(z.imag, k₂);

    dumpBlob(statement, idx₁, blob₁, k₁);
    dumpBlob(statement, idx₂, blob₂, k₂);
}

void Chunk::serialize(sqlite3_stmt * statement, int idx₀, int idx₁, int idx₂, int idx₃, int idx₄) {
    Bitfield<uint8_t> bitfield(0);

    bitfield.set(0, Math::isNeg(_pos.first.real));
    bitfield.set(1, Math::isNeg(_pos.first.imag));
    bitfield.set(2, Math::isNeg(_pos.second.real));
    bitfield.set(3, Math::isNeg(_pos.second.imag));

    sqlite3_bind_int(statement, idx₀, uint8_t(bitfield));
    dumpGaussian(statement, _pos.first,  idx₁, idx₂);
    dumpGaussian(statement, _pos.second, idx₃, idx₄);
}

void Chunk::load(ChunkOperator * generator, sqlite3 * engine) {
    sqlite3_stmt * statement = nullptr;

    if (_ready) return; _working = true;
    worker = std::async(std::launch::async, [statement, retval = 0, generator, engine, this]() mutable {
        _blob = new Blob();

        retval = sqlite3_prepare_v2(engine, loadcmd, -1, &statement, nullptr);
        if (retval != SQLITE_OK) { warn(engine); _needUnload = true; goto fin; }

        serialize(statement, 1, 2, 3, 4, 5);
        retval = sqlite3_step(statement);

        if (retval == SQLITE_ROW)
            memcpy(_blob, sqlite3_column_blob(statement, 0), sizeof(Blob));
        else { if (generator != nullptr) (*generator)(this); _dirty = true; }

        if (retval == SQLITE_ERROR) warn(engine);
        sqlite3_finalize(statement); requestRefresh();

        fin: _ready = true; _working = false;
    });
}

void Chunk::join() { worker.wait(); }

void Chunk::dump(sqlite3 * engine) {
    sqlite3_stmt * statement = nullptr;

    if (working()) return; _working = true;
    worker = std::async(std::launch::async, [statement, retval = 0, engine, this]() mutable {
        retval = sqlite3_prepare_v2(engine, insertcmd, -1, &statement, nullptr);
        if (retval != SQLITE_OK) { warn(engine); _working = false; return; }

        serialize(statement, 1, 2, 3, 4, 5);
        sqlite3_bind_blob(statement, 6, _blob, sizeof(Blob), SQLITE_TRANSIENT);

        retval = sqlite3_step(statement);

        if (retval != SQLITE_DONE) warn(engine); else _dirty = false;
        sqlite3_finalize(statement); _working = false;
    });
}

void Atlas::dump() {
    for (auto chunk : pool)
        if (chunk->dirty())
            chunk->dump(engine);
}
