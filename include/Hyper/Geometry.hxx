#pragma once

#include <optional>
#include <string>
#include <vector>

#include <future>
#include <chrono>

#include <GL/glew.h>
#include <sqlite3.h>

#include <Hyper/Fundamentals.hxx>
#include <Hyper/Shader.hxx>
#include <Hyper/Sheet.hxx>

#include <Math/Fuchsian.hxx>
#include <Math/AutD.hxx>

#include <Meta/DoubleBuffer.hxx>
#include <Meta/List.hxx>

using ℤi = Gaussian<Integer>;

namespace Tesselation {
    using namespace Fundamentals;

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
    using Grid = Array²<Gyrovector<Real>, sizeTileX + 1, sizeTileZ + 1>;

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

struct Cube { Texture top, bottom, left, right, front, back; };
struct NodeDef { std::string name; Cube cube; };

struct Node { NodeId id; };

class NodeRegistry {
private:
    NodeDef air; std::vector<NodeDef> table;

public:
    NodeRegistry();

    inline NodeId attach(const NodeDef & def)
    { table.push_back(def); return table.size() - 1; }

    inline NodeDef get(NodeId id) { return table[id]; }
    inline bool has(NodeId id) { return id < table.size(); }
};

struct Blob {
    Node data[Fundamentals::sizeTileX][Fundamentals::sizeTileY][Fundamentals::sizeTileZ];

    Blob() : data{} {}
} __attribute__((packed));

class Chunk; using ChunkOperator = void(Chunk *);

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

class Chunk {
private:
    Fuchsian<Integer> _isometry; Möbius<Real> _domain; Real _awayness; // used for drawing
    Gaussian²<Integer> _pos; // used for indexing, should be equal to `isometry.origin()`

    bool _working = false; std::future<void> worker;
    FaceShader::VAO faces; EdgeShader::VAO edges;

    DoubleBuffer<GLsizei, Fundamentals::sizeTileY> facesOffsetY, edgesLowerOffsetY, edgesUpperOffsetY;

    inline GLsizei facesLowerOffsetY(const int Y) const { return Y > 0 ? facesOffsetY(Y - 1) : 0; }
    inline GLsizei facesUpperOffsetY(const int Y) const { return facesOffsetY(Y);                 }

    bool _ready = false, _dirty = false, _needRefresh = false, _needUnload = false, needUpdateVAO = false;

    Blob * _blob = nullptr;
public:

    Chunk(const Fuchsian<Integer> & origin, const Fuchsian<Integer> & isometry);

    ~Chunk();

    void emitFaces(NodeRegistry &);
    void emitEdges(NodeRegistry &);

    void renderFaces(FaceShader *, int, int);
    void renderEdges(EdgeShader *, int, int);

    void updateMatrix(const Fuchsian<Integer> &);
    void refresh(NodeRegistry &);

    bool walkable(int, Real, int);

    void serialize(sqlite3_stmt *, int, int, int, int, int);
    void load(ChunkOperator *, sqlite3 *);
    void dump(sqlite3 *);
    void join();

    inline bool working() const { return _working; }

    inline constexpr bool ready()       { return _ready;       }
    inline constexpr bool dirty()       { return _dirty;       }
    inline constexpr bool needRefresh() { return _needRefresh; }
    inline constexpr bool needUnload()  { return _needUnload;  }

    inline constexpr void unload()         { _needUnload = true;  }
    inline constexpr auto requestRefresh() { _needRefresh = true; }

    inline constexpr auto awayness() const { return _awayness; }

    inline const auto isometry() const { return _isometry; }
    inline const auto domain()   const { return _domain;   }
    inline const auto pos()      const { return _pos;      }

    inline Blob * blob() { if (_blob != nullptr) _dirty = true; return _blob; }
    inline const Blob * blob() const { return _blob; }

    static inline int mod(int a, int b) { int r = a % b; return r < 0 ? r + b : r; }

    // TODO: make this to return a reference instead
    template<unsigned int mask = NONE> inline Node get(int X, int Y, int Z) const {
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

        return _blob->data[X][Y][Z];
    }

    inline void set(int X, int Y, int Z, const Node & node)
    { _dirty = true; _blob->data[X][Y][Z] = node; }

    static bool touch(const Gyrovector<Real> &, int, int);
    static std::pair<int, int> round(const Gyrovector<Real> &);

    static bool isInsideOfDomain(const Gyrovector<Real> &);
    static std::optional<size_t> matchNeighbour(const Gyrovector<Real> &);

    static inline Real clamp(Real x)
    { return Math::remainder<Real>(x, Fundamentals::sizeTileY); }
};

class Atlas {
private:
    sqlite3 * engine;

public:
    std::vector<Chunk *> pool;
    ChunkOperator * generator = nullptr;

    Atlas();
    ~Atlas();

    void connect(std::string &);
    void disconnect();

    void dump();

    Chunk * poll(const Fuchsian<Integer> & origin, const Fuchsian<Integer> & isometry);
    Chunk * lookup(const Gaussian²<Integer> &);

    void updateMatrix(const Fuchsian<Integer> &);
};

template<typename T> struct Bitfield {
    T value;

    Bitfield(T value) : value(value) {}
    inline operator T() { return value; }

    inline void set(size_t n, bool bit) { value = (value | (1 << n)) & ~(T(!bit) << n); }
    inline bool get(size_t n) const { return (value >> n) & 1; }
};
