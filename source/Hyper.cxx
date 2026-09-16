#include <cstdio>

#include <glm/ext/matrix_clip_space.hpp>
#include <glm/ext/matrix_transform.hpp>
#include <glm/mat4x4.hpp>

#include <Lua.hxx>

#include <Math/Gyrovector.hxx>
#include <Math/Fuchsian.hxx>

#include <Hyper/Config.hxx>
#include <Hyper/Shader.hxx>
#include <Hyper/Game.hxx>

void errorCallback(int error, const char * description) {
    fprintf(stderr, "Error: %s\n", description);
}

glm::mat4 view, projection;

DummyShader * dummyShader = nullptr;
FaceShader  * faceShader  = nullptr;
EdgeShader  * edgeShader  = nullptr;

DummyShader::VAO aimVao;

PBO<GLfloat, Action> pbo(GL_DEPTH_COMPONENT, 1, 1);

const auto origin = vec2(0.0f);

const auto white  = vec4(1.0f, 1.0f, 1.0f, 1.0f);
const auto black  = vec4(0.0f, 0.0f, 0.0f, 1.0f);

void drawAim(DummyShader::VAO & vao) {
    using namespace Game;

    vao.clear();

    auto wpixel = 1.0 / GLfloat(Window::width), hpixel = 1.0 / GLfloat(Window::height);

    vao.push(); vao.emit(vec3(-GLfloat(GUI::aimSize) * wpixel, 0, 0), white, origin, 1.0f);
    vao.push(); vao.emit(vec3(+GLfloat(GUI::aimSize) * wpixel, 0, 0), white, origin, 1.0f);
    vao.push(); vao.emit(vec3(0, -GLfloat(GUI::aimSize) * hpixel, 0), white, origin, 1.0f);
    vao.push(); vao.emit(vec3(0, +GLfloat(GUI::aimSize) * hpixel, 0), white, origin, 1.0f);

    vao.upload(GL_STATIC_DRAW);
}

void updateHotbar() {}

Real worldTileDiameter(const Real n) {
    using namespace Fundamentals;

    constexpr Gyrovector<Real> i(D½, 0), j(0, D½), k = Coadd(i, j);

    return (n * k).abs();
}

bool move(Entity & E, const Gyrovector<Real> & v, Real Δt) {
    constexpr Real Δtₘₐₓ = 1.0/5.0; bool P = false;

    while (Δt >= Δtₘₐₓ) { auto Q = E.move(v, Δtₘₐₓ); P = P || Q; Δt -= Δtₘₐₓ; }
    auto R = E.move(v, Δt); return P || R;
}

vec3 unproject(const glm::mat4 & view, const glm::mat4 & projection, const GLfloat depth) {
    auto v = glm::inverse(view) * glm::inverse(projection) * glm::vec4(0.0f, 0.0f, depth, 1.0f);
    return vec3(v.x / v.w, v.y / v.w, v.z / v.w);
}

vec3 trace(const glm::mat4 & view, const glm::mat4 & projection, const GLfloat zbuffer, const GLfloat y, bool forward) {
    using namespace Game;
    using namespace Render;

    auto v₀ = vec3(0.0f, y, 0.0f);

    auto v = standard->model.unapply(unproject(view, projection, 2.0f * zbuffer - 1.0f));
    auto dist₀ = glm::length(v - v₀);

    const GLfloat ε = standard->meter / 3.0f;
    GLfloat dist = dist₀ + (forward ? +ε : -ε);

    return (dist / dist₀) * (v - v₀) + v₀;
}

std::optional<std::pair<WorldTile *, Gyrovector<Real>>> getNeighbour(const Gyrovector<Real> & P) {
    using namespace Game;

    if (WorldTile::isInsideOfDomain(P)) {
        auto Q = player.tile()->cameraXZ().inverse().apply(P);
        return std::optional(std::pair(player.tile(), Q));
    }

    for (size_t k = 0; k < Tesselation::neighbours.size(); k++) {
        auto G = player.tile()->absoluteXZ() * Tesselation::neighbours[k];
        if (auto C = map.lookup(G.origin())) {
            auto Q = C->cameraXZ().inverse().apply(P);

            if (WorldTile::isInsideOfDomain(Q))
                return std::optional(std::pair(C, Q));
        }
    }

    return std::nullopt;
}

void setBlock(WorldTile * C, int X, Real y, int Z, NodeId id) {
    if (C == nullptr || !C->ready())
        return;

    if (Fundamentals::maxTileX < X || Fundamentals::maxTileZ < Z)
        return;

    int Y = std::floor(WorldTile::clamp(y));

    if (id != 0 && C->get(X, Y, Z).id != 0)
        return;

    C->set(X, Y, Z, {id});

    if (Game::player.stuck())
        C->set(X, Y, Z, {0});

    C->requestRefresh();
}

void click(const Aut𝔻<Real> & origin, const GLfloat zbuffer, const Action action) {
    const auto maxₕ = 5.0 * Tesselation::meter, maxᵥ = 4.0;
    const auto y₀ = Game::player.camera().climb + Game::player.eye;

    auto v = trace(view, projection, zbuffer, y₀, action == Action::Remove);
    auto P = Gyrovector(v.x, v.z);

    if (P.abs() <= maxₕ && fabs(v.y - y₀) <= maxᵥ) {
        if (auto ret = getNeighbour(origin.inverse().apply(P))) {
            auto [C, Q] = *ret; auto [i, k] = WorldTile::round(Q);

            if (action == Action::Place && Game::activeSlot < Game::hotbarSize) {
                auto id = Game::hotbar[Game::activeSlot];
                if (id != 0 && Game::Registry::node.has(id))
                    setBlock(C, i, v.y, k, id);
            }

            if (action == Action::Remove) setBlock(C, i, v.y, k, 0);
        }
    }
}

void pollNeighbours() {
    using namespace Game;

    map.updateMatrix(player.camera().position.absoluteXZ());

    for (size_t k = 0; k < Tesselation::neighbours.size(); k++) {
        auto G = player.tile()->absoluteXZ() * Tesselation::neighbours[k];
        map.poll(player.camera().position.absoluteXZ(), G);
    }

    /*for (size_t i = 0; i < Tesselation::neighbours.size(); i++)
        for (size_t j = 0; j < Tesselation::neighbours.size(); j++) {
            auto G = player.tile()->absoluteXZ() * Tesselation::neighbours[i] * Tesselation::neighbours[j];
            map.poll(player.camera().position.absoluteXZ(), G);
    }*/
}

template<ShaderSpec Spec>
inline void uploadMVP(ShaderProgram<Spec> * shader, Aut𝔻<Real> & cameraXZ, Real cameraY) {
    shader->uniform("view", view);
    shader->uniform("projection", projection);

    shader->uniform("cameraXZ.a", cameraXZ.a);
    shader->uniform("cameraXZ.b", cameraXZ.b);
    shader->uniform("cameraXZ.c", cameraXZ.c());
    shader->uniform("cameraXZ.d", cameraXZ.d());

    shader->uniform("cameraY", float(cameraY));
}

const double saveInterval = 1.0;

double globaltime = 0, saveTimer = 0;
void display(GLFWwindow * window, Config & config) {
    using namespace Game;

    auto dt = glfwGetTime() - globaltime;
    globaltime += dt; saveTimer += dt;

    using namespace std::complex_literals;

    auto dir = 0i;
    if (Keyboard::forward)  dir += +1i;
    if (Keyboard::backward) dir += -1i;
    if (Keyboard::left)     dir += +1;
    if (Keyboard::right)    dir += -1;

    if (dir != 0.0) dir /= std::abs(dir);

    auto n = std::polar(1.0, -player.camera().yaw);
    Gyrovector<Real> velocity(player.walkSpeed * dir * n);

    bool isTileChanged = move(player, velocity, dt);
    if (isTileChanged) pollNeighbours();

    for (auto it = map.pool.begin(); it != map.pool.end();) {
        auto tile = *it;

        if (!tile->ready()) { it++; continue; }

        if (tile->needRefresh())
            tile->refresh(Registry::node);

        if (Render::hmax < tile->cameraVerticalDistance())
            tile->unload();

        if (tile->needUnload() && !tile->dirty()) {
            delete tile; it = map.pool.erase(it);
        } else it++;
    }

    auto origin = player.camera().position.relativeXZ().inverse();

    if (Mouse::grabbed) {
        glfwGetCursorPos(window, &Mouse::xpos, &Mouse::ypos);
        glfwSetCursorPos(window, Window::width/2, Window::height/2);

        player.rotate(
            Mouse::speed * dt * (Window::width/2 - Mouse::xpos),
            Mouse::speed * dt * (Window::height/2 - Mouse::ypos),
            0.0f
        );
    }

    auto cameraY = player.camera().climb + player.eye;

    auto direction = player.camera().direction(), right = player.camera().right(), up = glm::cross(right, direction);
    auto eye = vec3(0.0f, -cameraY, 0.0f);

    view = glm::lookAt(vec3(0.0f), direction, up);
    view = glm::scale(view, vec3(1.0f, Render::standard->meter, 1.0f));
    view = glm::translate(view, eye);

    glClearColor(Render::background[0], Render::background[1], Render::background[2], Render::background[3]);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glBlendFunc(GL_ONE, GL_ZERO);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_MULTISAMPLE);

    const Real & vrd = config.camera.verticalRenderDistance;
    int Y₁ = std::floor(cameraY - vrd), Y₂ = std::floor(cameraY + vrd);

    faceShader->activate();
    uploadMVP(faceShader, origin, cameraY);

    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(1.0, 1.0);

    for (auto & tile : map.pool)
        if (tile->ready())
            tile->renderFaces(faceShader, Y₁, Y₂);

    glDisable(GL_POLYGON_OFFSET_FILL);

    edgeShader->activate();
    uploadMVP(edgeShader, origin, cameraY);

    for (auto & tile : map.pool)
        if (tile->ready())
            tile->renderEdges(edgeShader, Y₁, Y₂);

    if (auto value = pbo.read(Window::width/2 - 1, Window::height/2))
    { auto [zbuffer, action] = *value; click(origin, zbuffer, action); }

    glDisable(GL_DEPTH_TEST);
    glDisable(GL_MULTISAMPLE);

    dummyShader->activate();

    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glBlendFunc(GL_ONE_MINUS_DST_COLOR, GL_ZERO);
    aimVao.draw(GL_LINES);

    if (saveTimer >= saveInterval)
    { map.dump(); saveTimer = 0; }
}

void setupSheet() {
    using namespace Game;

    glEnable(GL_CULL_FACE);
}

void grabMouse(GLFWwindow * window) {
    using namespace Game;

    glfwSetCursorPos(window, Window::width/2, Window::height/2);
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
    if (glfwRawMouseMotionSupported())
        glfwSetInputMode(window, GLFW_RAW_MOUSE_MOTION, GLFW_TRUE);

    Mouse::grabbed = true;
}

void freeMouse(GLFWwindow * window) {
    using namespace Game;

    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
    glfwSetInputMode(window, GLFW_RAW_MOUSE_MOTION, GLFW_FALSE);

    Mouse::grabbed = false;
}

void mouseButtonCallback(GLFWwindow * window, int button, int action, int mods) {
    using namespace Game;

    if (Mouse::grabbed) {
        if (action == GLFW_PRESS) switch (button) {
            case GLFW_MOUSE_BUTTON_LEFT:  pbo.issue(Action::Remove); break;
            case GLFW_MOUSE_BUTTON_RIGHT: pbo.issue(Action::Place);  break;
        }
    } else if (Window::hovered && Window::focused) {
        if (action == GLFW_PRESS && button == GLFW_MOUSE_BUTTON_LEFT)
            grabMouse(window);
    }
}

void cursorEnterCallback(GLFWwindow * window, int entered) {
    using namespace Game;

    Window::hovered = entered;
    if (!Window::hovered) freeMouse(window);
}

void windowFocusCallback(GLFWwindow * window, int focused) {
    using namespace Game;

    Window::focused = focused;
    if (!Window::focused) freeMouse(window);
}

WorldTileData tileBuffer;

void copyTile() {
    using namespace Game;

    const auto & src = *player.tile();
    if (src.vxl() == nullptr) return;

    memcpy(&tileBuffer, src.vxl(), sizeof(WorldTileData));
}

void pasteTile() {
    using namespace Game;

    auto dest = player.tile()->vxl();
    if (dest == nullptr) return;

    memcpy(dest, &tileBuffer, sizeof(WorldTileData));
    player.tile()->requestRefresh();
}

void rotateTile() {
    using namespace Game;
    using namespace Fundamentals;

    auto src = player.tile()->vxl();
    if (src == nullptr) return;

    auto buf = new WorldTileData; memcpy(buf, src, sizeof(WorldTileData));

    static_assert(sizeTileX == sizeTileZ);

    for (int X = 0; X < sizeTileX; X++)
        for (int Y = 0; Y < sizeTileY; Y++)
            for (int Z = 0; Z < sizeTileZ; Z++)
                src->data[X][Y][Z] = buf->data[Z][Y][maxTileX - X];

    delete buf;

    player.tile()->requestRefresh();
}

const Real elevationRate = 3.0;

inline void pressLShift() { if (Game::player.flymode) Game::player.roc(-elevationRate); }
inline void releaseLShift() { if (Game::player.flymode) Game::player.roc(0); }

inline void pressSpace() {
    if (Game::player.flymode) Game::player.roc(elevationRate);
    else if (!Game::player.camera().flying) Game::player.jump();
}

inline void releaseSpace() { if (Game::player.flymode) Game::player.roc(0); }

inline void closeWindow(GLFWwindow * window) {
    glfwSetWindowShouldClose(window, GL_TRUE);
}

inline void hotbarSelect(size_t slot) {
    Game::activeSlot = slot;
    updateHotbar();
}

inline void returnToSpawn() {
    using namespace Game;

    player.teleport(Position(), 5);
    player.roc(0); pollNeighbours();
}

inline void toggleFlyMode() {
    using namespace Game;

    player.roc(0);
    player.flymode = !player.flymode;
}

inline void toggleNoclip() {
    using namespace Game;

    player.noclip = !player.noclip;
}

void keyboardCallback(GLFWwindow * window, int key, int scancode, int action, int mods) {
    using namespace Game;

    bool enabled = action != GLFW_RELEASE;

    if (action == GLFW_PRESS || action == GLFW_RELEASE) switch (key) {
        case GLFW_KEY_W:          Keyboard::forward  = enabled; break;
        case GLFW_KEY_S:          Keyboard::backward = enabled; break;
        case GLFW_KEY_A:          Keyboard::left     = enabled; break;
        case GLFW_KEY_D:          Keyboard::right    = enabled; break;
        case GLFW_KEY_SPACE:      Keyboard::space    = enabled; break;
        case GLFW_KEY_LEFT_SHIFT: Keyboard::lshift   = enabled; break;
    }

    if (action == GLFW_PRESS) switch (key) {
        case GLFW_KEY_ESCAPE:     closeWindow(window); break;
        case GLFW_KEY_O:          returnToSpawn();     break;
        case GLFW_KEY_H:          toggleNoclip();      break;
        case GLFW_KEY_K:          toggleFlyMode();     break;
        case GLFW_KEY_1:          hotbarSelect(0);     break;
        case GLFW_KEY_2:          hotbarSelect(1);     break;
        case GLFW_KEY_3:          hotbarSelect(2);     break;
        case GLFW_KEY_4:          hotbarSelect(3);     break;
        case GLFW_KEY_5:          hotbarSelect(4);     break;
        case GLFW_KEY_6:          hotbarSelect(5);     break;
        case GLFW_KEY_7:          hotbarSelect(6);     break;
        case GLFW_KEY_8:          hotbarSelect(7);     break;
        case GLFW_KEY_9:          hotbarSelect(8);     break;
        case GLFW_KEY_X:          rotateTile();        break;
        case GLFW_KEY_C:          copyTile();          break;
        case GLFW_KEY_V:          pasteTile();         break;
        case GLFW_KEY_BACKSLASH:  freeMouse(window);   break;
        case GLFW_KEY_SPACE:      pressSpace();        break;
        case GLFW_KEY_LEFT_SHIFT: pressLShift();       break;
    }

    if (action == GLFW_RELEASE) switch (key) {
        case GLFW_KEY_SPACE:      releaseSpace();  break;
        case GLFW_KEY_LEFT_SHIFT: releaseLShift(); break;
    }
}

void setupWindowSize(GLFWwindow * window, int width, int height) {
    using namespace Game;
    using namespace Render;

    Window::width  = width;
    Window::height = height;
    Window::aspect = Real(width) / Real(height);

    int frameBufferWidth, frameBufferHeight;
    glfwGetFramebufferSize(window, &frameBufferWidth, &frameBufferHeight);
    glViewport(0, 0, frameBufferWidth, frameBufferHeight);

    projection = glm::perspective(glm::radians(fov), Window::aspect, near, far);

    updateHotbar();
    drawAim(aimVao);
}

constexpr auto title = "Hypertest";
GLFWwindow * setupWindow(Config & config) {
    using namespace Game;

    glfwSetErrorCallback(errorCallback);
    if (!glfwInit()) throw std::runtime_error("glfwInit failure");

    if (config.window.msaa > 0) glfwWindowHint(GLFW_SAMPLES, config.window.msaa);

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    #ifdef __APPLE__
        glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    #endif

    Window::width  = config.window.width;
    Window::height = config.window.height;

    auto window = glfwCreateWindow(Window::width, Window::height, title, nullptr, nullptr);
    if (!window) throw std::runtime_error("glfwCreateWindow failure");

    grabMouse(window);

    glfwSetKeyCallback(window, keyboardCallback);
    glfwSetMouseButtonCallback(window, mouseButtonCallback);
    glfwSetWindowFocusCallback(window, windowFocusCallback);
    glfwSetCursorEnterCallback(window, cursorEnterCallback);
    glfwSetWindowSizeCallback(window, setupWindowSize);

    return window;
}

auto readText(const char * filepath) {
    size_t fsize; std::vector<char> retval;

    auto fd = std::fopen(filepath, "rb");
    if (fd == nullptr) goto error;

    if (std::fseek(fd, 0, SEEK_END) != 0) goto error;
    fsize = std::ftell(fd); if (fsize < 0) goto error;

    {
        std::rewind(fd);

        retval.resize(fsize + 1);
        std::fread(retval.data(), fsize, 1, fd);
        retval[fsize] = 0;

        std::fclose(fd);

        return retval;
    }

    error: fprintf(stderr, "Unable to read “%s”: %s\n", filepath, std::strerror(errno));
    return std::vector(1, '\0');
}

auto readModelShader(Model model) {
    switch (model) {
        case Poincaré:    return readText("shaders/Model/Poincare.glsl");
        case Klein:       return readText("shaders/Model/Klein.glsl");
        case Gans:        return readText("shaders/Model/Gans.glsl");
        case Equidistant: return readText("shaders/Model/Equidistant.glsl");
        case Lambert:     return readText("shaders/Model/Lambert.glsl");
        default:          return std::vector(1, '\0');
    }
}

void uploadShaders() {
    using namespace Game;

    auto ms = readModelShader(Render::standard->model);

    {
        delete faceShader;

        auto cs = readText("shaders/Voxel/Common.glsl");
        auto fs = readText("shaders/Voxel/FaceFragment.glsl");
        auto vs = readText("shaders/Voxel/FaceVertex.glsl");

        FragmentShader fragment(cs.data(), fs.data(), ms.data());
        VertexShader vertex(cs.data(), vs.data(), ms.data());

        faceShader = new FaceShader(fragment, vertex);
    }

    {
        delete edgeShader;

        auto cs = readText("shaders/Voxel/Common.glsl");
        auto fs = readText("shaders/Voxel/EdgeFragment.glsl");
        auto vs = readText("shaders/Voxel/EdgeVertex.glsl");

        FragmentShader fragment(cs.data(), fs.data(), ms.data());
        VertexShader vertex(cs.data(), vs.data(), ms.data());

        edgeShader = new EdgeShader(fragment, vertex);
    }

    {
        delete dummyShader;

        auto cs = readText("shaders/Dummy/Common.glsl");
        auto fs = readText("shaders/Dummy/Fragment.glsl");
        auto vs = readText("shaders/Dummy/Vertex.glsl");

        FragmentShader fragment(cs.data(), fs.data());
        VertexShader vertex(cs.data(), vs.data());

        dummyShader = new DummyShader(fragment, vertex);
    }
}

template<ShaderSpec Spec>
inline void uploadPrims(ShaderProgram<Spec> * shader, Config & config) {
    shader->uniform("vrd",         float(config.camera.verticalRenderDistance));
    shader->uniform("hrd",         float(config.camera.horizontalRenderDistance));
    shader->uniform("fog.enabled", config.fog.enabled);
    shader->uniform("fog.near",    config.fog.near);
    shader->uniform("fog.far",     config.fog.far);
    shader->uniform("fog.color",   config.fog.color);
}

void setupShaders(Config & config) {
    faceShader->activate(); uploadPrims(faceShader, config);
    edgeShader->activate(); uploadPrims(edgeShader, config);
}

void setupGL(GLFWwindow * window, Config & config) {
    using namespace Game;

    static Render::Standard standard(config.camera.model);
    Render::standard = &standard;

    glfwMakeContextCurrent(window);
    glewExperimental = GL_TRUE; glewInit();

    glEnable(GL_BLEND);

    uploadShaders();
    setupShaders(config);

    Render::fov  = config.camera.fov;
    Render::near = config.camera.near;
    Render::far  = config.camera.far;

    dummyShader->activate();

    aimVao.initialize();

    GUI::aimSize = config.gui.aimSize;
    setupWindowSize(window, Window::width, Window::height);

    pbo.initialize();
}

void buildFloor(WorldTile * C) {
    using namespace Fundamentals;

    /*for (int X = 0; X < sizeTileX; X++)
        for (int Z = 0; Z < sizeTileZ; Z++)
            C->set(X, 0, Z, {1});*/

    Node node = {1};

    for (int Y = 0; Y < sizeTileY; Y += 16) {
        for (int X = 0; X < sizeTileX; X++) for (int Z = 0; Z < sizeTileZ; Z++) {
            int H = abs(Z - sizeTileZ / 2) < 4 ? X : 0;

            C->set(X, Y + H, Z, node);
            if (Y + H < maxTileY) C->set(X, Y + H + 1, Z, node);
        }
    }

    for (int Y = 0; Y <= sizeTileY; Y++) {
        C->set(0,        Y, 0,        node);
        C->set(0,        Y, maxTileZ, node);
        C->set(maxTileX, Y, 0,        node);
        C->set(maxTileX, Y, maxTileZ, node);
    }
}

void setupGame(Config & config) {
    using namespace Tesselation;
    using namespace Game;

    map.mapgen = &buildFloor;

    Render::vmax = config.camera.verticalRenderDistance;
    Render::hmax = worldTileDiameter(config.camera.horizontalRenderDistance);

    map.poll(Tesselation::I, Tesselation::I);

    for (std::size_t k = 0; k < Tesselation::neighbours.size(); k++)
        map.poll(Tesselation::I, Tesselation::neighbours[k]);

    player.teleport(Position(), 4);
}

void cleanUp(GLFWwindow * window) {
    pbo.free();
    aimVao.free();

    delete dummyShader;
    delete faceShader;
    delete edgeShader;

    glfwDestroyWindow(window);
    glfwTerminate();
}

int main(int argc, char * argv[]) {
    using namespace Game;

    LuaJIT luajit;

    Config config(&luajit, "config.lua");

    auto window = setupWindow(config);
    setupGL(window, config);

    luajit.loadapi();

    for (int i = 1; i < argc; i++)
        luajit.go(argv[i]);

    map.connect(config.world);
    setupGame(config);
    setupSheet();

    updateHotbar();

    glfwSetTime(0);

    while (!glfwWindowShouldClose(window)) {
        display(window, config);
        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    map.disconnect();
    cleanUp(window);

    return 0;
}
