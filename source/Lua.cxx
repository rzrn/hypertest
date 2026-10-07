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

#include <libgen.h>
#include <string.h>
#include <stdio.h>

#include <Hyper/Game.hxx>
#include <Lua.hxx>

LuaJIT::LuaJIT() {
    vm = luaL_newstate();
    luaL_openlibs(vm);
}

LuaJIT::~LuaJIT() {
    lua_close(vm);
}

const char * error(const int errcode) {
    switch (errcode) {
        case LUA_ERRSYNTAX: return "LUA_ERRSYNTAX";
        case LUA_ERRMEM:    return "LUA_ERRSYNTAX";
        case LUA_ERRFILE:   return "LUA_ERRFILE";
        case LUA_ERRRUN:    return "LUA_ERRRUN";
        case LUA_ERRERR:    return "LUA_ERRERR";
        default:            return "LUA_UNKNOWN";
    }
}

inline LuaRef warning(lua_State * vm, int errcode) {
    fprintf(stderr, "[%s] Lua: %s\n", error(errcode), lua_tostring(vm, -1));
    lua_pop(vm, 1);

    lua_pushnil(vm);
    return LuaRef(vm);
}

LuaRef LuaJIT::require(const char * filename) {
    if (auto error = luaL_loadfile(vm, filename))
        return warning(vm, error);

    if (auto error = lua_pcall(vm, 0, 1, 0))
        return warning(vm, error);

    return LuaRef(vm);
}

static const char proxyname[] = "core";

LuaRef LuaJIT::go(const char * filename) {
    char * buff;

    lua_getglobal(vm, proxyname);

    lua_pushstring(vm, filename);
    lua_setfield(vm, -2, "filename");

    buff = strdup(filename); lua_pushstring(vm, dirname(buff));
    lua_setfield(vm, -2, "dirname"); free(buff);

    buff = strdup(filename); lua_pushstring(vm, basename(buff));
    lua_setfield(vm, -2, "basename"); free(buff);

    lua_pop(vm, 1);

    return require(filename);
}

namespace API {
    static int override(lua_State * vm) {
        using namespace Game;
        luaL_checktype(vm, 1, LUA_TTABLE);

        lua_getfield(vm, 1, "eye"); camera.eye = luaL_checknumber(vm, -1); lua_pop(vm, 1);
        lua_getfield(vm, 1, "height"); player.height = luaL_checknumber(vm, -1); lua_pop(vm, 1);
        lua_getfield(vm, 1, "gravity"); player.gravity = luaL_checknumber(vm, -1); lua_pop(vm, 1);
        lua_getfield(vm, 1, "jump"); player.jumpHeight = luaL_checknumber(vm, -1); lua_pop(vm, 1);
        lua_getfield(vm, 1, "walk"); player.walkingSpeed = luaL_checknumber(vm, -1) * Tesselation::meter; lua_pop(vm, 1);

        return 0;
    }
}

static const luaL_Reg externs[] = {
    {"override", API::override},
    {NULL,       NULL}
};

void LuaJIT::loadapi() {
    luaL_register(vm, proxyname, externs);

    lua_pop(vm, 1);
}
