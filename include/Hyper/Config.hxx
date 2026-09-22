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

#include <glm/vec4.hpp>
#include <GL/glew.h>

#include <Hyper/Fundamentals.hxx>
#include <Lua.hxx>

struct Config {
    std::string world = "world.sqlite3";

    struct {
        bool enabled = false;
        GLfloat near = 1.0, far = 5.0;

        glm::vec4 color = {1.0f, 1.0f, 1.0f, 1.0f};
    } fog;

    struct {
        int width  = 800;
        int height = 600;
        int msaa   = 0;
    } window;

    struct {
        Real horizontalRenderDistance = 10.0;
        Real verticalRenderDistance = 64.0;

        Real fov = 80.0, near = 1e-3, far = 150.0;
        Model model = {Gans};
    } camera;

    struct {
        GLfloat aimSize = 15.0;
    } gui;

    Config(LuaJIT *, const char *);
};