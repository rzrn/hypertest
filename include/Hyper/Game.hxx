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

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <Hyper/Sheet.hxx>
#include <Hyper/Physics.hxx>
#include <Hyper/Geometry.hxx>
#include <Hyper/Fundamentals.hxx>

enum class Action { Remove, Place };

namespace Game {
    namespace Registry {
        extern NodeRegistry node;
        extern Sheet sheet;
    }

    extern WorldMap map;
    extern Entity player;
    extern Camera camera;

    constexpr size_t hotbarSize = 9;
    extern NodeId hotbar[hotbarSize];
    extern size_t activeSlot;

    namespace Render {
        struct Standard {
            const Real meter;
            const Model model;

            Standard(const Model m) : meter(m.length(Tesselation::meter)), model(m) {}
        };

        extern unsigned int vmax; extern Real hmax;

        extern Real fov, near, far;
        extern Standard * standard;

        extern vec4 background;
    }

    namespace GUI {
        extern Real aimSize;
        extern Real aimTargetSize;
    }

    namespace Keyboard {
        extern bool forward;
        extern bool backward;
        extern bool left;
        extern bool right;
        extern bool space;
        extern bool lshift;
    }

    namespace Mouse {
        extern Real xpos, ypos, speed;
        extern bool grabbed;
    }

    namespace Window {
        extern bool hovered, focused;
        extern int width, height;
        extern Real aspect;
    }
}