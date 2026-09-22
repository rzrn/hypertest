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

#include <Hyper/Game.hxx>

namespace Game {

using namespace Fundamentals;

namespace Registry {
    Sheet sheet;
    NodeRegistry node;
}

WorldMap map;
Entity player(&map);
Camera camera;

NodeId hotbar[hotbarSize] = {0};
size_t activeSlot = 0;

namespace Render {
    unsigned int vmax; Real hmax;

    Real fov, near, far;
    Standard * standard;

    vec4 background = {1.0f, 1.0f, 1.0f, 1.0f};
}

namespace GUI {
    Real aimSize, aimTargetSize;
}

namespace Keyboard {
    bool forward  = false;
    bool backward = false;
    bool left     = false;
    bool right    = false;
    bool space    = false;
    bool lshift   = false;
}

namespace Mouse {
    bool grabbed = false;
    Real xpos, ypos, speed = 0.7;
}

namespace Window {
    bool hovered = true, focused = true;

    int width, height;
    Real aspect;
}

}