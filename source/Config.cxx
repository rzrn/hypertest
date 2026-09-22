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

#include <Hyper/Config.hxx>

Config::Config(LuaJIT * luajit, const char * filename) {
    using namespace Fundamentals;

    if (LuaTable config = luajit->require(filename)) {
        if (LuaString world_v = config.getitem("world"))
            world = world_v.decode();

        if (LuaTable window_v = config.getitem("window")) {
            if (LuaInteger width_v = window_v.getitem("width"))
                window.width = width_v.decode();

            if (LuaInteger height_v = window_v.getitem("height"))
                window.height = height_v.decode();

            if (LuaInteger msaa_v = window_v.getitem("msaa"))
                window.msaa = msaa_v.decode();
        }

        if (LuaTable camera_v = config.getitem("camera")) {
            if (LuaNumber vrd_v = camera_v.getitem("verticalRenderDistance"))
                camera.verticalRenderDistance = std::clamp(vrd_v.decode(), 8.0, sizeTileY / 2.0);

            if (LuaNumber hrd_v = camera_v.getitem("horizontalRenderDistance"))
                camera.horizontalRenderDistance = hrd_v.decode();

            if (LuaNumber fov_v = camera_v.getitem("fov"))
                camera.fov = fov_v.decode();

            if (LuaNumber near_v = camera_v.getitem("near"))
                camera.near = near_v.decode();

            if (LuaNumber far_v = camera_v.getitem("far"))
                camera.far = far_v.decode();

            if (LuaInteger model_v = camera_v.getitem("model"))
                camera.model = Model(model_v.decode());
        }

        if (LuaTable fog_v = config.getitem("fog")) {
            if (LuaBool enabled_v = fog_v.getitem("enabled"))
                fog.enabled = enabled_v.decode();

            if (LuaNumber near_v = fog_v.getitem("near"))
                fog.near = near_v.decode();

            if (LuaNumber far_v = fog_v.getitem("far"))
                fog.far = far_v.decode();

            if (LuaVec4 color_v = fog_v.getitem("color"))
                fog.color = color_v.decode();
        }

        if (LuaTable gui_v = config.getitem("gui")) {
            if (LuaNumber aim_v = gui_v.getitem("aimSize"))
                gui.aimSize = aim_v.decode();
        }
    }
}