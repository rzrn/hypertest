#pragma once

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

#include <vector>

#include <Hyper/Fundamentals.hxx>

class Texture {
    vec4 _lu, _ru, _rd, _ld;

public:
    inline Texture() {}

    inline Texture(const vec4 & rgba) :
    _lu(rgba), _ru(rgba), _rd(rgba), _ld(rgba) {}

    inline Texture(const vec4 & lu, const vec4 & ru, const vec4 & rd, const vec4 & ld) :
    _lu(lu), _ru(ru), _rd(rd), _ld(ld) {}

    inline constexpr auto & lu() { return _lu; }
    inline constexpr auto & ru() { return _ru; }
    inline constexpr auto & rd() { return _rd; }
    inline constexpr auto & ld() { return _ld; }
};

class Sheet {
private:
    std::vector<Texture> _textures;

public:

    template<typename... Ts> inline size_t attach(const Ts &... ts) {
        size_t index = _textures.size();
        _textures.emplace_back(ts...);

        return index;
    }

    inline auto size() const { return _textures.size(); }
    inline auto get(size_t idx) { return _textures[idx]; }
};