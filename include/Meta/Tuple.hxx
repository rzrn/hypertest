#pragma once

/*
    Copyright © 2023–2024 rzrn

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

// It is here simply because C++ standard does not require any specific order
// of fields in std::tuple, however we *need* such in VBO.
template<typename...> struct Tuple;

template<> struct Tuple<> {
    constexpr Tuple() {}
};

// Explicit overload is needed to avoid EBCO nuances.
template<typename T> struct Tuple<T> {
    T fin;

    constexpr Tuple() {}
    constexpr Tuple(const T & t) : fin(t) {}
};

template<typename T, typename... Ts> struct Tuple<T, Ts...> {
    T first; Tuple<Ts...> second;

    constexpr Tuple() {}
    constexpr Tuple(const T & t, Ts... ts) : first(t), second(Tuple<Ts...>(ts...)) {}
};

template<class... T> Tuple(T...) -> Tuple<T...>;