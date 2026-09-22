#pragma once

/*
    Copyright © 2022–2024 rzrn

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

template<typename T> struct Enumerator {};

template<typename T> concept Enumerable =
requires(T t) {
    { Enumerator<T>::cardinal   } -> std::convertible_to<size_t>;
    { Enumerator<T>::ordinal(t) } -> std::convertible_to<size_t>;
};

template<Enumerable T> constexpr size_t Card = Enumerator<T>::cardinal;

template<Enumerable T> inline constexpr size_t Ord(const T t)
{ return Enumerator<T>::ordinal(t); }

template<Enumerable T> inline constexpr size_t Ord²(const T t₁, const T t₂)
{ return Ord(std::pair(t₁, t₂)); }

template<Enumerable T> inline constexpr size_t Ord⁴(const T t₁, const T t₂, const T t₃, const T t₄)
{ return Ord(std::pair(std::pair(t₁, t₂), std::pair(t₃, t₄))); }

template<> struct Enumerator<bool> {
    static constexpr size_t cardinal = 2;
    static constexpr size_t ordinal(const bool P) { return P; }
};

template<Enumerable T, Enumerable U> struct Enumerator<std::pair<T, U>> {
    static constexpr size_t cardinal = Card<T> * Card<U>;

    static constexpr size_t ordinal(const std::pair<T, U> & w)
    { return Ord<T>(w.first) + Card<T> * Ord<U>(w.second); }
};

template<typename T> inline constexpr size_t Digit(T a, T b, T c, T d)
{ return Ord⁴(bool(a), bool(b), bool(c), bool(d)); }
