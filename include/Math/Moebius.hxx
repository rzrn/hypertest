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

#include <Math/Gyrovector.hxx>

// M = (az + b) / (cz + d)
template<typename T> struct Möbius {
    std::complex<T> a, b, c, d;

    constexpr Möbius() : a(1), b(0), c(0), d(1) {}
    constexpr Möbius(auto a, auto b, auto c, auto d) : a(a), b(b), c(c), d(d) {}

    constexpr std::complex<T> det() const { return a * d - b * c; }
    constexpr std::complex<T> tr()  const { return a + d; }
    constexpr std::complex<T> rot() const { return det() / Math::sqr(d); }

    constexpr inline Möbius<T> div(std::complex<T> k) const { return {a / k, b / k, c / k, d / k}; }
    constexpr inline void normalize() { auto σ = sqrt(det()); a /= σ; b /= σ; c /= σ; d /= σ; }

    constexpr Gyrovector<T> apply(const Gyrovector<T> & w) const
    { return (a * w.val + b) / (c * w.val + d); }

    constexpr inline auto origin() const { return Gyrovector<T>(b / d); }

    constexpr inline Möbius<T> inverse() const { return Möbius<T>(d, -b, -c, a); }

    constexpr static inline Möbius<T> identity() { return Möbius<T>(1, 0, 0, 1); }

    constexpr static Möbius<T> translate(const Gyrovector<T> & N)
    { return Möbius<T>(1, N.val, Math::conjc(N.val), 1); }

    friend std::ostream & operator<< (std::ostream & stream, const Möbius<T> & M)
    { return stream << "(" << M.a << ", " << M.b << ", " << M.c << ", " << M.d << ")"; }
};

template<typename T> Möbius<T> operator*(const Möbius<T> & A, const Möbius<T> & B) {
    return {
        A.a * B.a + A.b * B.c,
        A.a * B.b + A.b * B.d,
        A.c * B.a + A.d * B.c,
        A.c * B.b + A.d * B.d
    };
}
