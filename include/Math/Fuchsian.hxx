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

#include <Math/Gaussian.hxx>
#include <Math/Moebius.hxx>

template<EuclideanDomain T>
struct Fuchsian {
    Gaussian<T> a, b, c, d;

    constexpr Fuchsian() {}
    constexpr Fuchsian(auto a, auto b, auto c, auto d) : a(a), b(b), c(c), d(d) {}

    constexpr Gaussian<T> det() const { return a * d - b * c; }

    constexpr void normalize() {
        auto σ(a);

        σ = Gaussian<T>::hcf(b, σ);
        σ = Gaussian<T>::hcf(c, σ);
        σ = Gaussian<T>::hcf(d, σ);

        a.divexact(σ);
        b.divexact(σ);
        c.divexact(σ);
        d.divexact(σ);
    }

    template<typename U> constexpr inline Möbius<U> field() const {
        // See `source/Geometry.cxx` for “s” meaning.
        constexpr auto s = Math::sqrt<U>(6.0);

        return { a.template field<U>(),     b.template field<U>() / s,
                 c.template field<U>() * s, d.template field<U>() };
    }

    constexpr inline auto inverse() const { return Fuchsian<T>(d, -b, -c, a); }

    constexpr inline auto origin() const {
        if (b.isZero()) return std::pair(b, Gaussian<T>(Math::one<T>));

        auto α(b), β(d); auto σ = Gaussian<T>::hcf(α, β);
        α.divexact(σ); β.divexact(σ); α.normalize(β); return std::pair(α, β);
    }

    constexpr auto operator==(const Fuchsian<T> & G) const
    { return a == G.a && b == G.b && c == G.b && d == G.d; }

    friend std::ostream & operator<< (std::ostream & stream, const Fuchsian<T> & G)
    { return stream << "(" << G.a << ", " << G.b << ", " << G.c << ", " << G.d << ")"; }
};

template<typename T> constexpr Fuchsian<T> operator*(const Fuchsian<T> & A, const Fuchsian<T> & B) {
    return {
        A.a * B.a + A.b * B.c,
        A.a * B.b + A.b * B.d,
        A.c * B.a + A.d * B.c,
        A.c * B.b + A.d * B.d
    };
}

template<typename T> void operator*=(Fuchsian<T> & A, const Fuchsian<T> & B) { A = A * B; A.normalize(); }