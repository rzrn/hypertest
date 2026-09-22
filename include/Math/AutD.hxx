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

#include <Math/Moebius.hxx>

/*
    Matrix of Möbius transformation degrades drastically after successive multiplications
    due to the machine float precision limitations; in long term perspective
    it no longer maps 𝔻 into 𝔻, what results in weird graphical glitches.

    That’s why we use separate class to represent fundamental domain.
    Apart from accuracy issues, this class holds only two (instead of four)
    complex components, hence computations are done twice as fast.

    In general, automorphism from 𝔻 to 𝔻 can be represented as f(z) = exp(iφ)(z + z₀)/(zz₀* + 1) (*).
    https://en.wikipedia.org/wiki/M%C3%B6bius_transformation#Subgroups_of_the_M%C3%B6bius_group

    Multiplying this expression by exp(−iφ/2), we get:
        f(z) = exp(iφ)(z + z₀)/(zz₀* + 1)
             = (exp(iφ)z + exp(iφ)z₀)/(zb* + 1)
             = (exp(iφ/2)z + exp(iφ/2)z₀)/(zexp(−iφ/2)z₀* + exp(−iφ/2))
             = (exp(iφ/2)z + exp(iφ/2)z₀)/(z[exp(iφ/2)z₀]* + [exp(iφ/2)]*)

    Denote a = exp(iφ/2) and b = exp(iφ/2)z₀. Then:
        f(z) = (az + b)/(zb* + a*)

    It’s easy to see that product of two matrices of this type is again such matrix:
          [a₁, b₁; b₁*, a₁*] × [a₂, b₂; b₂*, a₂*]
        = [a₁a₂ + b₁b₂*, a₁b₂ + b₁a₂*; a₂b₁* + a₁*b₂*, b₂b₁* + a₁*a₂*]
        = [a₁a₂ + b₁b₂*, a₁b₂ + b₁a₂*; (a₁b₂ + b₁a₂*)*, (a₁a₂ + b₁b₂*)*]
        = [a₃, b₃; b₃*, a₃*],
    where a₃ = a₁a₂ + b₁b₂* and b₃ = a₁b₂ + b₁a₂*.

    These matrices map 𝔻 to 𝔻 iff |b| < |a|, while case |b| ≥ |a| corresponds to the jumping out of the Poincaré disk.
    (https://math.stackexchange.com/questions/3698149/m%c3%b6bius-transformation-that-maps-the-unit-circle-to-itself?rq=1)

    Of course, such a jump in our case is simply meaningless.

    It’s also important that if we divide such a matrix by a real number, we will again get a matrix
    of this form (since z = z* ↔ z ∈ ℝ). Its determinant is always real (aa* − bb* = |a|² − |b|²),
    so we can safely (without going out of class) normalize it as any other Möbius transformation.

    It also means that, conversely, we can (almost) always rearrange such matrix to the form (*).
    Indeed, dividing the matrix by |a| (if a ≠ 0, of course), we will normalize the first (a′) and fourth (a′*) components.
    Then:
        M = [a, b; b*, a*]
          ~ [a², ab; ab*, aa*]
          = [a², a²(b/a); a(ab/a)*, aa*]
          = [a², a²(b/a); aa*(b/a)*, aa*]
          = [a², a²(b/a); (b/a)*, 1] (because |a| = 1, then aa* = |a|² = 1)
          = [exp(iφ), exp(iφ)z₀; z₀*, 1],
    where exp(iφ) = a², z₀ = b/a and M₁ ~ M₂ ↔ ∃ z ∈ ℂ: M₁ = zM₂.
*/

// M = (az + b) / (zb* + a*)
template<typename T> struct Aut𝔻 {
    std::complex<T> a, b;

    constexpr inline std::complex<T> c() { return std::conj(b); }
    constexpr inline std::complex<T> d() { return std::conj(a); }

    constexpr Aut𝔻() : a(1), b(0) {}
    constexpr Aut𝔻(const Gyrovector<T> z₀) : a(1), b(z₀.val) {}
    constexpr Aut𝔻(const std::complex<T> z₀) : a(1), b(z₀) {}
    constexpr Aut𝔻(const std::complex<T> a, const std::complex<T> b) : a(a), b(b) {}

    constexpr Aut𝔻(const T φ, const std::complex<T> w)
    { auto ω = std::polar(1.0, φ / 2.0); a = ω; b = ω * w; }

    constexpr std::complex<T> det() const { return std::norm(a) - std::norm(b); }

    constexpr Gyrovector<T> apply(const Gyrovector<T> & w) const
    { return (a * w.val + b) / (std::conj(b) * w.val + std::conj(a)); }

    constexpr inline auto origin() const { return Gyrovector<T>(b / std::conj(a)); }

    constexpr inline Aut𝔻<T> inverse() const { return Aut𝔻<T>(std::conj(a), -b); }

    constexpr inline void normalize() { auto σ = sqrt(det()); a /= σ; b /= σ; }

    explicit operator Möbius<T>() const { return Möbius<T>(a, b, std::conj(b), std::conj(a)); }

    friend std::ostream & operator<< (std::ostream & stream, const Aut𝔻<T> & M)
    { return stream << "(" << M.a << ", " << M.b << ", " << std::conj(M.b) << ", " << std::conj(M.a) << ")"; }
};

template<typename T> Aut𝔻<T> operator*(const Aut𝔻<T> & A, const Aut𝔻<T> & B) {
    return {
        A.a * B.a + A.b * std::conj(B.b),
        A.a * B.b + A.b * std::conj(B.a)
    };
}