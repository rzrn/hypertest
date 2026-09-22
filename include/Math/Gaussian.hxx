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

#include <complex>

#include <Meta/Enumerable.hxx>
#include <Math/Euclidean.hxx>

template<EuclideanDomain T> struct Gaussian {
    T real, imag;

    constexpr Gaussian() : real(Math::zero<T>), imag(Math::zero<T>) {}
    constexpr Gaussian(const T & real) : real(real), imag(Math::zero<T>) {}
    constexpr Gaussian(const T & real, const T & imag) : real(real), imag(imag) {}

    constexpr inline auto operator-() const { return Gaussian(-real, -imag); }
    constexpr inline auto & operator+() const { return *this; }

    constexpr auto operator+(const Gaussian<T> & w) const
    { return Gaussian<T>(real + w.real, imag + w.imag); }

    constexpr void operator+=(const Gaussian<T> & w)
    { real += w.real; imag += w.imag; }

    constexpr auto operator-(const Gaussian<T> & w) const
    { return Gaussian<T>(real - w.real, imag - w.imag); }

    constexpr void operator-=(const Gaussian<T> & w)
    { real -= w.real; imag -= w.imag; }

    constexpr auto operator*(const Gaussian<T> & w) const
    { return Gaussian<T>(real * w.real - imag * w.imag, real * w.imag + imag * w.real); }

    constexpr void operator*=(const Gaussian<T> & w) {
        auto α(real), β(w.real);
        real = α * β - imag * w.imag;
        imag = α * w.imag + imag * β;
    }

    constexpr inline T norm() const { return real * real + imag * imag; }

    constexpr auto operator/(const T & k) const { return Gaussian<T>(real / k, imag / k); }
    constexpr void operator/=(const T & k) { real /= k; imag /= k; }

    constexpr auto operator/(const Gaussian<T> & w) const {
        auto N = w.norm();
        T x(real * w.real + imag * w.imag);
        T y(imag * w.real - real * w.imag);
        return Gaussian<T>(x / N, y / N);
    }

    constexpr inline void divexact(const T & k) {
        Math::divexact<T>(real, real, k);
        Math::divexact<T>(imag, imag, k);
    }

    constexpr void divexact(const Gaussian<T> & w) {
        auto N = w.norm();

        T x(real * w.real + imag * w.imag);
        T y(imag * w.real - real * w.imag);

        Math::divexact<T>(real, x, N);
        Math::divexact<T>(imag, y, N);
    }

    constexpr auto isZero() const { return Math::isZero(real) && Math::isZero(imag); }
    constexpr auto isUnit() const { return (Math::isUnit(real) && Math::isZero(imag))
                                        || (Math::isZero(real) && Math::isUnit(imag)); }

    constexpr void negate()  { real = -real; imag = -imag; }
    constexpr void mul2()    { Math::mul2(real); Math::mul2(imag); }
    constexpr void div2()    { Math::div2(real); Math::div2(imag); }
    constexpr void mulω()    { real -= imag; Math::mul2(imag); imag += real; }
    constexpr void mulnegi() { std::swap(real, imag); imag = -imag; }
    constexpr void muli()    { std::swap(real, imag); real = -real; }
    constexpr void divω()    { real += imag; Math::mul2(imag); imag -= real; div2(); }

    constexpr auto kind() const { return std::pair(Math::odd(real), Math::odd(imag)); }

    // https://www.researchgate.net/publication/269005874_A_Paper-and-Pencil_gcd_Algorithm_for_Gaussian_Integers
    // https://www.researchgate.net/publication/325472716_PERFORMANCE_OF_A_GCD_ALGO-RITHM_FOR_GAUSSIAN_INTEGERS
    constexpr static auto hcf(Gaussian<T> α, Gaussian<T> β) {
        Gaussian<T> δ(Math::one<T>);

        for (;;) {
            if (α == β || α == -β) return α * δ;
            if (α.isUnit() || β.isUnit()) return δ;
            if (α.isZero()) return β * δ;
            if (β.isZero()) return α * δ;

            switch (Ord(std::pair(α.kind(), β.kind()))) {
                case Digit(0, 0, 0, 0): α.div2(); β.div2(); δ.mul2(); break;
                case Digit(0, 0, 1, 1): α.div2(); β.divω(); δ.mulω(); break;
                case Digit(1, 1, 0, 0): α.divω(); β.div2(); δ.mulω(); break;
                case Digit(1, 1, 1, 1): α.divω(); β.divω(); δ.mulω(); break;

                case Digit(0, 1, 0, 0): case Digit(1, 0, 0, 0): β.div2(); break;
                case Digit(0, 0, 0, 1): case Digit(0, 0, 1, 0): α.div2(); break;
                case Digit(1, 0, 0, 1): case Digit(0, 1, 1, 0): β.muli(); break;
                case Digit(1, 1, 0, 1): case Digit(1, 1, 1, 0): α.divω(); break;
                case Digit(0, 1, 1, 1): case Digit(1, 0, 1, 1): β.divω(); break;

                case Digit(1, 0, 1, 0): case Digit(0, 1, 0, 1):
                α += β; β.mul2(); β -= α; α.div2(); β.div2(); break;
            }
        }
    }

    /* Given α = a + bi ∈ ℤ[i] and βᵢ ∈ ℤ[i] (1 ≤ i ≤ N), it multiplies all
       of them by u = ±1/±i so that Re(uα) > 0 and Im(uα) ≥ 0 if α ≠ 0 */
    template<std::same_as<Gaussian<T>>... Ts> constexpr void normalize(Ts &... ts) {
        using enum Math::Ordering;

        switch (Ord²(Math::compare(real, Math::zero<T>), Math::compare(imag, Math::zero<T>))) {
            /* (−a, −b) */ case Ord²(LT, LT): case Ord²(LT, EQ): negate();  (ts.negate(),  ...); break;
            /* (−b, +a) */ case Ord²(EQ, LT): case Ord²(GT, LT): muli();    (ts.muli(),    ...); break;
            /* (+b, −a) */ case Ord²(LT, GT): case Ord²(EQ, GT): mulnegi(); (ts.mulnegi(), ...); break;
            /* (+a, +b) */ case Ord²(GT, EQ): case Ord²(GT, GT): break;
        }
    }

    constexpr auto operator==(const Gaussian<T> & w) const
    { return Math::equal<T>(real, w.real) && Math::equal<T>(imag, w.imag); }

    constexpr auto operator!=(const Gaussian<T> & w) const
    { return Math::differ<T>(real, w.real) || Math::differ<T>(imag, w.imag); }

    template<EuclideanDomain U> auto transform() const { return Gaussian<U>(real, imag); }

    template<typename U> constexpr auto field() const
    { return std::complex<U>(Math::field<T, U>(real), Math::field<T, U>(imag)); }
};

template<EuclideanDomain T> using Gaussian² = std::pair<Gaussian<T>, Gaussian<T>>;