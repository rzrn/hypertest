#pragma once

/*
    Copyright © 2026 rzrn

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

#include <utility>

template<typename T, size_t size> class DoubleBuffer {
private:
    T buff1[size], buff2[size];
    T (&readBuff)[size], (&writeBuff)[size];

public:
    DoubleBuffer() : buff1{}, buff2{}, readBuff(buff1), writeBuff(buff2) {}

    inline void flip() { std::swap(readBuff, writeBuff); }

    T & operator[](size_t index) { return writeBuff[index]; }
    const T & operator()(size_t index) const { return readBuff[index]; }
};