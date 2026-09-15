#pragma once

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