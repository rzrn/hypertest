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

#include <stdexcept>
#include <bit>

#include <Meta/Basic.hxx>
#include <Hyper/HVXL.hxx>

template<typename T> constexpr T fnv1aBasis = 0xCBF29CE484222325;
template<typename T> constexpr T fnv1aPrime = 0x100000001B3;

template<typename T> inline void fnv1a(uint8_t value, T & hash) {
    // https://en.wikipedia.org/wiki/Fowler–Noll–Vo_hash_function

    hash = hash ^ value;
    hash = hash * fnv1aPrime<T>;
}

template<typename T> inline void fnv1a(const std::vector<uint8_t> & data, T & hash) {
    for (const auto value : data)
        fnv1a(value, hash);
}

template<typename T> inline void fmix64(T & hash) {
    /* MurmurHash3 finalizer: https://github.com/aappleby/smhasher/blob/master/src/MurmurHash3.cpp
       Tests suggest that it does not degrade the distribution but does improve the avalanche behaviour. */

    hash = hash ^ (hash >> 33);
    hash = hash * 0xFF51AFD7ED558CCD;

    hash = hash ^ (hash >> 33);
    hash = hash * 0xC4CEB9FE1A85EC53;

    hash = hash ^ (hash >> 33);
}

using Vector8u = std::vector<uint8_t>;

template<typename T> inline T readLE(FILE * fp) {
    if constexpr(std::same_as<T, uint8_t>) {
        uint8_t buf[1];

        if (std::fread(buf, sizeof(buf), 1, fp) <= 0)
            return 0;

        return buf[0];
    } else if constexpr(std::same_as<T, uint16_t>) {
        uint8_t buf[2];

        if (std::fread(buf, sizeof(buf), 1, fp) <= 0)
            return 0;

        return (uint16_t) buf[0] << 0
             | (uint16_t) buf[1] << 8;
    } else if constexpr(std::same_as<T, uint32_t>) {
        uint8_t buf[4];

        if (std::fread(buf, sizeof(buf), 1, fp) <= 0)
            return 0;

        return (uint32_t) buf[0] << 0
             | (uint32_t) buf[1] << 8
             | (uint32_t) buf[2] << 16
             | (uint32_t) buf[3] << 24;
    } else if constexpr(std::same_as<T, uint64_t>) {
        uint8_t buf[8];

        if (std::fread(buf, sizeof(buf), 1, fp) <= 0)
            return 0;

        return (uint64_t) buf[0] << 0
             | (uint64_t) buf[1] << 8
             | (uint64_t) buf[2] << 16
             | (uint64_t) buf[3] << 24
             | (uint64_t) buf[4] << 32
             | (uint64_t) buf[5] << 40
             | (uint64_t) buf[6] << 48
             | (uint64_t) buf[7] << 56;
    } else if constexpr(std::same_as<T, float>) {
        // I don’t know any platform that has C++20 and OpenGL but does not use IEEE 754.
        return std::bit_cast<float>(readLE<uint32_t>(fp));
    } else {
        static_assert(falsehood<T>, "unexpected type");
    }
}

template<typename T> size_t writeLE(FILE * fp, const T value) {
    if constexpr(std::same_as<T, uint8_t>) {
        uint8_t buf[] = {
            static_cast<uint8_t>(value)
        };

        return std::fwrite(buf, sizeof(buf), 1, fp);
    } else if constexpr(std::same_as<T, uint16_t>) {
        uint8_t buf[] = {
            static_cast<uint8_t>(value >> 0),
            static_cast<uint8_t>(value >> 8)
        };

        return std::fwrite(buf, sizeof(buf), 1, fp);
    } else if constexpr(std::same_as<T, uint32_t>) {
        uint8_t buf[] = {
            static_cast<uint8_t>(value >> 0),
            static_cast<uint8_t>(value >> 8),
            static_cast<uint8_t>(value >> 16),
            static_cast<uint8_t>(value >> 24)
        };

        return std::fwrite(buf, sizeof(buf), 1, fp);
    } else if constexpr(std::same_as<T, uint64_t>) {
        uint8_t buf[] = {
            static_cast<uint8_t>(value >> 0),
            static_cast<uint8_t>(value >> 8),
            static_cast<uint8_t>(value >> 16),
            static_cast<uint8_t>(value >> 24),
            static_cast<uint8_t>(value >> 32),
            static_cast<uint8_t>(value >> 40),
            static_cast<uint8_t>(value >> 48),
            static_cast<uint8_t>(value >> 56)
        };

        return std::fwrite(buf, sizeof(buf), 1, fp);
    } else if constexpr(std::same_as<T, float>) {
        return writeLE<uint32_t>(fp, std::bit_cast<uint32_t>(value));
    } else if constexpr(std::same_as<T, Vector8u>) {
        size_t retval = 0;

        retval += writeLE<uint64_t>(fp, value.size());
        retval += std::fwrite(value.data(), value.size(), 1, fp);

        return retval;
    } else {
        static_assert(falsehood<T>, "unexpected type");
    }
}

HVXL::HVXL() : mapname{"Untitled World"}, seed(0) {
    for (size_t i = 1; i <= 0xFF; i++) {
        // 215 web colors + 39 shades of gray

        if (i < 216) {
            size_t r = i % 6;
            size_t g = (i / 6) % 6;
            size_t b = (i / 36) % 6;

            color[i].r() = static_cast<float>(r) / 5.0F;
            color[i].g() = static_cast<float>(g) / 5.0F;
            color[i].b() = static_cast<float>(b) / 5.0F;
        } else {
            /* 215 is white, so we do not include white here.
               0 is reserved for air, so we do include black here. */
            float shade = static_cast<float>(i - 216) / static_cast<float>(256 - 216);

            color[i].r() = shade;
            color[i].g() = shade;
            color[i].b() = shade;
        }
    }
}

HVXL::~HVXL() { close(); }

void HVXL::open(const std::string & filepath) {
    // This creates the file if it does not exist.
    auto fpchk = std::fopen(filepath.c_str(), "ab");
    if (fpchk == nullptr) throw std::runtime_error("`fopen(..., \"ab\")` failed");

    std::fclose(fpchk);

    fp = std::fopen(filepath.c_str(), "rb+");
    if (fp == nullptr) throw std::runtime_error("`fopen(..., \"rb+\")` failed");

    if (offsetEOF() > 0)
        readMapHeader();
    else
        writeMapHeader();
}

void HVXL::close() {
    if (fp != nullptr)
        std::fclose(fp);
}

template<typename T> inline T HVXL::readAt(off_t offset) {
    if (fseeko(fp, offset, SEEK_SET) < 0)
        return 0;

    return readLE<T>(fp);
}

template<typename T> inline size_t HVXL::writeAt(off_t offset, const T value) {
    if (fseeko(fp, offset, SEEK_SET) < 0)
        return 0;

    return writeLE<T>(fp, value);
}

void HVXL::readMapHeader() {
    std::rewind(fp);

    if (readLE<uint64_t>(fp) != fmtver)
        throw std::runtime_error("`.hvxl` wrong version");

    if (std::fread(mapname, sizeof(mapname), 1, fp) <= 0)
        throw std::runtime_error("`.hvxl` too short");

    seed = readLE<uint32_t>(fp);

    for (size_t i = 1; i <= 0xFF; i++) {
        color[i].r() = readLE<float>(fp);
        color[i].g() = readLE<float>(fp);
        color[i].b() = readLE<float>(fp);
    }

    offsetHeader = ftello(fp);
}

size_t HVXL::writeMapHeader() {
    std::rewind(fp);

    size_t retval = 0;

    retval += writeLE<uint64_t>(fp, fmtver);
    retval += std::fwrite(mapname, sizeof(mapname), 1, fp);
    retval += writeLE<uint32_t>(fp, seed);

    for (size_t i = 1; i <= 0xFF; i++) {
        retval += writeLE<float>(fp, color[i].r());
        retval += writeLE<float>(fp, color[i].g());
        retval += writeLE<float>(fp, color[i].b());
    }

    offsetHeader = ftello(fp);
    retval += writeHashRecord();

    std::fflush(fp);

    return retval;
}

size_t HVXL::writeHashRecord() {
    if (fseeko(fp, sizeHashRecord, SEEK_CUR) < 0)
        return 0;

    return sizeHashRecord + writeLE<uint64_t>(fp, 0);
}

struct TileHeader {
    uint64_t hash; uint8_t signMask; Vector8u bufFstReal, bufFstImag, bufSndReal, bufSndImag;

    TileHeader(const Gaussian²<Integer> & absoluteOriginXZ) : signMask(0) {
        hash = fnv1aBasis<uint64_t>;

        bufFstReal.resize(Math::size(absoluteOriginXZ.first.real));
        bufFstImag.resize(Math::size(absoluteOriginXZ.first.imag));
        bufSndReal.resize(Math::size(absoluteOriginXZ.second.real));
        bufSndImag.resize(Math::size(absoluteOriginXZ.second.imag));

        Math::toBytes(absoluteOriginXZ.first.real,  bufFstReal.data());
        Math::toBytes(absoluteOriginXZ.first.imag,  bufFstImag.data());
        Math::toBytes(absoluteOriginXZ.second.real, bufSndReal.data());
        Math::toBytes(absoluteOriginXZ.second.imag, bufSndImag.data());

        if (Math::isNeg(absoluteOriginXZ.second.real))
            signMask |= static_cast<uint8_t>(1 << 0);

        if (Math::isNeg(absoluteOriginXZ.second.imag))
            signMask |= static_cast<uint8_t>(1 << 1);

        fnv1a<uint64_t>(signMask, hash);

        fnv1a<uint64_t>(bufFstReal, hash);
        fnv1a<uint64_t>(bufFstImag, hash);
        fnv1a<uint64_t>(bufSndReal, hash);
        fnv1a<uint64_t>(bufSndImag, hash);

        fmix64<uint64_t>(hash);
    }

    size_t write(FILE * fp) const {
        size_t retval = 0;

        retval += writeLE<uint64_t>(fp, hash);

        retval += writeLE<uint8_t>(fp, signMask);

        retval += writeLE<Vector8u>(fp, bufFstReal);
        retval += writeLE<Vector8u>(fp, bufFstImag);
        retval += writeLE<Vector8u>(fp, bufSndReal);
        retval += writeLE<Vector8u>(fp, bufSndImag);

        return retval;
    }

    static inline bool isEqualVector(FILE * fp, const Vector8u & vec) {
        if (readLE<uint64_t>(fp) != vec.size())
            return false;

        for (const auto value : vec)
            if (readLE<uint8_t>(fp) != value)
                return false;

        return true;
    }

    bool isEqualAt(FILE * fp, off_t offset) const {
        if (fseeko(fp, offset, SEEK_SET) < 0)
            return false;

        if (readLE<uint64_t>(fp) != hash)
            return false;

        if (readLE<uint8_t>(fp) != signMask)
            return false;

        if (!isEqualVector(fp, bufFstReal))
            return false;

        if (!isEqualVector(fp, bufFstImag))
            return false;

        if (!isEqualVector(fp, bufSndReal))
            return false;

        if (!isEqualVector(fp, bufSndImag))
            return false;

        return true;
    }
};

void HVXL::readTileData(uint8_t * buf, RGB3f & fog) {
    std::fread(buf, sizeTile, 1, fp);

    fog.r() = readLE<float>(fp);
    fog.g() = readLE<float>(fp);
    fog.b() = readLE<float>(fp);
}

size_t HVXL::writeTileData(const uint8_t * buf, RGB3f fog) {
    size_t retval = 0;

    if (buf == nullptr) {
        if (fseeko(fp, sizeTile, SEEK_CUR) < 0)
            return 0;

        retval += sizeTile;
    } else {
        retval += std::fwrite(buf, sizeTile, 1, fp);
    }

    retval += writeLE<float>(fp, fog.r());
    retval += writeLE<float>(fp, fog.g());
    retval += writeLE<float>(fp, fog.b());

    std::fflush(fp);

    return retval;
}

bool HVXL::findTileOffset(const Gaussian²<Integer> & absoluteOriginXZ, off_t & offsetOut, uint64_t & hash) {
    TileHeader thdr(absoluteOriginXZ); hash = thdr.hash;

    off_t offsetRecord = offsetHeader;

    for (;;) {
        static constexpr uint64_t hashMask = (1 << sizeHashInBits) - 1;
        off_t offsetTileRef = offsetRecord + (thdr.hash & hashMask) * sizeof(uint64_t);

        off_t offsetTile = readAt<uint64_t>(offsetTileRef);

        if (offsetTile == 0) {
            static RGB3f whiteRGB(1.0f);

            offsetTile = offsetEOF(); thdr.write(fp);

            offsetOut = ftello(fp); writeTileData(nullptr, whiteRGB);

            writeAt<uint64_t>(offsetTileRef, offsetTile);
            std::fflush(fp);

            return true;
        }

        if (thdr.isEqualAt(fp, offsetTile)) {
            offsetOut = ftello(fp); return false;
        }

        off_t offsetNextRecordRef = offsetRecord + sizeHashRecord;
        off_t offsetNextRecord = readAt<uint64_t>(offsetNextRecordRef);

        if (offsetNextRecord == 0) {
            offsetNextRecord = offsetEOF(); writeHashRecord();

            writeAt<uint64_t>(offsetNextRecordRef, offsetNextRecord);
        }

        offsetRecord = offsetNextRecord;
    }
}