#include "gf256.h"

GF256::GF256(uint16_t p_x) noexcept : p_x_(p_x) {}

uint8_t GF256::add(uint8_t a, uint8_t b) noexcept {
    return static_cast<uint8_t>(a ^ b);
}

uint8_t GF256::xtime(uint8_t b) const noexcept {
    const uint8_t hi      = static_cast<uint8_t>(b >> 7);
    const uint8_t mask    = static_cast<uint8_t>(0u - static_cast<unsigned>(hi));
    const uint8_t red     = modulus_byte();
    const uint8_t shifted = static_cast<uint8_t>(b << 1);
    return static_cast<uint8_t>(shifted ^ (mask & red));
}

uint8_t GF256::multiply(uint8_t a, uint8_t b) const noexcept {
    uint8_t result = 0;
    for (int i = 0; i < 8; ++i) {
        const uint8_t bit  = static_cast<uint8_t>(b & 1u);
        const uint8_t mask = static_cast<uint8_t>(0u - static_cast<unsigned>(bit));
        result ^= static_cast<uint8_t>(a & mask);
        a = xtime(a);
        b = static_cast<uint8_t>(b >> 1);
    }
    return result;
}

uint8_t GF256::inverse(uint8_t a) const noexcept {
    const uint8_t a2   = multiply(a, a);
    const uint8_t a4   = multiply(a2, a2);
    const uint8_t a8   = multiply(a4, a4);
    const uint8_t a16  = multiply(a8, a8);
    const uint8_t a32  = multiply(a16, a16);
    const uint8_t a64  = multiply(a32, a32);
    const uint8_t a128 = multiply(a64, a64);

    uint8_t res = a2;
    res = multiply(res, a4);
    res = multiply(res, a8);
    res = multiply(res, a16);
    res = multiply(res, a32);
    res = multiply(res, a64);
    res = multiply(res, a128);
    return res;
}

uint8_t GF256::affine(uint8_t b) noexcept {
    return static_cast<uint8_t>(
          b
        ^ rotl8(b, 1)
        ^ rotl8(b, 2)
        ^ rotl8(b, 3)
        ^ rotl8(b, 4)
        ^ 0x63);
}

uint8_t GF256::sbox(uint8_t b) const noexcept {
    return affine(inverse(b));
}

bool GF256::is_irreducible() const noexcept {
    if ((p_x_ & 0x100u) == 0) return false;

    auto poly_mod = [](uint16_t a, uint16_t b) -> uint16_t {
        if (b == 0) return a;
        auto deg = [](uint16_t v) -> int {
            int d = -1;
            while (v) { v >>= 1; ++d; }
            return d;
        };
        const int db = deg(b);
        while (a && deg(a) >= db) {
            const int da = deg(a);
            a ^= static_cast<uint16_t>(b << (da - db));
        }
        return a;
    };

    for (uint16_t d = 0x2; d <= 0x1F; ++d) {
        if (poly_mod(p_x_, d) == 0) return false;
    }
    return true;
}