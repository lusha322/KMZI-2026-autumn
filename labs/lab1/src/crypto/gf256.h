#pragma once
#include <cstdint>

class GF256 {
public:
    explicit GF256(uint16_t p_x = 0x11B) noexcept;

    static uint8_t add(uint8_t a, uint8_t b) noexcept;
    uint8_t multiply(uint8_t a, uint8_t b) const noexcept;
    uint8_t inverse(uint8_t a) const noexcept;
    uint8_t xtime(uint8_t b) const noexcept;

    static uint8_t affine(uint8_t b) noexcept;
    uint8_t sbox(uint8_t b) const noexcept;

    uint16_t modulus() const noexcept { return p_x_; }
    uint8_t  modulus_byte() const noexcept {
        return static_cast<uint8_t>(p_x_ & 0xFF);
    }
    bool is_irreducible() const noexcept;

private:
    uint16_t p_x_;
};

inline uint8_t rotl8(uint8_t x, unsigned n) noexcept {
    const unsigned m = n & 7u;
    return static_cast<uint8_t>((x << m) | (x >> ((8u - m) & 7u)));
}