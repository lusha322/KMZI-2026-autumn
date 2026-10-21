#pragma once
#include <cstdint>

class GaloisField {
public:
    explicit GaloisField(uint16_t poly) : poly_(poly) {}

    uint16_t poly() const { return poly_; }

    static uint8_t add(uint8_t a, uint8_t b) {
        return static_cast<uint8_t>(a ^ b);
    }

    static uint8_t xtime(uint8_t b, uint16_t p_x) {
        uint8_t mask = static_cast<uint8_t>(-(b >> 7));
        uint8_t reduce = static_cast<uint8_t>(p_x & 0xFF);
        return static_cast<uint8_t>((b << 1) ^ (mask & reduce));
    }

    static uint8_t multiply(uint8_t a, uint8_t b, uint16_t p_x) {
        uint8_t result = 0;
        for (int i = 0; i < 8; ++i) {
            uint8_t mask = static_cast<uint8_t>(-(b & 1));
            result ^= static_cast<uint8_t>(a & mask);
            a = xtime(a, p_x);
            b >>= 1;
        }
        return result;
    }

    static uint8_t inverse(uint8_t a, uint16_t p_x) {
        uint8_t power = multiply(a, a, p_x);
        uint8_t result = power;
        for (int i = 0; i < 6; ++i) {
            power = multiply(power, power, p_x);
            result = multiply(result, power, p_x);
        }
        return result;
    }

    static uint8_t rotl8(uint8_t b, int n) {
        return static_cast<uint8_t>((b << n) | (b >> (8 - n)));
    }

    static uint8_t affine(uint8_t b) {
        return static_cast<uint8_t>(b ^ rotl8(b, 1) ^ rotl8(b, 2) ^ rotl8(b, 3) ^ rotl8(b, 4) ^ 0x63);
    }

    uint8_t mul(uint8_t a, uint8_t b) const { return multiply(a, b, poly_); }
    uint8_t xt(uint8_t b) const { return xtime(b, poly_); }
    uint8_t inv(uint8_t a) const { return inverse(a, poly_); }
    uint8_t sub_byte(uint8_t a) const { return affine(inv(a)); }

    void generate_sbox(uint8_t table[256]) const {
        for (int i = 0; i < 256; ++i) {
            table[i] = sub_byte(static_cast<uint8_t>(i));
        }
    }

private:
    uint16_t poly_;
};
