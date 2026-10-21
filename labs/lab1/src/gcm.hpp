#pragma once
#include <cstddef>
#include <cstdint>
#include "aes128.hpp"

class GCM {
public:
    explicit GCM(const uint8_t key[16], uint16_t poly = 0x11B) : aes_(key, poly) {
        uint8_t zero[16] = {0};
        aes_.encrypt_block(zero, h_);
    }

    void encrypt(const uint8_t iv[12],
                 const uint8_t* aad, size_t aad_len,
                 const uint8_t* plain, size_t len,
                 uint8_t* cipher, uint8_t tag[16]) const {
        uint8_t j0[16];
        make_j0(iv, j0);
        ctr_crypt(j0, plain, len, cipher);
        compute_tag(j0, aad, aad_len, cipher, len, tag);
    }

    bool decrypt(const uint8_t iv[12],
                 const uint8_t* aad, size_t aad_len,
                 const uint8_t* cipher, size_t len,
                 const uint8_t tag[16], uint8_t* plain) const {
        uint8_t j0[16];
        make_j0(iv, j0);
        uint8_t expected[16];
        compute_tag(j0, aad, aad_len, cipher, len, expected);
        bool valid = tags_equal(expected, tag);
        ctr_crypt(j0, cipher, len, plain);
        uint8_t mask = static_cast<uint8_t>(-static_cast<int>(valid));
        for (size_t i = 0; i < len; ++i) {
            plain[i] &= mask;
        }
        return valid;
    }

    static bool tags_equal(const uint8_t a[16], const uint8_t b[16]) {
        uint8_t diff = 0;
        for (int i = 0; i < 16; ++i) {
            diff |= static_cast<uint8_t>(a[i] ^ b[i]);
        }
        return diff == 0;
    }

    static void gf128_multiply(const uint8_t x[16], const uint8_t y[16], uint8_t out[16]) {
        uint8_t z[16] = {0};
        uint8_t v[16];
        for (int i = 0; i < 16; ++i) {
            v[i] = y[i];
        }
        for (int i = 0; i < 128; ++i) {
            uint8_t bit = static_cast<uint8_t>((x[i / 8] >> (7 - (i % 8))) & 1);
            uint8_t mask = static_cast<uint8_t>(-bit);
            for (int j = 0; j < 16; ++j) {
                z[j] ^= static_cast<uint8_t>(v[j] & mask);
            }
            uint8_t lsb = static_cast<uint8_t>(v[15] & 1);
            uint8_t reduce = static_cast<uint8_t>(0xE1 & static_cast<uint8_t>(-lsb));
            for (int j = 15; j > 0; --j) {
                v[j] = static_cast<uint8_t>((v[j] >> 1) | (v[j - 1] << 7));
            }
            v[0] = static_cast<uint8_t>((v[0] >> 1) ^ reduce);
        }
        for (int i = 0; i < 16; ++i) {
            out[i] = z[i];
        }
    }

private:
    AES128 aes_;
    uint8_t h_[16];

    static void make_j0(const uint8_t iv[12], uint8_t j0[16]) {
        for (int i = 0; i < 12; ++i) {
            j0[i] = iv[i];
        }
        j0[12] = 0;
        j0[13] = 0;
        j0[14] = 0;
        j0[15] = 1;
    }

    static void inc32(uint8_t block[16]) {
        for (int i = 15; i >= 12; --i) {
            block[i] = static_cast<uint8_t>(block[i] + 1);
            if (block[i] != 0) {
                break;
            }
        }
    }

    void ctr_crypt(const uint8_t j0[16], const uint8_t* in, size_t len, uint8_t* out) const {
        uint8_t counter[16];
        for (int i = 0; i < 16; ++i) {
            counter[i] = j0[i];
        }
        size_t offset = 0;
        while (offset < len) {
            inc32(counter);
            uint8_t gamma[16];
            aes_.encrypt_block(counter, gamma);
            size_t chunk = len - offset < 16 ? len - offset : 16;
            for (size_t i = 0; i < chunk; ++i) {
                out[offset + i] = static_cast<uint8_t>(in[offset + i] ^ gamma[i]);
            }
            offset += chunk;
        }
    }

    void absorb(uint8_t y[16], const uint8_t* data, size_t len) const {
        size_t offset = 0;
        while (offset < len) {
            size_t chunk = len - offset < 16 ? len - offset : 16;
            for (size_t i = 0; i < chunk; ++i) {
                y[i] ^= data[offset + i];
            }
            gf128_multiply(y, h_, y);
            offset += chunk;
        }
    }

    void compute_tag(const uint8_t j0[16],
                     const uint8_t* aad, size_t aad_len,
                     const uint8_t* cipher, size_t len,
                     uint8_t tag[16]) const {
        uint8_t y[16] = {0};
        absorb(y, aad, aad_len);
        absorb(y, cipher, len);
        uint8_t lengths[16];
        uint64_t aad_bits = static_cast<uint64_t>(aad_len) * 8;
        uint64_t len_bits = static_cast<uint64_t>(len) * 8;
        for (int i = 0; i < 8; ++i) {
            lengths[i] = static_cast<uint8_t>(aad_bits >> (56 - 8 * i));
            lengths[8 + i] = static_cast<uint8_t>(len_bits >> (56 - 8 * i));
        }
        absorb(y, lengths, 16);
        uint8_t mask[16];
        aes_.encrypt_block(j0, mask);
        for (int i = 0; i < 16; ++i) {
            tag[i] = static_cast<uint8_t>(y[i] ^ mask[i]);
        }
    }
};
