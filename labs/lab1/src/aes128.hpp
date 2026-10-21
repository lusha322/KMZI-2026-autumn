#pragma once
#include <cstdint>
#include "galois_field.hpp"

class AES128 {
public:
    AES128(const uint8_t key[16], uint16_t poly = 0x11B) : gf_(poly) {
        expand_key(key);
    }

    void encrypt_block(const uint8_t in[16], uint8_t out[16]) const {
        uint8_t state[16];
        for (int i = 0; i < 16; ++i) {
            state[i] = in[i];
        }
        add_round_key(state, 0);
        for (int round = 1; round < 10; ++round) {
            sub_bytes(state);
            shift_rows(state);
            mix_columns(state);
            add_round_key(state, round);
        }
        sub_bytes(state);
        shift_rows(state);
        add_round_key(state, 10);
        for (int i = 0; i < 16; ++i) {
            out[i] = state[i];
        }
    }

    const GaloisField& field() const { return gf_; }

private:
    GaloisField gf_;
    uint8_t round_keys_[176];

    void expand_key(const uint8_t key[16]) {
        for (int i = 0; i < 16; ++i) {
            round_keys_[i] = key[i];
        }
        uint8_t rcon = 0x01;
        for (int i = 4; i < 44; ++i) {
            uint8_t temp[4];
            for (int j = 0; j < 4; ++j) {
                temp[j] = round_keys_[4 * (i - 1) + j];
            }
            if (i % 4 == 0) {
                uint8_t first = temp[0];
                temp[0] = gf_.sub_byte(temp[1]);
                temp[1] = gf_.sub_byte(temp[2]);
                temp[2] = gf_.sub_byte(temp[3]);
                temp[3] = gf_.sub_byte(first);
                temp[0] = GaloisField::add(temp[0], rcon);
                rcon = gf_.xt(rcon);
            }
            for (int j = 0; j < 4; ++j) {
                round_keys_[4 * i + j] = GaloisField::add(round_keys_[4 * (i - 4) + j], temp[j]);
            }
        }
    }

    void add_round_key(uint8_t state[16], int round) const {
        for (int i = 0; i < 16; ++i) {
            state[i] = GaloisField::add(state[i], round_keys_[16 * round + i]);
        }
    }

    void sub_bytes(uint8_t state[16]) const {
        for (int i = 0; i < 16; ++i) {
            state[i] = gf_.sub_byte(state[i]);
        }
    }

    static void shift_rows(uint8_t state[16]) {
        uint8_t tmp[16];
        for (int c = 0; c < 4; ++c) {
            for (int r = 0; r < 4; ++r) {
                tmp[r + 4 * c] = state[r + 4 * ((c + r) % 4)];
            }
        }
        for (int i = 0; i < 16; ++i) {
            state[i] = tmp[i];
        }
    }

    void mix_columns(uint8_t state[16]) const {
        static const uint8_t matrix[16] = {
            0x02, 0x03, 0x01, 0x01,
            0x01, 0x02, 0x03, 0x01,
            0x01, 0x01, 0x02, 0x03,
            0x03, 0x01, 0x01, 0x02
        };
        uint8_t tmp[16];
        for (int c = 0; c < 4; ++c) {
            for (int r = 0; r < 4; ++r) {
                uint8_t acc = 0;
                for (int k = 0; k < 4; ++k) {
                    acc = GaloisField::add(acc, gf_.mul(matrix[4 * r + k], state[k + 4 * c]));
                }
                tmp[r + 4 * c] = acc;
            }
        }
        for (int i = 0; i < 16; ++i) {
            state[i] = tmp[i];
        }
    }
};
