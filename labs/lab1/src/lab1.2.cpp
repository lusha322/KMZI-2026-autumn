#include "aes128.h"
#include <cstring>

namespace {
    constexpr uint8_t kRcon[11] = {
        0x00,
        0x01, 0x02, 0x04, 0x08, 0x10,
        0x20, 0x40, 0x80, 0x1B, 0x36
    };
}

AES128::AES128(const uint8_t key[kKeySize]) noexcept : gf_(0x11B) {
    expand_key(key);
}

void AES128::expand_key(const uint8_t key[kKeySize]) noexcept {
    uint8_t w[44][4];

    for (int i = 0; i < 4; ++i) {
        w[i][0] = key[4*i + 0];
        w[i][1] = key[4*i + 1];
        w[i][2] = key[4*i + 2];
        w[i][3] = key[4*i + 3];
    }

    for (int i = 4; i < 44; ++i) {
        uint8_t temp[4] = { w[i-1][0], w[i-1][1], w[i-1][2], w[i-1][3] };

        if (i % 4 == 0) {
            const uint8_t t0 = temp[0];
            temp[0] = temp[1];
            temp[1] = temp[2];
            temp[2] = temp[3];
            temp[3] = t0;

            temp[0] = sub_byte(temp[0]);
            temp[1] = sub_byte(temp[1]);
            temp[2] = sub_byte(temp[2]);
            temp[3] = sub_byte(temp[3]);

            temp[0] ^= kRcon[i / 4];
        }

        for (int j = 0; j < 4; ++j) {
            w[i][j] = static_cast<uint8_t>(w[i - 4][j] ^ temp[j]);
        }
    }

    for (int r = 0; r < 11; ++r) {
        for (int c = 0; c < 4; ++c) {
            round_keys_[r*16 + 4*c + 0] = w[4*r + c][0];
            round_keys_[r*16 + 4*c + 1] = w[4*r + c][1];
            round_keys_[r*16 + 4*c + 2] = w[4*r + c][2];
            round_keys_[r*16 + 4*c + 3] = w[4*r + c][3];
        }
    }
}

void AES128::sub_bytes(uint8_t s[4][4]) const noexcept {
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            s[i][j] = gf_.sbox(s[i][j]);
}

void AES128::inv_sub_bytes(uint8_t s[4][4]) const noexcept {
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            const uint8_t target = s[i][j];
            for (int b = 0; b < 256; ++b) {
                if (gf_.sbox(static_cast<uint8_t>(b)) == target) {
                    s[i][j] = static_cast<uint8_t>(b);
                    break;
                }
            }
        }
    }
}

void AES128::shift_rows(uint8_t s[4][4]) noexcept {
    uint8_t t;

    t = s[1][0]; s[1][0]=s[1][1]; s[1][1]=s[1][2]; s[1][2]=s[1][3]; s[1][3]=t;

    t = s[2][0]; s[2][0]=s[2][2]; s[2][2]=t;
    t = s[2][1]; s[2][1]=s[2][3]; s[2][3]=t;

    t = s[3][3]; s[3][3]=s[3][2]; s[3][2]=s[3][1]; s[3][1]=s[3][0]; s[3][0]=t;
}

void AES128::inv_shift_rows(uint8_t s[4][4]) noexcept {
    uint8_t t;

    t = s[1][3]; s[1][3]=s[1][2]; s[1][2]=s[1][1]; s[1][1]=s[1][0]; s[1][0]=t;

    t = s[2][0]; s[2][0]=s[2][2]; s[2][2]=t;
    t = s[2][1]; s[2][1]=s[2][3]; s[2][3]=t;

    t = s[3][0]; s[3][0]=s[3][1]; s[3][1]=s[3][2]; s[3][2]=s[3][3]; s[3][3]=t;
}

void AES128::mix_columns(uint8_t s[4][4]) const noexcept {
    for (int c = 0; c < 4; ++c) {
        const uint8_t a0 = s[0][c], a1 = s[1][c], a2 = s[2][c], a3 = s[3][c];
        const uint8_t t  = static_cast<uint8_t>(a0 ^ a1 ^ a2 ^ a3);

        s[0][c] = static_cast<uint8_t>(a0 ^ t ^ gf_.xtime(static_cast<uint8_t>(a0 ^ a1)));
        s[1][c] = static_cast<uint8_t>(a1 ^ t ^ gf_.xtime(static_cast<uint8_t>(a1 ^ a2)));
        s[2][c] = static_cast<uint8_t>(a2 ^ t ^ gf_.xtime(static_cast<uint8_t>(a2 ^ a3)));
        s[3][c] = static_cast<uint8_t>(a3 ^ t ^ gf_.xtime(static_cast<uint8_t>(a3 ^ a0)));
    }
}

void AES128::inv_mix_columns(uint8_t s[4][4]) const noexcept {
    for (int c = 0; c < 4; ++c) {
        const uint8_t a0 = s[0][c], a1 = s[1][c], a2 = s[2][c], a3 = s[3][c];
        s[0][c] = static_cast<uint8_t>(
            gf_.multiply(a0, 0x0E) ^ gf_.multiply(a1, 0x0B) ^
            gf_.multiply(a2, 0x0D) ^ gf_.multiply(a3, 0x09));
        s[1][c] = static_cast<uint8_t>(
            gf_.multiply(a0, 0x09) ^ gf_.multiply(a1, 0x0E) ^
            gf_.multiply(a2, 0x0B) ^ gf_.multiply(a3, 0x0D));
        s[2][c] = static_cast<uint8_t>(
            gf_.multiply(a0, 0x0D) ^ gf_.multiply(a1, 0x09) ^
            gf_.multiply(a2, 0x0E) ^ gf_.multiply(a3, 0x0B));
        s[3][c] = static_cast<uint8_t>(
            gf_.multiply(a0, 0x0B) ^ gf_.multiply(a1, 0x0D) ^
            gf_.multiply(a2, 0x09) ^ gf_.multiply(a3, 0x0E));
    }
}

void AES128::add_round_key(uint8_t s[4][4], const uint8_t rk[kBlockSize]) const noexcept {
    for (int c = 0; c < 4; ++c)
        for (int r = 0; r < 4; ++r)
            s[r][c] ^= rk[4*c + r];
}

void AES128::encrypt(const uint8_t in[kBlockSize], uint8_t out[kBlockSize]) const noexcept {
    uint8_t s[4][4];
    for (int c = 0; c < 4; ++c)
        for (int r = 0; r < 4; ++r)
            s[r][c] = in[4*c + r];

    add_round_key(s, round_keys_ + 0);

    for (int round = 1; round <= kRounds - 1; ++round) {
        sub_bytes(s);
        shift_rows(s);
        mix_columns(s);
        add_round_key(s, round_keys_ + round * kBlockSize);
    }

    sub_bytes(s);
    shift_rows(s);
    add_round_key(s, round_keys_ + kRounds * kBlockSize);

    for (int c = 0; c < 4; ++c)
        for (int r = 0; r < 4; ++r)
            out[4*c + r] = s[r][c];
}

void AES128::decrypt(const uint8_t in[kBlockSize], uint8_t out[kBlockSize]) const noexcept {
    uint8_t s[4][4];
    for (int c = 0; c < 4; ++c)
        for (int r = 0; r < 4; ++r)
            s[r][c] = in[4*c + r];

    add_round_key(s, round_keys_ + kRounds * kBlockSize);

    for (int round = kRounds - 1; round >= 1; --round) {
        inv_shift_rows(s);
        inv_sub_bytes(s);
        add_round_key(s, round_keys_ + round * kBlockSize);
        inv_mix_columns(s);
    }

    inv_shift_rows(s);
    inv_sub_bytes(s);
    add_round_key(s, round_keys_ + 0);

    for (int c = 0; c < 4; ++c)
        for (int r = 0; r < 4; ++r)
            out[4*c + r] = s[r][c];
}