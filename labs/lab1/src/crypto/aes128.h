#pragma once
#include <cstdint>
#include <cstddef>
#include "gf256.h"

class AES128 {
public:
    static constexpr size_t kBlockSize = 16;
    static constexpr size_t kKeySize   = 16;
    static constexpr int    kRounds    = 10;
    static constexpr size_t kRoundKeys = 44;

    explicit AES128(const uint8_t key[kKeySize]) noexcept;

    void encrypt(const uint8_t in[kBlockSize], uint8_t out[kBlockSize]) const noexcept;
    void decrypt(const uint8_t in[kBlockSize], uint8_t out[kBlockSize]) const noexcept;

    const GF256& gf() const noexcept { return gf_; }

private:
    void sub_bytes(uint8_t s[4][4]) const noexcept;
    void inv_sub_bytes(uint8_t s[4][4]) const noexcept;
    static void shift_rows(uint8_t s[4][4]) noexcept;
    static void inv_shift_rows(uint8_t s[4][4]) noexcept;
    void mix_columns(uint8_t s[4][4]) const noexcept;
    void inv_mix_columns(uint8_t s[4][4]) const noexcept;
    void add_round_key(uint8_t s[4][4], const uint8_t rk[kBlockSize]) const noexcept;
    void expand_key(const uint8_t key[kKeySize]) noexcept;
    uint8_t sub_byte(uint8_t b) const noexcept { return gf_.sbox(b); }

    GF256   gf_;
    uint8_t round_keys_[11 * kBlockSize];
};