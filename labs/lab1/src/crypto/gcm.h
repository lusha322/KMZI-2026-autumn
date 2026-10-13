#pragma once
#include <cstdint>
#include <cstddef>
#include <array>
#include "aes128.h"

class GCM {
public:
    static constexpr size_t kBlockSize = 16;
    static constexpr size_t kTagBytes  = 16;
    static constexpr size_t kNonceLen  = 12;

    using Tag = std::array<uint8_t, kTagBytes>;

    explicit GCM(const AES128& cipher) noexcept;

    void init() noexcept;

    void encrypt(const uint8_t* aad, size_t aad_len,
                 const uint8_t* pt,  size_t pt_len,
                 const uint8_t nonce[kNonceLen],
                 uint8_t* out_ct, Tag& out_tag) const noexcept;

    bool decrypt(const uint8_t* aad, size_t aad_len,
                 const uint8_t* ct,  size_t ct_len,
                 const uint8_t nonce[kNonceLen],
                 const Tag& tag,
                 uint8_t* out_pt) const noexcept;

    static bool ct_compare(const Tag& a, const Tag& b) noexcept;

private:
    void ghash(const uint8_t* aad, size_t aad_len,
               const uint8_t* ct,  size_t ct_len,
               uint8_t out_y[16]) const noexcept;

    static void gf128_mul(const uint8_t X[16], const uint8_t H[16],
                          uint8_t Z[16]) noexcept;

    static void xor_block(uint8_t Y[16], const uint8_t X[16]) noexcept;

    void ctr_block(const uint8_t CB[16], const uint8_t in[16], uint8_t out[16]) const noexcept;

    const AES128& aes_;
    uint8_t H_[16];
    uint8_t J0_[16];
};