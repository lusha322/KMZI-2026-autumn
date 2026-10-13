#include "gcm.h"
#include <cstring>

GCM::GCM(const AES128& cipher) noexcept : aes_(cipher) {
    std::memset(H_,  0, 16);
    std::memset(J0_, 0, 16);
}

void GCM::xor_block(uint8_t Y[16], const uint8_t X[16]) noexcept {
    for (int i = 0; i < 16; ++i)
        Y[i] ^= X[i];
}

void GCM::init() noexcept {
    const uint8_t zero[16] = {0};
    aes_.encrypt(zero, H_);
}

void GCM::gf128_mul(const uint8_t X[16], const uint8_t H[16],
                    uint8_t Z[16]) noexcept {
    uint8_t V[16];
    std::memcpy(V, H, 16);
    std::memset(Z, 0, 16);

    for (int i = 0; i < 128; ++i) {
        const uint8_t bit  = static_cast<uint8_t>((X[i >> 3] >> (7 - (i & 7))) & 1u);
        const uint8_t mask = static_cast<uint8_t>(0u - static_cast<unsigned>(bit));

        for (int j = 0; j < 16; ++j)
            Z[j] ^= static_cast<uint8_t>(V[j] & mask);

        const uint8_t lsb   = static_cast<uint8_t>(V[15] & 1u);
        const uint8_t rmask = static_cast<uint8_t>(0u - static_cast<unsigned>(lsb));

        uint8_t carry = 0;
        for (int j = 0; j < 16; ++j) {
            const uint8_t nxt = static_cast<uint8_t>(V[j] & 1u);
            V[j] = static_cast<uint8_t>((V[j] >> 1) | (carry << 7));
            carry = nxt;
        }

        V[0] ^= static_cast<uint8_t>(rmask & 0xE1);
    }
}

void GCM::ghash(const uint8_t* aad, size_t aad_len,
                const uint8_t* ct,  size_t ct_len,
                uint8_t out_y[16]) const noexcept {
    uint8_t Y[16] = {0};
    uint8_t block[16];

    size_t i = 0;
    while (i + 16 <= aad_len) {
        xor_block(Y, aad + i);
        gf128_mul(Y, H_, Y);
        i += 16;
    }
    if (i < aad_len) {
        std::memset(block, 0, 16);
        std::memcpy(block, aad + i, aad_len - i);
        xor_block(Y, block);
        gf128_mul(Y, H_, Y);
    }

    i = 0;
    while (i + 16 <= ct_len) {
        xor_block(Y, ct + i);
        gf128_mul(Y, H_, Y);
        i += 16;
    }
    if (i < ct_len) {
        std::memset(block, 0, 16);
        std::memcpy(block, ct + i, ct_len - i);
        xor_block(Y, block);
        gf128_mul(Y, H_, Y);
    }

    const uint64_t aad_bits = static_cast<uint64_t>(aad_len) * 8u;
    const uint64_t ct_bits  = static_cast<uint64_t>(ct_len)  * 8u;

    for (int k = 0; k < 8; ++k) {
        block[k]     = static_cast<uint8_t>(aad_bits >> (56 - 8*k));
        block[8 + k] = static_cast<uint8_t>(ct_bits  >> (56 - 8*k));
    }
    xor_block(Y, block);
    gf128_mul(Y, H_, Y);

    std::memcpy(out_y, Y, 16);
}

void GCM::ctr_block(const uint8_t CB[16], const uint8_t in[16], uint8_t out[16]) const noexcept {
    uint8_t ks[16];
    aes_.encrypt(CB, ks);
    for (int i = 0; i < 16; ++i)
        out[i] = static_cast<uint8_t>(in[i] ^ ks[i]);
}

void GCM::encrypt(const uint8_t* aad, size_t aad_len,
                  const uint8_t* pt,  size_t pt_len,
                  const uint8_t nonce[kNonceLen],
                  uint8_t* out_ct, Tag& out_tag) const noexcept {
    std::memset(J0_, 0, 16);
    std::memcpy(J0_, nonce, kNonceLen);
    J0_[15] = 0x01;

    uint8_t CB[16];
    std::memcpy(CB, J0_, 16);

    size_t i = 0;
    while (i + 16 <= pt_len) {
        uint32_t ctr = (static_cast<uint32_t>(CB[12]) << 24)
                     | (static_cast<uint32_t>(CB[13]) << 16)
                     | (static_cast<uint32_t>(CB[14]) << 8)
                     |  static_cast<uint32_t>(CB[15]);
        ++ctr;
        CB[12] = static_cast<uint8_t>(ctr >> 24);
        CB[13] = static_cast<uint8_t>(ctr >> 16);
        CB[14] = static_cast<uint8_t>(ctr >> 8);
        CB[15] = static_cast<uint8_t>(ctr);

        ctr_block(CB, pt + i, out_ct + i);
        i += 16;
    }
    if (i < pt_len) {
        uint32_t ctr = (static_cast<uint32_t>(CB[12]) << 24)
                     | (static_cast<uint32_t>(CB[13]) << 16)
                     | (static_cast<uint32_t>(CB[14]) << 8)
                     |  static_cast<uint32_t>(CB[15]);
        ++ctr;
        CB[12] = static_cast<uint8_t>(ctr >> 24);
        CB[13] = static_cast<uint8_t>(ctr >> 16);
        CB[14] = static_cast<uint8_t>(ctr >> 8);
        CB[15] = static_cast<uint8_t>(ctr);

        uint8_t ks[16];
        aes_.encrypt(CB, ks);
        const size_t rem = pt_len - i;
        for (size_t k = 0; k < rem; ++k)
            out_ct[i + k] = static_cast<uint8_t>(pt[i + k] ^ ks[k]);
    }

    uint8_t S[16];
    ghash(aad, aad_len, out_ct, pt_len, S);

    uint8_t ek_j0[16];
    aes_.encrypt(J0_, ek_j0);

    for (int k = 0; k < 16; ++k)
        out_tag[k] = static_cast<uint8_t>(S[k] ^ ek_j0[k]);
}

bool GCM::decrypt(const uint8_t* aad, size_t aad_len,
                  const uint8_t* ct,  size_t ct_len,
                  const uint8_t nonce[kNonceLen],
                  const Tag& tag,
                  uint8_t* out_pt) const noexcept {
    std::memset(J0_, 0, 16);
    std::memcpy(J0_, nonce, kNonceLen);
    J0_[15] = 0x01;

    uint8_t CB[16];
    std::memcpy(CB, J0_, 16);

    size_t i = 0;
    while (i + 16 <= ct_len) {
        uint32_t ctr = (static_cast<uint32_t>(CB[12]) << 24)
                     | (static_cast<uint32_t>(CB[13]) << 16)
                     | (static_cast<uint32_t>(CB[14]) << 8)
                     |  static_cast<uint32_t>(CB[15]);
        ++ctr;
        CB[12] = static_cast<uint8_t>(ctr >> 24);
        CB[13] = static_cast<uint8_t>(ctr >> 16);
        CB[14] = static_cast<uint8_t>(ctr >> 8);
        CB[15] = static_cast<uint8_t>(ctr);

        ctr_block(CB, ct + i, out_pt + i);
        i += 16;
    }
    if (i < ct_len) {
        uint32_t ctr = (static_cast<uint32_t>(CB[12]) << 24)
                     | (static_cast<uint32_t>(CB[13]) << 16)
                     | (static_cast<uint32_t>(CB[14]) << 8)
                     |  static_cast<uint32_t>(CB[15]);
        ++ctr;
        CB[12] = static_cast<uint8_t>(ctr >> 24);
        CB[13] = static_cast<uint8_t>(ctr >> 16);
        CB[14] = static_cast<uint8_t>(ctr >> 8);
        CB[15] = static_cast<uint8_t>(ctr);

        uint8_t ks[16];
        aes_.encrypt(CB, ks);
        const size_t rem = ct_len - i;
        for (size_t k = 0; k < rem; ++k)
            out_pt[i + k] = static_cast<uint8_t>(ct[i + k] ^ ks[k]);
    }

    uint8_t S[16];
    ghash(aad, aad_len, ct, ct_len, S);

    uint8_t ek_j0[16];
    aes_.encrypt(J0_, ek_j0);

    Tag computed;
    for (int k = 0; k < 16; ++k)
        computed[k] = static_cast<uint8_t>(S[k] ^ ek_j0[k]);

    return ct_compare(computed, tag);
}

bool GCM::ct_compare(const Tag& a, const Tag& b) noexcept {
    uint8_t acc = 0;
    for (size_t i = 0; i < kTagBytes; ++i) {
        acc = static_cast<uint8_t>(acc | (a[i] ^ b[i]));
    }
    return acc == 0;
}