#include "protocols.h"
#include "ct_secure.h"
#include <gmp.h>    // только арифметика больших чисел по модулю q
#include <cstdio>

// ---------------------------------------------------------------- SHA-256

static const uint32_t K256[64] = {
    0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
    0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
    0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
    0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
    0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
    0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
    0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
    0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2};

static inline uint32_t rotr(uint32_t x, int n) { return (x >> n) | (x << (32 - n)); }

static void sha256Block(uint32_t h[8], const uint8_t blk[64]) {
    uint32_t w[64];
    for (int i = 0; i < 16; ++i)
        w[i] = ((uint32_t)blk[4*i] << 24) | ((uint32_t)blk[4*i+1] << 16) | ((uint32_t)blk[4*i+2] << 8) | blk[4*i+3];
    for (int i = 16; i < 64; ++i) {
        uint32_t s0 = rotr(w[i-15], 7) ^ rotr(w[i-15], 18) ^ (w[i-15] >> 3);
        uint32_t s1 = rotr(w[i-2], 17) ^ rotr(w[i-2], 19) ^ (w[i-2] >> 10);
        w[i] = w[i-16] + s0 + w[i-7] + s1;
    }
    uint32_t a = h[0], b = h[1], c = h[2], d = h[3], e = h[4], f = h[5], g = h[6], hh = h[7];
    for (int i = 0; i < 64; ++i) {
        uint32_t S1 = rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25);
        uint32_t ch = (e & f) ^ (~e & g);
        uint32_t t1 = hh + S1 + ch + K256[i] + w[i];
        uint32_t S0 = rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22);
        uint32_t mj = (a & b) ^ (a & c) ^ (b & c);
        uint32_t t2 = S0 + mj;
        hh = g; g = f; f = e; e = d + t1; d = c; c = b; b = a; a = t1 + t2;
    }
    h[0]+=a; h[1]+=b; h[2]+=c; h[3]+=d; h[4]+=e; h[5]+=f; h[6]+=g; h[7]+=hh;
}

void sha256(const uint8_t* data, size_t len, uint8_t out[32]) {
    uint32_t h[8] = {0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19};
    size_t i = 0;
    for (; i + 64 <= len; i += 64) sha256Block(h, data + i);
    uint8_t tail[128] = {0};
    size_t rem = len - i;
    for (size_t j = 0; j < rem; ++j) tail[j] = data[i + j];
    tail[rem] = 0x80;
    size_t total = (rem + 9 <= 64) ? 64 : 128;
    uint64_t bits = (uint64_t)len * 8;
    for (int j = 0; j < 8; ++j) tail[total - 1 - j] = (uint8_t)(bits >> (8 * j));
    for (size_t j = 0; j < total; j += 64) sha256Block(h, tail + j);
    for (int j = 0; j < 8; ++j) {
        out[4*j] = (uint8_t)(h[j] >> 24); out[4*j+1] = (uint8_t)(h[j] >> 16);
        out[4*j+2] = (uint8_t)(h[j] >> 8); out[4*j+3] = (uint8_t)h[j];
    }
}

// ---------------------------------------------------------------- вспомогательные операции mod q (GMP)

static void toMpz(mpz_t z, const uint8_t b[32]) { mpz_import(z, 32, 1, 1, 1, 0, b); }

static void fromMpz(const mpz_t z, uint8_t out[32]) {
    for (int i = 0; i < 32; ++i) out[i] = 0;
    size_t cnt = 0;
    uint8_t tmp[32];
    mpz_export(tmp, &cnt, 1, 1, 1, 0, z);
    for (size_t i = 0; i < cnt && i < 32; ++i) out[32 - cnt + i] = tmp[i];
}

static void loadQ(mpz_t q) { toMpz(q, curve().q); }

bool scalarInRange(const uint8_t x[32]) {
    mpz_t z, q;
    mpz_init(z); mpz_init(q);
    toMpz(z, x); loadQ(q);
    bool ok = mpz_sgn(z) > 0 && mpz_cmp(z, q) < 0;
    mpz_clear(z); mpz_clear(q);
    return ok;
}

bool randomScalar(uint8_t out[32]) {
    FILE* f = std::fopen("/dev/urandom", "rb");
    if (!f) return false;
    bool ok = false;
    for (int attempt = 0; attempt < 1000 && !ok; ++attempt) {
        if (std::fread(out, 1, 32, f) != 32) break;
        ok = scalarInRange(out);
    }
    std::fclose(f);
    return ok;
}

bool keyGenPublic(const uint8_t d[32], Point* Q) {
    if (!scalarInRange(d)) return false;
    *Q = SecureMultiply(d, pointBase());
    return true;
}

bool ecdhSharedSecret(const uint8_t d[32], const Point& peer, uint8_t shared[32]) {
    if (!scalarInRange(d)) return false;
    if (peer.inf || !pointOnCurve(peer)) return false;
    Point K = SecureMultiply(d, peer);
    if (K.inf) return false;
    feToBytes(K.x, shared);
    return true;
}

// e mod q; возвращает false, если e = 0 mod q
static bool hashToScalar(const uint8_t hash[32], mpz_t e, const mpz_t q) {
    toMpz(e, hash);
    mpz_mod(e, e, q);
    return mpz_sgn(e) != 0;
}

bool stbSign(const uint8_t d[32], const uint8_t hash[32], const uint8_t k[32],
             uint8_t r[32], uint8_t s[32]) {
    if (!scalarInRange(d) || !scalarInRange(k)) return false;
    mpz_t q, e, dz, kz, rz, sz;
    mpz_init(q); mpz_init(e); mpz_init(dz); mpz_init(kz); mpz_init(rz); mpz_init(sz);
    loadQ(q);
    bool ok = hashToScalar(hash, e, q);
    if (ok) {
        Point C = SecureMultiply(k, pointBase());          // C = k*G, секретный k -> лестница Монтгомери
        uint8_t xb[32];
        feToBytes(C.x, xb);
        toMpz(rz, xb);
        mpz_mod(rz, rz, q);                                // r = x_c mod q
        toMpz(dz, d); toMpz(kz, k);
        mpz_mul(sz, rz, dz);                               // s = r*d + k*e mod q
        mpz_addmul(sz, kz, e);
        mpz_mod(sz, sz, q);
        ok = !C.inf && mpz_sgn(rz) != 0 && mpz_sgn(sz) != 0;
        if (ok) { fromMpz(rz, r); fromMpz(sz, s); }
    }
    mpz_clear(q); mpz_clear(e); mpz_clear(dz); mpz_clear(kz); mpz_clear(rz); mpz_clear(sz);
    return ok;
}

bool stbVerify(const Point& Q, const uint8_t hash[32], const uint8_t r[32], const uint8_t s[32]) {
    if (!scalarInRange(r) || !scalarInRange(s)) return false;
    if (Q.inf || !pointOnCurve(Q)) return false;
    mpz_t q, e, v, rz, sz, z1, z2;
    mpz_init(q); mpz_init(e); mpz_init(v); mpz_init(rz); mpz_init(sz); mpz_init(z1); mpz_init(z2);
    loadQ(q);
    bool result = false;
    if (hashToScalar(hash, e, q)) {
        toMpz(rz, r); toMpz(sz, s);
        mpz_invert(v, e, q);                               // v = e^-1 mod q
        mpz_mul(z1, sz, v); mpz_mod(z1, z1, q);            // z1 = s*v
        mpz_mul(z2, rz, v); mpz_neg(z2, z2); mpz_mod(z2, z2, q);   // z2 = -r*v
        uint8_t z1b[32], z2b[32], xb[32];
        fromMpz(z1, z1b); fromMpz(z2, z2b);
        // данные публичные: допускается обычный Double-and-Add
        Point A = pointMulDoubleAndAdd(z1b, pointBase());
        Point B = pointMulDoubleAndAdd(z2b, Q);
        Point C = pointAdd(A, B);
        if (!C.inf) {
            feToBytes(C.x, xb);
            mpz_t R; mpz_init(R);
            toMpz(R, xb); mpz_mod(R, R, q);                // R = x_c' mod q
            result = (mpz_cmp(R, rz) == 0);
            mpz_clear(R);
        }
    }
    mpz_clear(q); mpz_clear(e); mpz_clear(v); mpz_clear(rz); mpz_clear(sz); mpz_clear(z1); mpz_clear(z2);
    return result;
}
