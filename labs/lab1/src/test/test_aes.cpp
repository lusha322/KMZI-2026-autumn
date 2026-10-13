#include "../crypto/aes128.h"
#include <cstdio>
#include <cstring>
#include <cstdint>

static int g_passed = 0;
static int g_failed = 0;

#define CHECK(cond, msg)                                            \
    do {                                                            \
        if (cond) { ++g_passed; std::printf("  [OK]   %s\n", msg);  \
        } else    { ++g_failed; std::printf("  [FAIL] %s\n", msg); }\
    } while (0)

static void hex(const uint8_t* p, int n, char* out) {
    static const char* H = "0123456789ABCDEF";
    for (int i = 0; i < n; ++i) {
        out[2*i]   = H[p[i] >> 4];
        out[2*i+1] = H[p[i] & 0xF];
    }
    out[2*n] = 0;
}

static bool eq_hex(const uint8_t* buf, int n, const char* expect) {
    char tmp[64];
    hex(buf, n, tmp);
    return std::strcmp(tmp, expect) == 0;
}

static void test_fips_197_b() {
    std::printf("=== AES-128: FIPS-197 Прил. B ===\n");
    const uint8_t key[16] = {
        0x00,0x01,0x02,0x03,0x04,0x05,0x06,0x07,
        0x08,0x09,0x0A,0x0B,0x0C,0x0D,0x0E,0x0F
    };
    const uint8_t pt[16] = {
        0x00,0x11,0x22,0x33,0x44,0x55,0x66,0x77,
        0x88,0x99,0xAA,0xBB,0xCC,0xDD,0xEE,0xFF
    };
    AES128 aes(key);
    uint8_t ct[16];
    aes.encrypt(pt, ct);
    CHECK(eq_hex(ct, 16, "69C4E0D86A7B0430D8CDB78070B4C55A"),
          "E_K(001122...FF) = 69C4...C55A");
}

static void test_fips_197_c1() {
    std::printf("=== AES-128: FIPS-197 Прил. C.1 ===\n");
    const uint8_t key[16] = {0};
    const uint8_t pt[16]  = {0};
    AES128 aes(key);
    uint8_t ct[16];
    aes.encrypt(pt, ct);
    CHECK(eq_hex(ct, 16, "66E94BD4EF8A2C3B884CFA59CA342B2E"),
          "E_0(0) = 66E9...2B2E");
}

static void test_fips_197_c2() {
    std::printf("=== AES-128: FIPS-197 Прил. C.2 ===\n");
    const uint8_t key[16] = {0};
    const uint8_t pt[16] = {
        0x00,0x11,0x22,0x33,0x44,0x55,0x66,0x77,
        0x88,0x99,0xAA,0xBB,0xCC,0xDD,0xEE,0xFF
    };
    AES128 aes(key);
    uint8_t ct[16];
    aes.encrypt(pt, ct);
    CHECK(eq_hex(ct, 16, "69C4E0D86A7B0430D8CDB78070B4C55A"),
          "E_0(001122...FF) = 69C4...C55A");
}

static void test_roundtrip() {
    std::printf("=== AES-128: round-trip на 1000 блоках ===\n");
    uint8_t key[16] = {0};
    AES128 aes(key);

    uint32_t lcg = 0x12345678u;
    bool ok = true;
    for (int t = 0; t < 1000; ++t) {
        uint8_t pt[16], ct[16], back[16];
        for (int i = 0; i < 16; ++i) {
            lcg = lcg * 1103515245u + 12345u;
            pt[i] = static_cast<uint8_t>(lcg >> 16);
        }
        aes.encrypt(pt, ct);
        aes.decrypt(ct, back);
        if (std::memcmp(pt, back, 16) != 0) { ok = false; break; }
    }
    CHECK(ok, "decrypt(encrypt(P)) = P на 1000 блоках");
}

static void test_avalanche() {
    std::printf("=== AES-128: лавинный эффект ===\n");
    uint8_t key[16] = {0};
    AES128 aes(key);

    uint8_t pt1[16] = {0};
    uint8_t pt2[16] = {0};
    pt2[0] = 0x01;

    uint8_t c1[16], c2[16];
    aes.encrypt(pt1, c1);
    aes.encrypt(pt2, c2);

    int diff = 0;
    for (int i = 0; i < 16; ++i) {
        uint8_t x = static_cast<uint8_t>(c1[i] ^ c2[i]);
        while (x) { diff += (x & 1); x >>= 1; }
    }
    CHECK(diff >= 40 && diff <= 88,
          "изменение 1 бита PT → ~64 бита CT");
}

int run_test_aes() {
    g_passed = 0; g_failed = 0;
    test_fips_197_b();
    test_fips_197_c1();
    test_fips_197_c2();
    test_roundtrip();
    test_avalanche();
    std::printf("--- test_aes: %d passed, %d failed ---\n\n",
                g_passed, g_failed);
    return g_failed;
}