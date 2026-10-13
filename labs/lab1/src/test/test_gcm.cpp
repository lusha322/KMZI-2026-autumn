#include "../crypto/aes128.h"
#include "../crypto/gcm.h"
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
    char tmp[128];
    hex(buf, n, tmp);
    return std::strcmp(tmp, expect) == 0;
}

static void from_hex(const char* s, uint8_t* out, int n) {
    auto val = [](char c) -> int {
        if (c >= '0' && c <= '9') return c - '0';
        if (c >= 'A' && c <= 'F') return c - 'A' + 10;
        if (c >= 'a' && c <= 'f') return c - 'a' + 10;
        return 0;
    };
    for (int i = 0; i < n; ++i)
        out[i] = static_cast<uint8_t>((val(s[2*i]) << 4) | val(s[2*i+1]));
}

static void test_nist_tc1() {
    std::printf("=== GCM: NIST TC1 (пустой P/AAD) ===\n");
    const uint8_t key[16] = {0};
    const uint8_t nonce[12] = {0};
    AES128 aes(key);
    GCM gcm(aes);
    gcm.init();

    GCM::Tag tag;
    gcm.encrypt(nullptr, 0, nullptr, 0, nonce, nullptr, tag);

    CHECK(eq_hex(tag.data(), 16, "58E2FCCEFA7E3061367F1D57A4E7455A"),
          "T = 58E2...455A");
}

static void test_nist_tc2() {
    std::printf("=== GCM: NIST TC2 (P=0^128, IV=0) ===\n");
    const uint8_t key[16] = {0};
    const uint8_t nonce[12] = {0};
    uint8_t pt[16] = {0};

    AES128 aes(key);
    GCM gcm(aes);
    gcm.init();

    uint8_t ct[16];
    GCM::Tag tag;
    gcm.encrypt(nullptr, 0, pt, 16, nonce, ct, tag);

    CHECK(eq_hex(ct, 16, "0388DACE60B6A392F328C2B971B2FE78"),
          "C = 0388...FE78");
    CHECK(eq_hex(tag.data(), 16, "AB6E47D42CEC13BDF53A67B21257BDDF"),
          "T = AB6E...BDDF");
}

static void test_nist_tc3() {
    std::printf("=== GCM: NIST TC3 (60 байт P, без AAD) ===\n");
    uint8_t key[16], pt[60], nonce[12];
    from_hex("FEFFE9928665731C6D6A8F9467308308", key, 16);
    from_hex("CAFEBABEFACEDBADDECAF888", nonce, 12);
    from_hex("D9313225F88406E5A55909C5AFF5269A"
             "86A7A9531534F7DA2E4C303D8A318A72"
             "1C3C0C95956809532FCF0E2449A6B525"
             "B16AEDF5AA0DE657BA637B39", pt, 60);

    AES128 aes(key);
    GCM gcm(aes);
    gcm.init();

    uint8_t ct[60];
    GCM::Tag tag;
    gcm.encrypt(nullptr, 0, pt, 60, nonce, ct, tag);

    CHECK(eq_hex(ct, 60,
        "42831EC2217774244B7221B784D0D49C"
        "E3AA212F2C02A4E035C17E2329ACA12E"
        "21D514B25466931C7D8F6A5AAC84AA05"
        "1BA30B396A0AAC973D58E091"),
        "C = 4283...E091");
    CHECK(eq_hex(tag.data(), 16, "4D5C2AF327CD64A62CF35ABD2BA6FAB4"),
          "T = 4D5C...FAB4");
}

static void test_roundtrip_lengths() {
    std::printf("=== GCM: round-trip для разных длин ===\n");
    const uint8_t key[16] = {0};
    const uint8_t nonce[12] = {0};
    AES128 aes(key);
    GCM gcm(aes);
    gcm.init();

    const size_t lens[] = {0,1,15,16,17,31,32,33,64,100};
    bool ok = true;
    for (size_t L : lens) {
        uint8_t pt[100], ct[100], back[100];
        for (size_t i = 0; i < L; ++i)
            pt[i] = static_cast<uint8_t>(i * 7 + 3);

        GCM::Tag tag;
        gcm.encrypt(nullptr, 0, pt, L, nonce, ct, tag);

        bool ok_dec = gcm.decrypt(nullptr, 0, ct, L, nonce, tag, back);
        if (!ok_dec) { ok = false; break; }
        if (std::memcmp(pt, back, L) != 0) { ok = false; break; }
    }
    CHECK(ok, "round-trip для L ∈ {0,1,15,16,17,31,32,33,64,100}");
}

static void test_aad_only_affects_tag() {
    std::printf("=== GCM: AAD влияет только на тег ===\n");
    const uint8_t key[16] = {0};
    const uint8_t nonce[12] = {0};
    uint8_t pt[16] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
    uint8_t aad[8] = {0xFE,0xED,0xFA,0xCE,0xDE,0xAD,0xBE,0xEF};

    AES128 aes(key);
    GCM gcm(aes);
    gcm.init();

    uint8_t ct1[16], ct2[16];
    GCM::Tag t1, t2;
    gcm.encrypt(nullptr, 0, pt, 16, nonce, ct1, t1);
    gcm.encrypt(aad, 8, pt, 16, nonce, ct2, t2);

    CHECK(std::memcmp(ct1, ct2, 16) == 0, "CT не зависит от AAD");
    CHECK(std::memcmp(t1.data(), t2.data(), 16) != 0, "Tag зависит от AAD");
}

static void test_forgery_rejected() {
    std::printf("=== GCM: отказ при подделке ===\n");
    const uint8_t key[16] = {0};
    const uint8_t nonce[12] = {0};
    uint8_t pt[16] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};

    AES128 aes(key);
    GCM gcm(aes);
    gcm.init();

    uint8_t ct[16], back[16];
    GCM::Tag tag;
    gcm.encrypt(nullptr, 0, pt, 16, nonce, ct, tag);

    // Подмена CT
    uint8_t bad_ct[16];
    std::memcpy(bad_ct, ct, 16);
    bad_ct[0] ^= 0x01;
    CHECK(!gcm.decrypt(nullptr, 0, bad_ct, 16, nonce, tag, back),
          "подмена CT → отказ");

    // Подмена тега в каждом байте
    bool all_ok = true;
    for (int i = 0; i < 16; ++i) {
        GCM::Tag bad = tag;
        bad[i] ^= 0x01;
        if (gcm.decrypt(nullptr, 0, ct, 16, nonce, bad, back))
            all_ok = false;
    }
    CHECK(all_ok, "подмена любого байта тега → отказ (16 проверок)");

    // Нулевой тег
    GCM::Tag zero{};
    CHECK(!gcm.decrypt(nullptr, 0, ct, 16, nonce, zero, back),
          "нулевой тег → отказ");
}

static void test_ct_compare() {
    std::printf("=== GCM: Constant-time сравнение тегов ===\n");
    GCM::Tag a = {0}, b = {0};
    CHECK(GCM::ct_compare(a, b), "равные теги → true");

    for (int i = 0; i < 16; ++i) {
        GCM::Tag c = a;
        c[i] = 0x01;
        if (GCM::ct_compare(a, c)) { CHECK(false, "diff не пойман"); return; }
    }
    CHECK(true, "отличие в любом байте → false (16 проверок)");

    GCM::Tag d = a;
    d[0]  = 0xFF;
    GCM::Tag e = a;
    e[15] = 0xFF;
    CHECK(!GCM::ct_compare(a, d) && !GCM::ct_compare(a, e),
          "позиция различия не влияет на результат");
}

static void test_nonzero_tag() {
    std::printf("=== GCM: тег ненулевой на случайных данных ===\n");
    uint32_t lcg = 0xDEADBEEFu;
    bool all_nonzero = true;
    for (int t = 0; t < 64; ++t) {
        uint8_t key[16], nonce[12], pt[16];
        for (int i = 0; i < 16; ++i) {
            lcg = lcg * 1103515245u + 12345u;
            key[i] = static_cast<uint8_t>(lcg >> 16);
        }
        for (int i = 0; i < 12; ++i) {
            lcg = lcg * 1103515245u + 12345u;
            nonce[i] = static_cast<uint8_t>(lcg >> 16);
        }
        for (int i = 0; i < 16; ++i) {
            lcg = lcg * 1103515245u + 12345u;
            pt[i] = static_cast<uint8_t>(lcg >> 16);
        }

        AES128 aes(key);
        GCM gcm(aes);
        gcm.init();
        uint8_t ct[16];
        GCM::Tag tag;
        gcm.encrypt(nullptr, 0, pt, 16, nonce, ct, tag);

        bool zero = true;
        for (int i = 0; i < 16; ++i) if (tag[i]) { zero = false; break; }
        if (zero) { all_nonzero = false; break; }
    }
    CHECK(all_nonzero, "тег ≠ 0 на 64 случайных входах");
}

int run_test_gcm() {
    g_passed = 0; g_failed = 0;
    test_nist_tc1();
    test_nist_tc2();
    test_nist_tc3();
    test_roundtrip_lengths();
    test_aad_only_affects_tag();
    test_forgery_rejected();
    test_ct_compare();
    test_nonzero_tag();
    std::printf("--- test_gcm: %d passed, %d failed ---\n\n",
                g_passed, g_failed);
    return g_failed;
}