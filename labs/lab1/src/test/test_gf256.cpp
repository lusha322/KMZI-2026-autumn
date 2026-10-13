#include "../crypto/gf256.h"
#include <cstdio>
#include <cstdint>

static int g_passed = 0;
static int g_failed = 0;

#define CHECK(cond, msg)                                            \
    do {                                                            \
        if (cond) { ++g_passed; std::printf("  [OK]   %s\n", msg);  \
        } else    { ++g_failed; std::printf("  [FAIL] %s\n", msg); }\
    } while (0)

static void test_basic_axioms() {
    std::printf("=== GF(2^8) p(x)=0x11B: аксиомы поля ===\n");
    GF256 gf(0x11B);

    bool ok_mul_one = true, ok_add_zero = true, ok_comm = true;
    for (int a = 0; a < 256; ++a) {
        if (gf.multiply(static_cast<uint8_t>(a), 0x01) != a) ok_mul_one = false;
        if (GF256::add(static_cast<uint8_t>(a), 0x00) != a)  ok_add_zero = false;
    }
    for (int a = 0; a < 256; ++a)
        for (int b = 0; b < 256; ++b)
            if (gf.multiply(static_cast<uint8_t>(a), static_cast<uint8_t>(b))
             != gf.multiply(static_cast<uint8_t>(b), static_cast<uint8_t>(a)))
                ok_comm = false;

    CHECK(ok_mul_one,  "a * 1 = a для всех a");
    CHECK(ok_add_zero, "a + 0 = a для всех a");
    CHECK(ok_comm,     "a * b = b * a (все 65536 пар)");

    bool ok_distr = true;
    for (int a = 0; a < 256; a += 7)
        for (int b = 0; b < 256; b += 11)
            for (int c = 0; c < 256; c += 13) {
                const uint8_t A = static_cast<uint8_t>(a);
                const uint8_t B = static_cast<uint8_t>(b);
                const uint8_t C = static_cast<uint8_t>(c);
                const uint8_t lhs = gf.multiply(A, GF256::add(B, C));
                const uint8_t rhs = GF256::add(gf.multiply(A, B), gf.multiply(A, C));
                if (lhs != rhs) ok_distr = false;
            }
    CHECK(ok_distr, "a*(b+c) = a*b + a*c (выборка)");
}

static void test_xtime() {
    std::printf("=== GF(2^8): xtime(b) = b * 0x02 ===\n");
    GF256 gf(0x11B);

    bool ok = true;
    for (int b = 0; b < 256; ++b) {
        const uint8_t bb = static_cast<uint8_t>(b);
        if (gf.xtime(bb) != gf.multiply(bb, 0x02)) ok = false;
    }
    CHECK(ok, "xtime(b) = b * 0x02 для всех b");

    CHECK(gf.xtime(0x80) == 0x1B, "xtime(0x80) = 0x1B (переполнение)");
    CHECK(gf.xtime(0x57) == 0xAE, "xtime(0x57) = 0xAE");
    CHECK(gf.xtime(0xFF) == 0xE1, "xtime(0xFF) = 0xE1");
}

static void test_inverse() {
    std::printf("=== GF(2^8): обратный элемент (a^254) ===\n");
    GF256 gf(0x11B);

    CHECK(gf.inverse(0x00) == 0x00, "inverse(0x00) = 0x00 (соглашение)");
    CHECK(gf.inverse(0x01) == 0x01, "inverse(0x01) = 0x01");
    CHECK(gf.inverse(0x53) == 0xCA, "inverse(0x53) = 0xCA");

    bool ok = true;
    for (int a = 1; a < 256; ++a) {
        const uint8_t aa = static_cast<uint8_t>(a);
        if (gf.multiply(aa, gf.inverse(aa)) != 0x01) { ok = false; break; }
    }
    CHECK(ok, "inverse(a) * a = 1 для всех a != 0");
}

static void test_known_products() {
    std::printf("=== GF(2^8): эталонные произведения ===\n");
    GF256 gf(0x11B);

    CHECK(gf.multiply(0x57, 0x83) == 0xC1, "0x57 * 0x83 = 0xC1");
    CHECK(gf.multiply(0x57, 0x13) == 0xFE, "0x57 * 0x13 = 0xFE");
    CHECK(gf.multiply(0x00, 0xFF) == 0x00, "0x00 * 0xFF = 0x00");
    CHECK(gf.multiply(0x01, 0xFF) == 0xFF, "0x01 * 0xFF = 0xFF");
}

static void test_sbox() {
    std::printf("=== S-box (вычисление на лету, вариант 1) ===\n");
    GF256 gf(0x11B);

    CHECK(gf.sbox(0x00) == 0x63, "S(0x00) = 0x63");
    CHECK(gf.sbox(0x01) == 0x7C, "S(0x01) = 0x7C");
    CHECK(gf.sbox(0x53) == 0xED, "S(0x53) = 0xED");
    CHECK(gf.sbox(0xFF) == 0x16, "S(0xFF) = 0x16");
    CHECK(gf.sbox(0x10) == 0xCA, "S(0x10) = 0xCA");

    // Биекция
    bool seen[256] = {false};
    bool biject = true;
    for (int b = 0; b < 256; ++b) {
        const uint8_t v = gf.sbox(static_cast<uint8_t>(b));
        if (seen[v]) { biject = false; break; }
        seen[v] = true;
    }
    CHECK(biject, "S-box — биекция на {0..255}");

    // Нет неподвижных точек
    bool no_fixed = true;
    for (int b = 0; b < 256; ++b)
        if (gf.sbox(static_cast<uint8_t>(b)) == static_cast<uint8_t>(b))
            no_fixed = false;
    CHECK(no_fixed, "нет неподвижных точек S(b) = b");
}

static void test_irreducibility() {
    std::printf("=== Проверка неприводимости p(x) ===\n");
    CHECK(GF256(0x11B).is_irreducible(), "0x11B неприводим (AES)");
    CHECK(GF256(0x12B).is_irreducible(), "0x12B неприводим");
    CHECK(GF256(0x14D).is_irreducible(), "0x14D неприводим");
    CHECK(!GF256(0x100).is_irreducible(), "0x100 = x^8 приводим");
    CHECK(!GF256(0x101).is_irreducible(), "0x101 = x^8+1 приводим");
}

int run_test_gf256() {
    g_passed = 0; g_failed = 0;
    test_basic_axioms();
    test_xtime();
    test_inverse();
    test_known_products();
    test_sbox();
    test_irreducibility();
    std::printf("--- test_gf256: %d passed, %d failed ---\n\n",
                g_passed, g_failed);
    return g_failed;
}