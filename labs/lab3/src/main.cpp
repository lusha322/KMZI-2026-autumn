#include <cstdio>
#include <cstring>
#include <chrono>
#include <gmp.h>     // в main - только для сверки арифметики поля с эталоном
#include "ec_algebra.h"
#include "ct_secure.h"
#include "protocols.h"

static int g_failed = 0;

static void check(const char* name, bool ok) {
    std::printf("%-76s %s\n", name, ok ? "PASS" : "FAIL");
    if (!ok) ++g_failed;
}

static void hexToBytes(const char* hex, uint8_t* out, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        unsigned v = 0;
        std::sscanf(hex + 2 * i, "%2x", &v);
        out[i] = (uint8_t)v;
    }
}

static void printHex(const char* label, const uint8_t* d, size_t n) {
    std::printf("%s: ", label);
    for (size_t i = 0; i < n; ++i) std::printf("%02x", d[i]);
    std::printf("\n");
}

static void printPoint(const char* label, const Point& P) {
    std::printf("%s\n", label);
    if (P.inf) { std::printf("  O (бесконечно удалённая точка)\n"); return; }
    uint8_t b[32];
    feToBytes(P.x, b); printHex("  x", b, 32);
    feToBytes(P.y, b); printHex("  y", b, 32);
}

static void feToMpz(mpz_t z, const Fe& a) {
    uint8_t b[32]; feToBytes(a, b);
    mpz_import(z, 32, 1, 1, 1, 0, b);
}

// псевдослучайный Fe для тестов арифметики
static uint64_t g_rng = 0x1234567890abcdefULL;
static uint64_t rnd64() { g_rng ^= g_rng << 13; g_rng ^= g_rng >> 7; g_rng ^= g_rng << 17; return g_rng; }

static void fieldTests() {
    std::printf("\n--- Арифметика GF(p), сверка с GMP ---\n");
    mpz_t p, x, y, r, ref;
    mpz_init(p); mpz_init(x); mpz_init(y); mpz_init(r); mpz_init(ref);
    mpz_set_str(p, "ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff43", 16);
    bool okMul = true, okAdd = true, okSub = true, okInv = true;
    for (int it = 0; it < 300; ++it) {
        Fe a, b;   // сырые 256-битные значения, затем приводятся к [0, p) через GMP
        for (int i = 0; i < 4; ++i) { a.v[i] = rnd64(); b.v[i] = rnd64(); }
        if (it < 20) {                          // граничные значения: близко к p
            for (int i = 1; i < 4; ++i) a.v[i] = ~0ULL;
            a.v[0] = ~0ULL - (rnd64() % 190);
        }
        feToMpz(x, a); feToMpz(y, b);
        mpz_mod(x, x, p); mpz_mod(y, y, p);     // канонические представители
        Fe ca, cb;
        // построим канонические a, b из x, y
        char hx[65], hy[65];
        gmp_snprintf(hx, sizeof hx, "%064Zx", x);
        gmp_snprintf(hy, sizeof hy, "%064Zx", y);
        feFromHex(hx, &ca); feFromHex(hy, &cb);
        feToMpz(r, feMul(ca, cb)); mpz_mul(ref, x, y); mpz_mod(ref, ref, p);
        okMul &= mpz_cmp(r, ref) == 0;
        feToMpz(r, feAdd(ca, cb)); mpz_add(ref, x, y); mpz_mod(ref, ref, p);
        okAdd &= mpz_cmp(r, ref) == 0;
        feToMpz(r, feSub(ca, cb)); mpz_sub(ref, x, y); mpz_mod(ref, ref, p);
        okSub &= mpz_cmp(r, ref) == 0;
        if (mpz_sgn(x) != 0) {
            Fe inv = feInvert(ca);
            okInv &= feEqual(feMul(inv, ca), feOne()) == 1;
        }
    }
    check("feMul совпадает с GMP (300 пар, включая значения около p)", okMul);
    check("feAdd совпадает с GMP", okAdd);
    check("feSub совпадает с GMP", okSub);
    check("feInvert: a * a^-1 = 1", okInv);
    check("feInvert(0) = 0 (безопасное вырожденное значение)", feIsZero(feInvert(feZero())) == 1);
    mpz_clear(p); mpz_clear(x); mpz_clear(y); mpz_clear(r); mpz_clear(ref);
}

static void curveTests() {
    std::printf("\n--- Параметры кривой и групповые законы ---\n");
    const CurveParams& c = curve();
    Point G = pointBase();
    check("G лежит на кривой y^2 = x^3 + ax + b", pointOnCurve(G));
    check("4a^3 + 27b^2 != 0 (кривая невырождена)", [&]{
        Fe a3 = feMul(feSquare(c.a), c.a);
        Fe t = feAdd(feMul(feFromU64(4), a3), feMul(feFromU64(27), feSquare(c.b)));
        return feIsZero(t) == 0; }());
    Point qG = SecureMultiply(c.q, G);
    check("q * G = O (порядок базовой точки)", qG.inf == 1);
    uint8_t qm1[32]; std::memcpy(qm1, c.q, 32); qm1[31] -= 1;
    check("(q-1) * G = -G", pointEqual(SecureMultiply(qm1, G), pointNeg(G)) == 1);

    Point O = pointInfinity();
    Point G2 = pointDouble(G);
    Point G3 = pointAdd(G2, G);
    check("2G лежит на кривой", pointOnCurve(G2));
    check("G + O = G, O + G = G", pointEqual(pointAdd(G, O), G) && pointEqual(pointAdd(O, G), G));
    check("G + (-G) = O", pointAdd(G, pointNeg(G)).inf == 1);
    check("G + G = 2G (удвоение внутри полной формулы)", pointEqual(pointAdd(G, G), G2) == 1);
    check("(G + 2G) = (2G + G), коммутативность", pointEqual(G3, pointAdd(G, G2)) == 1);
    check("(G + 2G) + 3G = G + (2G + 3G), ассоциативность",
          pointEqual(pointAdd(G3, pointAdd(G2, G)), pointAdd(G, pointAdd(G2, G3))) == 1);
    check("O + O = O", pointAdd(O, O).inf == 1);

    uint8_t s2[32] = {0}, s3[32] = {0}, s0[32] = {0}, s1[32] = {0};
    s2[31] = 2; s3[31] = 3; s1[31] = 1;
    check("Лестница: 0*G = O", SecureMultiply(s0, G).inf == 1);
    check("Лестница: 1*G = G", pointEqual(SecureMultiply(s1, G), G) == 1);
    check("Лестница: 2*G = G + G", pointEqual(SecureMultiply(s2, G), G2) == 1);
    check("Лестница: 3*G = 2G + G", pointEqual(SecureMultiply(s3, G), G3) == 1);

    bool ok = true;
    for (int i = 0; i < 20; ++i) {
        uint8_t k[32];
        randomScalar(k);
        Point A = SecureMultiply(k, G), B = pointMulDoubleAndAdd(k, G);
        ok &= pointEqual(A, B) == 1 && pointOnCurve(A);
    }
    check("Лестница = Double-and-Add для 20 случайных скаляров", ok);
    {
        // (a+b)G = aG + bG
        uint8_t a[32] = {0}, b[32] = {0}, ab[32] = {0};
        a[31] = 0x37; b[30] = 0x01; b[31] = 0x05;
        unsigned carry = 0;
        for (int i = 31; i >= 0; --i) { unsigned t = a[i] + b[i] + carry; ab[i] = (uint8_t)t; carry = t >> 8; }
        check("(a+b)*G = a*G + b*G", pointEqual(SecureMultiply(ab, G),
              pointAdd(SecureMultiply(a, G), SecureMultiply(b, G))) == 1);
    }
}

static void cswapTests() {
    std::printf("\n--- Условный обмен точек (CSWAP) ---\n");
    Point A = pointBase(), B = pointDouble(A), A0 = A, B0 = B;
    ctCswapPoint(0, &A, &B);
    check("CSWAP(b=0): точки не изменились", pointEqual(A, A0) && pointEqual(B, B0));
    ctCswapPoint(1, &A, &B);
    check("CSWAP(b=1): точки поменялись местами", pointEqual(A, B0) && pointEqual(B, A0));
}

static void spaTrace() {
    std::printf("\n--- Профиль операций лестницы (SPA): число сложений/удвоений ---\n");
    Point G = pointBase();
    struct { const char* name; const char* hex; } sc[] = {
        {"d = 1 (вес Хэмминга 1)  ", "0000000000000000000000000000000000000000000000000000000000000001"},
        {"d = 2^255 + 1           ", "8000000000000000000000000000000000000000000000000000000000000001"},
        {"d = q-1 (плотный)       ", "ffffffffffffffffffffffffffffffffd95c8ed60dfb4dfc7e5abf99263d6606"},
        {"d = 0xAAAA... (чередов.)", "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"},
    };
    bool same = true; uint64_t a0 = 0, d0 = 0;
    for (int i = 0; i < 4; ++i) {
        uint8_t k[32]; hexToBytes(sc[i].hex, k, 32);
        g_addCalls = g_dblCalls = 0;
        SecureMultiply(k, G);
        std::printf("  Лестница        %s: сложений %llu, удвоений %llu\n", sc[i].name,
                    (unsigned long long)g_addCalls, (unsigned long long)g_dblCalls);
        if (i == 0) { a0 = g_addCalls; d0 = g_dblCalls; }
        same &= (g_addCalls == a0 && g_dblCalls == d0);
    }
    for (int i = 0; i < 2; ++i) {
        uint8_t k[32]; hexToBytes(sc[i == 0 ? 0 : 2].hex, k, 32);
        g_addCalls = g_dblCalls = 0;
        pointMulDoubleAndAdd(k, G);
        std::printf("  Double-and-Add  %s: сложений %llu, удвоений %llu\n", sc[i == 0 ? 0 : 2].name,
                    (unsigned long long)g_addCalls, (unsigned long long)g_dblCalls);
    }
    check("Лестница: профиль операций не зависит от скаляра (256 + 256)", same && a0 == 256 && d0 == 256);
}

static void timingTest() {
    std::printf("\n--- Эмпирическое время скалярного умножения (индикативно) ---\n");
    Point G = pointBase();
    uint8_t lo[32] = {0}, hi[32];
    lo[31] = 1;
    hexToBytes("ffffffffffffffffffffffffffffffffd95c8ed60dfb4dfc7e5abf99263d6606", hi, 32);
    const int N = 15;
    double t[4];
    for (int round = 0; round < 2; ++round) {          // раунд 0 - прогрев
        for (int m = 0; m < 4; ++m) {
            const uint8_t* k = (m % 2 == 0) ? lo : hi;
            auto t0 = std::chrono::steady_clock::now();
            for (int i = 0; i < N; ++i) {
                Point R = (m < 2) ? SecureMultiply(k, G) : pointMulDoubleAndAdd(k, G);
                if (R.inf == 2) std::printf("x");
            }
            auto t1 = std::chrono::steady_clock::now();
            t[m] = std::chrono::duration<double, std::milli>(t1 - t0).count() / N;
        }
    }
    std::printf("  Лестница Монтгомери, d = 1:        %.2f мс\n", t[0]);
    std::printf("  Лестница Монтгомери, d = q-1:      %.2f мс\n", t[1]);
    std::printf("  Double-and-Add,      d = 1:        %.2f мс\n", t[2]);
    std::printf("  Double-and-Add,      d = q-1:      %.2f мс\n", t[3]);
    double dl = t[1] > t[0] ? t[1] / t[0] : t[0] / t[1];
    std::printf("  Отношение времён (лестница): %.3f; (Double-and-Add): %.3f\n", dl, t[3] > t[2] ? t[3]/t[2] : t[2]/t[3]);
}

struct Msg { const char* name; const char* msgHex; const char* hashHex; };

static void ecdhAndSignatures() {
    std::printf("\n--- SHA-256: проверка хэш-функции ---\n");
    {
        uint8_t abc[3] = {'a', 'b', 'c'}, h[32], ex[32];
        sha256(abc, 3, h);
        hexToBytes("ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad", ex, 32);
        check("SHA-256(\"abc\") совпадает с эталоном", std::memcmp(h, ex, 32) == 0);
    }

    std::printf("\n--- ECDH: тестовые ключи (hex) ---\n");
    uint8_t dA[32], dB[32];
    hexToBytes("1f2e3d4c5b6a79881726354453627180918a9bacbdcedfe0f1021324354657a8", dA, 32);
    hexToBytes("0a1b2c3d4e5f60718293a4b5c6d7e8f90123456789abcdef0fedcba987654321", dB, 32);
    Point QA, QB;
    check("dA в диапазоне (0, q), dB в диапазоне (0, q)", scalarInRange(dA) && scalarInRange(dB));
    keyGenPublic(dA, &QA); keyGenPublic(dB, &QB);
    printHex("dA", dA, 32); printPoint("QA = dA*G", QA);
    printHex("dB", dB, 32); printPoint("QB = dB*G", QB);
    check("QA, QB лежат на кривой", pointOnCurve(QA) && pointOnCurve(QB));
    uint8_t KA[32], KB[32];
    bool o1 = ecdhSharedSecret(dA, QB, KA), o2 = ecdhSharedSecret(dB, QA, KB);
    printHex("Shared secret Алисы (x(dA*QB))", KA, 32);
    printHex("Shared secret Боба  (x(dB*QA))", KB, 32);
    check("ECDH: Алиса и Боб получили одинаковый секрет", o1 && o2 && ctCompareBytes(KA, KB, 32) == 1);

    bool okRnd = true;
    for (int i = 0; i < 5; ++i) {
        uint8_t a[32], b[32], k1[32], k2[32]; Point qa, qb;
        randomScalar(a); randomScalar(b);
        keyGenPublic(a, &qa); keyGenPublic(b, &qb);
        okRnd &= ecdhSharedSecret(a, qb, k1) && ecdhSharedSecret(b, qa, k2) && ctCompareBytes(k1, k2, 32);
    }
    check("ECDH с эфемерными ключами /dev/urandom (5 запусков)", okRnd);
    Point bad = QB; bad.y = feAdd(bad.y, feOne());
    uint8_t tmp[32];
    check("ECDH отвергает точку не на кривой", !ecdhSharedSecret(dA, bad, tmp));

    std::printf("\n--- ЭЦП СТБ 34.101.45: цикл «подпись - проверка» ---\n");
    static const Msg msgs[] = {
        {"M1", "616263", "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"},
        {"M2", "5354422033342e3130312e3435206c616233206d657373616765202331",
               "6719d14d3312ff4e7295e50d2440ecbd78435c1b52217e3da8970c91ac36a004"},
        {"M3", "45434448202b207369676e6174757265207465737420766563746f72",
               "63b1d20bc37de4bae6adf99f863c55209f975431ed910fa15c55c44eba700720"},
    };
    static const char* kHex[3] = {
        "3c4d5e6f708192a3b4c5d6e7f8091a2b3c4d5e6f708192a3b4c5d6e7f8091a2b",
        "0f1e2d3c4b5a69788796a5b4c3d2e1f00f1e2d3c4b5a69788796a5b4c3d2e1f0",
        "5a5a5a5a5a5a5a5a5a5a5a5a5a5a5a5a1234567890abcdef1234567890abcdef"};
    for (int i = 0; i < 3; ++i) {
        uint8_t msg[64], h[32], hex[32], k[32], r[32], s[32];
        size_t mlen = std::strlen(msgs[i].msgHex) / 2;
        hexToBytes(msgs[i].msgHex, msg, mlen);
        sha256(msg, mlen, h);
        hexToBytes(msgs[i].hashHex, hex, 32);
        hexToBytes(kHex[i], k, 32);
        std::printf("\nСообщение %s (hex): %s\n", msgs[i].name, msgs[i].msgHex);
        printHex("  e = SHA-256(M)", h, 32);
        check("  хэш совпадает с заданной в коде hex-строкой", std::memcmp(h, hex, 32) == 0);
        printHex("  k", k, 32);
        bool sg = stbSign(dA, h, k, r, s);
        check("  подпись выработана (r != 0, s != 0)", sg);
        printHex("  r", r, 32); printHex("  s", s, 32);
        check("  проверка корректной подписи: ACCEPT", stbVerify(QA, h, r, s));
        uint8_t s2[32]; std::memcpy(s2, s, 32); s2[31] ^= 0x01;
        check("  искажён 1 байт (младший бит) s: REJECT", !stbVerify(QA, h, r, s2));
        std::memcpy(s2, s, 32); s2[0] ^= 0x80;
        check("  искажён 1 байт (старший бит) s: REJECT", !stbVerify(QA, h, r, s2));
        uint8_t r2[32]; std::memcpy(r2, r, 32); r2[15] ^= 0x10;
        check("  искажён байт r: REJECT", !stbVerify(QA, h, r2, s));
        uint8_t h2[32]; std::memcpy(h2, h, 32); h2[7] ^= 0x01;
        check("  искажён хэш сообщения: REJECT", !stbVerify(QA, h2, r, s));
        check("  чужой открытый ключ QB: REJECT", !stbVerify(QB, h, r, s));
        uint8_t zero[32] = {0};
        check("  r = 0 или s = 0: REJECT", !stbVerify(QA, h, zero, s) && !stbVerify(QA, h, r, zero));
    }
    bool okR = true;
    for (int i = 0; i < 5; ++i) {
        uint8_t d[32], k[32], h[32], r[32], s[32]; Point Q;
        randomScalar(d); randomScalar(k); randomScalar(h);
        keyGenPublic(d, &Q);
        okR &= stbSign(d, h, k, r, s) && stbVerify(Q, h, r, s);
    }
    check("Подпись/проверка со случайными d, k, e (5 запусков)", okR);
}

int main() {
    std::printf("Лабораторная работа №3, вариант 1: СТБ 34.101.45, Вейерштрасс / аффинные, лестница Монтгомери\n");
    const CurveParams& c = curve();
    std::printf("\nПараметры кривой (bign-curve256v1):\n");
    uint8_t b[32];
    feToBytes(c.a, b); printHex("  a ", b, 32);
    feToBytes(c.b, b); printHex("  b ", b, 32);
    printHex("  q ", c.q, 32);
    feToBytes(c.Gx, b); printHex("  Gx", b, 32);
    feToBytes(c.Gy, b); printHex("  Gy", b, 32);
    fieldTests();
    curveTests();
    cswapTests();
    spaTrace();
    ecdhAndSignatures();
    timingTest();
    std::printf("\nИтог: %s\n", g_failed == 0 ? "ВСЕ ТЕСТЫ ПРОЙДЕНЫ" : "ЕСТЬ ОШИБКИ");
    return g_failed == 0 ? 0 : 1;
}
