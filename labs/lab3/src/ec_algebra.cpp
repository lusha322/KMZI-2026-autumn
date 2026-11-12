#include "ec_algebra.h"

typedef unsigned __int128 u128;

// p = 2^256 - 189
static const uint64_t P_LIMBS[4] = {0xffffffffffffff43ULL, 0xffffffffffffffffULL,
                                    0xffffffffffffffffULL, 0xffffffffffffffffULL};
static const uint64_t EXP_P_MINUS_2[4] = {0xffffffffffffff41ULL, 0xffffffffffffffffULL,
                                          0xffffffffffffffffULL, 0xffffffffffffffffULL};
static const uint64_t C189 = 189;   // 2^256 mod p

uint64_t g_addCalls = 0;
uint64_t g_dblCalls = 0;

// ---------------------------------------------------------------- поле

Fe feZero() { Fe r = {{0, 0, 0, 0}}; return r; }
Fe feOne()  { Fe r = {{1, 0, 0, 0}}; return r; }
Fe feFromU64(uint64_t x) { Fe r = {{x, 0, 0, 0}}; return r; }

Fe feSelect(uint64_t bit, const Fe& a, const Fe& b) {
    uint64_t mask = 0 - (bit & 1);
    Fe r;
    for (int i = 0; i < 4; ++i) r.v[i] = (a.v[i] & mask) | (b.v[i] & ~mask);
    return r;
}

// Constant-time приведение: если a >= P, вычесть P (маска вместо while/if).
static Fe condSubP(const Fe& a) {
    Fe t;
    uint64_t borrow = 0;
    for (int i = 0; i < 4; ++i) {
        u128 d = (u128)a.v[i] - P_LIMBS[i] - borrow;
        t.v[i] = (uint64_t)d;
        borrow = (uint64_t)(d >> 64) & 1;
    }
    return feSelect(borrow, a, t);   // borrow = 1 <=> a < P
}

Fe feAdd(const Fe& a, const Fe& b) {
    Fe s;
    uint64_t c = 0;
    for (int i = 0; i < 4; ++i) {
        u128 t = (u128)a.v[i] + b.v[i] + c;
        s.v[i] = (uint64_t)t;
        c = (uint64_t)(t >> 64);
    }
    // 2^256 = 189 (mod p): при переносе добавляем 189 (всегда выполняем сложение, масштабируя на c)
    uint64_t add = C189 * c;
    uint64_t carry = add;
    for (int i = 0; i < 4; ++i) {
        u128 t = (u128)s.v[i] + carry;
        s.v[i] = (uint64_t)t;
        carry = (uint64_t)(t >> 64);
    }
    return condSubP(s);
}

Fe feSub(const Fe& a, const Fe& b) {
    Fe t;
    uint64_t borrow = 0;
    for (int i = 0; i < 4; ++i) {
        u128 d = (u128)a.v[i] - b.v[i] - borrow;
        t.v[i] = (uint64_t)d;
        borrow = (uint64_t)(d >> 64) & 1;
    }
    uint64_t mask = 0 - borrow;
    uint64_t c = 0;
    for (int i = 0; i < 4; ++i) {
        u128 s = (u128)t.v[i] + (P_LIMBS[i] & mask) + c;
        t.v[i] = (uint64_t)s;
        c = (uint64_t)(s >> 64);
    }
    return t;
}

Fe feNeg(const Fe& a) { return feSub(feZero(), a); }

// r += carry * 189, возвращает перенос из старшего слова
static uint64_t foldInto(uint64_t r[4], uint64_t carry) {
    u128 t = (u128)carry * C189;
    uint64_t c = 0;
    for (int i = 0; i < 4; ++i) {
        u128 s = (u128)r[i] + (uint64_t)t + c;
        r[i] = (uint64_t)s;
        c = (uint64_t)(s >> 64);
        t >>= 64;
    }
    return c;
}

Fe feMul(const Fe& a, const Fe& b) {
    uint64_t t[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    for (int i = 0; i < 4; ++i) {
        uint64_t carry = 0;
        for (int j = 0; j < 4; ++j) {
            u128 cur = (u128)a.v[i] * b.v[j] + t[i + j] + carry;
            t[i + j] = (uint64_t)cur;
            carry = (uint64_t)(cur >> 64);
        }
        t[i + 4] = carry;
    }
    // lo + hi * 189
    uint64_t r[4] = {t[0], t[1], t[2], t[3]};
    uint64_t carry = 0, top;
    for (int i = 0; i < 4; ++i) {
        u128 cur = (u128)t[4 + i] * C189 + r[i] + carry;
        r[i] = (uint64_t)cur;
        carry = (uint64_t)(cur >> 64);
    }
    top = carry;                       // < 2^8
    top = foldInto(r, top);            // перенос 0 или 1
    top = foldInto(r, top);
    foldInto(r, top);
    Fe res = {{r[0], r[1], r[2], r[3]}};
    return condSubP(res);
}

Fe feSquare(const Fe& a) { return feMul(a, a); }

uint64_t feEqual(const Fe& a, const Fe& b) {
    uint64_t acc = 0;
    for (int i = 0; i < 4; ++i) acc |= a.v[i] ^ b.v[i];
    return 1 ^ ((acc | (0 - acc)) >> 63);
}

uint64_t feIsZero(const Fe& a) { return feEqual(a, feZero()); }

// Фиксированная публичная экспонента p - 2: последовательность операций не зависит от a.
Fe feInvert(const Fe& a) {
    Fe result = feOne();
    for (int i = 255; i >= 0; --i) {
        result = feSquare(result);
        Fe m = feMul(result, a);
        uint64_t bit = (EXP_P_MINUS_2[i / 64] >> (i % 64)) & 1;
        result = feSelect(bit, m, result);
    }
    return result;
}

void feToBytes(const Fe& a, uint8_t out[32]) {
    for (int i = 0; i < 32; ++i) out[31 - i] = (uint8_t)(a.v[i / 8] >> (8 * (i % 8)));
}

static int hexVal(char c) {
    if (c >= '0' && c <= '9') return c - '0';
    if (c >= 'a' && c <= 'f') return c - 'a' + 10;
    if (c >= 'A' && c <= 'F') return c - 'A' + 10;
    return -1;
}

bool feFromHex(const char* hex, Fe* out) {
    uint8_t be[32];
    for (int i = 0; i < 32; ++i) {
        int h = hexVal(hex[2 * i]), l = hexVal(hex[2 * i + 1]);
        if (h < 0 || l < 0) return false;
        be[i] = (uint8_t)(h * 16 + l);
    }
    for (int w = 0; w < 4; ++w) {
        uint64_t x = 0;
        for (int k = 0; k < 8; ++k) x = (x << 8) | be[24 - 8 * w + k];
        out->v[w] = x;
    }
    // значение должно быть < p
    Fe t = *out;
    uint64_t borrow = 0;
    for (int i = 0; i < 4; ++i) {
        u128 d = (u128)t.v[i] - P_LIMBS[i] - borrow;
        borrow = (uint64_t)(d >> 64) & 1;
    }
    return borrow == 1;
}

// ---------------------------------------------------------------- параметры кривой (СТБ 34.101.45, уровень 128)

static CurveParams makeCurve() {
    CurveParams c;
    feFromHex("ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff40", &c.a);   // p - 3
    feFromHex("77ce6c1515f3a8edd2c13aabe4d8fbbe4cf55069978b9253b22e7d6bd69c03f1", &c.b);
    feFromHex("0000000000000000000000000000000000000000000000000000000000000000", &c.Gx);
    feFromHex("6bf7fc3cfb16d69f5ce4c9a351d6835d78913966c408f6521e29cf1804516a93", &c.Gy);
    static const char* QH = "ffffffffffffffffffffffffffffffffd95c8ed60dfb4dfc7e5abf99263d6607";
    for (int i = 0; i < 32; ++i)
        c.q[i] = (uint8_t)(hexVal(QH[2 * i]) * 16 + hexVal(QH[2 * i + 1]));
    return c;
}

const CurveParams& curve() {
    static const CurveParams c = makeCurve();
    return c;
}

// ---------------------------------------------------------------- группа точек

Point pointInfinity() { Point P; P.x = feZero(); P.y = feZero(); P.inf = 1; return P; }
Point pointFromAffine(const Fe& x, const Fe& y) { Point P; P.x = x; P.y = y; P.inf = 0; return P; }
Point pointBase() { return pointFromAffine(curve().Gx, curve().Gy); }

static Point pointSelect(uint64_t bit, const Point& A, const Point& B) {
    Point r;
    r.x = feSelect(bit, A.x, B.x);
    r.y = feSelect(bit, A.y, B.y);
    r.inf = (A.inf & (0 - (bit & 1))) | (B.inf & ~(0 - (bit & 1)));
    return r;
}

Point pointNeg(const Point& P) {
    Point r = P;
    r.y = feNeg(P.y);
    return r;
}

// Полная формула. Все ветви вычисляются всегда, результат выбирается масками.
Point pointAdd(const Point& P, const Point& Q) {
    ++g_addCalls;
    Fe dx = feSub(Q.x, P.x);
    uint64_t sameX = feIsZero(dx);
    uint64_t sameY = feEqual(P.y, Q.y);
    uint64_t dbl = sameX & sameY;                 // P == Q
    uint64_t yZero = feIsZero(P.y);

    // числитель/знаменатель наклона для обоих случаев
    Fe numAdd = feSub(Q.y, P.y);
    Fe denAdd = dx;
    Fe x1sq = feSquare(P.x);
    Fe numDbl = feAdd(feAdd(feAdd(x1sq, x1sq), x1sq), curve().a);   // 3x^2 + a
    Fe denDbl = feAdd(P.y, P.y);                                      // 2y

    Fe num = feSelect(dbl, numDbl, numAdd);
    Fe den = feSelect(dbl, denDbl, denAdd);
    Fe s = feMul(num, feInvert(den));             // feInvert(0) = 0 -> безопасно для вырожденных случаев

    Fe x3 = feSub(feSub(feSquare(s), P.x), Q.x);
    Fe y3 = feSub(feMul(s, feSub(P.x, x3)), P.y);

    Point R;
    R.x = x3; R.y = y3;
    // результат O: P = -Q (x равны, y различны) или удвоение точки с y = 0
    R.inf = (sameX & (1 ^ sameY)) | (dbl & yZero);
    R = pointSelect(P.inf, Q, R);                 // O + Q = Q
    R = pointSelect(Q.inf, P, R);                 // P + O = P (если обе O - результат O, т.к. P = O)
    return R;
}

Point pointDouble(const Point& P) {
    ++g_dblCalls;
    --g_addCalls;                                 // внутри pointAdd увеличится ещё раз
    return pointAdd(P, P);
}

uint64_t pointEqual(const Point& P, const Point& Q) {
    uint64_t bothInf = P.inf & Q.inf;
    uint64_t noInf = (1 ^ P.inf) & (1 ^ Q.inf);
    return bothInf | (noInf & feEqual(P.x, Q.x) & feEqual(P.y, Q.y));
}

bool pointOnCurve(const Point& P) {
    if (P.inf) return true;
    Fe lhs = feSquare(P.y);
    Fe rhs = feAdd(feAdd(feMul(feSquare(P.x), P.x), feMul(curve().a, P.x)), curve().b);
    return feEqual(lhs, rhs) == 1;
}

Point pointMulDoubleAndAdd(const uint8_t scalar[32], const Point& P) {
    Point R = pointInfinity();
    for (int i = 255; i >= 0; --i) {
        R = pointDouble(R);
        if ((scalar[31 - i / 8] >> (i % 8)) & 1) R = pointAdd(R, P);   // ветвление по биту - небезопасно для секретов
    }
    return R;
}
