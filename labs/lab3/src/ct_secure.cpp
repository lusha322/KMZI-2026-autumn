#include "ct_secure.h"

void ctCswapPoint(uint64_t b, Point* A, Point* B) {
    uint64_t mask = 0 - (b & 1);
    for (int i = 0; i < 4; ++i) {
        uint64_t d = mask & (A->x.v[i] ^ B->x.v[i]);
        A->x.v[i] ^= d; B->x.v[i] ^= d;
        d = mask & (A->y.v[i] ^ B->y.v[i]);
        A->y.v[i] ^= d; B->y.v[i] ^= d;
    }
    uint64_t d = mask & (A->inf ^ B->inf);
    A->inf ^= d; B->inf ^= d;
}

uint64_t ctCompareBytes(const uint8_t* a, const uint8_t* b, size_t len) {
    uint64_t acc = 0;
    for (size_t i = 0; i < len; ++i) acc |= (uint64_t)(a[i] ^ b[i]);
    return 1 ^ ((acc | (0 - acc)) >> 63);
}

uint64_t SecureCompare(const Point& A, const Point& B) { return pointEqual(A, B); }

// Инвариант: R1 - R0 = P. Бит 0: R1 = R0 + R1, R0 = 2R0. Бит 1: R0 = R0 + R1, R1 = 2R1.
// Оба случая сводятся к одному коду: swap(R0,R1) по биту; S = R0 + R1; R0 = 2R0; swap обратно.
Point SecureMultiply(const uint8_t scalar[32], const Point& P) {
    Point R0 = pointInfinity();
    Point R1 = P;
    for (int i = 255; i >= 0; --i) {
        uint64_t bit = (scalar[31 - i / 8] >> (i % 8)) & 1;
        ctCswapPoint(bit, &R0, &R1);
        R1 = pointAdd(R0, R1);
        R0 = pointDouble(R0);
        ctCswapPoint(bit, &R0, &R1);
    }
    return R0;
}
