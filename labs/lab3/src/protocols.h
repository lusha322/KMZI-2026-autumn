#ifndef PROTOCOLS_H
#define PROTOCOLS_H

// Модуль 2. SHA-256, ECDH и ЭЦП СТБ 34.101.45 (схема Шнорра-типа: s = r*d + k*e mod q).
// Все скаляры, хэши и координаты - 32 байта big-endian.

#include <cstdint>
#include <cstddef>
#include "ec_algebra.h"

// ---- SHA-256 (собственная реализация) ----
void sha256(const uint8_t* data, size_t len, uint8_t out[32]);

// ---- скаляры ----
// Проверка 0 < x < q
bool scalarInRange(const uint8_t x[32]);
// Случайный скаляр 0 < x < q (источник - /dev/urandom, rejection sampling)
bool randomScalar(uint8_t out[32]);

// ---- ключи ----
bool keyGenPublic(const uint8_t d[32], Point* Q);               // Q = d*G (SecureMultiply)

// ---- ECDH ----
// K = d * Qpeer, в качестве общего секрета - координата x. false, если Qpeer не на кривой / равна O.
bool ecdhSharedSecret(const uint8_t d[32], const Point& peer, uint8_t shared[32]);

// ---- ЭЦП СТБ 34.101.45 ----
// e - хэш сообщения (32 байта), k - эфемерный ключ. false, если r = 0, s = 0 или e = 0 mod q.
bool stbSign(const uint8_t d[32], const uint8_t hash[32], const uint8_t k[32],
             uint8_t r[32], uint8_t s[32]);
bool stbVerify(const Point& Q, const uint8_t hash[32], const uint8_t r[32], const uint8_t s[32]);

#endif
