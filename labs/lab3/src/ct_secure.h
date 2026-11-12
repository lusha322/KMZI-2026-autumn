#ifndef CT_SECURE_H
#define CT_SECURE_H

// Модуль 3. Constant-time защита (вариант 1): лестница Монтгомери для скалярного умножения,
// условный обмен точек без ветвлений, безопасное сравнение.

#include <cstdint>
#include "ec_algebra.h"

// Условный обмен точек: если b = 1, A и B меняются местами (маска, без if)
void ctCswapPoint(uint64_t b, Point* A, Point* B);

// Сравнение байтовых массивов, XOR-аккумулятор без раннего выхода. 1 при равенстве.
uint64_t ctCompareBytes(const uint8_t* a, const uint8_t* b, size_t len);

// Сравнение точек (SecureCompare)
uint64_t SecureCompare(const Point& A, const Point& B);

// Лестница Монтгомери: на каждый из 256 бит скаляра ровно одно сложение и одно удвоение.
// scalar - 32 байта big-endian.
Point SecureMultiply(const uint8_t scalar[32], const Point& P);

#endif
