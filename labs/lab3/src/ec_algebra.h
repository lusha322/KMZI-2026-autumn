#ifndef EC_ALGEBRA_H
#define EC_ALGEBRA_H

// Модуль 1. Алгебра эллиптической кривой в короткой форме Вейерштрасса y^2 = x^3 + a x + b над GF(p)
// (СТБ 34.101.45, bign-curve256v1, p = 2^256 - 189). Аффинные координаты.
// Элемент поля - 4 беззнаковых 64-битных слова (little-endian) в каноническом виде [0, p).
// Динамическая память не используется, все структуры - непрерывные блоки фиксированного размера.

#include <cstdint>
#include <cstddef>

struct Fe {
    uint64_t v[4];
};

// Аффинная точка. inf = 1 для бесконечно удалённой точки O (тогда x, y не имеют значения).
struct Point {
    Fe x, y;
    uint64_t inf;
};

// ---- арифметика поля GF(p) (без ветвлений по значениям) ----
Fe feZero();
Fe feOne();
Fe feFromU64(uint64_t x);
Fe feAdd(const Fe& a, const Fe& b);
Fe feSub(const Fe& a, const Fe& b);
Fe feNeg(const Fe& a);
Fe feMul(const Fe& a, const Fe& b);
Fe feSquare(const Fe& a);
Fe feSelect(uint64_t bit, const Fe& a, const Fe& b);   // bit ? a : b, без ветвления
uint64_t feEqual(const Fe& a, const Fe& b);            // 1 если равны, XOR-аккумулятор без раннего выхода
uint64_t feIsZero(const Fe& a);
Fe feInvert(const Fe& a);                              // a^(p-2), малая теорема Ферма (feInvert(0) = 0)
void feToBytes(const Fe& a, uint8_t out[32]);          // big-endian
bool feFromHex(const char* hex64, Fe* out);            // big-endian hex, 64 символа

// ---- параметры кривой (заданы hex-строками) ----
struct CurveParams {
    Fe a, b;
    Fe Gx, Gy;
    uint8_t q[32];           // порядок подгруппы (big-endian)
};
const CurveParams& curve();

// ---- группа точек ----
Point pointInfinity();
Point pointFromAffine(const Fe& x, const Fe& y);
Point pointBase();
Point pointNeg(const Point& P);
// Полная (branchless) формула сложения: учитывает P+Q, P+P, P+(-P), P+O, O+Q без if по значениям.
Point pointAdd(const Point& P, const Point& Q);
Point pointDouble(const Point& P);                    // = pointAdd(P, P)
uint64_t pointEqual(const Point& P, const Point& Q);  // 1 при равенстве, без раннего выхода
bool pointOnCurve(const Point& P);

// Классическое скалярное умножение Double-and-Add (с ветвлением по биту, ТОЛЬКО для публичных скаляров)
// Скаляр - 32 байта big-endian.
Point pointMulDoubleAndAdd(const uint8_t scalar[32], const Point& P);

// Счётчики вызовов (для демонстрации одинаковой последовательности операций, SPA)
extern uint64_t g_addCalls;
extern uint64_t g_dblCalls;

#endif
