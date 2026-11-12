Министерство образования Республики Беларусь

Учреждение образования
«Брестский Государственный технический университет»

Кафедра ИИТ

**Лабораторная работа №3**

По дисциплине «КМЗИ»

Тема: «Арифметика точечных групп эллиптических кривых, протокол ECDH и алгоритмы ЭЦП»

> **Выполнил:** Студент 3 курса, группы ИИ-27/24 Глаголева Е.В.
>
> **Проверил:** Хацкевич А.С.

Брест 2026

---

**Цель работы:** программная реализация защищённого математического ядра для работы с точечной группой эллиптической кривой в форме Вейерштрасса над $\mathbb{GF}(p)$, протокола ECDH и ЭЦП по СТБ 34.101.45 с защитой скалярного умножения от атак по побочным каналам.

**Индивидуальный вариант — №1:**

| Параметр | Значение |
|---|---|
| Алгоритм ЭЦП | СТБ 34.101.45 (РБ) |
| Тип кривой / координаты | Вейерштрасс / аффинные |
| Задача Constant-time | Лестница Монтгомери (Montgomery Ladder) для скалярного умножения |

**Язык и средства:** C++17. Собственная реализация арифметики $\mathbb{GF}(p)$ (4×64 бит, без динамической памяти), SHA-256 и всех протоколов. Библиотека GMP используется **исключительно** для арифметики больших чисел по модулю порядка $q$ (скаляры подписи) и в тестах как эталон для сверки арифметики поля. Парсинг PEM/DER/ASN.1 не реализуется; все параметры, сообщения и хэши заданы hex-строками в коде.

**Параметры кривой** (bign-curve256v1, СТБ 34.101.45, уровень 128):

- $p = 2^{256} - 189$ = `ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff43`
- $a = p - 3$ = `ffff…ff40`
- $b$ = `77ce6c1515f3a8edd2c13aabe4d8fbbe4cf55069978b9253b22e7d6bd69c03f1`
- $q$ = `ffffffffffffffffffffffffffffffffd95c8ed60dfb4dfc7e5abf99263d6607`
- $G = (0,\ $ `6bf7fc3cfb16d69f5ce4c9a351d6835d78913966c408f6521e29cf1804516a93` $)$

Корректность параметров подтверждается программно: $G$ лежит на кривой и $q\cdot G=\mathcal{O}$ (см. раздел 4).

---

## 1. Математическая часть

### 1.1. Закон сложения в аффинных координатах

Кривая: $y^2 = x^3 + ax + b$, $4a^3+27b^2\not\equiv 0 \pmod p$. Группа — аффинные точки и $\mathcal{O}$ (нейтральный элемент). Три точки кривой, лежащие на одной прямой, в сумме дают $\mathcal{O}$: $P_1+P_2+P_3=\mathcal{O}$.

**Сложение $P(x_1,y_1)\ne \pm Q(x_2,y_2)$.** Прямая через $P$ и $Q$: $y = s x + c$, $s=\dfrac{y_2-y_1}{x_2-x_1}$. Подстановка в уравнение кривой даёт кубическое уравнение
$$x^3 - s^2x^2 + (a-2sc)x + (b-c^2)=0,$$
сумма корней которого равна $s^2$ (теорема Виета). Корни — абсциссы $x_1, x_2$ и третьей точки пересечения $x_3'$, значит $x_1+x_2+x_3'=s^2$, откуда
$$x_3 = s^2 - x_1 - x_2.$$
Третья точка пересечения $(x_3, s x_3 + c)$, а искомая сумма — её отражение относительно оси $x$: $y_3 = -(s x_3 + c)$. Так как $c = y_1 - s x_1$:
$$y_3 = s(x_1 - x_3) - y_1.$$

**Удвоение $P(x_1,y_1)$, $y_1\neq 0$.** Вместо секущей берётся касательная. Дифференцируя $y^2=x^3+ax+b$: $2y\,dy=(3x^2+a)\,dx$, откуда
$$s=\frac{3x_1^2+a}{2y_1}.$$
Дальше рассуждение то же: касательная пересекает кривую в $P$ (дважды) и в третьей точке, поэтому $2x_1+x_3'=s^2$,
$$x_3=s^2-2x_1,\qquad y_3=s(x_1-x_3)-y_1.$$

**Особые случаи.** $P+\mathcal{O}=P$; $P+(-P)=\mathcal{O}$ при $x_1=x_2,\ y_1=-y_2$ (вертикальная прямая); $2P=\mathcal{O}$ при $y_1=0$.

### 1.2. Единая (полная) реализация формулы без ветвлений

Аффинные формулы требуют ветвления: «$P=Q$?», «$P=-Q$?», «$P=\mathcal{O}$?». Если эти ветвления зависят от секрета, они образуют канал утечки. Поэтому `pointAdd` вычисляет **всё одновременно** и выбирает результат масками:

1. $dx=x_2-x_1$; $\text{sameX}=[dx=0]$; $\text{sameY}=[y_1=y_2]$; $\text{dbl}=\text{sameX}\wedge\text{sameY}$.
2. Числитель и знаменатель наклона выбираются маской: $\text{num}=\text{dbl}\ ?\ 3x_1^2+a : y_2-y_1$, $\text{den}=\text{dbl}\ ?\ 2y_1 : dx$.
3. $s=\text{num}\cdot \text{den}^{-1}$, где обращение — по Ферма $\text{den}^{p-2}$ (при $\text{den}=0$ получается $0$ без исключений).
4. $x_3, y_3$ по формулам 1.1.
5. Признак результата $\mathcal{O}$: $(\text{sameX}\wedge\neg\text{sameY})\vee(\text{dbl}\wedge[y_1=0])$.
6. Если $P=\mathcal{O}$, результат $Q$; если $Q=\mathcal{O}$ — результат $P$ (масочный выбор).

Удвоение — это $\text{pointAdd}(P,P)$; отдельной ветки в коде нет.

### 1.3. ЭЦП СТБ 34.101.45 (по заданию)

Пусть $d\in(0,q)$, $Q=dG$, $e=\text{Hash}(M)\bmod q$, $e\ne 0$.

**Подпись.** $k\in(0,q)$; $C=kG$; $r=x_C \bmod q$ ($r\neq 0$); $s=(r d + k e)\bmod q$ ($s\ne 0$). Подпись — $(r,s)$.

**Проверка.** $0<r,s<q$; $v=e^{-1}\bmod q$; $z_1=sv$, $z_2=-rv \pmod q$; $C'=z_1G+z_2Q$; $R=x_{C'}\bmod q$; подпись верна при $R=r$.

**Доказательство корректности ($R=r$).** Подставим $s = rd+ke$ и $Q = dG$:
$$z_1 G+z_2 Q = (sv)G+(-rv)(dG) = v\,(s - rd)\,G = e^{-1}(rd+ke-rd)\,G = e^{-1}\cdot ke\cdot G = kG = C.$$
Скаляры приводятся по модулю $q$, а $qG=\mathcal{O}$, поэтому приведение не меняет точку. Следовательно $C'=C$, $x_{C'}=x_C$ и $R = x_{C'} \bmod q = x_C \bmod q = r$. Условие $e\not\equiv 0$ необходимо для существования $e^{-1}$, поэтому подпись при $e\equiv 0$ не вырабатывается.

*Замечание о хэш-функции.* Стандарт использует belt-hash; в работе (по заданию «работа с хэш-функциями») применена собственная реализация SHA-256, хэши тестовых сообщений дополнительно заданы hex-строками и сверяются с вычисленными. Схема подписи от выбора 256-битной хэш-функции не зависит.

### 1.4. ECDH

$Q_A=d_AG$, $Q_B=d_BG$; $K=d_AQ_B=d_BQ_A=(d_Ad_B)G$; общий секрет — координата $x_K$. Перед вычислением проверяется $0<d<q$, $Q_{peer}\ne\mathcal{O}$ и $Q_{peer}\in E$ (защита от атаки неверной точкой).

---

## 2. Архитектура

| Файл | Модуль | Содержимое |
|---|---|---|
| `ec_algebra.h/.cpp` | 1. Алгебра ЭК | `Fe` (4×uint64), арифметика GF(p) по $p=2^{256}-189$, `Point` (аффинная, флаг `inf`), `pointAdd`, `pointDouble`, `pointEqual`, `pointOnCurve`, публичный Double-and-Add |
| `ct_secure.h/.cpp` | 3. Constant-time | `ctCswapPoint`, `ctCompareBytes`, `SecureCompare`, **`SecureMultiply` (лестница Монтгомери)** |
| `protocols.h/.cpp` | 2. ECDH и ЭЦП | SHA-256, `keyGenPublic`, `ecdhSharedSecret`, `stbSign`, `stbVerify`, арифметика по модулю $q$ (GMP) |
| `main.cpp` | — | тесты и демонстрация |

Все структуры — непрерывные блоки фиксированного размера на стеке/в статической памяти; `std::vector`, `new`, `malloc` и вложенные контейнеры не используются. Секретные операции (генерация открытого ключа, $kG$ при подписи, ECDH) идут через `SecureMultiply`; проверка подписи работает с публичными данными и использует обычный Double-and-Add.

Приведение по модулю в поле выполняется без цикла `while (x >= p)`: $2^{256}\equiv 189\pmod p$, после чего одно условное вычитание $p$ маской (`condSubP`).

---

## 3. Анализ безопасности (Constant-time)

### 3.1. Проблема классического Double-and-Add

```cpp
for (i = 255; i >= 0; --i) {
    R = pointDouble(R);
    if ((scalar[31 - i/8] >> (i%8)) & 1) R = pointAdd(R, P);   // ветвление по биту секрета
}
```
Сложение выполняется только для единичных битов. Время и энергопотребление раскрывают вес Хэмминга и положение единиц в $d$ (timing / SPA). Эта функция в коде помечена «ТОЛЬКО для публичных скаляров» и применяется лишь при проверке подписи и как эталон в тестах.

### 3.2. Лестница Монтгомери

Поддерживаются $R_0, R_1$ с инвариантом $R_1-R_0=P$. По биту $b$ секрета: при $b=0$: $R_1=R_0+R_1,\ R_0=2R_0$; при $b=1$: $R_0=R_0+R_1,\ R_1=2R_1$. Чтобы убрать `if`, оба случая сведены к одному коду с условным обменом:

```cpp
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
```
Проверка эквивалентности: при $b=1$ после обмена $R_0'=R_1,\ R_1'=R_0$; затем $R_1'=R_0'+R_1'=R_0+R_1$ и $R_0'=2R_1$; после обратного обмена $R_0=R_0+R_1$, $R_1=2R_1$ — ровно ветка $b=1$. При $b=0$ обмены тождественны — ветка $b=0$. Цикл всегда выполняет все 256 итераций (старшие нулевые биты обрабатываются корректно благодаря $R_0=\mathcal{O}$ и полной формуле), в каждой — ровно одно сложение и одно удвоение.

### 3.3. Отсутствие ветвлений, зависящих от секрета

1. **Обмен точек (CSWAP)** — маска вместо `if`:
```cpp
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
```
2. **Сложение/удвоение** — единая функция, случаи $P=Q$, $P=-Q$, $P=\mathcal{O}$, $Q=\mathcal{O}$ различаются масками (`feSelect`, `pointSelect`), а не `if` (раздел 1.2). Промежуточные точки лестницы зависят от секретного префикса $d$, поэтому именно полнота формулы без ветвлений здесь существенна: на $R_0=\mathcal{O}$ (первые итерации при ведущих нулях) ветвь выбирается по данным, но не по коду.
3. **Обращение** — $a^{p-2}$ по фиксированной публичной экспоненте, 256 квадратов и 256 умножений независимо от $a$; умножение выполняется всегда, результат выбирается маской:
```cpp
for (int i = 255; i >= 0; --i) {
    result = feSquare(result);
    Fe m = feMul(result, a);
    result = feSelect((EXP_P_MINUS_2[i/64] >> (i%64)) & 1, m, result);
}
```
4. **Приведение по модулю** — `condSubP`: вычитание $p$ всегда, результат выбирается по знаку заёма (`feSelect(borrow, a, t)`), цикла `while` нет.
5. **Сравнения** (`feEqual`, `pointEqual`, `ctCompareBytes`) — XOR-аккумулятор без раннего выхода.
6. **Доступ к памяти** — индексы массивов зависят только от публичных счётчиков цикла, таблиц, индексируемых секретом, нет.

### 3.4. Эмпирическая проверка

*Профиль операций (SPA-модель).* Счётчики вызовов (см. `run_output.txt`):

| Скаляр | Лестница: сложений / удвоений | Double-and-Add: сложений / удвоений |
|---|---|---|
| $d=1$ | 256 / 256 | 1 / 256 |
| $d=2^{255}+1$ | 256 / 256 | — |
| $d=q-1$ | 256 / 256 | 201 / 256 |
| $d=\texttt{0xAAAA…}$ | 256 / 256 | — |

У лестницы последовательность операций не зависит от $d$, у Double-and-Add число и расположение сложений раскрывает биты $d$.

*Время (индикативно, Apple Silicon, `-O2`, среднее по 15 запускам).*

| | $d=1$ | $d=q-1$ | отношение |
|---|---|---|---|
| Лестница Монтгомери | 5.28 мс | 5.28 мс | 1.001 |
| Double-and-Add | 2.66 мс | 4.73 мс | 1.78 |

Время лестницы практически не зависит от веса ключа; время Double-and-Add — зависит почти вдвое.

### 3.5. Ограничения

- Аффинные формулы требуют обращения на каждый шаг (по заданию варианта), поэтому лестница заметно медленнее проективных реализаций (≈5 мс на умножение).
- Арифметика по модулю $q$ при выработке подписи ($s=rd+ke \bmod q$) выполняется GMP (`mpz_mul`, `mpz_mod`), он не гарантирует constant-time. Это допустимо по условию работы (GMP — только для больших чисел), но в промышленной реализации потребовалась бы собственная constant-time арифметика по $q$.
- Гарантии на уровне исходного кода: компилятор теоретически может заменить маски ветвлениями; для строгого вывода нужен анализ ассемблера (в `-O2` для AArch64 `feSelect`/`ctCswapPoint` компилируются в `csel`/логические инструкции).
- Генерация $k$ в `randomScalar` использует rejection sampling по `/dev/urandom`; отбракованные значения не связаны с итоговым ключом.
- В качестве хэша использована SHA-256, а не belt-hash.

---

## 4. Результаты тестирования

Полный вывод программы — `doc/run_output.txt`. Все проверки пройдены (`ВСЕ ТЕСТЫ ПРОЙДЕНЫ`).

**Корректность поля и группы.** `feMul/feAdd/feSub` совпадают с GMP на 300 случайных и граничных парах; $a\cdot a^{-1}=1$; $G\in E$; $q\cdot G=\mathcal{O}$ (подтверждает параметры $p, a, b, G, q$); $(q-1)G=-G$; коммутативность и ассоциативность сложения; $P+(-P)=\mathcal{O}$; $P+P=2P$; лестница совпадает с Double-and-Add на 20 случайных скалярах; $(a+b)G=aG+bG$.

**ECDH** (тестовые ключи в hex):

```
dA = 1f2e3d4c5b6a79881726354453627180918a9bacbdcedfe0f1021324354657a8
QA.x = 8a5808285b33a94f04a0aae018d7450a79cadbcdb7937709e3b86f84bc39773b
QA.y = ae4a4dc6e78730df0d630dbbdb00ce434772f91edca9ba528439e01b715914a3
dB = 0a1b2c3d4e5f60718293a4b5c6d7e8f90123456789abcdef0fedcba987654321
QB.x = 0979523bd51cc56af8c6056b517ccb25d74995af4bf8a209894085acd10b882e
QB.y = 3639b04367c9e895f03e1acb2c94e071fc685a646ebc3351570b5b6563bbcc76
Shared (Алиса) = 500b69c1c47fc6a60490419972205ae33368b2e53d1f7a89e5ab0ed5de1f3634
Shared (Боб)   = 500b69c1c47fc6a60490419972205ae33368b2e53d1f7a89e5ab0ed5de1f3634
```
Секреты совпадают. Дополнительно: 5 запусков с эфемерными ключами из `/dev/urandom`; точка не на кривой отвергается.

**ЭЦП** — три сообщения, подпись ключом $d_A$ (пример для M1 = `616263`):

```
e = ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad
k = 3c4d5e6f708192a3b4c5d6e7f8091a2b3c4d5e6f708192a3b4c5d6e7f8091a2b
r = cea0130fbdf89458729d5c9866977e57c88abd7dde16cb10f136fe21df30fc12
s = 4eab6f1952c9882473f3e7131fd155d9a320ca3b6f0ddd0feac0e96c6d8fefcc
```
Для каждого сообщения: корректная подпись — **ACCEPT**; искажение одного байта $s$ (младший и старший бит) — **REJECT**; искажение $r$, хэша, чужой ключ, $r=0$ или $s=0$ — **REJECT**. Ещё 5 циклов со случайными $d,k,e$ прошли успешно.

**Независимая проверка.** Значения $Q_A$, общего секрета, $r$ и $s$ дополнительно воспроизведены отдельным скриптом на Python (простая аффинная реализация) — результаты совпали с выводом программы.

**Сборка и запуск:**
```
cmake -S src -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build && ./build/main
```
Требуется GMP (`brew install gmp`).

---

## 5. Код программы

### `ec_algebra.h`

```cpp
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
```

### `ec_algebra.cpp`

```cpp
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
```

### `ct_secure.h`

```cpp
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
```

### `ct_secure.cpp`

```cpp
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
```

### `protocols.h`

```cpp
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
```

### `protocols.cpp`

```cpp
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
```

### `main.cpp`

```cpp
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
```

### `CMakeLists.txt`

```cmake
cmake_minimum_required(VERSION 3.10)
project(lab3_stb34101_45_variant1 CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# GMP используется только для арифметики больших чисел по модулю q (скаляры подписи)
find_path(GMP_INCLUDE_DIR gmp.h HINTS /opt/homebrew/include /usr/local/include /usr/include)
find_library(GMP_LIBRARY gmp HINTS /opt/homebrew/lib /usr/local/lib /usr/lib)
include_directories(${GMP_INCLUDE_DIR})

add_executable(main
    ./ec_algebra.cpp
    ./ct_secure.cpp
    ./protocols.cpp
    ./main.cpp
)
target_link_libraries(main ${GMP_LIBRARY})
```
