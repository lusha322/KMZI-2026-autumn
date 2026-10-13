#include "../crypto/aes128.h"
#include "../crypto/gcm.h"
#include <cstdio>
#include <cstring>
#include <cstdint>

static void hex(const uint8_t* p, int n, char* out) {
    static const char* H = "0123456789ABCDEF";
    for (int i = 0; i < n; ++i) {
        out[2*i]   = H[p[i] >> 4];
        out[2*i+1] = H[p[i] & 0xF];
    }
    out[2*n] = 0;
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

int run_demo() {
    std::printf("====================================================================\n");
    std::printf("  ЛАБОРАТОРНАЯ РАБОТА №1. ВАРИАНТ 1\n");
    std::printf("====================================================================\n");
    std::printf("  Алгоритм          : AES-128\n");
    std::printf("  Неприводимый p(x) : x^8 + x^4 + x^3 + x + 1  (0x11B)\n");
    std::printf("  Полином GHASH     : x^128 + x^7 + x^2 + x + 1\n");
    std::printf("  Constant-time     : S-box вычисляется через a^254 на лету\n");
    std::printf("--------------------------------------------------------------------\n\n");

    // [1] Инициализация GF(2^8)
    std::printf("[1] Инициализация калькулятора Галуа GF(2^8), p(x) = 0x11B\n");
    GF256 gf(0x11B);
    std::printf("    Проверка неприводимости p(x) : %s\n",
                gf.is_irreducible() ? "OK" : "FAIL");
    std::printf("    inverse(0x00) = 0x%02X  (соглашение)\n", gf.inverse(0x00));
    std::printf("    inverse(0x01) = 0x%02X\n", gf.inverse(0x01));
    std::printf("    inverse(0x53) = 0x%02X\n", gf.inverse(0x53));
    std::printf("    S-box[0x00] = 0x%02X   S-box[0x01] = 0x%02X\n",
                gf.sbox(0x00), gf.sbox(0x01));
    std::printf("    S-box[0x53] = 0x%02X   S-box[0xFF] = 0x%02X\n",
                gf.sbox(0x53), gf.sbox(0xFF));
    std::printf("\n");

    // [2] Ключ AES-128
    std::printf("[2] Установка 128-битного ключа AES-128\n");
    uint8_t key[16];
    from_hex("000102030405060708090A0B0C0D0E0F", key, 16);
    char kbuf[40];
    hex(key, 16, kbuf);
    std::printf("  K                : %s\n", kbuf);
    AES128 aes(key);
    std::printf("  Развёртка в 44 слова (11 раундовых ключей) выполнена\n\n");

    // [3] GCM init
    std::printf("[3] Инициализация AEAD-режима GCM\n");
    GCM gcm(aes);
    gcm.init();
    std::printf("    Вспомогательный ключ H = E_K(0^128) вычислен\n\n");

    // [4] Входные данные
    std::printf("[4] Входные данные\n");
    uint8_t nonce[12];
    from_hex("CAFEBABEFACEDBADDECAF888", nonce, 12);
    char nbuf[32];
    hex(nonce, 12, nbuf);
    std::printf("  Nonce (96 бит)   : %s\n", nbuf);

    uint8_t aad[8];
    from_hex("FEEDFACEDEADBEEF", aad, 8);
    char abuf[24];
    hex(aad, 8, abuf);
    std::printf("  AAD              : %s\n", abuf);

    const char* pt_hex = "00112233445566778899AABBCCDDEEFF";
    uint8_t pt[16];
    from_hex(pt_hex, pt, 16);
    std::printf("  Открытый текст   : %s\n", pt_hex);
    std::printf("    Длина PT         : 16 байт\n\n");

    // [5] AEAD-шифрование
    std::printf("[5] AEAD-шифрование (GCM)\n");
    std::printf("--------------------------------------------------------------------\n");
    std::printf("  ВЫХОДНЫЕ ДАННЫЕ\n");
    std::printf("--------------------------------------------------------------------\n");

    uint8_t ct[16];
    GCM::Tag tag;
    gcm.encrypt(aad, 8, pt, 16, nonce, ct, tag);

    char cbuf[40], tbuf[40];
    hex(ct, 16, cbuf);
    hex(tag.data(), 16, tbuf);
    std::printf("  Шифротекст (C)   : %s  (16 байт)\n", cbuf);
    std::printf("  Тег целостности  : %s\n", tbuf);
    std::printf("--------------------------------------------------------------------\n\n");

    // [6] Расшифрование
    std::printf("[6] AEAD-расшифрование с корректным тегом\n");
    uint8_t back[16];
    bool ok = gcm.decrypt(aad, 8, ct, 16, nonce, tag, back);
    std::printf("    Верификация тега : %s\n", ok ? "ПРОЙДЕНА" : "ПРОВАЛЕНА");
    char pbuf[40];
    hex(back, 16, pbuf);
    std::printf("  Восстановленный PT: %s\n", pbuf);
    std::printf("    Совпадение с PT  : %s\n\n",
                std::memcmp(pt, back, 16) == 0 ? "ДА" : "НЕТ");

    // [7] Подделка
    std::printf("[7] Демонстрация защиты от подделки\n");
    GCM::Tag forged = tag;
    forged[0] ^= 0x01;
    char fbuf[40];
    hex(forged.data(), 16, fbuf);
    std::printf("    Подменённый тег  : %s\n", fbuf);
    uint8_t dummy[16];
    bool ok_forged = gcm.decrypt(aad, 8, ct, 16, nonce, forged, dummy);
    std::printf("    Результат проверки : %s\n",
                ok_forged ? "FAIL — подделка принята" : "OK — подделка отклонена");
    std::printf("    Открытый текст не выдан: %s\n\n",
                ok_forged ? "НЕТ" : "ДА");

    // [8] Constant-time сравнение
    std::printf("[8] Constant-time сравнение тегов (XOR-аккумулятор)\n");
    bool self_eq = GCM::ct_compare(tag, tag);
    bool diff_eq = GCM::ct_compare(tag, forged);
    std::printf("    ct_compare(tag, tag)          = %s (совпали)\n",
                self_eq ? "true" : "false");
    std::printf("    ct_compare(tag, forged_tag)   = %s (отличаются)\n",
                diff_eq ? "true" : "false");
    std::printf("    Early-exit отсутствует: проверяются все 16 байт\n");

    return 0;
}