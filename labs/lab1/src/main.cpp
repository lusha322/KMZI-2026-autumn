#include <cstdint>
#include <cstdio>
#include <cstring>
#include "galois_field.hpp"
#include "aes128.hpp"
#include "gcm.hpp"

static int passed = 0;
static int failed = 0;

static void check(const char* name, bool ok) {
    std::printf("[%s] %s\n", ok ? "OK  " : "FAIL", name);
    if (ok) {
        ++passed;
    } else {
        ++failed;
    }
}

static int hex_digit(char c) {
    if (c >= '0' && c <= '9') return c - '0';
    if (c >= 'a' && c <= 'f') return c - 'a' + 10;
    return c - 'A' + 10;
}

static size_t from_hex(const char* hex, uint8_t* out) {
    size_t n = std::strlen(hex) / 2;
    for (size_t i = 0; i < n; ++i) {
        out[i] = static_cast<uint8_t>(hex_digit(hex[2 * i]) * 16 + hex_digit(hex[2 * i + 1]));
    }
    return n;
}

static bool equals_hex(const uint8_t* data, size_t len, const char* hex) {
    uint8_t expected[256];
    size_t n = from_hex(hex, expected);
    return n == len && std::memcmp(data, expected, len) == 0;
}

static void print_hex(const char* label, const uint8_t* data, size_t len) {
    std::printf("%s", label);
    for (size_t i = 0; i < len; ++i) {
        std::printf("%02x", data[i]);
    }
    std::printf("\n");
}

static void test_field() {
    const uint16_t p = 0x11B;
    check("multiply(0x57, 0x83) = 0xC1", GaloisField::multiply(0x57, 0x83, p) == 0xC1);
    check("multiply(0x80, 0x02) = 0x1B (overflow reduction)", GaloisField::multiply(0x80, 0x02, p) == 0x1B);
    check("multiply(0xFF, 0xFF) = 0x13", GaloisField::multiply(0xFF, 0xFF, p) == 0x13);
    check("add(0x57, 0x83) = 0xD4", GaloisField::add(0x57, 0x83) == 0xD4);
    check("inverse(0x00) = 0x00", GaloisField::inverse(0x00, p) == 0x00);
    check("inverse(0x01) = 0x01", GaloisField::inverse(0x01, p) == 0x01);
    check("inverse(0x53) = 0xCA", GaloisField::inverse(0x53, p) == 0xCA);
    bool all_ok = true;
    for (int a = 1; a < 256; ++a) {
        uint8_t inv = GaloisField::inverse(static_cast<uint8_t>(a), p);
        if (GaloisField::multiply(static_cast<uint8_t>(a), inv, p) != 1) {
            all_ok = false;
        }
    }
    check("a * a^-1 = 1 for all a from 1 to 255", all_ok);
}

static void test_sbox() {
    GaloisField gf(0x11B);
    uint8_t sbox[256];
    gf.generate_sbox(sbox);
    check("S-box[0x00] = 0x63", sbox[0x00] == 0x63);
    check("S-box[0x01] = 0x7C", sbox[0x01] == 0x7C);
    check("S-box[0x53] = 0xED", sbox[0x53] == 0xED);
    check("S-box[0xFF] = 0x16", sbox[0xFF] == 0x16);
    bool on_the_fly = true;
    bool permutation = true;
    bool seen[256] = {false};
    for (int i = 0; i < 256; ++i) {
        if (sbox[i] != gf.sub_byte(static_cast<uint8_t>(i))) {
            on_the_fly = false;
        }
        if (seen[sbox[i]]) {
            permutation = false;
        }
        seen[sbox[i]] = true;
    }
    check("table S-box equals on-the-fly a^254 computation", on_the_fly);
    check("S-box is a permutation of 0..255", permutation);
}

static void test_aes() {
    uint8_t key[16], plain[16], cipher[16];
    from_hex("000102030405060708090a0b0c0d0e0f", key);
    from_hex("00112233445566778899aabbccddeeff", plain);
    AES128 aes(key);
    aes.encrypt_block(plain, cipher);
    print_hex("AES-128 ciphertext: ", cipher, 16);
    check("AES-128 FIPS-197 Appendix C.1", equals_hex(cipher, 16, "69c4e0d86a7b0430d8cdb78070b4c55a"));
}

static void run_gcm_case(const char* name, const char* key_hex, const char* iv_hex,
                         const char* aad_hex, const char* plain_hex,
                         const char* cipher_hex, const char* tag_hex) {
    uint8_t key[16], iv[12], aad[64], plain[128], cipher[128], tag[16], back[128];
    from_hex(key_hex, key);
    from_hex(iv_hex, iv);
    size_t aad_len = from_hex(aad_hex, aad);
    size_t len = from_hex(plain_hex, plain);
    GCM gcm(key);
    gcm.encrypt(iv, aad, aad_len, plain, len, cipher, tag);
    char label[96];
    std::snprintf(label, sizeof(label), "%s: encryption + tag", name);
    check(label, equals_hex(cipher, len, cipher_hex) && equals_hex(tag, 16, tag_hex));
    bool valid = gcm.decrypt(iv, aad, aad_len, cipher, len, tag, back);
    std::snprintf(label, sizeof(label), "%s: decryption restores plaintext", name);
    check(label, valid && std::memcmp(back, plain, len) == 0);
}

static void test_gcm() {
    const char* k1 = "00000000000000000000000000000000";
    const char* iv1 = "000000000000000000000000";
    run_gcm_case("GCM case 1", k1, iv1, "", "", "", "58e2fccefa7e3061367f1d57a4e7455a");
    run_gcm_case("GCM case 2", k1, iv1, "", "00000000000000000000000000000000",
                 "0388dace60b6a392f328c2b971b2fe78", "ab6e47d42cec13bdf53a67b21257bddf");
    const char* k3 = "feffe9928665731c6d6a8f9467308308";
    const char* iv3 = "cafebabefacedbaddecaf888";
    run_gcm_case("GCM case 3", k3, iv3, "",
                 "d9313225f88406e5a55909c5aff5269a86a7a9531534f7da2e4c303d8a318a72"
                 "1c3c0c95956809532fcf0e2449a6b525b16aedf5aa0de657ba637b391aafd255",
                 "42831ec2217774244b7221b784d0d49ce3aa212f2c02a4e035c17e2329aca12e"
                 "21d514b25466931c7d8f6a5aac84aa051ba30b396a0aac973d58e091473f5985",
                 "4d5c2af327cd64a62cf35abd2ba6fab4");
    run_gcm_case("GCM case 4", k3, iv3, "feedfacedeadbeeffeedfacedeadbeefabaddad2",
                 "d9313225f88406e5a55909c5aff5269a86a7a9531534f7da2e4c303d8a318a72"
                 "1c3c0c95956809532fcf0e2449a6b525b16aedf5aa0de657ba637b39",
                 "42831ec2217774244b7221b784d0d49ce3aa212f2c02a4e035c17e2329aca12e"
                 "21d514b25466931c7d8f6a5aac84aa051ba30b396a0aac973d58e091",
                 "5bc94fbc3221a5db94fae95ae7121a47");
}

static void test_tamper() {
    uint8_t key[16], iv[12], aad[20], plain[60], cipher[60], tag[16], back[60];
    from_hex("feffe9928665731c6d6a8f9467308308", key);
    from_hex("cafebabefacedbaddecaf888", iv);
    from_hex("feedfacedeadbeeffeedfacedeadbeefabaddad2", aad);
    for (int i = 0; i < 60; ++i) {
        plain[i] = static_cast<uint8_t>(i);
    }
    GCM gcm(key);
    gcm.encrypt(iv, aad, 20, plain, 60, cipher, tag);

    uint8_t bad_tag[16];
    std::memcpy(bad_tag, tag, 16);
    bad_tag[15] ^= 0x01;
    check("tampered last tag byte is rejected", !gcm.decrypt(iv, aad, 20, cipher, 60, bad_tag, back));

    std::memcpy(bad_tag, tag, 16);
    bad_tag[0] ^= 0x80;
    check("tampered first tag byte is rejected", !gcm.decrypt(iv, aad, 20, cipher, 60, bad_tag, back));

    uint8_t bad_cipher[60];
    std::memcpy(bad_cipher, cipher, 60);
    bad_cipher[30] ^= 0x01;
    check("tampered ciphertext is rejected", !gcm.decrypt(iv, aad, 20, bad_cipher, 60, tag, back));

    uint8_t bad_aad[20];
    std::memcpy(bad_aad, aad, 20);
    bad_aad[5] ^= 0x01;
    check("tampered AAD is rejected", !gcm.decrypt(iv, bad_aad, 20, cipher, 60, tag, back));

    bool zeroed = true;
    gcm.decrypt(iv, aad, 20, cipher, 60, bad_tag, back);
    for (int i = 0; i < 60; ++i) {
        if (back[i] != 0) {
            zeroed = false;
        }
    }
    check("plaintext is wiped when the tag is invalid", zeroed);
    check("valid message is accepted", gcm.decrypt(iv, aad, 20, cipher, 60, tag, back));
}

int main() {
    test_field();
    test_sbox();
    test_aes();
    test_gcm();
    test_tamper();
    std::printf("\nPassed: %d, failed: %d\n", passed, failed);
    return failed == 0 ? 0 : 1;
}
