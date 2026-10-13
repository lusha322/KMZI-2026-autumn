#include <cstdio>

int run_test_gf256();
int run_test_aes();
int run_test_gcm();
int run_demo();

int main() {
    int failed = 0;
    failed += run_test_gf256();
    failed += run_test_aes();
    failed += run_test_gcm();

    std::printf("=================================================\n");
    std::printf(" ИТОГО: %s\n", (failed == 0) ? "все тесты пройдены" : "есть ошибки");
    std::printf("=================================================\n\n");

    run_demo();

    std::printf("\n=================================================\n");
    std::printf(" РЕЗУЛЬТАТ: %s\n", (failed == 0) ? "OK" : "FAILED");
    std::printf("=================================================\n");
    return failed == 0 ? 0 : 1;
}