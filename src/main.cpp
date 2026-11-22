#include "run.hpp"

#include <cstdio>

int main(int argc, char** argv) {
    const int rc = tpms::run({argv + 1, argv + argc}, [](const std::string& text) {
        std::FILE* f = text.rfind("error:", 0) == 0 ? stderr : stdout;
        std::fputs(text.c_str(), f);
        std::fputc('\n', f);
    });
    if (rc == 2) std::fputs("Run 'tpms-forge --help' for usage.\n", stderr);
    return rc;
}
