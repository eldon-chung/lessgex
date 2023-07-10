#define PCRE2_CODE_UNIT_WIDTH 8

#include "lessgex.h"

#include <pcre2.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

#include <chrono>
#include <optional>
#include <string>

#include <iomanip>

#include <iostream>

// you need to promise to null terminate these
std::optional<std::chrono::duration<double, std::nano>>
pcre_bench(std::string const &in_pattern, std::string const &in_subject) {
    pcre2_code *re;
    PCRE2_SPTR pattern =
        (PCRE2_SPTR)in_pattern.c_str(); /* PCRE2_SPTR is a pointer to
                            unsigned code units of */

    PCRE2_SPTR subject =
        (PCRE2_SPTR)in_subject.c_str(); /* the appropriate width (in this
                                case, 8 bits). */
    // PCRE2_SPTR name_table;

    // uint32_t option_bits;
    // uint32_t namecount;
    // uint32_t name_entry_size;
    // uint32_t newline;

    PCRE2_SIZE erroroffset;
    PCRE2_SIZE *ovector;
    PCRE2_SIZE subject_length = (PCRE2_SIZE)in_subject.length();

    pcre2_match_data *match_data;

    int errornumber;
    re = pcre2_compile(
        pattern,               /* the pattern */
        PCRE2_ZERO_TERMINATED, /* indicates pattern is zero-terminated */
        0,                     /* default options */
        &errornumber,          /* for error number */
        &erroroffset,          /* for error offset */
        NULL);                 /* use default compile context */

    /* Compilation failed: print the error message and exit. */

    if (re == NULL) {
        PCRE2_UCHAR buffer[256];
        pcre2_get_error_message(errornumber, buffer, sizeof(buffer));
        printf("PCRE2 compilation failed at offset %d: %s\n", (int)erroroffset,
               buffer);
        return std::nullopt;
    }

    match_data = pcre2_match_data_create_from_pattern(re, NULL);

    // start timer just before matcher
    auto t0 = std::chrono::high_resolution_clock::now();
    auto rc = pcre2_match(re,             /* the compiled pattern */
                          subject,        /* the subject string */
                          subject_length, /* the length of the subject */
                          0,              /* start at offset 0 in the subject */
                          0,              /* default options */
                          match_data,     /* block for storing the result */
                          NULL);          /* use default match context */
    auto t1 = std::chrono::high_resolution_clock::now();

    // you probably need this stuff
    ovector = pcre2_get_ovector_pointer(match_data);
    // printf("Match succeeded at offset %d\n", (int)ovector[0]);

    // for (int i = 0; i < rc; i++) {
    //     PCRE2_SPTR substring_start = subject + ovector[2 * i];
    //     PCRE2_SIZE substring_length = ovector[2 * i + 1] - ovector[2 * i];
    //     printf("%2d: %.*s\n", i, (int)substring_length,
    //            (char *)substring_start);
    // }

    pcre2_match_data_free(match_data); /* Release the memory that was used */
    pcre2_code_free(re);

    return t1 - t0;
}

std::optional<std::chrono::duration<double, std::nano>>
my_bench(std::string const &in_pattern, std::string const &in_subject) {
    Parser p{in_pattern};
    if (!p.parse()) {
        return std::nullopt;
    }
    Matcher m = p.get_compiled_matcher();
    auto t0 = std::chrono::high_resolution_clock::now();
    m.greedy_match_view(in_subject);
    auto t1 = std::chrono::high_resolution_clock::now();

    return t1 - t0;
}

int main(int argc, char **argv) {

    if (argc < 3) {
        std::cout << "provide both a pattern and a subject, and optionally a "
                     "number of iters"
                  << std::endl;
    }
    std::string pattern{argv[1]};
    std::string subject{argv[2]};
    uint32_t total_iter = 0;
    if (argc >= 4) {
        char *ignore;
        total_iter = strtoumax(argv[3], &ignore, 10);
    }

    total_iter = (total_iter == 0) ? 1000'000 : total_iter;

    double pcre_sum = 0;
    double our_sum = 0;
    for (size_t iter = 0; iter < total_iter; ++iter) {
        auto duration_in_nanos = pcre_bench(pattern, subject);
        if (!duration_in_nanos) {
            std::cerr << "pcre failed to compile pattern " << pattern
                      << std::endl;
            return 0;
        } else {
            pcre_sum += duration_in_nanos->count();
        }

        duration_in_nanos = my_bench(pattern, subject);
        if (!duration_in_nanos) {
            std::cerr << "we failed to compile pattern " << pattern
                      << std::endl;
            return 0;
        } else {
            our_sum += duration_in_nanos->count();
        }
        std::cout << iter << "        \r";
    }
    std::cout << std::endl;

    std::cout << "pattern: " << pattern << std::endl;
    std::cout << "subject: " << subject << std::endl;
    std::cout << "total_iter: " << total_iter << std::endl;
    std::cout << "pcre total time: " << pcre_sum << std::endl;
    std::cout << "our total time: " << our_sum << std::endl;
}
