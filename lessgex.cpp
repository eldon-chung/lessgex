#include "lessgex.h"

#include <iostream>

void parse_test_empty1() {
    bool expected = true;
    Parser p("");

    bool result = p.parse();
    std::cout << __func__ << ": "
              << ((result == expected) ? "passed" : "failed") << std::endl;
}

void parse_test_empty2() {
    bool expected = true;
    Parser p("()");

    bool result = p.parse();
    std::cout << __func__ << ": "
              << ((result == expected) ? "passed" : "failed") << std::endl;
}

void parse_test_empty3() {
    bool expected = true;
    Parser p("()()");

    bool result = p.parse();
    std::cout << __func__ << ": "
              << ((result == expected) ? "passed" : "failed") << std::endl;
}

void parse_test_empty4() {
    bool expected = false;
    Parser p("(");

    bool result = p.parse();
    std::cout << __func__ << ": "
              << ((result == expected) ? "passed" : "failed") << std::endl;
}

void parse_test_empty5() {
    bool expected = false;
    Parser p(")");

    bool result = p.parse();
    std::cout << __func__ << ": "
              << ((result == expected) ? "passed" : "failed") << std::endl;
}

void parse_basic_char1() {
    bool expected = true;

    Parser p("a");

    bool result = p.parse();
    std::cout << __func__ << ": "
              << ((result == expected) ? "passed" : "failed") << std::endl;
}

void parse_basic_char2() {
    bool expected = true;

    Parser p("z");

    bool result = p.parse();
    std::cout << __func__ << ": "
              << ((result == expected) ? "passed" : "failed") << std::endl;
}

void parse_basic_char3() {
    bool expected = true;

    Parser p(" ");

    bool result = p.parse();
    std::cout << __func__ << ": "
              << ((result == expected) ? "passed" : "failed") << std::endl;
}

void parse_basic_char4() {
    bool expected = true;

    Parser p("\t");

    bool result = p.parse();
    std::cout << __func__ << ": "
              << ((result == expected) ? "passed" : "failed") << std::endl;
}

void parse_basic_char5() {
    bool expected = true;

    Parser p("0");

    bool result = p.parse();
    std::cout << __func__ << ": "
              << ((result == expected) ? "passed" : "failed") << std::endl;
}

void parse_basic_char6() {
    bool expected = true;

    Parser p("9");

    bool result = p.parse();
    std::cout << __func__ << ": "
              << ((result == expected) ? "passed" : "failed") << std::endl;
}

void parse_basic_char7() {
    bool expected = true;

    // Parser p("ab10()$\t ");
    Parser p("ab1");

    bool result = p.parse();
    std::cout << __func__ << ": "
              << ((result == expected) ? "passed" : "failed") << std::endl;
}

void parse_basic_char8() {
    bool expected = true;

    // Parser p("ab10()$\t ");
    Parser p("(a)b((c)d)e");

    bool result = p.parse();
    std::cout << __func__ << ": "
              << ((result == expected) ? "passed" : "failed") << std::endl;
}

void parse_basic_set1() {
    bool expected = true;

    // Parser p("ab10()$\t ");
    Parser p(".");

    bool result = p.parse();
    std::cout << __func__ << ": "
              << ((result == expected) ? "passed" : "failed") << std::endl;
}

void parse_basic_set2() {
    bool expected = true;

    // Parser p("ab10()$\t ");
    Parser p("..");

    bool result = p.parse();
    std::cout << __func__ << ": "
              << ((result == expected) ? "passed" : "failed") << std::endl;
}

void parse_basic_set3() {
    bool expected = true;

    // Parser p("ab10()$\t ");
    Parser p("\\s");

    bool result = p.parse();
    std::cout << __func__ << ": "
              << ((result == expected) ? "passed" : "failed") << std::endl;
}

void parse_basic_op1() {
    bool expected = true;

    // Parser p("ab10()$\t ");
    Parser p("a|b");

    bool result = p.parse();
    std::cout << __func__ << ": "
              << ((result == expected) ? "passed" : "failed") << std::endl;
}

void parse_basic_op2() {
    bool expected = true;

    // Parser p("ab10()$\t ");
    Parser p("a*");

    bool result = p.parse();
    std::cout << __func__ << ": "
              << ((result == expected) ? "passed" : "failed") << std::endl;
}

void parse_basic_op3() {
    bool expected = true;

    // Parser p("ab10()$\t ");
    Parser p("ab*");

    bool result = p.parse();
    std::cout << __func__ << ": "
              << ((result == expected) ? "passed" : "failed") << std::endl;
}

void parse_basic_op4() {
    bool expected = true;

    // Parser p("ab10()$\t ");
    Parser p("(a)+");

    bool result = p.parse();
    std::cout << __func__ << ": "
              << ((result == expected) ? "passed" : "failed") << std::endl;
}

void parse_basic_op5() {
    bool expected = true;

    // Parser p("ab10()$\t ");
    Parser p("ab+");

    bool result = p.parse();
    std::cout << __func__ << ": "
              << ((result == expected) ? "passed" : "failed") << std::endl;
}

void parse_basic_op6() {
    bool expected = true;

    // Parser p("ab10()$\t ");
    Parser p("ab?");

    bool result = p.parse();
    std::cout << __func__ << ": "
              << ((result == expected) ? "passed" : "failed") << std::endl;
}

void parse_basic_op7() {
    bool expected = true;

    // Parser p("ab10()$\t ");
    Parser p("(ab)?");

    bool result = p.parse();
    std::cout << __func__ << ": "
              << ((result == expected) ? "passed" : "failed") << std::endl;
}

void parse_basic_op8() {
    bool expected = false;

    // Parser p("ab10()$\t ");
    Parser p("(ab)??");

    bool result = p.parse();
    std::cout << __func__ << ": "
              << ((result == expected) ? "passed" : "failed") << std::endl;
}

void parse_basic_op9() {
    bool expected = false;

    // Parser p("ab10()$\t ");
    Parser p("a?+b?");

    bool result = p.parse();
    std::cout << __func__ << ": "
              << ((result == expected) ? "passed" : "failed") << std::endl;
}

void parse_basic_op10() {
    bool expected = false;

    // Parser p("ab10()$\t ");
    Parser p("*");

    bool result = p.parse();
    std::cout << __func__ << ": "
              << ((result == expected) ? "passed" : "failed") << std::endl;
}

void parse_set_combi1() {
    bool expected = false;

    // Parser p("ab10()$\t ");
    Parser p("[]");

    bool result = p.parse();
    std::cout << __func__ << ": "
              << ((result == expected) ? "passed" : "failed") << std::endl;
}

void parse_set_combi2() {
    bool expected = false;

    // Parser p("ab10()$\t ");
    Parser p("[^]");

    bool result = p.parse();
    std::cout << __func__ << ": "
              << ((result == expected) ? "passed" : "failed") << std::endl;
}

void parse_set_combi3() {
    bool expected = true;

    // Parser p("ab10()$\t ");
    Parser p("[a]");

    bool result = p.parse();
    std::cout << __func__ << ": "
              << ((result == expected) ? "passed" : "failed") << std::endl;
}

void parse_set_combi4() {
    bool expected = true;

    // Parser p("ab10()$\t ");
    Parser p("[a-]");

    bool result = p.parse();
    std::cout << __func__ << ": "
              << ((result == expected) ? "passed" : "failed") << std::endl;
}

void parse_set_combi5() {
    bool expected = true;

    // Parser p("ab10()$\t ");
    Parser p("[a-b]");

    bool result = p.parse();
    std::cout << __func__ << ": "
              << ((result == expected) ? "passed" : "failed") << std::endl;
}

void parse_set_combi6() {
    bool expected = true;

    // Parser p("ab10()$\t ");
    Parser p("[a-zf-h]");

    bool result = p.parse();
    std::cout << __func__ << ": "
              << ((result == expected) ? "passed" : "failed") << std::endl;
}

void parse_set_combi7() {
    bool expected = true;

    Parser p("[a-z9-5]");

    bool result = p.parse();
    std::cout << __func__ << ": "
              << ((result == expected) ? "passed" : "failed") << std::endl;
}

void parse_word_boundary1() {
    bool expected = true;
    Parser p("\\ba");

    bool result = p.parse();
    std::cout << __func__ << ": "
              << ((result == expected) ? "passed" : "failed") << std::endl;
}

void parse_word_boundary2() {
    bool expected = true;
    Parser p("\\Ba");

    bool result = p.parse();
    std::cout << __func__ << ": "
              << ((result == expected) ? "passed" : "failed") << std::endl;
}

void parse_word_boundary3() {
    bool expected = true;
    Parser p("\\B\\ba");

    bool result = p.parse();
    std::cout << __func__ << ": "
              << ((result == expected) ? "passed" : "failed") << std::endl;
}

void parse_word_boundary4() {
    bool expected = true;
    Parser p(".+\btg");

    bool result = p.parse();
    std::cout << __func__ << ": "
              << ((result == expected) ? "passed" : "failed") << std::endl;
}

void parse_arb_1() {
    bool expected = true;
    Parser p("^.+\\btg&[a-]$");
    // Parser p("^.+\\btg&[a-]$");

    bool result = p.parse();
    std::cout << __func__ << ": "
              << ((result == expected) ? "passed" : "failed") << std::endl;
}

void parse_arb_2() {
    bool expected = true;
    Parser p("\\b^");
    // Parser p("^.+\\btg&[a-]$");

    bool result = p.parse();
    std::cout << __func__ << ": "
              << ((result == expected) ? "passed" : "failed") << std::endl;
}

void parse_arb_3() {
    bool expected = true;
    Parser p("^\\b");
    // Parser p("^.+\\btg&[a-]$");

    bool result = p.parse();
    std::cout << __func__ << ": "
              << ((result == expected) ? "passed" : "failed") << std::endl;
}

void test_table(std::string_view regex_pattern) {
    std::cout << "compiling:[" << regex_pattern << "]" << std::endl;
    Parser p{regex_pattern};
    p.parse();
    std::cout << p.get_compiled_table() << std::endl;
}

void test_match(std::string_view regex_pattern, std::string_view subject,
                std::optional<Matcher::Result> expected_result) {

    // std::cout << "compiling:[" << regex_pattern << "]" << std::endl;
    Parser p{regex_pattern};
    if (!p.parse()) {
        std::cerr << "failed to parse this" << std::endl;
        return;
    }
    // std::cout << "getting matcher" << std::endl;
    Matcher m = p.get_compiled_matcher();

    // std::cout << "calling matcher" << std::endl;
    auto res = m.greedy_match_view(subject);

    // std::cout << "pattern: [" << regex_pattern << "] subject: [" << subject
    //           << "]" << std::endl;

    if (expected_result != res) {
        std::cout << "expected: [" << expected_result << "] result: [" << res
                  << "]" << std::endl;
    } else {
        std::cout << "passed" << std::endl;
    }
}

int main(void) {
    std::cout << std::boolalpha;

    Parser p(".*.*");
    bool parse_result = p.parse();
    std::cout << "parsed" << std::endl;
    if (!parse_result) {
        std::cout << "failed to parse for some reason" << std::endl;
        return 0;
    } else {
        std::cout << "parsed, continuing" << std::endl;
    }

    // Matcher matcher = Matcher(std::move(tb));
    Matcher matcher = p.get_compiled_matcher();
    std::cerr << "compiled" << std::endl;
    // std::cerr << matcher << std::endl;
    // for (size_t idx = 0; idx < 1000'000'0; ++idx) {
    auto maybe_res =
        matcher.greedy_match_view("aababaaaaabababaaaababbbbabbbbabbbababab"
                                  "bbababbaabbbabbaabbabbababab");
    // }
    std::cout << maybe_res << std::endl;

    return 0;
}