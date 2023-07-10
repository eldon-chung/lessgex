#pragma once

#include <assert.h>

#include <stdint.h>
#include <stdlib.h>

#include <algorithm>
#include <array>
#include <iostream>
#include <iterator>
#include <limits>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <ranges>

// Stolen from boost
template <typename... Ts> static size_t combine(size_t seed, Ts &&...vs) {
    auto xorshift = [](size_t x, int i) -> size_t { return x ^ (x >> i); };

    auto distribute = [&xorshift](uint64_t n) -> uint64_t {
        uint64_t p = 0x5555555555555555ull;   // pattern of alternating 0 and 1
        uint64_t c = 17316035218449499591ull; // random uneven integer constant;
        return c * xorshift(p * xorshift(n, 32), 32);
    };

    return (std::rotl(seed, std::numeric_limits<size_t>::digits / 3) ^ ... ^
            distribute(std::hash<Ts>()(vs)));
}

struct TransitionState {
    static const uint8_t NONE = 0x00;
    static const uint8_t WB = 0x01;
    static const uint8_t NWB = 0x02;
    static inline size_t global_id = 0;

    size_t id;
    uint8_t front_entered_modifier;

    TransitionState() : id(global_id++), front_entered_modifier(NONE) {}

    TransitionState(size_t _id, uint8_t _front_entered_modifier)
        : id(_id), front_entered_modifier(_front_entered_modifier) {}
    TransitionState(uint8_t _front_entered_modifier)
        : id(global_id++), front_entered_modifier(_front_entered_modifier) {}

    bool operator==(TransitionState const &other) const {
        return id == other.id;
    };
    void set_word_boundary() { front_entered_modifier |= WB; }
    void set_non_word_boundary() { front_entered_modifier |= NWB; }

    uint8_t get_front_modifier() const { return front_entered_modifier; }
    void set_front_modifier(uint8_t val) { front_entered_modifier = val; }
    void add_front_modifier(uint8_t val) { front_entered_modifier |= val; }
};

std::ostream &operator<<(std::ostream &os, TransitionState const &ts) {
    os << "id:" << ts.id
       << " front_entered_modifier: " << (int)ts.front_entered_modifier;
    return os;
}

template <> struct std::hash<TransitionState> {
    std::size_t operator()(const TransitionState &k) const {
        return std::hash<size_t>()(k.id);
    }
};

template <> struct std::hash<std::pair<TransitionState, char>> {
    std::size_t operator()(const std::pair<TransitionState, char> &pair) const {
        return combine(69420, std::hash<TransitionState>()(pair.first),
                       std::hash<char>()(pair.second));
    }
};

using StateChar = std::pair<TransitionState, char>;

struct TransitionTable {
    std::unordered_map<StateChar, std::vector<TransitionState>> table;

    std::unordered_set<TransitionState> starting_states;
    std::unordered_set<TransitionState> accepting_states;
};

std::ostream &operator<<(std::ostream &os, StateChar const &sc) {
    os << "state id: " << sc.first << " char: " << sc.second;
    return os;
}

std::ostream &operator<<(std::ostream &os, TransitionTable const &tb) {
    std::cout << "starting_states ====================" << std::endl;
    for (auto const &s : tb.starting_states) {
        std::cout << s << std::endl;
    }
    std::cout << "starting_states end ====================" << std::endl;

    std::cout << "accepting_states ====================" << std::endl;
    for (auto const &s : tb.accepting_states) {
        std::cout << s << std::endl;
    }
    std::cout << "accepting_states end ====================" << std::endl;

    os << "transition table ====================" << std::endl;
    for (auto const &transition : tb.table) {
        os << "state char pair:{" << transition.first << "}" << std::endl;
        ;
        os << "next states: " << std::endl;
        for (auto const &next_states : transition.second) {
            os << next_states << std::endl;
        }
        os << "next states end" << std::endl;
    }
    os << "transition table end ====================" << std::endl;

    return os;
}

/*
    right now i'd rather take in a string view so we dont make needless copies
    but this means the lifetime of our matcher is tied to the string view
*/
class Matcher {

  public:
    struct Result {
        size_t start;
        size_t length;

        bool operator==(Result const &other) const = default;
    };

    struct ResultHash {
        size_t operator()(Result const &r) const {
            return combine(69420, std::hash<size_t>()(r.start),
                           std::hash<size_t>()(r.length));
        }
    };

  private:
    TransitionTable tb;

    // intermediate state?
    struct State {
        size_t offset;
        TransitionState nfa_state;
        bool operator==(State const &other) const = default;
    };

    struct StateHash {
        size_t operator()(State const &r) const {
            return combine(69420, std::hash<size_t>()(r.offset),
                           std::hash<TransitionState>()(r.nfa_state));
        }
    };

  public:
    Matcher(TransitionTable _tb) : tb(std::move(_tb)) {}

    // returns an offset into the subject of the earliest matching
    // result
    // note: perhaps we'll eventually want something more stateful
    std::optional<Result> greedy_match_view(std::string_view str) {
        std::unordered_set<Result, ResultHash>
            result_list; // stores active results
        std::unordered_set<State, StateHash> active_states;
        std::unordered_set<State, StateHash> next_states;

        auto process_accepting_states = [&](size_t s_idx) {
            for (State state : active_states) {
                if (!tb.accepting_states.contains(state.nfa_state)) {
                    continue;
                }
                // annoying but do the lookahead test here i guess
                if (state.nfa_state.get_front_modifier() ==
                    TransitionState::NONE) {
                    result_list.insert({state.offset, s_idx - state.offset});
                } else if (state.nfa_state.get_front_modifier() ==
                           TransitionState::WB) {
                    char lookahead = (s_idx < str.size()) ? str[s_idx] : 3;
                    if (is_word_char(str[s_idx]) != is_word_char(lookahead)) {
                        result_list.insert(
                            {state.offset, s_idx - state.offset});
                    }
                } else if (state.nfa_state.get_front_modifier() ==
                           TransitionState::NWB) {
                    char lookahead = (s_idx < str.size()) ? str[s_idx] : 3;
                    if (is_word_char(str[s_idx]) == is_word_char(lookahead)) {
                        result_list.insert(
                            {state.offset, s_idx - state.offset});
                    }
                } else {
                    continue;
                }
            }
        };

        // DEBUG:
        // {
        //     std::cerr << "about to start matching, accepting states:"
        //               << std::endl;
        //     std::for_each(tb.accepting_states.begin(),
        //                   tb.accepting_states.end(),
        //                   [&](TransitionState const &s) { std::cout << s; });
        //     std::cerr << "accepting states end" << std::endl;
        // }

        // run one round with the start of text char
        // add new starting states from this offset;
        std::transform(tb.starting_states.begin(), tb.starting_states.end(),
                       std::inserter(active_states, active_states.end()),
                       [](TransitionState const &nfa_state) -> State {
                           return {0, nfa_state};
                       });

        // get one transition for each of the states
        // add new starting states from this offset;
        transition_char(active_states, next_states, 2, 2);
        char prev_char = 2;
        std::swap(active_states, next_states);
        next_states.clear();

        for (size_t s_idx = 0; s_idx < str.size(); ++s_idx) {

            // add new starting states from this offset;
            std::transform(tb.starting_states.begin(), tb.starting_states.end(),
                           std::inserter(active_states, active_states.end()),
                           [s_idx](TransitionState const &nfa_state) -> State {
                               return {s_idx, nfa_state};
                           });

            process_accepting_states(s_idx);
            // DEBUG:
            // {
            //     std::cerr << "======================================="
            //               << std::endl;
            //     std::cerr << "at idx:" << s_idx << std::endl;
            //     std::cerr << "active states: " << std::endl;
            //     std::for_each(active_states.begin(), active_states.end(),
            //                   [&](State const &s) {
            //                       std::cout << "[[[[nfa_state: " <<
            //                       s.nfa_state
            //                                 << " offset: " << s.offset
            //                                 << "]]]]";
            //                   });
            //     std::cout << std::endl;
            //     std::cerr << "end active states" << std::endl;
            // }

            // get one transition for each of the states
            // add new starting states from this offset;
            transition_char(active_states, next_states, str[s_idx], prev_char);

            // DEBUG
            // {
            //     std::cerr << "next states: " << std::endl;
            //     std::for_each(next_states.begin(), next_states.end(),
            //                   [&](State const &s) {
            //                       std::cout << "[[[[nfa_state: " <<
            //                       s.nfa_state
            //                                 << " offset: " << s.offset
            //                                 << "]]]]";
            //                   });
            //     std::cout << std::endl;
            //     std::cerr << "end next states" << std::endl;
            // }

            // update prev char
            prev_char = str[s_idx];

            // swap the two states
            std::swap(active_states, next_states);
            next_states.clear();
        }

        // last round here
        process_accepting_states(str.size());

        next_states.clear();
        transition_char(active_states, next_states, 3, prev_char);
        std::swap(active_states, next_states);

        // last round here with EOF
        process_accepting_states(str.size());

        // take the earliest result and return it
        if (result_list.empty()) {
            return {};
        }
        Result res = *result_list.begin();

        for (Result r : result_list) {
            if (r.start < res.start) {
                res = r;
            }

            if (r.start == res.start && r.length > res.length) {
                res = r;
            }
        }
        return res;
    }

    void
    progress_states_into(State state, char ch,
                         std::unordered_set<State, StateHash> &next_states) {

        auto it = tb.table.find({state.nfa_state, ch});
        if (it == tb.table.end()) {
            return;
        }
        auto const &new_states = it->second;

        std::transform(new_states.begin(), new_states.end(),
                       std::inserter(next_states, next_states.end()),
                       [state](TransitionState const &ts) {
                           return State{state.offset, ts};
                       });
    }

    bool is_word_char(char ch) {
        return ch != ' ' && ch != '\t' && ch != 3 && ch != 2;
    }

    void
    transition_char(std::unordered_set<State, StateHash> const &active_states,
                    std::unordered_set<State, StateHash> &next_states, char ch,
                    char prev_char) {
        // get one transition for each of the states
        for (State s : active_states) {
            auto modifier = s.nfa_state.get_front_modifier();
            switch (modifier) {
            case TransitionState::NONE:
                progress_states_into(s, ch, next_states);
                break;
            case TransitionState::WB:
                if (is_word_char(prev_char) != is_word_char(ch)) {
                    progress_states_into(s, ch, next_states);
                }
                break;

            case TransitionState::NWB:
                if (is_word_char(prev_char) == is_word_char(ch)) {
                    progress_states_into(s, ch, next_states);
                }
                break;

            default:
                // this means both were switched on so we do nothing;
                break;
            }
        }
    }
};

struct MatcherBuilder;
std::ostream &operator<<(std::ostream &os, MatcherBuilder const &tb);

struct MatcherBuilder {

    std::vector<TransitionState> starting_states;
    std::vector<TransitionState> accepting_states;
    std::unordered_map<StateChar, std::vector<TransitionState>> table;
    std::unordered_map<TransitionState, std::vector<StateChar>> reverse_table;

    MatcherBuilder()
        : starting_states(1, TransitionState()), accepting_states{
                                                     starting_states.front()} {}

    MatcherBuilder(MatcherBuilder const &to_copy) {
        std::unordered_map<TransitionState, TransitionState> old_to_new;

        // will bump the numbers by a lot but who cares
        for (const auto &p : to_copy.table) {
            old_to_new[p.first.first] = TransitionState();
            for (const auto &t : p.second) {
                old_to_new[t] = TransitionState(t.get_front_modifier());
            }
        }
        for (const auto &p : to_copy.starting_states) {
            old_to_new[p] = TransitionState(p.get_front_modifier());
        }
        for (const auto &p : to_copy.accepting_states) {
            old_to_new[p] = TransitionState(p.get_front_modifier());
        }

        starting_states.reserve(to_copy.starting_states.size());
        for (auto const &tc_s : to_copy.starting_states) {
            starting_states.push_back(old_to_new.at(tc_s));
        }
        accepting_states.reserve(to_copy.accepting_states.size());
        for (auto const &tc_s : to_copy.accepting_states) {
            accepting_states.push_back(old_to_new.at(tc_s));
        }

        std::vector<TransitionState> new_next;
        for (auto const &t_p : to_copy.table) {
            // prep the next vec first
            new_next.clear();
            new_next.reserve(t_p.second.size());
            for (auto const &next : t_p.second) {
                new_next.push_back(old_to_new.at(next));
            }
            table.insert({{old_to_new.at(t_p.first.first), t_p.first.second},
                          std::move(new_next)});
        }

        std::vector<StateChar> new_next_reverse;
        for (auto const &r_p : to_copy.reverse_table) {
            new_next_reverse.clear();
            new_next_reverse.reserve(r_p.second.size());
            for (auto const &sc : r_p.second) {
                new_next_reverse.push_back(
                    StateChar{old_to_new.at(sc.first), sc.second});
            }

            reverse_table[old_to_new.at(r_p.first)].insert(
                reverse_table[old_to_new.at(r_p.first)].end(),
                new_next_reverse.begin(), new_next_reverse.end());
        }
    }

    friend void swap(MatcherBuilder &mb1, MatcherBuilder &mb2) {
        using std::swap;
        swap(mb1.accepting_states, mb2.accepting_states);
        swap(mb1.starting_states, mb2.starting_states);
        swap(mb1.table, mb2.table);
        swap(mb1.reverse_table, mb2.reverse_table);
    }

    MatcherBuilder &operator=(MatcherBuilder const &to_copy) {
        MatcherBuilder temp{to_copy};
        swap(*this, temp);
        return *this;
    }

    MatcherBuilder(char c)
        : starting_states(1, TransitionState()),
          accepting_states(1, TransitionState()),
          table{{{starting_states.front(), c}, {accepting_states.front()}}},
          reverse_table{
              {accepting_states.front(), {{starting_states.front(), c}}}} {}

    MatcherBuilder(char lb, char rb)
        : starting_states(1, TransitionState()),
          accepting_states(1, TransitionState()) {
        assert(lb < rb);
        TransitionState const &s_state = starting_states.front();
        TransitionState const &e_state = accepting_states.front();
        for (char c = lb; c < rb; ++c) {
            table[{s_state, c}].push_back(e_state);
            reverse_table[e_state].push_back({s_state, c});
        }
    }

    // returns true if there is only a starting state and
    // we have no transitions
    bool fresh_state() const {
        return starting_states.size() == 1 && accepting_states.size() == 1 &&
               table.empty();
    }

    MatcherBuilder operator+(MatcherBuilder const &other) const {
        MatcherBuilder to_ret = *this;
        to_ret += other;
        return to_ret;
    }

    MatcherBuilder operator|(MatcherBuilder const &other) {
        MatcherBuilder to_ret = *this;
        to_ret |= other;
        return to_ret;
    }

    MatcherBuilder &append_char_range(char begin, char end) {
        assert(begin < end);
        MatcherBuilder mb{begin};
        for (char c = begin + 1; c < end; ++c) {
            mb.table.at({mb.starting_states.front(), c})
                .push_back(mb.accepting_states.front());

            mb.reverse_table.at(mb.accepting_states.front())
                .push_back({mb.accepting_states.front(), c});
        }

        return *this;
    }

    MatcherBuilder &operator+=(char c) {
        if (fresh_state()) {
            assert(starting_states.size() == 1);
            TransitionState new_state = TransitionState();
            table[{starting_states.front(), c}].push_back(new_state);
            reverse_table[new_state].push_back({starting_states.front(), c});
            // remove the old accepting state
            accepting_states.clear();
            accepting_states.push_back(new_state);
        } else {
            // this is so inefficient. fix this later
            *this += MatcherBuilder(c);
        }

        return *this;
    }

    MatcherBuilder &operator+=(MatcherBuilder const &mb) {
        // might need to deal with initially empty states that
        // have 0 transitions;
        if (fresh_state()) {
            // front modify next dfa
            uint8_t old_boundary = starting_states.front().get_front_modifier();
            *this = mb;

            std::for_each(starting_states.begin(), starting_states.end(),
                          [old_boundary](TransitionState &s) {
                              s.add_front_modifier(old_boundary);
                          });
        } else {
            add_parallel_transitions(accepting_states, mb.starting_states);
            // append the tables too you dolt
            table.insert(mb.table.begin(), mb.table.end());
            reverse_table.insert(mb.reverse_table.begin(),
                                 mb.reverse_table.end());
            accepting_states = mb.accepting_states;
        }
        return *this;
    }

    MatcherBuilder &operator|=(MatcherBuilder const &mb) {
        if (fresh_state()) {
            *this = mb;
        } else {
            accepting_states.insert(accepting_states.end(),
                                    mb.accepting_states.begin(),
                                    mb.accepting_states.end());
            starting_states.insert(starting_states.end(),
                                   mb.starting_states.begin(),
                                   mb.starting_states.end());
            table.insert(mb.table.begin(), mb.table.end());
            reverse_table.insert(mb.reverse_table.begin(),
                                 mb.reverse_table.end());
        }
        return *this;
    }

    MatcherBuilder &question_modify() {
        accepting_states.insert(accepting_states.end(), starting_states.begin(),
                                starting_states.end());
        return *this;
    }

    MatcherBuilder &star_modify() {
        add_parallel_transitions(accepting_states, starting_states);
        accepting_states.insert(accepting_states.end(), starting_states.begin(),
                                starting_states.end());
        return *this;
    }

    MatcherBuilder &plus_modify() {
        MatcherBuilder mb = *this;
        mb.star_modify();
        *this += mb;

        return *this;
    }

    MatcherBuilder &boundary_modify(uint8_t modifier) {
        // this feels like a hotfix but this
        // seems like the only place that it'll happen so
        std::for_each(
            starting_states.begin(), starting_states.end(),
            [=](TransitionState &s) { s.add_front_modifier(modifier); });

        std::unordered_set<TransitionState> starting_set{
            starting_states.begin(), starting_states.end()};
        for (auto &s : accepting_states) {
            if (starting_set.contains(s)) {
                s.add_front_modifier(modifier);
            }
        }

        return *this;
    }

  private:
    //   updates the table
    void
    add_parallel_transitions(std::vector<TransitionState> const &old_target,
                             std::vector<TransitionState> const &new_target) {
        for (TransitionState s : old_target) {
            // assert(reverse_table.contains(s));
            for (auto [prev_s, ch] : reverse_table[s]) {
                table[{prev_s, ch}].insert(table[{prev_s, ch}].end(),
                                           new_target.begin(),
                                           new_target.end());

                for (TransitionState new_s : new_target) {
                    reverse_table[new_s].push_back({prev_s, ch});
                }
            }
        }
    }
};

struct Token {
    enum class Type {
        CHAR,
        PLUS,
        STAR,
        QUESTION,
        ESCAPED,
        PAREN_OPEN,
        PAREN_CLOSE,
        SET_OPEN,
        SET_CLOSE,
        DASH,
        DOLLAR,
        CARET,
        END,
    };

    Type type;
    char c;
};

struct Parser {
    // handles escaped chars
    struct EChar {
        char ch;
        bool escaped = false;

        EChar(char _ch) : ch(_ch) {}
        EChar(char _ch, bool _escaped) : ch(_ch), escaped(_escaped) {}

        // should never be true to a char if escpaed
        bool operator==(char c) const { return !escaped && ch == c; }
        bool is_escaped() const { return escaped; }
        char get_ch() const { return ch; }
    };

    std::string_view sv;
    std::optional<MatcherBuilder> result;

    Parser(std::string_view _sv) : sv(_sv) {}

    // needs to be called after a successful parse
    TransitionTable get_compiled_table() {
        if (!result) {
            throw std::runtime_error("there isnt anything parsed");
        }
        TransitionTable to_ret = TransitionTable{
            result->table,
            std::unordered_set<TransitionState>{result->starting_states.begin(),
                                                result->starting_states.end()},
            std::unordered_set<TransitionState>{
                result->accepting_states.begin(),
                result->accepting_states.end()}};

        return to_ret;
    }

    Matcher get_compiled_matcher() { return Matcher(get_compiled_table()); }

    // skips over escaped chars too
    void sv_advance() {
        if (sv.size() >= 2 && sv[0] == '\\') {
            sv = sv.substr(2, std::string::npos);
        } else if (!sv.empty()) {
            sv = sv.substr(1, std::string::npos);
        }
    }

    // does not advance sv
    EChar sv_front() const {
        if (sv.empty()) {
            return {0};
        }

        if (sv.size() >= 2 && sv[0] == '\\') {
            // try getting another one
            return {sv[1], true};
        } else {
            return {sv.front()};
        }
    }

    bool sv_empty() const { return sv.empty(); }

    std::vector<std::pair<char, char>>
    negate_set_ranges(std::vector<std::pair<char, char>> &s_ranges) {
        // we're going to assume input is nice-ish and maybe catch issues if
        // need be
        std::sort(s_ranges.begin(), s_ranges.end(),
                  [](std::pair<char, char> pair1, std::pair<char, char> pair2) {
                      return pair1.first < pair2.first;
                  });

        std::vector<std::pair<char, char>> to_ret;

        size_t idx = 0;
        if (auto [s, e] = s_ranges.front(); s != '\t') {
            assert(s_ranges.front().first > '\t');
            to_ret.push_back({'\t', '\t' + 1});
        } else {
            ++idx;
        }

        char left_end = ' ';
        while (idx < s_ranges.size()) {
            char right_end = s_ranges[idx].first;
            assert(left_end <= right_end);
            if (left_end < right_end) {
                to_ret.push_back({left_end, right_end});
            }
            left_end = s_ranges[idx].second;
            ++idx;
        }

        if (left_end < 127) {
            to_ret.push_back({left_end, 127});
        }

        return to_ret;
    }

    std::vector<std::pair<char, char>> parse_predefined_set_into_ranges() {
        EChar curr_char = sv_front();
        if (!curr_char.is_escaped()) {
            return {};
        }

        switch (curr_char.get_ch()) {
        case 's':
            sv_advance();
            return {{'\t', '\t' + 1}, {' ', ' ' + 1}};

        case 'S':
            sv_advance();
            return {{33, 127}};

        case 'd':
            sv_advance();
            return {{'0', '9' + 1}};
        case 'D':
            sv_advance();
            return {{'\t', '\t' + 1}, {' ', '0'}, {'9' + 1, 127}};
        case 'w':
            sv_advance();
            return {{'A', 'Z' + 1}, {'a', 'z' + 1}};
        case 'W':
            sv_advance();
            return {
                {'\t', '\t' + 1}, {' ', 'A'}, {'Z' + 1, 'a'}, {'z' + 1, 127}};

        default:
            return {};
        }
    }

    // modifies the given mb with post-modifiers, if any
    void parse_post_modifiers(MatcherBuilder &mb) {
        EChar curr_char = sv_front();
        if (curr_char.is_escaped()) {
            return;
        }

        switch (curr_char.get_ch()) {
        case '?':
            mb.question_modify();
            sv_advance();
            break;
        case '+':
            mb.plus_modify();
            sv_advance();
            break;
        case '*':
            mb.star_modify();
            sv_advance();
            break;
        default:
            break;
        }
    }

    // because -sighs-
    void tidy_ranges(std::vector<std::pair<char, char>> &ranges) {
        std::vector<std::pair<char, char>> to_ret;

        std::for_each(ranges.begin(), ranges.end(), [](auto &p1) {
            if (p1.first > p1.second) {
                auto temp = p1.first;
                p1.first = p1.second;
                p1.second = temp;
            }
        });

        std::sort(ranges.begin(), ranges.end(),
                  [](auto p1, auto p2) { return p1.first < p2.first; });

        auto [left_end, right_end] = ranges.front();
        for (auto ran : ranges) {
            if (right_end >= ran.first) {
                right_end = (right_end > ran.second) ? right_end : ran.second;
                continue;
            }

            to_ret.push_back({left_end, right_end});
            left_end = ran.first;
            right_end = ran.second;
        }

        to_ret.push_back({left_end, right_end});

        ranges = to_ret;
    }

    // enters parsing mode for sets [ ]
    std::optional<MatcherBuilder> parse_set() {
        // quit the moment we see a ']'
        bool neg = (sv_front() == '^');
        if (neg) {
            sv_advance();
        }

        // don't accept empty sets like this
        if (sv_front() == ']') {
            return {};
        }

        // prepare the char ranges needed
        std::vector<std::pair<char, char>> ranges_to_create;
        while (!sv_empty() && sv_front() != ']') {
            EChar current_char = sv_front();
            auto predef_range = parse_predefined_set_into_ranges();
            // case 1: it's a predefined set
            if (!predef_range.empty()) {
                ranges_to_create.insert(ranges_to_create.end(),
                                        predef_range.begin(),
                                        predef_range.end());
                continue;
            }

            // case 2: it MIGHT be a range or just a bunch of chars
            sv_advance();
            EChar lookahead = sv_front();
            if (lookahead != '-') {
                // save lookahead for the next iter
                ranges_to_create.push_back(
                    {current_char.get_ch(), current_char.get_ch() + 1});
                continue;
            }

            // else we saw a '-'; then:
            sv_advance();
            EChar next = sv_front();
            // case 3.1: we didnt see a ']'
            if (next != ']') {
                sv_advance();
                ranges_to_create.push_back(
                    {current_char.get_ch(), next.get_ch() + 1});
                continue;
            }

            // case 3.2: we saw a ']', then we just add the current char and -
            // and break
            ranges_to_create.push_back(
                {current_char.get_ch(), current_char.get_ch() + 1});
            ranges_to_create.push_back({'-', '-' + 1});
            break;

            // kind of weird this ends up being the break case
        }
        tidy_ranges(ranges_to_create);
        if (neg) {
            ranges_to_create = negate_set_ranges(ranges_to_create);
        }

        MatcherBuilder mb;
        for (auto r : ranges_to_create) {
            assert(r.first < r.second);
            mb |= MatcherBuilder(r.first, r.second);
        }
        return mb;
    }

    std::optional<MatcherBuilder> parse_char() {
        MatcherBuilder char_seq;
        EChar curr_char = sv_front();
        // case char set
        if (curr_char.is_escaped()) {

            if (curr_char.get_ch() == 's' || curr_char.get_ch() == 'S' ||
                curr_char.get_ch() == 'd' || curr_char.get_ch() == 'D' ||
                curr_char.get_ch() == 'w' || curr_char.get_ch() == 'W') {
                // note this advances the sv
                auto set_ranges = parse_predefined_set_into_ranges();
                for (auto ran : set_ranges) {
                    char_seq |= MatcherBuilder(ran.first, ran.second);
                }
            } else if (curr_char.get_ch() == 'b' || curr_char.get_ch() == 'B') {
                return std::nullopt;
            } else {
                char_seq += curr_char.get_ch();
                sv_advance();
            }

        } else if (curr_char == '.') {
            // this is the match all pattern
            (char_seq += '\t') |= MatcherBuilder(32, 127);
            sv_advance();

        } else if (curr_char == '(' || curr_char == ')' || curr_char == '[' ||
                   curr_char == ']' || curr_char == '^' || curr_char == '$' ||
                   curr_char == '+' || curr_char == '?' || curr_char == '*' ||
                   curr_char == '|') {
            return std::nullopt;
        } else {
            char_seq += curr_char.get_ch();
            sv_advance();
        }

        // now check for post modifier
        parse_post_modifiers(char_seq);
        return char_seq;
    }

    bool parse() {
        MatcherBuilder accum_b;
        bool res = parse_help(accum_b);
        if (res) {
            result = accum_b;
        }
        // success is if we successfully parsed something and the entire string
        // was consumed
        return (bool)res && sv.empty();
    }

    bool parse_help(MatcherBuilder &accum_b) {

        if (sv_empty()) {
            return true;
        }

        if (sv_front() == ')' || sv_front() == ']') {
            return true;
        }

        // handles the current expression
        MatcherBuilder mb;
        // handle bol
        while (sv_front() == '^') {
            // deal with this at some point
            mb += (char)2; // let 2 be the line beginning of text
            sv_advance();
        }

        // handle border modifier
        while (sv_front().is_escaped() &&
               (sv_front().get_ch() == 'b' || sv_front().get_ch() == 'B')) {
            mb.boundary_modify((sv_front().get_ch() == 'b')
                                   ? TransitionState::WB
                                   : TransitionState::NWB);
            sv_advance();
        }

        if (sv_front() == '(') {
            sv_advance(); // move past the '('
            // parse the sub expression
            bool valid = parse_help(mb);
            if (!valid) {
                return false;
            }

            // expect the rparen
            if (sv_front() != ')') {
                return false;
            }
            sv_advance(); // move past the ')'

            // handle post modifiers
            parse_post_modifiers(mb);

            // handle eol
            if (sv_front() == '$') {
                mb += (char)3; // let 3 be the ending of text
                sv_advance();
            }
            accum_b += mb;
            // next parse
            auto next_parse = parse_help(accum_b);
            return next_parse;
        }

        if (sv_front() == '[') {
            sv_advance(); // move past the '['
            auto maybe_set = parse_set();

            if (!maybe_set) {
                return false;
            }

            // expect the rparen
            if (sv_front() != ']') {
                return false;
            }
            // move past the ']'
            sv_advance();
            // handle post modifiers
            parse_post_modifiers(*maybe_set);

            // handle eol
            if (sv_front() == '$') {
                *maybe_set += (char)3; // let 3 be the ending of text
                sv_advance();
            }
            accum_b += (mb += *maybe_set);
            return true;
        }

        // need to handle charsets,
        // normal chars, certain escaped chars
        // dot

        while (!sv_empty()) {
            auto maybe_char = parse_char();
            if (!maybe_char) {
                break;
            }
            mb += *maybe_char;
        }

        if (sv_front() == '|') {
            sv_advance();
            bool res = parse_help(mb);
            if (!res) {
                return false;
            } else {
                accum_b |= mb;
                return true;
            }
        }

        // handle eol
        if (sv_front() == '$') {
            mb += (char)3; // let 3 be the ending of text
            sv_advance();
        }

        accum_b += mb;
        return parse_help(accum_b);
    }
};

std::ostream &operator<<(std::ostream &os, MatcherBuilder const &tb) {
    os << "starting states: " << std::endl;
    for (const auto &s : tb.starting_states) {
        os << s << std::endl;
    }
    os << "========================" << std::endl;
    os << "accepting states: " << std::endl;
    for (const auto &s : tb.accepting_states) {
        os << s << std::endl;
    }
    os << "========================" << std::endl;
    os << "transition table: " << std::endl;
    for (const auto &s : tb.table) {
        os << "from pair :" << s.first << std::endl;
        for (const auto &n : s.second) {
            os << n << std::endl;
        }
        os << "------------------" << std::endl;
    }
    os << "========================" << std::endl;
    return os;
}

std::ostream &operator<<(std::ostream &os, Matcher::Result const &res) {
    os << "start: " << res.start << " length: " << res.length;
    return os;
}

template <typename T>
std::ostream &operator<<(std::ostream &os, std::optional<T> const &opt) {
    if (!opt) {
        os << "<empty>";
    } else {
        os << *opt;
    }
    return os;
}