#pragma once

#include <assert.h>

#include <stdint.h>
#include <stdlib.h>

#include <algorithm>
#include <array>
#include <bit>
#include <cstring>
#include <deque>
#include <iostream>
#include <iterator>
#include <limits>
#include <optional>
#include <set>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <variant>
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

#define NONE 0x00
#define WB 0x01
#define NWB 0x02

struct TransitionState {
    static inline size_t global_id = 0;

    size_t id;
    uint8_t front_entered_modifier;

    TransitionState() : id(global_id++), front_entered_modifier(NONE) {}

    TransitionState(size_t _id, uint8_t _front_entered_modifier)
        : id(_id), front_entered_modifier(_front_entered_modifier) {}
    TransitionState(uint8_t _front_entered_modifier)
        : id(global_id++), front_entered_modifier(_front_entered_modifier) {}

    TransitionState(TransitionState const &other) = default;
    TransitionState &operator=(TransitionState const &other) = default;

    operator size_t() const { return id; }

    bool operator<(TransitionState const &other) const { return id < other.id; }

    bool operator==(TransitionState const &other) const {
        return id == other.id;
    };

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

// template <> struct std::hash<std::pair<TransitionState, char>> {
//     std::size_t operator()(const std::pair<TransitionState, char> &pair)
//     const {
//         return combine(69420, std::hash<TransitionState>()(pair.first),
//                        std::hash<char>()(pair.second));
//     }
// };

template <typename U, typename T> struct std::hash<std::pair<U, T>> {
    std::size_t operator()(const std::pair<U, T> &pair) const {
        return combine(69420, std::hash<U>()(pair.first),
                       std::hash<T>()(pair.second));
    }
};

using StateChar = std::pair<TransitionState, char>;

struct TransitionTable {
    // std::vector<std::array<TransitionState, 256>> table;
    // std::vector<TransitionState> starting_states;
    // std::vector<bool> accepting_states;

    using State = size_t;
    using Row = std::array<std::vector<State>, 256>;

    size_t num_states;
    std::vector<Row> trans_rows;
    std::vector<State> starting_states;
    std::vector<bool> accepting_states;
    std::vector<uint8_t> front_entered_modifiers;
    std::vector<State> acc_state_list;
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
    for (size_t state_id = 0; state_id < tb.trans_rows.size(); ++state_id) {
        os << "-------> state id: " << state_id << std::endl;
        os << "-------> front entered mod: "
           << (size_t)tb.front_entered_modifiers[state_id] << std::endl;

        for (size_t char_val = 0; char_val < 128; ++char_val) {
            if (tb.trans_rows[state_id][char_val].empty()) {
                continue;
            }
            for (auto const &next_st : tb.trans_rows[state_id][char_val]) {
                os << "{targ: " << next_st << " ; transition char : ["
                   << (char)char_val << "] ; modifier: "
                   << (size_t)tb.front_entered_modifiers[next_st] << " }"
                   << std::endl;
            }
        }
    }
    os << "transition table end ====================" << std::endl;

    os << "acc list ====================" << std::endl;
    for (auto const &s : tb.acc_state_list) {
        std::cout << s << std::endl;
    }
    os << "acc list end ====================" << std::endl;

    return os;
}

/*
    right now i'd rather take in a string view so we dont make needless
   copies but this means the lifetime of our matcher is tied to the string
   view
*/
class Matcher {

    using TState = TransitionTable::State;
    using TRow = TransitionTable::Row;

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

    static bool is_word_char(char ch) {
        return ch != ' ' && ch != '\t' && ch != 3 && ch != 2;
    }

    // intermediate state?
    struct State {
        size_t offset;
        TransitionTable::State nfa_state;

        State(size_t _offset, TState _state)
            : offset(_offset), nfa_state(_state) {}

        State(State &&other) = default;
        State(State const &other) = default;
        State &operator=(State const &other) = default;
        State &operator=(State &&other) = default;

        bool operator==(State const &other) const = default;

        // using this just to define it here cause im lazy
        friend std::ostream &operator<<(std::ostream &os, State s) {
            os << "{offset: " << s.offset << " nfa_state: " << s.nfa_state
               << "}" << std::endl;
            return os;
        }
    };

    struct StateHash {
        size_t operator()(State const &r) const {
            return combine(69420, std::hash<size_t>()(r.offset),
                           std::hash<TState>()(r.nfa_state));
        }
    };

    // mostly screw around with this to change how the
    // states are updated?
    struct StateList {
        static const size_t NOFFSET = (size_t)-1;
        std::vector<size_t> active_state_to_offset;
        std::vector<size_t> buffer_state_to_offset;

        StateList(size_t fixed_size)
            : active_state_to_offset(fixed_size, NOFFSET),
              buffer_state_to_offset(fixed_size, NOFFSET) {}

        void insert_state(size_t state_idx, size_t offset) {
            if (active_state_to_offset[state_idx] == NOFFSET ||
                offset < active_state_to_offset[state_idx]) {
                active_state_to_offset[state_idx] = offset;
            }
        }

        void insert_buffer_state(size_t state_idx, size_t offset) {
            if (buffer_state_to_offset[state_idx] == NOFFSET ||
                offset < buffer_state_to_offset[state_idx]) {
                buffer_state_to_offset[state_idx] = offset;
            }
        }

        void progress_states(TransitionTable const &tran_tab, char prev_char,
                             char curr_char) {
            for (size_t idx = 0; idx < buffer_state_to_offset.size(); ++idx) {
                buffer_state_to_offset[idx] = NOFFSET;
            }
            for (size_t state_idx = 0;
                 state_idx < active_state_to_offset.size(); ++state_idx) {

                if (active_state_to_offset[state_idx] == NOFFSET) {
                    continue;
                }

                uint8_t mod = tran_tab.front_entered_modifiers[state_idx];
                for (auto next_state :
                     tran_tab.trans_rows[state_idx][curr_char]) {

                    if (mod == NONE) {
                        insert_buffer_state(next_state,
                                            active_state_to_offset[state_idx]);
                        continue;
                    }

                    if (mod == (WB | NWB)) {
                        continue;
                    }

                    bool is_word_boundary =
                        (is_word_char(prev_char) != is_word_char(curr_char));

                    if ((mod == WB && is_word_boundary) ||
                        (mod == NWB && !is_word_boundary)) {
                        insert_buffer_state(next_state,
                                            active_state_to_offset[state_idx]);
                    }
                }
            }
            std::swap(active_state_to_offset, buffer_state_to_offset);
        }

        size_t &operator[](size_t state_idx) {
            return active_state_to_offset[state_idx];
        }
    };

  public:
    Matcher(TransitionTable _tb) : tb(std::move(_tb)) {}

    std::optional<Result> greedy_match_view(std::string_view str) {

        // should we write our own contains for the states?
        // idea: for small enough stuff we just std::vector<bool>
        StateList state_list(tb.num_states);
        std::optional<Result> result;

        auto update_result = [&](size_t offset, size_t length) {
            if (!result) {
                result = Result{offset, length};
            } else if (offset < result->start ||
                       (offset == result->start && length > result->length)) {
                result = Result{offset, length};
            }
        };

        auto process_accepts = [&, update_result](size_t final_offset,
                                                  char curr_char,
                                                  char next_char) {
            for (auto const &acc_state : tb.acc_state_list) {
                size_t offset = state_list[acc_state];

                if (offset == StateList::NOFFSET) {
                    continue;
                }

                if (offset == final_offset) {
                    continue;
                }

                uint8_t mod = tb.front_entered_modifiers[acc_state];
                if (mod == NONE) {
                    update_result(offset, final_offset - offset);
                    continue;
                }

                if (mod == (WB | NWB)) {
                    continue;
                }

                bool word_boundary =
                    (is_word_char(curr_char) != is_word_char(next_char));
                if (mod == WB && word_boundary) {
                    update_result(offset, final_offset - offset);
                } else if (mod == NWB && !word_boundary) {
                    update_result(offset, final_offset - offset);
                }
            }
        };

        auto insert_states = [&](size_t offset, char curr_char,
                                 char prev_char) {
            for (auto const &state : tb.starting_states) {
                uint8_t mod = tb.front_entered_modifiers[state];
                if (mod == NONE) {
                    state_list.insert_state(state, offset);
                    continue;
                }

                if (mod == (WB | NWB)) {
                    continue;
                }

                bool word_boundary =
                    (is_word_char(curr_char) != is_word_char(prev_char));
                if ((mod == WB && word_boundary) ||
                    (mod == NWB && !word_boundary)) {
                    state_list.insert_state(state, offset);
                }
            }
        };

        // feed the line start char
        char prev_char = 2;
        // start the initial list
        insert_states(0, 2, 2);
        // feed the bol character
        state_list.progress_states(tb, 2, 2);

        for (size_t idx = 0; idx < str.length(); ++idx) {
            if (!result) {
                insert_states(idx, str[idx], prev_char);
            } else {
                // do we want to do pruning?
            }
            state_list.progress_states(tb, prev_char, str[idx]);
            // technically it's after idx?
            process_accepts(idx + 1, str[idx],
                            (idx + 1 >= str.length()) ? 3 : str[idx + 1]);

            prev_char = str[idx];
        }

        state_list.progress_states(tb, prev_char, 3);
        process_accepts(str.length(), str.back(), 3);
        return result;
    }

    friend std::ostream &operator<<(std::ostream &os, Matcher const &m) {
        os << m.tb << std::endl;
        return os;
    }
};

struct MatcherBuilder;
std::ostream &operator<<(std::ostream &os, MatcherBuilder const &tb);

void prune_states(MatcherBuilder &mb);
std::unordered_map<TransitionState, TransitionState>
resequence_states(MatcherBuilder &mb);
TransitionTable compile_transition_table(MatcherBuilder &mb);

struct MatcherBuilder {

    std::vector<TransitionState> starting_states;
    std::vector<TransitionState> accepting_states;
    std::unordered_map<StateChar, std::vector<TransitionState>> table;
    std::unordered_map<TransitionState, std::vector<StateChar>> reverse_table;

    MatcherBuilder()
        : starting_states(1, TransitionState()),
          accepting_states{starting_states.front()} {}

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

    MatcherBuilder get_matcher_builder() {
        if (!result) {
            throw std::runtime_error("there isnt anything parsed");
        }
        return *result;
    }

    // needs to be called after a successful parse
    TransitionTable get_compiled_table() {
        if (!result) {
            throw std::runtime_error("there isnt anything parsed");
        }

        prune_states(*result);
        return compile_transition_table(*result);
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
            mb.boundary_modify((sv_front().get_ch() == 'b') ? WB : NWB);
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
            return parse_help(accum_b);
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
            accum_b += mb;
            sv_advance();
            MatcherBuilder sub_expr;
            bool res = parse_help(sub_expr);
            if (!res) {
                return false;
            } else {
                accum_b |= sub_expr;
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

struct TSCMP {
    bool operator()(TransitionState const &a, TransitionState const &b) const {
        return a.id < b.id;
    }
};

// consumes the matcher builder
TransitionTable compile_transition_table(MatcherBuilder &mb) {
    // takes a resequenced table and creates a vec of TransitionTable::row
    std::unordered_map<TransitionState, TransitionState> old_to_new =
        resequence_states(mb);

    std::vector<TransitionTable::Row> final_table(old_to_new.size());
    std::vector<uint8_t> final_modifiers(old_to_new.size());

    for (auto const &[sc_pair, next_row] : mb.table) {
        // should just be our first time setting it
        assert(final_table.at(sc_pair.first)[sc_pair.second].empty());
        final_table[sc_pair.first][sc_pair.second].reserve(next_row.size());

        for (auto const &next_st : next_row) {
            final_table[sc_pair.first][sc_pair.second].push_back(next_st);
        }
    }

    for (auto const &otn : old_to_new) {
        final_modifiers[otn.second.id] = otn.second.get_front_modifier();
    }

    std::vector<bool> final_accepts(old_to_new.size(), false);
    for (auto const &acc_st : mb.accepting_states) {
        final_accepts[acc_st] = true;
    }

    std::vector<TransitionTable::State> final_starts;
    final_starts.reserve(mb.starting_states.size());
    for (auto const &start_st : mb.starting_states) {
        final_starts.push_back(start_st.id);
    }

    std::vector<TransitionTable::State> acc_state_list;
    acc_state_list.reserve(mb.accepting_states.size());
    for (auto const &acc_st : mb.accepting_states) {
        acc_state_list.push_back(acc_st.id);
    }

    return {old_to_new.size(), final_table,     final_starts,
            final_accepts,     final_modifiers, acc_state_list};
}

// Given a matcher builder, resequences the states
// so that it's sequential
std::unordered_map<TransitionState, TransitionState>
resequence_states(MatcherBuilder &mb) {
    std::unordered_map<TransitionState, TransitionState> old_to_new;

    size_t running_id = 0;
    auto add_old_to_new_entry = [&](TransitionState ts) {
        auto it = old_to_new.find(ts);
        if (it == old_to_new.end()) {
            old_to_new.insert({ts, {running_id++, ts.get_front_modifier()}});
        }
    };

    // we first need to collect all the states, go through
    // starting, accepting and the table itself

    for (TransitionState ts : mb.starting_states) {
        add_old_to_new_entry(ts);
    }

    for (TransitionState ts : mb.accepting_states) {
        add_old_to_new_entry(ts);
    }

    for (auto &row : mb.table) {
        add_old_to_new_entry(row.first.first);
        for (auto &next : row.second) {
            add_old_to_new_entry(next);
        }
    }

    // now we in place transform we can in place transform for
    // the other states at least
    for (TransitionState &ts : mb.starting_states) {
        ts = old_to_new.at(ts);
    }

    for (TransitionState &ts : mb.accepting_states) {
        ts = old_to_new.at(ts);
    }

    std::unordered_map<StateChar, std::vector<TransitionState>> new_table;
    std::vector<TransitionState> new_row;
    for (auto transition_row : mb.table) {
        new_row.clear();
        new_row.reserve(transition_row.second.size());
        for (auto next : transition_row.second) {
            new_row.push_back(old_to_new.at(next));
        }
        new_table[{old_to_new.at(transition_row.first.first),
                   transition_row.first.second}] = new_row;
    }

    mb.table = std::move(new_table);
    return old_to_new;

    // Do we ignore the reverse table? feels like we don't need
    // it after this step
}

void prune_states(MatcherBuilder &mb) {
    std::unordered_set<std::pair<StateChar, TransitionState>> all_edges;
    for (auto table_pair : mb.table) {
        for (auto targ_state : table_pair.second) {
            all_edges.insert({table_pair.first, targ_state});
        }
    }

    // create a transition state queue
    std::deque<TransitionState> ts_queue{mb.accepting_states.begin(),
                                         mb.accepting_states.end()};

    std::unordered_set<TransitionState> visited;

    while (!ts_queue.empty()) {
        TransitionState curr_ts = ts_queue.front();
        ts_queue.pop_front();

        // if we have been visited we can trust
        // our edges have been removed
        if (visited.contains(curr_ts)) {
            continue;
        }

        // mark all outgoing edges as important
        if (!mb.reverse_table.contains(curr_ts)) {
            continue;
        }
        for (auto anc_pair : mb.reverse_table.at(curr_ts)) {
            all_edges.erase({anc_pair, curr_ts});
            ts_queue.push_back(anc_pair.first);
        }
    }

    auto remove_from_vec = [](std::vector<TransitionState> &vec,
                              TransitionState targ_state) -> bool {
        auto it = vec.begin();
        while (it != vec.end()) {
            if (*it == targ_state) {
                vec.erase(it);
                return true;
            }
            ++it;
        }
        return false;
        // throw std::runtime_error("I need it removed and it's not here");
    };

    // now we delete all edges
    for (auto edge : all_edges) {
        // each edge should only be removed once
        assert(remove_from_vec(mb.table.at(edge.first), edge.second));
        if (mb.table.at(edge.first).empty()) {
            mb.table.erase(edge.first);
        }
    }
}
