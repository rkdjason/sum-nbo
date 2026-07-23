#include <boost/multiprecision/cpp_int.hpp>
#include <algorithm>
#include <climits>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using boost::multiprecision::cpp_int;

static const std::string N_HEX =
    "e505004fb5d34eb712d48ff4bbe8d27fc388133c6c0e734001061c0ee0a4edc6"
    "37c04fe8dd376185de8ba04d0ccdbabb93ab7c371b88d92e865eec42b028c61d"
    "d7004ebf2ebb5d69d0a09142be5c9de4da16e514eea318172ecda6cd192073eb"
    "afb1e02d522ec05334590ea6d75960c4937bf64f9700db177a4aa3da6aae6807"
    "e5e32c0d0e428a0db68d299f20c235d84ef459b0cf11828659c31663c9ea8204"
    "4b28152c89a9c36c3ec4303bd36664fd77fb02c58340bdae21120326d83fc017"
    "34bc90048dec9fe35f08c8fdc523abf84a91ec430f49567237c3153a2035ff62"
    "5613b6dc3e6cb14d50e18b8a79b25d678465b3ad02f5b7d818a1e2d635a0baf1";
static const std::string MASK_HEX =
    "fffffffffffffffffffffffff0ffffffffffffffffffffffc00000000000fffe"
    "0000000000000000000003ffe00000000000000000ffffffffffffffffffff"
    "fffffffffffffffffffffff000000000000003ffe00000000000000000001ff"
    "fffffffffffffffffffffffffc3fffffffffffffffffffffffffffffffffffff";
static const std::string PMASK_HEX =
    "ffa360d46885c534d186538170633fafc2c0548a2e24a2c1c0000000000039e2"
    "0000000000000000000000a52000000000000000003e2de4c436d2ca740a6246"
    "99e1a1af94045c63261323c000000000000003bba00000000000000000000e5"
    "0b0bc2461fcbac0726360c2c0809450a9a892cbf1d98ceee48827591ccc593c9";

cpp_int from_hex(const std::string &s) {
    cpp_int x = 0;
    for (char c : s) {
        int v;
        if ('0' <= c && c <= '9') v = c - '0';
        else if ('a' <= c && c <= 'f') v = c - 'a' + 10;
        else if ('A' <= c && c <= 'F') v = c - 'A' + 10;
        else continue;
        x <<= 4;
        x += v;
    }
    return x;
}

int bit(const cpp_int &x, int i) {
    return static_cast<int>((x >> i) & 1);
}

cpp_int low_mask(int k) {
    return (cpp_int(1) << k) - 1;
}

cpp_int mod_pos(cpp_int x, const cpp_int &m) {
    x %= m;
    if (x < 0) x += m;
    return x;
}

cpp_int inv_mod_pow2(const cpp_int &a, int k) {
    if ((a & 1) == 0) throw std::runtime_error("inverse requires odd input");
    cpp_int x = 1;
    int bits = 1;
    while (bits < k) {
        int next = std::min(2 * bits, k);
        cpp_int mod = cpp_int(1) << next;
        x = mod_pos(x * (2 - a * x), mod);
        bits = next;
    }
    return x;
}

cpp_int ceil_div(const cpp_int &a, const cpp_int &b) {
    return (a + b - 1) / b;
}

class CNF {
  public:
    static constexpr int FALSE_LIT = 0;
    static constexpr int TRUE_LIT = INT_MAX;

    explicit CNF(const std::string &tmp_path) : tmp_path_(tmp_path), out_(tmp_path) {
        if (!out_) throw std::runtime_error("cannot open temporary clause file");
    }

    int new_var() { return ++vars_; }
    long long clauses() const { return clauses_; }
    int vars() const { return vars_; }

    static bool is_true(int a) { return a == TRUE_LIT; }
    static bool is_false(int a) { return a == FALSE_LIT; }
    static int lnot(int a) {
        if (a == TRUE_LIT) return FALSE_LIT;
        if (a == FALSE_LIT) return TRUE_LIT;
        return -a;
    }

    void add_clause(std::vector<int> xs) {
        std::vector<int> v;
        v.reserve(xs.size());
        for (int x : xs) {
            if (is_true(x)) return;
            if (is_false(x)) continue;
            v.push_back(x);
        }
        std::sort(v.begin(), v.end(), [](int a, int b) {
            if (std::abs(a) != std::abs(b)) return std::abs(a) < std::abs(b);
            return a < b;
        });
        std::vector<int> u;
        for (int x : v) {
            if (!u.empty() && u.back() == x) continue;
            if (!u.empty() && u.back() == -x) return;
            u.push_back(x);
        }
        for (int x : u) out_ << x << ' ';
        out_ << "0\n";
        ++clauses_;
    }

    void force(int a, bool value) { add_clause({value ? a : lnot(a)}); }

    int mk_and(int a, int b) {
        if (is_false(a) || is_false(b)) return FALSE_LIT;
        if (is_true(a)) return b;
        if (is_true(b)) return a;
        if (a == b) return a;
        if (a == -b) return FALSE_LIT;
        int z = new_var();
        add_clause({-z, a});
        add_clause({-z, b});
        add_clause({z, -a, -b});
        return z;
    }

    int mk_or(int a, int b) {
        if (is_true(a) || is_true(b)) return TRUE_LIT;
        if (is_false(a)) return b;
        if (is_false(b)) return a;
        if (a == b) return a;
        if (a == -b) return TRUE_LIT;
        int z = new_var();
        add_clause({z, -a});
        add_clause({z, -b});
        add_clause({-z, a, b});
        return z;
    }

    int mk_xor2(int a, int b) {
        if (is_false(a)) return b;
        if (is_false(b)) return a;
        if (is_true(a)) return lnot(b);
        if (is_true(b)) return lnot(a);
        if (a == b) return FALSE_LIT;
        if (a == -b) return TRUE_LIT;
        int z = new_var();
        add_clause({a, b, -z});
        add_clause({-a, -b, -z});
        add_clause({a, -b, z});
        add_clause({-a, b, z});
        return z;
    }

    int mk_xnor2(int a, int b) { return lnot(mk_xor2(a, b)); }

    std::pair<int,int> full_adder(int a, int b, int c) {
        std::vector<int> in = {a,b,c};
        int ones = 0;
        std::vector<int> v;
        for (int x : in) {
            if (is_true(x)) ++ones;
            else if (!is_false(x)) v.push_back(x);
        }
        if (v.empty()) {
            return {ones & 1 ? TRUE_LIT : FALSE_LIT, ones >= 2 ? TRUE_LIT : FALSE_LIT};
        }
        if (v.size() == 1) {
            if (ones == 0) return {v[0], FALSE_LIT};
            if (ones == 1) return {lnot(v[0]), v[0]};
            return {v[0], TRUE_LIT};
        }
        if (v.size() == 2) {
            if (v[0] == v[1]) {
                if (ones == 0) return {FALSE_LIT, v[0]};
                return {TRUE_LIT, v[0]};
            }
            if (v[0] == -v[1]) {
                if (ones == 0) return {TRUE_LIT, FALSE_LIT};
                return {FALSE_LIT, TRUE_LIT};
            }
            if (ones == 0) return {mk_xor2(v[0],v[1]), mk_and(v[0],v[1])};
            return {mk_xnor2(v[0],v[1]), mk_or(v[0],v[1])};
        }
        // Three nonconstant inputs.
        if (v[0] == v[1]) return {v[2], v[0]};
        if (v[0] == v[2]) return {v[1], v[0]};
        if (v[1] == v[2]) return {v[0], v[1]};
        if (v[0] == -v[1]) return {lnot(v[2]), v[2]};
        if (v[0] == -v[2]) return {lnot(v[1]), v[1]};
        if (v[1] == -v[2]) return {lnot(v[0]), v[0]};

        int s = new_var();
        for (int mask = 0; mask < 8; ++mask) {
            int va = (mask >> 0) & 1;
            int vb = (mask >> 1) & 1;
            int vc = (mask >> 2) & 1;
            int parity = va ^ vb ^ vc;
            int bad_s = parity ^ 1;
            add_clause({va ? -v[0] : v[0], vb ? -v[1] : v[1], vc ? -v[2] : v[2], bad_s ? -s : s});
        }
        int carry = new_var();
        add_clause({-v[0], -v[1], carry});
        add_clause({-v[0], -v[2], carry});
        add_clause({-v[1], -v[2], carry});
        add_clause({v[0], v[1], -carry});
        add_clause({v[0], v[2], -carry});
        add_clause({v[1], v[2], -carry});
        return {s, carry};
    }

    void finish(const std::string &path) {
        out_.flush();
        out_.close();
        std::ofstream final(path, std::ios::binary);
        final << "p cnf " << vars_ << ' ' << clauses_ << "\n";
        std::ifstream tmp(tmp_path_, std::ios::binary);
        final << tmp.rdbuf();
        final.close();
    }

  private:
    std::string tmp_path_;
    std::ofstream out_;
    int vars_ = 0;
    long long clauses_ = 0;
};

int main(int argc, char **argv) {
    if (argc != 5) {
        std::cerr << "usage: sat_gen2 D output.cnf output.meta temp.clauses\n";
        return 2;
    }
    int d = std::stoi(argv[1]);
    if (d < 0 || d > 15) throw std::runtime_error("D out of range");
    const std::string cnf_path = argv[2];
    const std::string meta_path = argv[3];
    const std::string tmp_path = argv[4];

    const cpp_int N = from_hex(N_HEX);
    const cpp_int MASK = from_hex(MASK_HEX);
    const cpp_int PMASK = from_hex(PMASK_HEX);
    const cpp_int ALL1024 = low_mask(1024);
    const cpp_int D_MASK = cpp_int(15) << 920;

    std::vector<cpp_int> qlow(16);
    const cpp_int M265 = cpp_int(1) << 265;
    for (int a = 0; a < 16; ++a) {
        cpp_int plow = (PMASK | (cpp_int(a) << 150)) & (M265 - 1);
        qlow[a] = (N * inv_mod_pow2(plow, 265)) & (M265 - 1);
    }

    cpp_int remaining_unknown = (ALL1024 ^ MASK) & (ALL1024 ^ D_MASK);
    cpp_int pmin = PMASK | (cpp_int(d) << 920);
    cpp_int pmax = pmin | remaining_unknown;
    cpp_int qmin = ceil_div(N, pmax);
    cpp_int qmax = N / pmin;
    if (qmin > qmax) throw std::runtime_error("bad q interval");

    std::vector<int> q_fixed(1024, -1);
    // Low bits common to all 16 low-hole candidates.
    for (int i = 0; i < 265; ++i) {
        int v = bit(qlow[0], i);
        bool same = true;
        for (int a = 1; a < 16; ++a) if (bit(qlow[a], i) != v) same = false;
        if (same) q_fixed[i] = v;
    }
    // Common prefix implied by the exact p interval after fixing d.
    int common_high = 0;
    for (int i = 1023; i >= 0; --i) {
        if (bit(qmin, i) != bit(qmax, i)) break;
        q_fixed[i] = bit(qmin, i);
        ++common_high;
    }

    CNF cnf(tmp_path);
    std::vector<int> sel(16);
    for (int a = 0; a < 16; ++a) sel[a] = cnf.new_var();

    std::vector<int> p(1024), q(1024);
    for (int i = 0; i < 1024; ++i) {
        if (920 <= i && i < 924) p[i] = ((d >> (i - 920)) & 1) ? CNF::TRUE_LIT : CNF::FALSE_LIT;
        else if (bit(MASK, i)) p[i] = bit(PMASK, i) ? CNF::TRUE_LIT : CNF::FALSE_LIT;
        else p[i] = cnf.new_var();
    }
    for (int i = 0; i < 1024; ++i) {
        if (q_fixed[i] == 0) q[i] = CNF::FALSE_LIT;
        else if (q_fixed[i] == 1) q[i] = CNF::TRUE_LIT;
        else q[i] = cnf.new_var();
    }

    // Exactly one selector, with selector values tied bidirectionally to p[150..153].
    cnf.add_clause(sel);
    for (int a = 0; a < 16; ++a) {
        std::vector<int> reverse_clause = {sel[a]};
        for (int j = 0; j < 4; ++j) {
            bool expected = (a >> j) & 1;
            cnf.add_clause({-sel[a], expected ? p[150+j] : CNF::lnot(p[150+j])});
            reverse_clause.push_back(expected ? CNF::lnot(p[150+j]) : p[150+j]);
        }
        cnf.add_clause(reverse_clause);
    }

    // Once a is selected, p modulo 2^265 is fixed, hence q modulo 2^265 is exact.
    for (int a = 0; a < 16; ++a) {
        for (int i = 0; i < 265; ++i) {
            int expected = bit(qlow[a], i);
            cnf.add_clause({-sel[a], expected ? q[i] : CNF::lnot(q[i])});
        }
    }

    // Build a fully simplified Wallace-tree multiplier.
    constexpr int COLS = 2053;
    std::vector<std::vector<int>> columns(COLS);
    long long products = 0;
    for (int i = 0; i < 1024; ++i) {
        if (CNF::is_false(p[i])) continue;
        for (int j = 0; j < 1024; ++j) {
            if (CNF::is_false(q[j])) continue;
            int t = cnf.mk_and(p[i], q[j]);
            if (!CNF::is_false(t)) columns[i+j].push_back(t);
            ++products;
        }
    }

    long long adders = 0;
    for (int k = 0; k < COLS - 1; ++k) {
        auto &v = columns[k];
        while (v.size() > 2) {
            int a = v.back(); v.pop_back();
            int b = v.back(); v.pop_back();
            int c = v.back(); v.pop_back();
            auto [s, carry] = cnf.full_adder(a,b,c);
            if (!CNF::is_false(s)) v.push_back(s);
            if (!CNF::is_false(carry)) columns[k+1].push_back(carry);
            ++adders;
        }
    }

    int carry = CNF::FALSE_LIT;
    for (int k = 0; k < COLS - 1; ++k) {
        int a = columns[k].size() > 0 ? columns[k][0] : CNF::FALSE_LIT;
        int b = columns[k].size() > 1 ? columns[k][1] : CNF::FALSE_LIT;
        auto [s, next] = cnf.full_adder(a,b,carry);
        bool target = k < 2048 ? bit(N,k) : false;
        cnf.force(s, target);
        carry = next;
        ++adders;
    }
    cnf.force(carry, false);
    for (int x : columns[COLS-1]) cnf.force(x, false);

    cnf.finish(cnf_path);

    std::ofstream meta(meta_path);
    meta << "D " << d << "\n";
    meta << "COMMON_HIGH " << common_high << "\n";
    meta << "QMIN " << qmin << "\nQMAX " << qmax << "\n";
    meta << "P";
    for (int x : p) {
        if (CNF::is_true(x)) meta << " -1";
        else if (CNF::is_false(x)) meta << " 0";
        else meta << ' ' << x;
    }
    meta << "\nQ";
    for (int x : q) {
        if (CNF::is_true(x)) meta << " -1";
        else if (CNF::is_false(x)) meta << " 0";
        else meta << ' ' << x;
    }
    meta << "\nSEL";
    for (int x : sel) meta << ' ' << x;
    meta << "\n";
    meta.close();

    std::cerr << "D=" << d << " common_high=" << common_high
              << " products=" << products << " adders=" << adders
              << " vars=" << cnf.vars() << " clauses=" << cnf.clauses() << "\n";
    return 0;
}
