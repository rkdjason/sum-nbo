#define main sat_gen2_renamed_main
#include "sat_gen2.cpp"
#undef main

#include <map>

static std::vector<int> weighted_sum(CNF &cnf, const std::vector<int> &unit_terms,
                                     const std::vector<int> &binary_carry) {
    std::vector<std::vector<int>> levels(24);
    for (int x : unit_terms) if (!CNF::is_false(x)) levels[0].push_back(x);
    for (size_t i = 0; i < binary_carry.size(); ++i) {
        if (i >= levels.size()) levels.resize(i + 8);
        if (!CNF::is_false(binary_carry[i])) levels[i].push_back(binary_carry[i]);
    }
    for (size_t r = 0; r < levels.size(); ++r) {
        auto &v = levels[r];
        while (v.size() >= 3) {
            int a = v.back(); v.pop_back();
            int b = v.back(); v.pop_back();
            int c = v.back(); v.pop_back();
            auto [s, carry] = cnf.full_adder(a,b,c);
            if (!CNF::is_false(s)) v.push_back(s);
            if (!CNF::is_false(carry)) {
                if (r + 1 >= levels.size()) levels.resize(levels.size() + 8);
                levels[r+1].push_back(carry);
            }
        }
        if (v.size() == 2) {
            int a = v.back(); v.pop_back();
            int b = v.back(); v.pop_back();
            auto [s, carry] = cnf.full_adder(a,b,CNF::FALSE_LIT);
            if (!CNF::is_false(s)) v.push_back(s);
            if (!CNF::is_false(carry)) {
                if (r + 1 >= levels.size()) levels.resize(levels.size() + 8);
                levels[r+1].push_back(carry);
            }
        }
    }
    std::vector<int> result;
    result.reserve(levels.size());
    for (auto &v : levels) result.push_back(v.empty() ? CNF::FALSE_LIT : v[0]);
    while (!result.empty() && CNF::is_false(result.back())) result.pop_back();
    return result;
}

static std::vector<int> add_one(CNF &cnf, const std::vector<int> &x, int one) {
    std::vector<int> out;
    out.reserve(x.size() + 1);
    int carry = one;
    for (int b : x) {
        if (CNF::is_false(carry)) {
            out.push_back(b);
        } else {
            out.push_back(cnf.mk_xor2(b, carry));
            carry = cnf.mk_and(b, carry);
        }
    }
    if (!CNF::is_false(carry)) out.push_back(carry);
    while (!out.empty() && CNF::is_false(out.back())) out.pop_back();
    return out;
}

int main(int argc, char **argv) {
    if (argc != 5) {
        std::cerr << "usage: sat_dir D output.cnf output.meta temp.clauses\n";
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

    std::vector<int> q_fixed(1024, -1);
    for (int i = 0; i < 265; ++i) {
        int v = bit(qlow[0], i);
        bool same = true;
        for (int a = 1; a < 16; ++a) if (bit(qlow[a], i) != v) same = false;
        if (same) q_fixed[i] = v;
    }
    int common_high = 0;
    for (int i = 1023; i >= 0; --i) {
        if (bit(qmin,i) != bit(qmax,i)) break;
        q_fixed[i] = bit(qmin,i);
        ++common_high;
    }

    CNF cnf(tmp_path);
    std::vector<int> sel(16);
    for (int a = 0; a < 16; ++a) sel[a] = cnf.new_var();

    std::vector<int> p(1024), q(1024, CNF::FALSE_LIT);
    for (int i = 0; i < 1024; ++i) {
        if (920 <= i && i < 924) p[i] = ((d >> (i-920)) & 1) ? CNF::TRUE_LIT : CNF::FALSE_LIT;
        else if (bit(MASK,i)) p[i] = bit(PMASK,i) ? CNF::TRUE_LIT : CNF::FALSE_LIT;
        else p[i] = cnf.new_var();
    }

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

    std::vector<int> carry;
    long long product_gates = 0;
    long long column_terms = 0;

    // Low half: derive q_k functionally from p and all previously derived q bits.
    for (int k = 0; k < 1024; ++k) {
        std::vector<int> terms;
        terms.reserve(k + 1);
        for (int i = 1; i <= k && i < 1024; ++i) {
            int j = k - i;
            int t = cnf.mk_and(p[i], q[j]);
            if (!CNF::is_false(t)) terms.push_back(t);
            ++product_gates;
        }
        column_terms += terms.size();
        std::vector<int> total = weighted_sum(cnf, terms, carry);
        int t0 = total.empty() ? CNF::FALSE_LIT : total[0];
        bool nbit = bit(N,k);
        int computed_q = nbit ? CNF::lnot(t0) : t0;
        if (q_fixed[k] >= 0) {
            cnf.force(computed_q, q_fixed[k]);
            q[k] = q_fixed[k] ? CNF::TRUE_LIT : CNF::FALSE_LIT;
        } else {
            int qb = cnf.new_var();
            cnf.add_clause({-qb, computed_q});
            cnf.add_clause({qb, CNF::lnot(computed_q)});
            q[k] = qb;
        }
        int low_carry = nbit ? CNF::FALSE_LIT : t0;
        std::vector<int> upper;
        if (total.size() > 1) upper.assign(total.begin()+1, total.end());
        carry = add_one(cnf, upper, low_carry);
    }

    // Exact low q bits selected by the four-bit low hole.
    for (int a = 0; a < 16; ++a) {
        for (int i = 0; i < 265; ++i) {
            bool expected = bit(qlow[a],i);
            cnf.add_clause({-sel[a], expected ? q[i] : CNF::lnot(q[i])});
        }
    }

    // High half: all q bits now exist; enforce the remaining 1024 product bits.
    for (int k = 1024; k < 2048; ++k) {
        std::vector<int> terms;
        int imin = std::max(0, k - 1023);
        int imax = std::min(1023, k);
        terms.reserve(imax - imin + 1);
        for (int i = imin; i <= imax; ++i) {
            int j = k - i;
            if (j < 0 || j >= 1024) continue;
            int t = cnf.mk_and(p[i],q[j]);
            if (!CNF::is_false(t)) terms.push_back(t);
            ++product_gates;
        }
        column_terms += terms.size();
        std::vector<int> total = weighted_sum(cnf, terms, carry);
        int outbit = total.empty() ? CNF::FALSE_LIT : total[0];
        cnf.force(outbit, bit(N,k));
        carry.clear();
        if (total.size() > 1) carry.assign(total.begin()+1,total.end());
    }
    for (int b : carry) cnf.force(b, false);

    cnf.finish(cnf_path);
    std::ofstream meta(meta_path);
    meta << "D " << d << "\nCOMMON_HIGH " << common_high << "\n";
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
              << " product_gates=" << product_gates << " terms=" << column_terms
              << " vars=" << cnf.vars() << " clauses=" << cnf.clauses() << "\n";
    return 0;
}
