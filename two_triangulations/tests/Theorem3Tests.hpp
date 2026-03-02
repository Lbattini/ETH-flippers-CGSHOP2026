// Theorem3Tests.hpp
#pragma once
#include <vector>
#include <array>
#include <stdexcept>
#include <string>
#include <iostream>
#include <unordered_set>
#include <algorithm>
#include <cmath>

#include "../include/CutTriangulation.hpp"       // your file (namespaces & helpers)
#include "../include/Theorem2Transformer.hpp"    // polygon engine used internally
#include "../include/Theorem3Transformer.hpp"    // point-set driver we just wrote

namespace theorem3_tests {

// ---------- tiny test helpers ----------
    [[noreturn]] inline void fail(const std::string& msg) { throw std::runtime_error(msg); }
    inline void assert_true(bool cond, const std::string& msg) { if (!cond) fail(msg); }

    using Point = ::Point;   // from CutTriangulation.hpp
    using Edge  = ::Edge;

// normalized edge key
    static inline long long ekey(int a,int b){ if(a>b) std::swap(a,b); return ( (long long)a<<32 ) ^ (unsigned long long)b; }

    static inline std::unordered_set<long long> edge_set_of(const std::vector<Edge>& E){
        std::unordered_set<long long> S; S.reserve(E.size()*2);
        for (auto e: E) if (e.u!=e.v) S.insert(ekey(e.u,e.v));
        return S;
    }

    static inline void require_triangulation_faces_extractable(const std::vector<Point>& P,
                                                               const std::vector<Edge>&  E,
                                                               const char* name)
    {
        auto F = triangulation_cut::extractFacesAsTriangles(P, E);
        // Basic sanity: #triangles should be >= 1 for n>=3; we won’t over-police here.
        if (P.size() >= 3) {
            assert_true(!F.triangles.empty(), std::string(name) + ": could not extract any triangular faces.");
        }
    }

// ---------- point-set builders ----------
// Make bi-chain point set (x-monotone):
//  lower y=0: L0..L_{p-1} (x=0..p-1)
//  upper y=1: U0..U_{m-1} (x=0..m-1)
// Returns indices arrays U (left->right) and L (left->right).
    static void make_bichain_points(int m_upper, int p_lower,
                                    std::vector<Point>& P,
                                    std::vector<int>& U,
                                    std::vector<int>& L)
    {
        assert_true(m_upper>=2 && p_lower>=2, "Need at least 2 per chain.");
        P.clear(); U.clear(); L.clear();
        P.reserve(m_upper + p_lower);

        for (int i=0;i<p_lower;++i) { P.push_back(Point{double(i), 0.0}); L.push_back((int)P.size()-1); }
        for (int j=0;j<m_upper;++j){ P.push_back(Point{double(j), 1.0}); U.push_back((int)P.size()-1); }
    }

// Triangulation “Strip A” (edges):
//   core quads (L_{i-1}, L_i, U_i, U_{i-1}) use diagonal (L_i, U_{i-1})
// plus caps.
    static void build_strip_A_edges(const std::vector<int>& U,
                                    const std::vector<int>& L,
                                    std::vector<Edge>& E)
    {
        E.clear();
        auto add=[&](int a,int b){ if(a!=b) E.push_back({std::min(a,b),std::max(a,b)}); };
        const int m=(int)U.size(), p=(int)L.size(), k=std::min(m,p);

        // boundary polygon: lower chain + upper chain + hull closes implicitly via edges we add below
        for (int i=1;i<p;++i) add(L[i-1], L[i]);
        for (int j=1;j<m;++j) add(U[j-1], U[j]);
        add(L.front(), U.front());   // left cap boundary
        add(L.back(),  U.back());    // right cap boundary

        // core strips -> two triangles per cell
        for (int i=1;i<k;++i) {
            int L_im1=L[i-1], L_i=L[i], U_im1=U[i-1], U_i=U[i];
            // Triangles: (L_{i-1}, L_i, U_{i-1}) and (L_i, U_i, U_{i-1})
            add(L_im1, L_i); add(L_i, U_im1); add(U_im1, L_im1);
            add(L_i, U_i);   add(U_i, U_im1); add(U_im1, L_i);
        }
        // top cap (m>p): fan from L_{p-1}
        if (m>p) {
            int base=L[p-1];
            for (int j=p;j<m;++j) { add(base, U[j-1]); add(U[j-1], U[j]); add(U[j], base); }
        }
        // bottom cap (p>m): fan from U_{m-1}
        if (p>m) {
            int base=U[m-1];
            for (int i=m;i<p;++i) { add(L[i-1], L[i]); add(L[i], base); add(base, L[i-1]); }
        }

        // dedup
        auto S=edge_set_of(E);
        std::vector<Edge> ded; ded.reserve(S.size());
        for (auto k : S) ded.push_back({int(k>>32), int(k & 0xffffffff)});
        E.swap(ded);
    }

// Triangulation “Strip B” (edges):
//   core quads use diagonal (L_{i-1}, U_i)
    static void build_strip_B_edges(const std::vector<int>& U,
                                    const std::vector<int>& L,
                                    std::vector<Edge>& E)
    {
        E.clear();
        auto add=[&](int a,int b){ if(a!=b) E.push_back({std::min(a,b),std::max(a,b)}); };
        const int m=(int)U.size(), p=(int)L.size(), k=std::min(m,p);

        for (int i=1;i<p;++i) add(L[i-1], L[i]);
        for (int j=1;j<m;++j) add(U[j-1], U[j]);
        add(L.front(), U.front());
        add(L.back(),  U.back());

        for (int i=1;i<k;++i) {
            int L_im1=L[i-1], L_i=L[i], U_im1=U[i-1], U_i=U[i];
            // Triangles: (L_{i-1}, U_i, U_{i-1}) and (L_{i-1}, L_i, U_i)
            add(L_im1, U_i); add(U_i, U_im1); add(U_im1, L_im1);
            add(L_im1, L_i); add(L_i, U_i);   add(U_i, L_im1);
        }
        if (m>p) {
            int base=L[p-1];
            for (int j=p;j<m;++j) { add(base, U[j-1]); add(U[j-1], U[j]); add(U[j], base); }
        }
        if (p>m) {
            int base=U[m-1];
            for (int i=m;i<p;++i) { add(L[i-1], L[i]); add(L[i], base); add(base, L[i-1]); }
        }

        auto S=edge_set_of(E);
        std::vector<Edge> ded; ded.reserve(S.size());
        for (auto k : S) ded.push_back({int(k>>32), int(k & 0xffffffff)});
        E.swap(ded);
    }

// ---------- Tests ----------

    static void test_theorem3_balanced_6x6() {
        std::cout << "[TEST] Theorem 3 — balanced point set 6x6: Strip A → Strip B\n";

        std::vector<Point> P; std::vector<int> U, L;
        make_bichain_points(6,6, P, U, L);

        std::vector<Edge> E_start, E_target;
        build_strip_A_edges(U, L, E_start);
        build_strip_B_edges(U, L, E_target);

        require_triangulation_faces_extractable(P, E_start, "start");
        require_triangulation_faces_extractable(P, E_target, "target");

        auto stats = theorem3_pset::transform(P, E_start, E_target, /*dbg_out=*/nullptr);

        auto A = edge_set_of(E_start);
        auto B = edge_set_of(E_target);
        assert_true(A == B, "Edge sets differ after Theorem 3 (6x6).");
        std::cout << "  - depth=" << stats.depth
                  << ", nodes=" << stats.nodes
                  << ", prop1_phases=" << stats.prop1_phases
                  << ", lemma4_calls=" << stats.lemma4_calls
                  << ", lemma5_calls=" << stats.lemma5_calls
                  << ", flips_total=" << stats.flips_total << " ✅\n";
    }

    static void test_theorem3_unbalanced_7x4() {
        std::cout << "[TEST] Theorem 3 — unbalanced point set 7x4: Strip A → Strip B\n";

        std::vector<Point> P; std::vector<int> U, L;
        make_bichain_points(7,4, P, U, L);

        std::vector<Edge> E_start, E_target;
        build_strip_A_edges(U, L, E_start);
        build_strip_B_edges(U, L, E_target);

        require_triangulation_faces_extractable(P, E_start, "start");
        require_triangulation_faces_extractable(P, E_target, "target");

        auto stats = theorem3_pset::transform(P, E_start, E_target, /*dbg_out=*/nullptr);

        auto A = edge_set_of(E_start);
        auto B = edge_set_of(E_target);
        assert_true(A == B, "Edge sets differ after Theorem 3 (7x4).");
        std::cout << "  - depth=" << stats.depth
                  << ", nodes=" << stats.nodes
                  << ", prop1_phases=" << stats.prop1_phases
                  << ", lemma4_calls=" << stats.lemma4_calls
                  << ", lemma5_calls=" << stats.lemma5_calls
                  << ", flips_total=" << stats.flips_total << " ✅\n";
    }

    static void test_theorem3_jittered_6x5() {
        std::cout << "[TEST] Theorem 3 — jittered 6x5: Strip A → Strip B\n";

        std::vector<Point> P; std::vector<int> U, L;
        make_bichain_points(6,5, P, U, L);

        // add small vertical jitter to the upper chain to avoid accidental collinearities
        for (size_t i=0;i<U.size();++i) {
            P[ U[i] ].y += ((i%2)? +1e-9 : -1e-9);
        }

        std::vector<Edge> E_start, E_target;
        build_strip_A_edges(U, L, E_start);
        build_strip_B_edges(U, L, E_target);

        require_triangulation_faces_extractable(P, E_start, "start");
        require_triangulation_faces_extractable(P, E_target, "target");

        auto stats = theorem3_pset::transform(P, E_start, E_target, /*dbg_out=*/&std::cout);

        auto A = edge_set_of(E_start);
        auto B = edge_set_of(E_target);
        assert_true(A == B, "Edge sets differ after Theorem 3 (jittered 6x5).");
        std::cout << "  - depth=" << stats.depth
                  << ", nodes=" << stats.nodes
                  << ", prop1_phases=" << stats.prop1_phases
                  << ", lemma4_calls=" << stats.lemma4_calls
                  << ", lemma5_calls=" << stats.lemma5_calls
                  << ", flips_total=" << stats.flips_total << " ✅\n";
    }

    inline void run_all_theorem3_tests() {
        std::cout << "===== Running Theorem 3 (Point-Set) Tests =====\n\n";
        test_theorem3_balanced_6x6();
        std::cout << "\n";
        test_theorem3_unbalanced_7x4();
        std::cout << "\n";
        test_theorem3_jittered_6x5();
        std::cout << "\nAll Theorem 3 tests passed ✅\n";
    }

} // namespace theorem3_tests
