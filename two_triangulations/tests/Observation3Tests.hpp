// Observation3Tests.hpp
#pragma once
#include <vector>
#include <array>
#include <stdexcept>
#include <string>
#include <iostream>
#include <algorithm>

#include "../include/TriCommon.hpp"
#include "../include/Proposition1Inserter.hpp"
#include "../include/Theorem2Transformer.hpp"

namespace observation3_tests {

    [[noreturn]] inline void fail(const std::string& msg) { throw std::runtime_error(msg); }
    inline void assert_true(bool cond, const std::string& msg) { if (!cond) fail(msg); }

    using Point = cgshop::Point;
    using Face  = std::array<int,3>;
    using EdgeKey = cgshop::EdgeKey;
    using EdgeKeyHash = cgshop::EdgeKeyHash;

// ========== helpers (bi-chain + ε self-touch) ==========

// Standard bi-chain boundary (CCW): L0..L_{p-1}, U_{m-1}..U0
    static void make_bichain_polygon_with_boundary(
            int m_upper, int p_lower,
            std::vector<Point>& pts,
            std::vector<int>& upper_lr,
            std::vector<int>& lower_lr,
            std::vector<int>& boundary_ccw)
    {
        assert_true(m_upper >= 2 && p_lower >= 2, "Need at least 2 vertices per chain.");
        pts.clear(); upper_lr.clear(); lower_lr.clear(); boundary_ccw.clear();
        pts.reserve(m_upper + p_lower);

        for (int i = 0; i < p_lower; ++i) {
            pts.push_back(Point{double(i), 0.0});          // lower y=0
            lower_lr.push_back((int)pts.size()-1);
        }
        for (int j = 0; j < m_upper; ++j) {
            pts.push_back(Point{double(j), 1.0});          // upper y=1
            upper_lr.push_back((int)pts.size()-1);
        }
        for (int v : lower_lr) boundary_ccw.push_back(v);
        for (int k = (int)upper_lr.size()-1; k >= 0; --k) boundary_ccw.push_back(upper_lr[k]);
    }

// ε-regularized external touch: nudge the rightmost upper vertex almost onto L_last
// keeping interior simply connected and visibility unchanged.
    static void apply_epsilon_touch_right(std::vector<Point>& pts,
                                          int upper_last, int lower_last,
                                          double eps = 1e-12)
    {
        // Move U_last to (x_L_last, eps) — arbitrarily close to L_last=(x_L_last, 0).
        pts[ upper_last ].x = pts[ lower_last ].x;
        pts[ upper_last ].y = eps;
    }

// Strip A triangulation (diagonal (L_i, U_{i-1}) in each quad)
    static void build_strip_A(const std::vector<int>& U_lr,
                              const std::vector<int>& L_lr,
                              std::vector<Face>& faces)
    {
        faces.clear();
        const int m = (int)U_lr.size();
        const int p = (int)L_lr.size();
        const int k = std::min(m, p);

        for (int i = 1; i < k; ++i) {
            int L_im1 = L_lr[i-1], L_i = L_lr[i];
            int U_im1 = U_lr[i-1], U_i = U_lr[i];
            faces.push_back({L_im1, L_i, U_im1});
            faces.push_back({L_i, U_i, U_im1});
        }
        if (m > p) {
            int base = L_lr[p-1];
            for (int j = p; j < m; ++j) faces.push_back({base, U_lr[j-1], U_lr[j]});
        } else if (p > m) {
            int base = U_lr[m-1];
            for (int i = m; i < p; ++i) faces.push_back({L_lr[i-1], L_lr[i], base});
        }

        const int n = m + p;
        assert_true((int)faces.size() == n - 2, "Strip A produced wrong face count.");
    }

// Strip B triangulation (diagonal (L_{i-1}, U_i) in each quad)
    static void build_strip_B(const std::vector<int>& U_lr,
                              const std::vector<int>& L_lr,
                              std::vector<Face>& faces)
    {
        faces.clear();
        const int m = (int)U_lr.size();
        const int p = (int)L_lr.size();
        const int k = std::min(m, p);

        for (int i = 1; i < k; ++i) {
            int L_im1 = L_lr[i-1], L_i = L_lr[i];
            int U_im1 = U_lr[i-1], U_i = U_lr[i];
            faces.push_back({L_im1, U_i, U_im1});
            faces.push_back({L_im1, L_i, U_i});
        }
        if (m > p) {
            int base = L_lr[p-1];
            for (int j = p; j < m; ++j) faces.push_back({base, U_lr[j-1], U_lr[j]});
        } else if (p > m) {
            int base = U_lr[m-1];
            for (int i = m; i < p; ++i) faces.push_back({L_lr[i-1], L_lr[i], base});
        }

        const int n = m + p;
        assert_true((int)faces.size() == n - 2, "Strip B produced wrong face count.");
    }

    static std::unordered_set<EdgeKey, EdgeKeyHash>
    edges_from_faces(const std::vector<Face>& faces)
    {
        std::unordered_set<EdgeKey, EdgeKeyHash> E;
        E.reserve(faces.size()*3);
        for (const auto& T : faces) {
            E.insert(EdgeKey(T[0],T[1]));
            E.insert(EdgeKey(T[1],T[2]));
            E.insert(EdgeKey(T[2],T[0]));
        }
        return E;
    }

    static void assert_same_edge_sets(const std::vector<Face>& A,
                                      const std::vector<Face>& Bf,
                                      const std::string& msg_if_not)
    {
        auto EA = edges_from_faces(A);
        auto EB = edges_from_faces(Bf);
        if (EA.size() != EB.size()) {
            fail(msg_if_not + " (edge count mismatch)");
        }
        for (const auto& e : EA) if (!EB.count(e)) {
                fail(msg_if_not + " (missing edge " + std::to_string(e.a) + "-" + std::to_string(e.b) + ")");
            }
    }

// ========== TESTS ==========

// Proposition 1 on ε-self-touching polygon (6x5)
    static void test_prop1_on_epsilon_self_touching() {
        std::cout << "[TEST] Observation 3 — Proposition 1 on ε-self-touching (6x5)\n";

        std::vector<Point> pts; std::vector<int> U, L, boundary; std::vector<Face> faces;
        make_bichain_polygon_with_boundary(6, 5, pts, U, L, boundary);

        // ε-touch at the right tip: bring U_last almost onto L_last
        apply_epsilon_touch_right(pts, /*upper_last=*/U.back(), /*lower_last=*/L.back(), /*eps=*/1e-12);

        // Start from strip A
        build_strip_A(U, L, faces);
        const int u = U.back();      // rightmost upper (near-touch to L_last)
        const int v = L.front();     // leftmost lower

        cgshop::Proposition1Inserter ins(pts, boundary, faces);
        // ins.enable_debug(std::cout);
        auto st = ins.insert(u, v);

        assert_true(st.inserted, "Prop.1 failed to insert target diagonal on ε-self-touching polygon.");
        // Face count invariant
        assert_true((int)faces.size() == (int)pts.size() - 2, "Face count changed.");

        std::cout << "  - phases=" << st.phases
                  << ", lemma4_calls=" << st.lemma4_calls
                  << ", lemma5_calls=" << st.lemma5_calls
                  << ", flips_total=" << st.flips_total << " ✅\n";
    }

// Theorem 2 A→B on ε-self-touching polygon (5x5)
    static void test_theorem2_on_epsilon_self_touching() {
        std::cout << "[TEST] Observation 3 — Theorem 2 A→B on ε-self-touching (5x5)\n";

        std::vector<Point> pts; std::vector<int> U, L, boundary;
        make_bichain_polygon_with_boundary(5, 5, pts, U, L, boundary);

        // ε-touch at the right tip
        apply_epsilon_touch_right(pts, U.back(), L.back(), 1e-12);

        std::vector<Face> F_start, F_target;
        build_strip_A(U, L, F_start);
        build_strip_B(U, L, F_target);

        cgshop::Theorem2Transformer T2(pts, boundary, F_start, F_target);
        // T2.enable_debug(std::cout);
        auto stats = T2.transform();

        assert_true(stats.finished, "Theorem 2 did not finish on ε-self-touching polygon.");
        assert_same_edge_sets(F_start, F_target, "Edge sets differ after transform (ε-self-touching).");
        assert_true((int)F_start.size() == (int)pts.size() - 2, "Face count changed.");

        std::cout << "  - nodes=" << stats.recursion_nodes
                  << ", prop1_phases=" << stats.prop1_phases
                  << ", lemma4_calls=" << stats.lemma4_calls
                  << ", lemma5_calls=" << stats.lemma5_calls
                  << ", flips_total=" << stats.flips_total << " ✅\n";
    }

    inline void run_all_observation3_tests() {
        std::cout << "===== Running Observation 3 (ε self-touching) Tests =====\n\n";
        test_prop1_on_epsilon_self_touching();
        std::cout << "\n";
        test_theorem2_on_epsilon_self_touching();
        std::cout << "\nAll Observation 3 tests passed ✅\n";
    }

} // namespace observation3_tests
