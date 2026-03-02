// Theorem2Tests.hpp
#pragma once
#include <vector>
#include <array>
#include <stdexcept>
#include <string>
#include <iostream>
#include <unordered_set>
#include <algorithm>

#include "../include/TriCommon.hpp"
#include "../include/Theorem2Transformer.hpp"

namespace theorem2_tests {

// ---------- tiny test helpers ----------
    [[noreturn]] inline void fail(const std::string& msg) { throw std::runtime_error(msg); }
    inline void assert_true(bool cond, const std::string& msg) { if (!cond) fail(msg); }

    using Point = cgshop::Point;
    using Face  = std::array<int,3>;
    using EdgeKey = cgshop::EdgeKey;
    using EdgeKeyHash = cgshop::EdgeKeyHash;

// Build a bi-chain, x-monotone polygon with boundary (CCW):
//   lower chain y=0: L0..L_{p-1}  (left->right)
//   upper chain y=1: U0..U_{m-1}  (left->right)
// CCW boundary we return: L0..L_{p-1}, U_{m-1}..U0
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

        // Lower y=0
        for (int i = 0; i < p_lower; ++i) {
            pts.push_back(Point{double(i), 0.0});
            lower_lr.push_back((int)pts.size()-1);
        }
        // Upper y=1
        for (int j = 0; j < m_upper; ++j) {
            pts.push_back(Point{double(j), 1.0});
            upper_lr.push_back((int)pts.size()-1);
        }
        // Boundary CCW: lower forward, upper reversed
        for (int v : lower_lr) boundary_ccw.push_back(v);
        for (int k = (int)upper_lr.size()-1; k >= 0; --k) boundary_ccw.push_back(upper_lr[k]);
    }

// Strip A: for each i=1..min-1, triangulate quad (L_{i-1},L_i,U_i,U_{i-1}) with diagonal (L_i,U_{i-1})
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

// Strip B (the “other diagonal”): choose (L_{i-1},U_i) inside each quad
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

// Build normalized edge set from faces
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

// Compare full edge sets (boundary + diagonals)
    static void assert_same_edge_sets(const std::vector<Face>& A,
                                      const std::vector<Face>& B,
                                      const std::string& msg_if_not)
    {
        auto EA = edges_from_faces(A);
        auto EB = edges_from_faces(B);
        if (EA.size() != EB.size()) {
            fail(msg_if_not + " (different edge counts: " + std::to_string(EA.size()) + " vs " + std::to_string(EB.size()) + ")");
        }
        for (const auto& e : EA) {
            if (!EB.count(e)) {
                fail(msg_if_not + " (missing edge " + std::to_string(e.a) + "-" + std::to_string(e.b) + ")");
            }
        }
    }

// ------------- TESTS -------------

    static void test_theorem2_balanced_6x6() {
        std::cout << "[TEST] Theorem 2 — balanced 6x6: strip A → strip B\n";

        std::vector<Point> pts; std::vector<int> U, L, boundary;
        make_bichain_polygon_with_boundary(6, 6, pts, U, L, boundary);

        std::vector<Face> F_start, F_target;
        build_strip_A(U, L, F_start);
        build_strip_B(U, L, F_target);

        // Run transformer
        cgshop::Theorem2Transformer T2(pts, boundary, F_start, F_target);
        // T2.enable_debug(std::cout);
        auto stats = T2.transform();

        assert_true(stats.finished, "Theorem 2 transformer did not finish (6x6).");
        // Strong check: identical edge sets
        assert_same_edge_sets(F_start, F_target, "Edge sets differ after transform (6x6).");

        const int n = (int)pts.size();
        assert_true((int)F_start.size() == n - 2, "Face count changed (6x6).");

        std::cout << "  - nodes=" << stats.recursion_nodes
                  << ", prop1_phases=" << stats.prop1_phases
                  << ", lemma4_calls=" << stats.lemma4_calls
                  << ", lemma5_calls=" << stats.lemma5_calls
                  << ", flips_total=" << stats.flips_total << " ✅\n";
    }

    static void test_theorem2_unbalanced_7x4() {
        std::cout << "[TEST] Theorem 2 — unbalanced 7x4: strip A → strip B\n";

        std::vector<Point> pts; std::vector<int> U, L, boundary;
        make_bichain_polygon_with_boundary(7, 4, pts, U, L, boundary);

        std::vector<Face> F_start, F_target;
        build_strip_A(U, L, F_start);
        build_strip_B(U, L, F_target);

        cgshop::Theorem2Transformer T2(pts, boundary, F_start, F_target);
        auto stats = T2.transform();

        assert_true(stats.finished, "Theorem 2 transformer did not finish (7x4).");
        assert_same_edge_sets(F_start, F_target, "Edge sets differ after transform (7x4).");

        const int n = (int)pts.size();
        assert_true((int)F_start.size() == n - 2, "Face count changed (7x4).");

        std::cout << "  - nodes=" << stats.recursion_nodes
                  << ", prop1_phases=" << stats.prop1_phases
                  << ", lemma4_calls=" << stats.lemma4_calls
                  << ", lemma5_calls=" << stats.lemma5_calls
                  << ", flips_total=" << stats.flips_total << " ✅\n";
    }

    static void test_theorem2_jittered_6x5() {
        std::cout << "[TEST] Theorem 2 — 6x5 with slight jitter: strip A → strip B\n";

        std::vector<Point> pts; std::vector<int> U, L, boundary;
        make_bichain_polygon_with_boundary(6, 5, pts, U, L, boundary);

        // slight vertical jitter on upper chain (keeps monotone, avoids strict-collinearity)
        for (size_t i = 0; i < U.size(); ++i) {
            pts[ U[i] ].y += ((i%2) ? +1e-9 : -1e-9);
        }

        std::vector<Face> F_start, F_target;
        build_strip_A(U, L, F_start);
        build_strip_B(U, L, F_target);

        cgshop::Theorem2Transformer T2(pts, boundary, F_start, F_target);
        auto stats = T2.transform();

        assert_true(stats.finished, "Theorem 2 transformer did not finish (6x5 jittered).");
        assert_same_edge_sets(F_start, F_target, "Edge sets differ after transform (6x5 jittered).");

        const int n = (int)pts.size();
        assert_true((int)F_start.size() == n - 2, "Face count changed (6x5 jittered).");

        std::cout << "  - nodes=" << stats.recursion_nodes
                  << ", prop1_phases=" << stats.prop1_phases
                  << ", lemma4_calls=" << stats.lemma4_calls
                  << ", lemma5_calls=" << stats.lemma5_calls
                  << ", flips_total=" << stats.flips_total << " ✅\n";
    }

    inline void run_all_theorem2_tests() {
        std::cout << "===== Running Theorem 2 Transformer Tests =====\n\n";
        test_theorem2_balanced_6x6();
        std::cout << "\n";
        test_theorem2_unbalanced_7x4();
        std::cout << "\n";
        test_theorem2_jittered_6x5();
        std::cout << "\nAll Theorem 2 tests passed ✅\n";
    }

} // namespace theorem2_tests
