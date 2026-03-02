// Proposition1Tests.hpp
#pragma once
#include <vector>
#include <array>
#include <stdexcept>
#include <string>
#include <iostream>
#include <algorithm>

// Project headers
#include "../include/TriCommon.hpp"
#include "../include/Proposition1Inserter.hpp"

namespace proposition1_tests {

// ---------- tiny test helpers ----------
    [[noreturn]] inline void fail(const std::string& msg) { throw std::runtime_error(msg); }
    inline void assert_true(bool cond, const std::string& msg) { if (!cond) fail(msg); }

    using Point = cgshop::Point;
    using Face  = std::array<int,3>;

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

// Simple strip triangulation over the bi-chain polygon.
// For each i=1..min-1, triangulate quad (L_{i-1},L_i,U_i,U_{i-1}) with diagonal (L_i,U_{i-1}),
// then cap on the longer chain side.
    static void build_strip_triangulation(const std::vector<int>& U_lr,
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
        assert_true((int)faces.size() == n - 2, "Strip triangulation produced wrong face count.");
    }

// Utility: does (a,b) exist in faces? (unordered)
    static bool triangulation_has_edge(const std::vector<Face>& faces, int a, int b) {
        if (a > b) std::swap(a,b);
        for (const auto& T : faces) {
            int x=T[0], y=T[1], z=T[2];
            if ((std::min(x,y)==a && std::max(x,y)==b) ||
                (std::min(y,z)==a && std::max(y,z)==b) ||
                (std::min(z,x)==a && std::max(z,x)==b)) return true;
        }
        return false;
    }

// ------------- TESTS -------------

// Balanced case: 6x6 — target (u=U5, v=L0)
    static void test_prop1_balanced_6x6() {
        std::cout << "[TEST] Proposition 1 on balanced 6x6 (target: U_last — L0)\n";
        std::vector<Point> pts; std::vector<int> U, L, boundary; std::vector<Face> faces;
        make_bichain_polygon_with_boundary(6, 6, pts, U, L, boundary);
        build_strip_triangulation(U, L, faces);

        const int u = U.back();
        const int v = L.front();

        assert_true(!triangulation_has_edge(faces, u, v), "Target diagonal already present.");

        cgshop::Proposition1Inserter ins(pts, boundary, faces);
        // ins.enable_debug(std::cout); // uncomment for verbose phase logs
        auto st = ins.insert(u, v);

        assert_true(st.inserted, "Failed to insert (u,v) in 6x6.");
        assert_true(triangulation_has_edge(faces, u, v), "Edge (u,v) not found after insertion.");
        assert_true((int)faces.size() == (int)pts.size() - 2, "Face count changed.");

        std::cout << "  - phases=" << st.phases
                  << ", lemma4_calls=" << st.lemma4_calls
                  << ", lemma5_calls=" << st.lemma5_calls
                  << ", flips_total=" << st.flips_total << " ✅\n";
    }

// Unbalanced: 7x4 — target (u=U6, v=L0)
    static void test_prop1_unbalanced_7x4() {
        std::cout << "[TEST] Proposition 1 on unbalanced 7x4 (target: U_last — L0)\n";
        std::vector<Point> pts; std::vector<int> U, L, boundary; std::vector<Face> faces;
        make_bichain_polygon_with_boundary(7, 4, pts, U, L, boundary);
        build_strip_triangulation(U, L, faces);

        const int u = U.back();
        const int v = L.front();

        assert_true(!triangulation_has_edge(faces, u, v), "Target diagonal already present.");

        cgshop::Proposition1Inserter ins(pts, boundary, faces);
        auto st = ins.insert(u, v);

        assert_true(st.inserted, "Failed to insert (u,v) in 7x4.");
        assert_true(triangulation_has_edge(faces, u, v), "Edge (u,v) not found after insertion.");
        assert_true((int)faces.size() == (int)pts.size() - 2, "Face count changed.");

        std::cout << "  - phases=" << st.phases
                  << ", lemma4_calls=" << st.lemma4_calls
                  << ", lemma5_calls=" << st.lemma5_calls
                  << ", flips_total=" << st.flips_total << " ✅\n";
    }

// Near-coincident right upper (Observation 2-style): 5x5 with tiny nudge, target (U_last, L0)
    static void test_prop1_near_coincident_right_5x5() {
        std::cout << "[TEST] Proposition 1 on 5x5 with near-coincident right upper (target: U_last — L0)\n";
        std::vector<Point> pts; std::vector<int> U, L, boundary; std::vector<Face> faces;
        make_bichain_polygon_with_boundary(5, 5, pts, U, L, boundary);

        // Nudge rightmost upper slightly to the right to simulate coincidence geometry
        pts[ U.back() ].x += 1e-9;

        build_strip_triangulation(U, L, faces);

        const int u = U.back();
        const int v = L.front();

        cgshop::Proposition1Inserter ins(pts, boundary, faces);
        auto st = ins.insert(u, v);

        assert_true(st.inserted, "Failed to insert (u,v) in near-coincident 5x5.");
        assert_true(triangulation_has_edge(faces, u, v), "Edge (u,v) not found after insertion.");
        assert_true((int)faces.size() == (int)pts.size() - 2, "Face count changed.");

        std::cout << "  - phases=" << st.phases
                  << ", lemma4_calls=" << st.lemma4_calls
                  << ", lemma5_calls=" << st.lemma5_calls
                  << ", flips_total=" << st.flips_total << " ✅\n";
    }

// Light jitter test: 6x5 with small random-ish perturbations (deterministic values)
    static void test_prop1_with_jitter_6x5() {
        std::cout << "[TEST] Proposition 1 on 6x5 with small coordinate jitter\n";
        std::vector<Point> pts; std::vector<int> U, L, boundary; std::vector<Face> faces;
        make_bichain_polygon_with_boundary(6, 5, pts, U, L, boundary);

        // Add slight vertical jitter to upper chain (keeps x-monotone)
        for (size_t i = 0; i < U.size(); ++i) {
            double j = ( (int)i % 2 ? +1 : -1 ) * 1e-9;
            pts[ U[i] ].y += j;
        }

        build_strip_triangulation(U, L, faces);

        const int u = U.back();
        const int v = L.front();

        cgshop::Proposition1Inserter ins(pts, boundary, faces);
        auto st = ins.insert(u, v);

        assert_true(st.inserted, "Failed to insert (u,v) in jittered 6x5.");
        assert_true(triangulation_has_edge(faces, u, v), "Edge (u,v) not found after insertion.");
        assert_true((int)faces.size() == (int)pts.size() - 2, "Face count changed.");

        std::cout << "  - phases=" << st.phases
                  << ", lemma4_calls=" << st.lemma4_calls
                  << ", lemma5_calls=" << st.lemma5_calls
                  << ", flips_total=" << st.flips_total << " ✅\n";
    }

    inline void run_all_proposition1_tests() {
        std::cout << "===== Running Proposition 1 Inserter Tests =====\n\n";
        test_prop1_balanced_6x6();
        std::cout << "\n";
        test_prop1_unbalanced_7x4();
        std::cout << "\n";
        test_prop1_near_coincident_right_5x5();
        std::cout << "\n";
        test_prop1_with_jitter_6x5();
        std::cout << "\nAll Proposition 1 tests passed ✅\n";
    }

} // namespace proposition1_tests
