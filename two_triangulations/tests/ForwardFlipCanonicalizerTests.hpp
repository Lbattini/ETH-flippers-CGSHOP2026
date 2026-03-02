// ForwardFlipTests.hpp
#pragma once
#include <vector>
#include <array>
#include <stdexcept>
#include <string>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include "../include/ForwardFlipCanonicalizer.hpp"

namespace forward_flip_tests {

// ----------------- tiny test helpers -----------------
    [[noreturn]] inline void fail(const std::string& msg) {
        throw std::runtime_error(msg);
    }
    inline void assert_true(bool cond, const std::string& msg) {
        if (!cond) fail(msg);
    }

    using Point = cgshop::Point;
    using Face  = std::array<int,3>;

// Build a bi-chain, x-monotone polygon: lower y=0 (L0..L_{p-1}), upper y=1 (U0..U_{m-1})
// CCW boundary order: L0, L1, ..., L_{p-1}, U_{m-1}, ..., U0
    static void make_bichain_polygon(int m_upper, int p_lower,
                                     std::vector<Point>& pts,
                                     std::vector<int>& upperChain,
                                     std::vector<int>& lowerChain)
    {
        assert_true(m_upper >= 1 && p_lower >= 1, "Need at least 1 vertex per chain.");
        assert_true(m_upper + p_lower >= 3, "Polygon needs at least 3 vertices.");
        pts.clear(); upperChain.clear(); lowerChain.clear();
        pts.reserve(m_upper + p_lower);

        // Lower chain (left -> right), y = 0
        for (int i = 0; i < p_lower; ++i) {
            pts.push_back(Point{double(i), 0.0});
            lowerChain.push_back((int)pts.size()-1);
        }
        // Upper chain (left -> right), y = 1
        std::vector<int> tempU;
        tempU.reserve(m_upper);
        for (int j = 0; j < m_upper; ++j) {
            pts.push_back(Point{double(j), 1.0});
            tempU.push_back((int)pts.size()-1);
        }
        upperChain = tempU; // left->right
    }

// Create a simple strip triangulation
    static void build_strip_triangulation(const std::vector<int>& U, const std::vector<int>& L,
                                          std::vector<Face>& faces)
    {
        faces.clear();
        const int m = (int)U.size();
        const int p = (int)L.size();
        const int k = std::min(m, p);

        for (int i = 1; i < k; ++i) {
            int L_im1 = L[i-1], L_i = L[i];
            int U_im1 = U[i-1], U_i = U[i];
            faces.push_back({L_im1, L_i, U_im1});
            faces.push_back({L_i, U_i, U_im1});
        }
        if (m > p) {
            int base = L[p-1];
            for (int j = p; j < m; ++j)
                faces.push_back({base, U[j-1], U[j]});
        }
        if (p > m) {
            int base = U[m-1];
            for (int i = m; i < p; ++i)
                faces.push_back({L[i-1], L[i], base});
        }
        const int n = m + p;
        assert_true((int)faces.size() == n - 2, "Strip triangulation produced wrong face count.");
    }

// Extract “top” upper vertex for each lower base (L_{i-1}, L_i)
    static std::vector<int> compute_lower_base_tops(const std::vector<Face>& faces,
                                                    const std::vector<int>& U, const std::vector<int>& L)
    {
        std::vector<int> rankU, rankL;
        int nMax = 0;
        for (int v : U) nMax = std::max(nMax, v);
        for (int v : L) nMax = std::max(nMax, v);
        rankU.assign(nMax+1, -1);
        rankL.assign(nMax+1, -1);
        for (int i=0;i<(int)U.size();++i) rankU[U[i]] = i;
        for (int i=0;i<(int)L.size();++i) rankL[L[i]] = i;

        struct EdgeKey { int a,b; bool operator==(const EdgeKey& o) const { return a==o.a && b==o.b; } };
        struct EdgeHash { size_t operator()(const EdgeKey& e) const { return (size_t)e.a*1315423911u ^ (size_t)e.b; } };
        auto norm = [](int a,int b){ return (a<b)? EdgeKey{a,b}:EdgeKey{b,a}; };

        std::unordered_map<EdgeKey, std::vector<int>, EdgeHash> base2faces;
        base2faces.reserve(L.size()*3);

        for (int f=0; f<(int)faces.size(); ++f) {
            auto tri = faces[f];
            int a=tri[0], b=tri[1], c=tri[2];
            if (rankL[a]>=0 && rankL[b]>=0) base2faces[norm(a,b)].push_back(f);
            if (rankL[b]>=0 && rankL[c]>=0) base2faces[norm(b,c)].push_back(f);
            if (rankL[c]>=0 && rankL[a]>=0) base2faces[norm(c,a)].push_back(f);
        }

        std::vector<int> tops; tops.resize(L.size()>0? (L.size()-1):0, -1);
        for (int i=1; i<(int)L.size(); ++i) {
            auto ek = norm(L[i-1], L[i]);
            auto it = base2faces.find(ek);
            if (it == base2faces.end()) fail("Missing lower base in faces.");
            int found = -1;
            for (int f : it->second) {
                auto tri = faces[f];
                for (int t=0;t<3;++t) {
                    int v = tri[t];
                    if (v != L[i-1] && v != L[i] && v < (int)rankU.size() && rankU[v] >= 0) {
                        found = v; break;
                    }
                }
                if (found >= 0) break;
            }
            if (found < 0) fail("Could not find upper top for a lower base.");
            tops[i-1] = found;
        }
        return tops;
    }

// ===== TESTS =====

    static void test_forward_flip_on_balanced_grid() {
        std::cout << "[TEST] Forward flip canonicalization on balanced 4x4 bichain\n";
        std::vector<Point> pts; std::vector<int> U, L;
        make_bichain_polygon(4,4, pts, U, L);


        std::vector<Face> faces;
        build_strip_triangulation(U, L, faces);

        cgshop::ForwardFlipCanonicalizer canon(pts, U, L, faces);
        auto stats = canon.run();

        assert_true((int)faces.size() == (int)(pts.size() - 2), "Face count changed.");
        auto tops = compute_lower_base_tops(faces, U, L);
        for (int v : tops)
            if (v != U.back()) fail("Not all lower-base tops reached rightmost upper vertex.");

        std::cout << "  - rounds=" << stats.rounds << ", flips=" << stats.flips_total << " ✅\n";
        std::cout << "Flips rounds:\n";
        for (auto& i : stats.flips){
            for (auto& j : i)
                std::cout << j.e << '\n';
            std::cout << "----------\n";
        }
    }

    static void test_forward_flip_on_unbalanced_5x3() {
        std::cout << "[TEST] Forward flip canonicalization on unbalanced 5x3 bichain\n";
        std::vector<Point> pts; std::vector<int> U, L;
        make_bichain_polygon(5,3, pts, U, L);

        std::vector<Face> faces;
        build_strip_triangulation(U, L, faces);

        cgshop::ForwardFlipCanonicalizer canon(pts, U, L, faces);
        auto stats = canon.run();

        auto tops = compute_lower_base_tops(faces, U, L);
        for (int v : tops)
            if (v != U.back()) fail("Not all lower-base tops reached rightmost upper vertex.");

        std::cout << "  - rounds=" << stats.rounds << ", flips=" << stats.flips_total << " ✅\n";
    }

    static void test_idempotence_against_extra_round() {
        std::cout << "[TEST] Idempotence of canonicalization\n";
        std::vector<Point> pts; std::vector<int> U, L;
        make_bichain_polygon(4,4, pts, U, L);
        std::vector<Face> faces;
        build_strip_triangulation(U, L, faces);

        cgshop::ForwardFlipCanonicalizer canon1(pts, U, L, faces);
        canon1.run();
        cgshop::ForwardFlipCanonicalizer canon2(pts, U, L, faces);
        auto stats2 = canon2.run();

        assert_true(stats2.flips_total == 0, "Canonical triangulation should be fixed point.");
        std::cout << "  - flips on second run: " << stats2.flips_total << " ✅\n";
    }

    static void test_forward_flip_three_convex_case() {
        std::cout << "[TEST] Forward flip — three-convex (U size = 1)\n";
        std::vector<Point> pts; std::vector<int> U, L;

        // U ima 1 vrh (apeks), L ima 5 vrhova — „levak”.
        // Sve baze su na donjem lancu; „top” je uvek jedini gornji vrh => već kanonski.
        make_bichain_polygon(/*m_upper=*/1, /*p_lower=*/5, pts, U, L);

        // Napravi fan triangulaciju od apeksa preko donjeg lanca:
        std::vector<Face> faces;
        faces.clear();
        for (int i = 1; i < (int)L.size(); ++i) {
            // trougao (apeks, L_{i-1}, L_i)
            faces.push_back({U[0], L[i-1], L[i]});
        }
        // sanity
        assert_true((int)faces.size() == (int)pts.size() - 2, "Face count mismatch.");

        cgshop::ForwardFlipCanonicalizer canon(pts, U, L, faces);
        auto stats = canon.run();

        // Ne bi trebalo da bude flipova; već je kanon.
        assert_true(stats.flips_total == 0, "Three-convex case should be already canonical.");
        std::cout << "  - flips=" << stats.flips_total << " (expected 0) ✅\n";
    }


    inline void run_all_forward_flip_tests() {
        std::cout << "===== Running Forward Flip Tests =====\n\n";
        test_forward_flip_on_balanced_grid();
        std::cout << "\n";
        test_forward_flip_on_unbalanced_5x3();
        std::cout << "\n";
        test_idempotence_against_extra_round();
        std::cout << "\n";
        test_forward_flip_three_convex_case();
        std::cout << "\nAll forward flip tests passed ✅\n";
    }

} // namespace forward_flip_tests
