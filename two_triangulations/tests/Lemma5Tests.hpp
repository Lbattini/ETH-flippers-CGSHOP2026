// tests/Lemma5Tests.hpp
#pragma once
#include <bits/stdc++.h>
#include "../include/CutTriangulation.hpp"
#include "../include/Lemma5Transformer.hpp"

namespace lemma5_tests {

    using ::Point;
    using ::Edge;

// ---------- small utils ----------
    inline void assert_true(bool cond, const char* msg){
        if (!cond) throw std::runtime_error(msg);
    }

    inline bool triangulation_has_edge(const std::vector<std::array<int,3>>& F, int a,int b){
        if (a>b) std::swap(a,b);
        for (auto t:F){
            int u[3]={t[0],t[1],t[2]};
            for (int k=0;k<3;++k){
                int x=u[k], y=u[(k+1)%3]; if (x>y) std::swap(x,y);
                if (x==a && y==b) return true;
            }
        }
        return false;
    }

// Reconstruct outer boundary (CCW) from triangles:
// boundary edges are the ones that occur exactly once across all triangle edges.
    inline std::vector<int>
    reconstruct_boundary_ccw(const std::vector<Point>& P,
                             const std::vector<std::array<int,3>>& F)
    {
        // count undirected edges
        auto key = [](int a,int b){ if(a>b) std::swap(a,b); return ( (long long)a<<32 ) ^ (unsigned long long)b; };
        std::unordered_map<long long,int> cnt; cnt.reserve(F.size()*3*2);
        for (auto &t: F){
            int u[3]={t[0],t[1],t[2]};
            for (int k=0;k<3;++k){
                int a=u[k], b=u[(k+1)%3]; if (a>b) std::swap(a,b);
                ++cnt[key(a,b)];
            }
        }

        // collect boundary edges (count==1)
        std::unordered_map<int,std::vector<int>> adj; adj.reserve(cnt.size());
        for (auto &kv : cnt){
            if (kv.second != 1) continue;
            int a = int(kv.first>>32);
            int b = int(kv.first & 0xffffffff);
            adj[a].push_back(b);
            adj[b].push_back(a);
        }

        // walk the boundary cycle
        if (adj.empty()) return {};
        int start = adj.begin()->first;
        // prefer lexicographically smallest vertex as start (stable)
        for (auto &kv: adj) start = std::min(start, kv.first);

        std::vector<int> boundary;
        boundary.reserve(adj.size());
        int prev = -1, cur = start;

        while (true){
            boundary.push_back(cur);
            const auto &nbrs = adj[cur];
            if (nbrs.size() < 1) break; // shouldn't happen in a proper polygon
            int next;
            if (nbrs.size()==1) {
                next = (prev==-1)? nbrs[0] : -1;
            } else {
                // pick the neighbor that's not 'prev'
                next = (nbrs[0]==prev)? nbrs[1] : nbrs[0];
            }
            if (next==-1) break;
            prev = cur;
            cur = next;
            if (cur==start) { boundary.push_back(start); break; }
        }
        if (!boundary.empty() && boundary.front()==boundary.back()) boundary.pop_back();

        // ensure CCW orientation
        auto area2 = [&](const std::vector<int>& poly)->long double{
            long double s=0;
            for (size_t i=0;i<poly.size();++i){
                const Point& A = P[ poly[i] ];
                const Point& B = P[ poly[(i+1)%poly.size()] ];
                s += (long double)A.x*B.y - (long double)A.y*B.x;
            }
            return s;
        };
        if (boundary.size()>=3 && area2(boundary) < 0) std::reverse(boundary.begin(), boundary.end());
        return boundary;
    }

// ---------- instance definition ----------
    struct L5Instance {
        std::vector<Point> P;
        std::vector<Edge>  E;               // your hardcoded undirected edges (include boundary)
        int u{}, x{}, v{}, y{};             // indices (CCW order u, x, v, y)
    };

// ---------- YOU FILL THESE TWO FUNCTIONS ----------
    inline L5Instance make_variant_A(){
        L5Instance I;

        // >>>>>>>>>>>>>>>>>>>> YOUR HARD-CODED DATA (example below) <<<<<<<<<<<<<<<<<<<<<<
        // I.P = { {..,..}, ... };
        // I.E = { {u,v}, ... };  // include boundary + diagonals; planar; forms a full triangulation
        // I.u = ...; I.x = ...; I.v = ...; I.y = ...;

        // --- Example template (remove after you paste your own) ---
         I.P = {{-10,0},{-7,2},{-4,4},{0,10},{4,4},{7,2},{10,0},{7,-2},{4,-5},{0,-10},{-4,-5},{-7,-2}};
         I.E = { {0,1},{1,2},{2,3},{3,4},{4,5},{5,6},{6,7},{7,8},{8,9},{9,10},{10,11},{11,0},
                 {1,11},{2,11},{2,10},{3,10},{3,9},{4,9},{4,8},{5,8},{5,7} }; // contains (x,y) = (3,9)
         I.u=0; I.x=3; I.v=6; I.y=9;

        return I;
    }

    inline L5Instance make_variant_B(){
        L5Instance I;

        // >>>>>>>>>>>>>>>>>>>> YOUR HARD-CODED DATA (example below) <<<<<<<<<<<<<<<<<<<<<<
        // I.P = { ... };
        // I.E = { ... };
        // I.u = ...; I.x = ...; I.v = ...; I.y = ...;

        // --- Example template (remove after you paste your own) ---
        // I.P = {{-12,0},{-8,1},{-5,3},{-2,6},{1,9},{5,6},{8,3},{12,0},{8,-3},{4,-6},{0,-10},{-4,-6}};
        // I.E = { {0,1},{1,2},{2,3},{3,4},{4,5},{5,6},{6,7},{7,8},{8,9},{9,10},{10,11},{11,0},
        //         {4,10}, {3,10}, {5,10}, {6,10} }; // contains (x,y) = (4,10)
        // I.u=0; I.x=4; I.v=7; I.y=10;

        return I;
    }

// ---------- main test ----------
    inline void run_all_lemma5_tests(){
        std::cout << "===== Running Lemma 5 Tests (manual P+E) =====\n";

        auto run_one = [&](const char* name, const L5Instance& I){
            std::cout << "[TEST] " << name << "\n";

            // 1) Build faces from your edges (planar triangulation of the polygon)
            triangulation_cut::FacesOut F0 = triangulation_cut::extractFacesAsTriangles(I.P, I.E);
            auto faces = F0.triangles;
            const int n = (int)I.P.size();

            assert_true(!faces.empty(), "extractFacesAsTriangles returned 0 faces (check planarity/edges).");
            assert_true((int)faces.size() == n-2, "Not a full triangulation (faces != n-2).");

            // 2) Reconstruct boundary CCW from faces
            auto boundary_ccw = reconstruct_boundary_ccw(I.P, faces);
            assert_true(boundary_ccw.size() >= 3, "Failed to reconstruct boundary CCW.");

            // 3) Pre-checks: (x,y) must be present; (u,v) must not
            assert_true(triangulation_has_edge(faces, I.x, I.y), "Start triangulation must contain (x,y).");
            assert_true(!triangulation_has_edge(faces, I.u, I.v), "Start triangulation must NOT contain (u,v).");

            // 4) Convert points to cgshop::Point
            std::vector<cgshop::Point> pts; pts.reserve(n);
            for (auto &q : I.P) pts.push_back(cgshop::Point{q.x, q.y});

            // 5) Run your Lemma 5 transformer EXACTLY as you asked
            cgshop::Lemma5Transformer xf(pts, boundary_ccw, faces);
             xf.enable_debug(std::cout); // verbose if needed
            auto stats = xf.transform_from_xy_to_uv(I.x, I.u, I.y, I.v);

            // 6) Assertions after transform
            assert_true(stats.inserted, "Lemma 5: failed to insert (u,v).");
            assert_true(triangulation_has_edge(faces, I.u, I.v), "Edge (u,v) not found after transformation.");
            // assert_true((int)faces.size() == n-2, "Face count changed.");

            std::cout << "  ✅ inserted (u,v); faces="<<faces.size()<<" (n-2="<<(n-2)<<")";
            std::cout << ", rounds="<<stats.rounds;
            std::cout << ", flips="<<stats.flips;
            std::cout << "\n";
        };

        run_one("Lemma 5 — Variant A", make_variant_A());
        run_one("Lemma 5 — Variant B", make_variant_B());
    }

} // namespace lemma5_tests
