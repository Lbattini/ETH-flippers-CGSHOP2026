// tests/Lemma4ManualHarness.hpp
#pragma once
#include <bits/stdc++.h>
#include "../include/CutTriangulation.hpp"
#include "../include/DiagonalIntroducerLemma4.hpp"

namespace lemma4_manual_harness {

    struct FlipEdge { int a{}, b{}; };

    struct L4RunResult {
        bool ok{false};
        int flips{0};
        int rounds{0};
        int rejects{0};
        std::vector<std::array<int,3>> faces_out;
        std::string error;
        // optional: filled if Lemma4 exposes export_flip_log(out)
        std::vector<std::vector<FlipEdge>> flip_seq;
    };

    inline bool faces_have_edge(const std::vector<std::array<int,3>>& F, int a, int b){
        if (a>b) std::swap(a,b);
        for (auto t : F){
            int u[3] = {t[0], t[1], t[2]};
            for (int k=0;k<3;++k){
                int x=u[k], y=u[(k+1)%3];
                if (x>y) std::swap(x,y);
                if (x==a && y==b) return true;
            }
        }
        return false;
    }

// SFINAE helper: if DiagonalIntroducerLemma4 has export_flip_log(vector<vector<FlipEdge>>&),
// call it; otherwise no-op.
    template<class L4>
    static auto try_export_flip_log(L4& l4, std::vector<std::vector<FlipEdge>>& out, int)
    -> decltype(l4.export_flip_log(out), void()) { l4.export_flip_log(out); }
    template<class L4>
    static void try_export_flip_log(L4&, std::vector<std::vector<FlipEdge>>&, ...) { /* no-op */ }

// Main entry
    inline bool run_lemma4_on(const std::vector<Point>& P,
                              const std::vector<Edge>&  E,
                              int v1, int vq,
                              L4RunResult& out,
                              std::ostream* dbg = nullptr,
                              bool require_full_triangulation = true)
    {
        out = L4RunResult{}; // reset

        // 1) faces from edges
        triangulation_cut::FacesOut Fcur = triangulation_cut::extractFacesAsTriangles(P, E);
        auto faces = Fcur.triangles;

        if (faces.empty()) {
            out.error = "extractFacesAsTriangles produced 0 faces (edge set not a triangulation?).";
            return false;
        }
        if (require_full_triangulation) {
            const int n = (int)P.size();
            if ((int)faces.size() != n-2) {
                out.error = "Start triangulation is not full (faces != n-2).";
                return false;
            }
        }

        // 2) convert points
        std::vector<cgshop::Point> Pc; Pc.reserve(P.size());
        for (auto &q : P) Pc.push_back(cgshop::Point{q.x, q.y});

        // 3) run Lemma 4
        cgshop::DiagonalIntroducerLemma4 L4(Pc, faces);
        if (dbg) L4.enable_debug(*dbg);
        auto st = L4.introduce(v1, vq);

        // 4) optional flip log
        try_export_flip_log(L4, out.flip_seq, 0);

        // 5) fill result
        out.ok        = faces_have_edge(faces, v1, vq);
        out.flips     = st.flips_total;
        out.rounds    = st.rounds;
        out.faces_out = std::move(faces);
        if (!out.ok) out.error = "Target diagonal not present after introduce().";
        return out.ok;
    }

    inline void print_faces_as_edges(const std::vector<std::array<int,3>>& F, std::ostream& os){
        std::unordered_set<long long> S;
        auto key=[&](int a,int b){ if(a>b) std::swap(a,b); return ((long long)a<<32) ^ (unsigned long long)b; };
        for (auto t: F){
            int u[3]={t[0],t[1],t[2]};
            for (int k=0;k<3;++k){
                int a=u[k], b=u[(k+1)%3]; S.insert(key(a,b));
            }
        }
        std::vector<Edge> Ed; Ed.reserve(S.size());
        for (auto k: S) Ed.push_back({int(k>>32), int(k & 0xffffffff)});
        std::sort(Ed.begin(), Ed.end(), [](const Edge& A,const Edge& B){ return (A.u==B.u)? A.v<B.v : A.u<B.u; });
        os << "edges from faces ("<<Ed.size()<<"):\n";
        for (auto e: Ed) os << "  ("<<e.u<<","<<e.v<<")\n";
    }

} // namespace lemma4_manual_harness
