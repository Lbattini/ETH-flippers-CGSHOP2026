// tests/DiagonalIntroducerLemma4ProperTests.hpp
#pragma once
#include <bits/stdc++.h>
#include "../include/DiagonalIntroducerLemma4.hpp"
#include "../include/CutTriangulation.hpp"

namespace lemma4_proper_tests {

    using ::Point;
    using ::Edge;

// -------------------- geometry helpers --------------------
    inline long double cross(const Point& A, const Point& B, const Point& C){
        return (long double)(B.x-A.x)*(C.y-A.y) - (long double)(B.y-A.y)*(C.x-A.x);
    }
    inline int orient_ccw(const std::vector<Point>& P, int a,int b,int c){
        long double s = cross(P[a],P[b],P[c]);
        return (s>0) - (s<0);
    }
    inline bool seg_strict_intersect(Point A,Point B,Point C,Point D){
        auto d = [](Point A,Point B,Point C){
            long double s = (long double)(B.x-A.x)*(C.y-A.y) - (long double)(B.y-A.y)*(C.x-A.x);
            return (s>0) - (s<0);
        };
        int o1=d(A,B,C), o2=d(A,B,D), o3=d(C,D,A), o4=d(C,D,B);
        return (o1*o2<0 && o3*o4<0);
    }
    inline bool visible_segment_inside(const std::vector<Point>& P, int a, int b){
        int n=(int)P.size();
        for (int i=0;i<n;++i){
            int j=(i+1)%n;
            if (i==a || i==b || j==a || j==b) continue;
            if (seg_strict_intersect(P[a],P[b], P[i],P[j])) return false;
        }
        return true;
    }

// point in triangle (strict barycentric)
    inline bool point_in_tri_strict(const Point& A,const Point& B,const Point& C,const Point& X){
        long double c1 = cross(A,B,X);
        long double c2 = cross(B,C,X);
        long double c3 = cross(C,A,X);
        bool pos = (c1>0 && c2>0 && c3>0);
        bool neg = (c1<0 && c2<0 && c3<0);
        return pos || neg;
    }

// Ear clipping that forbids creating the diagonal (forbid_u, forbid_v)
    inline std::vector<std::array<int,3>>
    ear_clip_triangulate_no_edge(const std::vector<Point>& P, int forbid_u, int forbid_v)
    {
        const int n=(int)P.size();
        std::vector<int> V(n);
        std::iota(V.begin(), V.end(), 0);

        auto next_idx = [&](int i){ return (i+1)%V.size(); };
        auto prev_idx = [&](int i){ return (i-1+V.size())%V.size(); };

        std::vector<std::array<int,3>> T; T.reserve(n-2);

        auto is_forbidden = [&](int a,int b){
            if (a>b) std::swap(a,b);
            int u=std::min(forbid_u,forbid_v), v=std::max(forbid_u,forbid_v);
            return (a==u && b==v);
        };

        // polygon is CCW; accept an ear at i only if triangle prev-i-next is convex,
        // no other vertex inside, and clipping doesn't create the forbidden diagonal.
        while ((int)V.size() > 3){
            bool clipped=false;

            for (int i=0;i<(int)V.size();++i){
                int i0 = prev_idx(i), i1 = i, i2 = next_idx(i);
                int a=V[i0], b=V[i1], c=V[i2];

                // convex?
                if (orient_ccw(P, a,b,c) <= 0) continue;

                // would ear create forbidden diagonal? (the diagonal is between a and c)
                if (is_forbidden(a,c)) continue;

                // any vertex strictly inside tri?
                bool any_inside=false;
                for (int k=0;k<(int)V.size();++k){
                    if (k==i0 || k==i1 || k==i2) continue;
                    if (point_in_tri_strict(P[a],P[b],P[c], P[V[k]])) { any_inside=true; break; }
                }
                if (any_inside) continue;

                // it's an ear: add triangle and remove b

                std:: cout << a << ' ' << b << ' ' << c << '\n';
                T.push_back({a,b,c});
                V.erase(V.begin()+i1);
                clipped=true;
                break;
            }

            if (!clipped) {
                // degenerate or numeric trouble; fallback: allow any ear (still avoid forbidden if possible)
                for (int i=0;i<(int)V.size();++i){
                    int i0 = prev_idx(i), i1 = i, i2 = next_idx(i);
                    int a=V[i0], b=V[i1], c=V[i2];
                    if (is_forbidden(a,c)) continue;
                    std::cout << a << ' ' << b << ' ' << c << '\n';
                    T.push_back({a,b,c});
                    V.erase(V.begin()+i1);
                    clipped=true;
                    break;
                }
                if (!clipped) break; // give up
            }
        }
        if (V.size()==3) {
            T.push_back({V[0], V[1], V[2]});
            std::cout << V[0] << ' ' << V[1] << ' ' << V[2] << '\n';
        }
        return T;
    }

// -------------------- test instance --------------------
    struct L4Instance {
        std::vector<Point> P;   // CCW polygon
        std::vector<std::array<int,3>> F0; // initial triangulation as triangles (NO (v1,vq))
        int v1{}, vp{}, vq{}, vqp1{}, vn{};
    };

    inline L4Instance make_instance_strict_no_target() {
        std::vector<Point> P;

        // v_q on baseline
        int vq = 0;                       P.push_back({0.0, 0.0});  // 0

        // upper chain with a tall peak v_p
        P.push_back({1.0, 0.6});          // 1
        P.push_back({2.0, 1.3});          // 2
        int vp = (int)P.size();           P.push_back({3.0, 2.3});  // 3 = v_p
        P.push_back({4.0, 1.1});          // 4
        P.push_back({5.0, 0.4});          // 5

        // v_1 on baseline, right
        int v1 = (int)P.size();           P.push_back({6.0, 0.0});  // 6 = v_1

        // bottom chain with sharp dip after v1 and convex up near v_q
        int vqp1 = (int)P.size();         P.push_back({5.2, -1.0}); // 7 = v_{q+1}
        P.push_back({4.4, -1.5});         // 8
        P.push_back({3.6, -1.6});         // 9
        P.push_back({2.6, -1.2});         //10
        P.push_back({1.6, -1.35});        //11
        int vn = (int)P.size();           P.push_back({0.4, -0.6}); //12 = v_n

        // Ear-clip triangulate while forbidding (v1,vq)
        auto F0 = ear_clip_triangulate_no_edge(P, v1, vq);

        return {P, F0, v1, vp, vq, vqp1, vn};
    }

// -------------------- precheck & test --------------------
    inline bool check_preconditions(const L4Instance& I, std::ostream& os){
        const auto &P = I.P;
        int n=(int)P.size();
        std::vector<int> convexIdx;
        for (int i=0;i<n;++i){
            int a=(i-1+n)%n, b=i, c=(i+1)%n;
            if (orient_ccw(P,a,b,c) > 0) convexIdx.push_back(i);
        }
        bool vis = visible_segment_inside(P, I.v1, I.vq);
        os << "Lemma4 precheck: convex="<<convexIdx.size()
           <<", visible(v1,vq)="<<(vis?"yes":"no")<<"\n  convex at: ";
        for (int id : convexIdx) os << id << " ";
        os << "\n";

        return (convexIdx.size() >= 5) && vis;
    }

    inline bool faces_have_edge(const std::vector<std::array<int,3>>& F, int a,int b){
        if (a>b) std::swap(a,b);
        for (auto t: F){
            int u[3]={t[0],t[1],t[2]};
            for (int k=0;k<3;++k){
                int x=u[k], y=u[(k+1)%3]; if (x>y) std::swap(x,y);
                if (x==a && y==b) return true;
            }
        }
        return false;
    }

    inline void run_proper_lemma4_test() {
        std::cout << "===== Lemma 4 (proper) — must insert v1–vq (no-op forbidden) =====\n";
        auto inst = make_instance_strict_no_target();

        if (!check_preconditions(inst, std::cout))
            throw std::runtime_error("Constructed instance does not meet Lemma 4 preconditions.");

        const int n = (int)inst.P.size();
        if ((int)inst.F0.size() != n-2)
            throw std::runtime_error("Start triangulation is not full (faces != n-2).");

        if (faces_have_edge(inst.F0, inst.v1, inst.vq))
            throw std::runtime_error("Start triangulation unexpectedly has (v1,vq).");

        // Convert to cgshop types
        std::vector<cgshop::Point> Pcg; Pcg.reserve(n);
        for (auto &q: inst.P) Pcg.push_back(cgshop::Point{q.x, q.y});
        std::vector<std::array<int,3>> faces = inst.F0;

        // Run Lemma 4
        cgshop::DiagonalIntroducerLemma4 L4(Pcg, faces);
        L4.enable_debug(std::cout); // comment if too verbose
        auto st = L4.introduce(inst.v1, inst.vq);

        std::cout << "Lemma4 stats: flips=" << st.flips_total
                  << ", rounds=" << st.rounds << "\n";

        if (!faces_have_edge(faces, inst.v1, inst.vq))
            throw std::runtime_error("Resulting triangulation does not contain (v1,vq).");

        if (st.flips_total == 0)
            throw std::runtime_error("Lemma 4 reported 0 flips, but target was absent initially.");

        std::cout << "✅ Lemma 4 succeeded with flips and inserted ("<<inst.v1<<","<<inst.vq<<").\n";
    }

} // namespace lemma4_proper_tests
