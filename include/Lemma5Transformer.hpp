// Lemma5Transformer.hpp
#pragma once
#include <vector>
#include <array>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <ostream>
#include <cassert>
#include "./TriCommon.hpp"
#include "./DiagonalIntroducerLemma4.hpp"

namespace cgshop {

    class Lemma5Transformer {
    public:
        using Face = std::array<int,3>;

        struct Stats {
            bool inserted = false;
            int  flips    = 0;
            int  rounds   = 0;
            int  steps    = 0;  // number of Lemma 4 *iterations* (including zero-flip skips)
        };

        Lemma5Transformer(const std::vector<Point>& pts,
                          const std::vector<int>& boundary_ccw,
                          std::vector<Face>& faces)
                : P_(pts), B_(boundary_ccw), F_(faces)
        {
            pos_.reserve(B_.size());
            build_adjacency_from_faces(F_, pts.size(), adj_);
            for (int i = 0; i < (int)B_.size(); ++i) pos_[B_[i]] = i;
        }

        void enable_debug(std::ostream& os) { dbg_ = &os; debug_ = true; }

        Stats transform_from_xy_to_uv(int x, int u, int y, int v) {
            Stats out;


            if (has_edge(u, v)) { out.inserted = true; return out; }

            if (!has_edge(x, y)) {out.inserted = false; return out; }

            std::cout << "Valid conditions for Lemma 5\n";
            int n = P_.size();
            bool swapped = 0;
            if (degree_of_vertex(y) > 3){
                //good
            } else {
                swapped = 1;
                std::swap(x, y);
            }

            int ix = pos_[x];  // position of x in boundary
            int iu = pos_[u];  // position of u in boundary
            ix++;
            ix %= n;

            int x1 = -1;

            for (int k = ix; k != iu; k = (k + 1) % n) {
                int v = B_[k];
                if (has_edge(v, y)){
                    x1 = v;
                    break;
                }
            }


            if (x1 == -1){
                std::swap(u, v);
                swapped = 1-swapped;

                ix = pos_[x];  // position of x in boundary
                iu = pos_[u];  // position of u in boundary

                ix = ix + n-1;
                ix %= n;

                for (int k = ix; k != iu; k = (k - 1 + n) % n) {
                    int v = B_[k];
                    if (has_edge(v, y)){
                        x1 = v;
                        break;
                    }
                }
            }


            if (x1 == -1){
                out.inserted = 0;
                std::cout << "FAILED ON FINDING X1\n";
                return out;
            }

            std::cout << "x y u v x1\n";
            std::cout << x << ' ' << y << ' ' << u << ' ' << v << ' ' << x1 << '\n';

            int w = find_tangent_to_chain(x1, x, v);

            std::cout << "w\n";
            std::cout << w << '\n';

            int add = 1;
            if (swapped)
                add = n-1;

            int w1 = B_[(pos_[w] + add ) %n ];
            std::cout << "w1=" << w1 << '\n';

            add = n-1;
            if (swapped)
                add=1;

            int st = (pos_[w] + add) %n, y1 = -1;
            for (st; st != pos_[w1]; st = (st +add) % n){
                if (has_edge(B_[st], w) && has_edge(B_[st], w1)){
                    std::cout << "y1 " << B_[st] << '\n';
                    y1 = B_[st];
                    break;
                }
            }

            if (y1 == -1){
                std::cout << "Picking y1 failed\n";
                return out;
            }

            add = 1;
            if (swapped)
                add = n-1;

            std::vector<int> poly;
            for (int i = pos_[x]; i != pos_[x1]; i=(i+add)%n){
                poly.push_back(B_[i]);
            }
            poly.push_back(B_[pos_[x1]]);
            for (int i = pos_[y]; i != pos_[y1]; i=(i+add)%n){
                poly.push_back(B_[i]);
            }
            poly.push_back(B_[pos_[y1]]);
            for (int i = pos_[w]; i != pos_[x]; i=(i+add)%n){
                poly.push_back(B_[i]);
            }

            std::cout << "POLY\n";
            for (int i : poly)
                std::cout << i << ' ';
            std::cout << '\n';

            callLemma4(poly, x1, w, y1);

            //!!!!!! Maybetangent should work cw and ccw?
            int z = find_tangent_to_chain(y1, y, u);
            std::cout << "z=" << z <<'\n';

            //!!!!!! consider above
            int z1 = B_[(pos_[z] + 1) % n];
            add = 1;
            if (swapped)
                add=n-1;

            st = (pos_[z] + add) %n;
            int t = -1;
            for (st; st != pos_[z1]; st = (st +add) % n){
                if (has_edge(B_[st], z) && has_edge(B_[st], z1)){
                    std::cout << "t " << B_[st] << '\n';
                    t = B_[st];
                    break;
                }
            }

            if (t == -1){
                std::cout << "Picking t failed\n";
                return out;
            }

            poly.clear();

            add = 1;
            if (swapped)
                add = n-1;



            int g = -1;
            add = 1;
            if (swapped)
                add = n - 1;

            int iy = pos_[y];
            int iv = pos_[v];
            for (int i = iy; i != iv; i = (i + add) % n){
                if (i == iy) continue;
                std::cout << B_[i] << ' ';
                if (has_edge(B_[i], B_[x])){
                    g = i;
                    break;
                }
            }

            if (g == -1){
                std::cout << "FAILED ON PICIKING " << g << '\n';
                return out;
            }

            /*
             - t is y1
             - w is z
             - w1 is z1
             - x is y
             - x1 is ??
             */

            for (int i = pos_[y]; i != pos_[y1]; i=(i+add)%n){
                poly.push_back(B_[i]);
            }
            poly.push_back(B_[pos_[y1]);
            for (int i = pos_[x]; i != pos_[t]; i=(i+add)%n){
                poly.push_back(B_[i]);
            }
            poly.push_back(B_[pos_[t]]);
            for (int i = pos_[z]; i != pos_[y]; i=(i+add)%n){
                poly.push_back(B_[i]);
            }

            std::cout << "POLY\n";
            for (int i : poly)
                std::cout << i << ' ';
            std::cout << '\n';

            callLemma4(poly, y1, z,t);

            out.inserted = has_edge(u, v);

            return out;
        }

    private:
        const std::vector<Point>& P_;
        const std::vector<int>&   B_;   // CCW boundary
        std::vector<Face>&        F_;

        std::unordered_map<int,int> pos_;
        std::vector<std::unordered_set<int>> adj_;

        bool        debug_ = false;
        std::ostream* dbg_ = nullptr;

        static void build_adjacency_from_faces(const std::vector<Face>& F,
                                               size_t n_vertices,
                                               std::vector<std::unordered_set<int>>& adj)
        {
            adj.assign(n_vertices, {});
            for (const auto& T : F) {
                int a=T[0], b=T[1], c=T[2];
                adj[a].insert(b); adj[b].insert(a);
                adj[b].insert(c); adj[c].insert(b);
                adj[c].insert(a); adj[a].insert(c);
            }
        }


        int degree_of_vertex(int v)
        {
           // build_adjacency_from_faces(F, n_vertices, adj);
            return (int)adj_[v].size();
        }


        static std::vector<int> neighbors_of_vertex(int v,
                                                    const std::vector<Face>& F,
                                                    size_t n_vertices)
        {
            std::vector<std::unordered_set<int>> adj;
            //build_adjacency_from_faces(F, n_vertices, adj);
            return std::vector<int>(adj[v].begin(), adj[v].end());
        }

        inline void callLemma4(std::vector<int>& poly, int x1, int w, int y1){

            auto poly1 = poly;
            sort(poly1.begin(), poly1.end());

            std::vector<std::array<int, 3>> newF, oldF;
            for (auto f : F_){
                bool ok =1;
                for (int i = 0; i < 3; ++i)
                    if (!std::binary_search(poly1.begin(), poly1.end(), f[i])){
                        ok = 0;
                        break;
                    }
                if (ok) {
                    newF.push_back(f);
                } else
                    oldF.push_back(f);
            }

            std::vector<Point> pt;
            int vq = -1, v1 = -1, vy = -1;
            std::vector<int> lbl(P_.size());
            for (int i : poly){
                pt.push_back(P_[i]);
                lbl[P_[i].lbl] = (int)pt.size() -1;
                if (P_[i].lbl == x1)
                    v1 = (int)pt.size() - 1;
                if (P_[i].lbl == w)
                    vq= (int)pt.size() - 1;
                if (P_[i].lbl == y1)
                    vy = (int)pt.size() - 1;
            }

            for (auto& f : newF){
                for (int i = 0; i < 3; ++i)
                    f[i] = lbl[f[i]];
            }

            cgshop::DiagonalIntroducerLemma4 L4(pt, newF);
            L4.enable_debug(std::cout);

            auto st = L4.introduce(v1, vq);

            for (auto& f : newF){
                for (int i = 0; i < 3; ++i)
                    f[i] = pt[f[i]].lbl;
            }

            F_.clear();
            for (auto& f : newF)
                F_.push_back(f);
            for (auto& f : oldF)
                F_.push_back(f);


            if (!has_edge(x1, w)){
                std::cout << "L4 Failed first\n";
            } else {
                std::cout << "L4 Succeeded first\n";
            }

            std::cout << x1 << ' ' << y1 << '\n';

            if (!has_edge(x1,y1)){

                for (auto& f : newF){
                    for (int i = 0; i < 3; ++i)
                        f[i] = lbl[f[i]];
                }
                cgshop::DiagonalIntroducerLemma4 LN4(pt, newF);
                LN4.enable_debug(std::cout);
                st = LN4.introduce(v1,vy);
                for (auto& f : newF){
                    for (int i = 0; i < 3; ++i)
                        f[i] = pt[f[i]].lbl;
                }
                F_.clear();
                for (auto& f : newF)
                    F_.push_back(f);
                for (auto& f : oldF)
                    F_.push_back(f);
            }

            if (!has_edge(x1, w)){
                std::cout << "L4 Failed Redo\n";
            }

            if (!has_edge(y1, x1)){
                std::cout << "L4 Failed Second\n";
            }

            std::cout << "L4 Passed\n";
        }


        bool has_edge(int a, int b) const {
            if (a > b) std::swap(a, b);
            for (const auto& T : F_) {
                int e0a = std::min(T[0], T[1]), e0b = std::max(T[0], T[1]);
                int e1a = std::min(T[1], T[2]), e1b = std::max(T[1], T[2]);
                int e2a = std::min(T[2], T[0]), e2b = std::max(T[2], T[0]);
                if ((e0a==a && e0b==b) || (e1a==a && e1b==b) || (e2a==a && e2b==b)) return true;
            }
            return false;
        }

        inline int d_ccw(int i, int j)
        {
            int n = (int)B_.size();
            int di = pos_.at(i), dj = pos_.at(j);
            int d = dj - di; if (d < 0) d += n; return d;
        }

        // iterate CCW from index ia to ib inclusive
        static inline std::vector<int> chain_ccw(const std::vector<int>& B,
                                                 const std::unordered_map<int,int>& pos,
                                                 int a, int b)
        {
            int n = (int)B.size();
            int ia = pos.at(a), ib = pos.at(b);
            std::vector<int> C;
            for (int k = ia; ; k = (k+1) % n) {
                C.push_back(B[k]);
                if (k == ib) break;
            }
            return C;
        }


        static inline bool on_open_segment(const Point& a, const Point& b, const Point& c) {
            if (orient2d(a,b,c) != 0) return false;
            double minx = std::min(a.x,b.x), maxx = std::max(a.x,b.x);
            double miny = std::min(a.y,b.y), maxy = std::max(a.y,b.y);
            return (c.x > minx && c.x < maxx && c.y > miny && c.y < maxy);
        }

        static inline bool proper_intersection(const Point& a,const Point& b,
                                               const Point& c,const Point& d)
        {
            int o1 = orient2d(a,b,c), o2 = orient2d(a,b,d);
            int o3 = orient2d(c,d,a), o4 = orient2d(c,d,b);
            if (o1==0 && on_open_segment(a,b,c)) return true;
            if (o2==0 && on_open_segment(a,b,d)) return true;
            if (o3==0 && on_open_segment(c,d,a)) return true;
            if (o4==0 && on_open_segment(c,d,b)) return true;
            return (o1*o2 < 0) && (o3*o4 < 0);
        }

// visibility of segment (x1,w) against polygon boundary (CCW)
// allows touching at endpoints, disallows any other touch
         inline bool visible_on_boundary(int x1, int w)
        {
            int n = (int)B_.size();
            for (int i = 0; i < n; ++i) {
                int u = B_[i], v = B_[(i+1)%n];
                // skip edges incident to x1 or w
                if (u==x1 || v==x1 || u==w || v==w) continue;
                if (proper_intersection(P_[x1], P_[w], P_[u], P_[v])) return false;
            }
            return true;
        }

        // local tangency test at w:
        // for CCW polygon, polygon interior is to the left of each oriented boundary edge.
        // A supporting tangent from x1 to the chain at w means both incident edges at w
        // lie on the same side of the line (x1 -> w). We use "right-or-collinear" w.r.t. that line.
        inline bool tangent_wedge_ok(int x1, int w, int w_prev, int w_next)
        {
            // sign of oriented area for triples (x1, w, neighbor)
            int s_prev = orient2d(P_[x1], P_[w], P_[w_prev]);
            int s_next = orient2d(P_[x1], P_[w], P_[w_next]);
            // require both neighbors not to be strictly on the left of line x1->w
            // i.e., polygon near w is on or to the right side of the supporting line.
            return (s_prev <= 0) && (s_next <= 0);
        }


        int find_tangent_to_chain(int x1, int x, int v)
        {
            // build positions for quick wrap-around
            std::unordered_map<int,int> pos;
            pos.reserve(B_.size());
            for (int i = 0; i < (int)B_.size(); ++i) pos[B_[i]] = i;

            auto C = chain_ccw(B_, pos, x, v); // includes both endpoints
            int m = (int)C.size();
            if (m <= 1) return -1;

            int best = -1;
            for (int j = 0; j < m; ++j) {
                int w = C[j];
                if (w == x1) continue; // skip trivial

                // neighbors of w along the chain C (use chain neighbors, not full boundary neighbors)
                int w_prev = C[(j-1+m)%m];
                int w_next = C[(j+1)%m];

                // segment must be visible
                if (!visible_on_boundary(x1, w)) continue;

                // local tangency condition at w
                if (!tangent_wedge_ok(x1, w, w_prev, w_next)) continue;

                best = w; // keep the farthest along the chain because we iterate CCW from x to v
            }
            return best; // -1 if none
        }

    };

} // namespace cgshop
