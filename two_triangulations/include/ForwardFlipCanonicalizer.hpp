#pragma once
#include <vector>
#include <array>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <limits>
#include <cassert>
#include <utility>
#include <cmath>

/*
    ForwardFlipCanonicalizer.hpp

    Model:
      - Polygon Q has vertices indexed 0..n-1 in CCW boundary order.
      - Q is “bi-chain”: its boundary decomposes into two chains from a left splitter to a right splitter:
          * upperChain:  u0, u1, ..., um   (left-to-right order along the upper boundary)
          * lowerChain:  l0, l1, ..., lp   (left-to-right order along the lower boundary)
        (The splitter/convex hull details don’t matter; just ensure each vertex index is in exactly one chain and
         the order along each chain is “left to right”.)
      - Triangulation is given as CCW faces: vector<array<int,3>>; all vertices are from the polygon.

    What this class does:
      - It repeatedly performs a PARALLEL round of flips. In each round, it selects all flippable diagonals (li, ut)
        such that flipping (li, ut) strictly increases the upper-chain rank of the top of the lower-base triangle
        Δ_i = (l_{i-1}, l_i, ut). Those flips are then filtered to be pairwise flip-independent and applied in parallel.
      - It terminates when no such flip exists, at which point the triangulation is the canonical T* of Lemma 3.

*/

#include "./TriCommon.hpp"

namespace cgshop {


    class ForwardFlipCanonicalizer {
    public:
        using Face = std::array<int,3>;
        struct FlipCand { EdgeKey e; int fLeft=-1, fRight=-1; /* incident faces */ };


        struct Stats {
            int rounds = 0;           // number of forward-flip rounds
            int flips_total = 0;      // total number of single-edge flips
            std::vector<std::vector<FlipCand>> flips;
        };

        ForwardFlipCanonicalizer(const std::vector<Point>& pts,
                                 const std::vector<int>& upperChain,
                                 const std::vector<int>& lowerChain,
                                 std::vector<Face>& faces /* in/out */)
                : P(pts), U(upperChain), L(lowerChain), F(faces)
        {
            const int n = static_cast<int>(P.size());
            rankU.assign(n, -1);
            rankL.assign(n, -1);
            for (int i = 0; i < (int)U.size(); ++i) rankU[U[i]] = i;
            for (int i = 0; i < (int)L.size(); ++i) rankL[L[i]] = i;
            buildIncidence();
            // Precompute base-pairs for lower chain (li-1, li), i = 1..|L|-1
            basePairs.reserve(L.size() ? L.size()-1 : 0);
            for (int i = 1; i < (int)L.size(); ++i)
                basePairs.emplace_back(L[i-1], L[i]);
        }

        // Run the forward-flip canonicalization. Mutates F in-place.
        Stats run() {
            Stats s;
            // Sanity: each lower base must have exactly one incident triangle with 3rd vertex on upper chain.
            // We recompute tops each round, so just assert existence on the fly inside the loop.

            while (true) {
                // 1) Compute current top vertex u_t(i) for every lower-base pair (l_{i-1}, l_i)
                std::vector<int> topU; topU.resize(basePairs.size(), -1);
                if (!computeCurrentTops(topU)) {
                    // Malformed triangulation for the given chains (no triangle over some base).
                    // We fail fast; in your pipeline, you might throw or return partial stats.
                    break;
                }

                // 2) Collect ALL flippable edge candidates (li, u_t(i)) that strictly advance top rightwards.
                //    Also collect the pair of adjacent faces for each such edge (for independence filtering).
                std::vector<FlipCand> cands;
                cands.reserve(basePairs.size());

                for (int i = 0; i < (int)basePairs.size(); ++i) {
                    int l_im1 = basePairs[i].first;
                    int l_i   = basePairs[i].second;
                    int u = topU[i];
                    if (u < 0) continue; // missing (shouldn’t happen if computeCurrentTops passed)
                    // Candidate is edge (l_i, u). Check if flippable and if flip increases upper rank.
                    EdgeKey e(l_i, u);
                    auto it = edge2faces.find(e);
                    if (it == edge2faces.end()) continue; // boundary or not present
                    const auto& inc = it->second;
                    if (inc.size() != 2) continue; // not flippable (must be interior)

                    int fA = inc[0], fB = inc[1];
                    int otherUpper = otherUpperVertexAcross(l_i, u, fA, fB);
                    if (otherUpper < 0) continue; // the opposite tri didn't have an upper vertex (degenerate for our two-chain model)
                    if (rankU[u] < rankU[otherUpper]) {
                        // Check convexity of quadrilateral around e, i.e., ensure flip is geometrically valid.
                        if (isFlippableConvex(l_i, u, fA, fB)) {
                            cands.push_back({e, fA, fB});
                        }
                    }
                }

                if (cands.empty()) break; // done — reached T*

                // 3) Filter to get a flip-independent subset: no two candidates share a face.
                //    Greedy is enough (the theoretical construction guarantees uniqueness, but we stay safe).
                std::unordered_set<int> usedFaces;
                std::vector<FlipCand> chosen;
                chosen.reserve(cands.size());
                for (const auto& c : cands) {
                    if (!faceValid(c.fLeft) || !faceValid(c.fRight)) continue;
                    if (usedFaces.count(c.fLeft) || usedFaces.count(c.fRight)) continue;
                    chosen.push_back(c);
                    usedFaces.insert(c.fLeft);
                    usedFaces.insert(c.fRight);
                }

                if (chosen.empty()) break; // no independent flips possible (should rarely happen)
                // 4) Apply flips in parallel: for each chosen edge, replace with the other diagonal in its quad.
                for (const auto& c : chosen) {
                    doFlip(c.e, c.fLeft, c.fRight);
                }
                s.rounds += 1;
                s.flips_total += static_cast<int>(chosen.size());
                s.flips.push_back(chosen);
                // 5) Rebuild incidence (simplest & robust). For performance you can local-update instead.
                rebuildIncidence();
            }

            return s;
        }

    private:
        const std::vector<Point>& P;
        const std::vector<int>& U; // upper chain, left -> right
        const std::vector<int>& L; // lower chain, left -> right
        std::vector<Face>& F;      // triangulation faces (CCW), mutated

        std::vector<int> rankU, rankL; // -1 if not on that chain
        std::vector<std::pair<int,int>> basePairs;

        // Edge -> incident faces indices
        std::unordered_map<EdgeKey, std::vector<int>, EdgeKeyHash> edge2faces;

        // ----------------- geometry helpers -----------------
        static double orient2d(const Point& a, const Point& b, const Point& c) {
            // robust enough for our usage; tune EPS if needed
            long double x1 = b.x - a.x;
            long double y1 = b.y - a.y;
            long double x2 = c.x - a.x;
            long double y2 = c.y - a.y;
            long double v = x1*y2 - y1*x2;
            return static_cast<double>(v);
        }
        static bool ccw(const Point& a, const Point& b, const Point& c, double eps = 0) {
            return orient2d(a,b,c) > eps;
        }

        // Return vertices of face f in CCW order
        inline std::array<int,3> faceVerts(int f) const {
            return F[f];
        }
        inline bool faceValid(int f) const {
            return f >= 0 && f < (int)F.size();
        }

        // Make sure e=(x,y) is an edge of face f, and return the third vertex (opposite)
        static int thirdVertexOf(const Face& T, int x, int y) {
            int c = -1;
            int cnt = 0;
            for (int k=0;k<3;++k) {
                int v = T[k];
                if (v != x && v != y) { c = v; ++cnt; }
            }
            return (cnt==1) ? c : -1;
        }

        void buildIncidence() {
            edge2faces.clear();
            edge2faces.reserve(F.size()*3);
            for (int i=0;i<(int)F.size();++i) {
                const auto& t = F[i];
                addEdgeInc(EdgeKey(t[0], t[1]), i);
                addEdgeInc(EdgeKey(t[1], t[2]), i);
                addEdgeInc(EdgeKey(t[2], t[0]), i);
            }
        }
        void rebuildIncidence() { buildIncidence(); }

        void addEdgeInc(const EdgeKey& e, int fIdx) {
            auto& v = edge2faces[e];
            if (v.empty() || (v.size()==1 && v[0]!=fIdx)) v.push_back(fIdx);
            // if >2, triangulation invalid; we ignore silently.
        }

        // Compute current "top" for each lower base (l_{i-1}, l_i): find the unique face with those two and a 3rd on upper chain.
        bool computeCurrentTops(std::vector<int>& topU) const {
            std::unordered_map<EdgeKey, std::vector<int>, EdgeKeyHash> base2faces;
            base2faces.reserve(basePairs.size()*2+7);

            // Gather faces by lower edges (unordered key)
            for (int f=0; f<(int)F.size(); ++f) {
                const auto& T = F[f];
                for (int e=0;e<3;++e) {
                    int a = T[e], b = T[(e+1)%3], c = T[(e+2)%3];
                    if (rankL[a] >= 0 && rankL[b] >= 0) {
                        base2faces[EdgeKey(a,b)].push_back(f);
                    }
                    if (rankL[b] >= 0 && rankL[c] >= 0) {
                        base2faces[EdgeKey(b,c)].push_back(f);
                    }
                    if (rankL[c] >= 0 && rankL[a] >= 0) {
                        base2faces[EdgeKey(c,a)].push_back(f);
                    }
                }
            }
            // For each base edge (li-1, li), locate its (single) face whose 3rd vertex is on upper chain; that third is the current top.
            for (int i = 0; i < (int)basePairs.size(); ++i) {
                EdgeKey base(basePairs[i].first, basePairs[i].second);
                auto it = base2faces.find(base);
                if (it == base2faces.end()) return false;
                int foundTop = -1;
                for (int f : it->second) {
                    const auto& T = F[f];
                    int c = thirdVertexOf(T, base.a, base.b);
                    if (c < 0) continue;
                    if (rankU[c] >= 0) {
                        foundTop = c; break;
                    }
                }
                if (foundTop < 0) return false;
                topU[i] = foundTop;
            }
            return true;
        }

        // Given interior edge (li, u) with incident faces fA, fB, return the OTHER upper-chain vertex across that edge (if any).
        // If neither adjacent triangle has an upper vertex != u, returns -1.
        int otherUpperVertexAcross(int li, int u, int fA, int fB) const {
            const int cand[2] = { fA, fB };
            for (int k=0;k<2;++k) {
                int f = cand[k];
                const auto& T = F[f];
                int c = thirdVertexOf(T, li, u);
                if (c >= 0 && c != li && c != u && rankU[c] >= 0) return c;
            }
            return -1;
        }

        // Check if flipping diagonal (li, u) between faces fA and fB is geometrically valid (forms convex quad).
        bool isFlippableConvex(int li, int u, int fA, int fB) const {
            if (!faceValid(fA) || !faceValid(fB)) return false;
            int a = li, b = u;
            int x = thirdVertexOf(F[fA], a, b);
            int y = thirdVertexOf(F[fB], a, b);
            if (x < 0 || y < 0 || x == y || x == a || x == b || y == a || y == b) return false;

            // Convexity test: quad (x, a, y, b) or (x, b, y, a) must be convex. Determine correct cyclic order.
            // We know faces are CCW; pick order so that (a,b) is the diagonal between triangles (a,b,x) and (a,b,y).
            // We test both plausible orders and accept if either yields convexity.
            return quadConvex(x,a,y,b) || quadConvex(x,b,y,a);
        }

        bool quadConvex(int q0, int q1, int q2, int q3) const {
            const Point& A = P[q0];
            const Point& B = P[q1];
            const Point& C = P[q2];
            const Point& D = P[q3];
            // All consecutive turns must be CCW for strict convexity.
            const double eps = 0.0;
            return ccw(A,B,C,eps) && ccw(B,C,D,eps) && ccw(C,D,A,eps) && ccw(D,A,B,eps);
        }

        // Perform the flip of diagonal e=(li,u) between faces fA and fB, replacing it with the other diagonal (x,y),
        // where x = thirdVertex(F[fA], li, u) and y = thirdVertex(F[fB], li, u). Updates faces and edge incidence maps.
        void doFlip(const EdgeKey& e, int fA, int fB) {
            assert(faceValid(fA) && faceValid(fB));
            int a = e.a, b = e.b;
            int x = thirdVertexOf(F[fA], a, b);
            int y = thirdVertexOf(F[fB], a, b);
            if (x < 0 || y < 0 || x == y) return; // safety

            // Faces before flip: (a,b,x) and (a,b,y) in some CCW order.
            // Faces after flip: (x,y,a) and (y,x,b) (choose orientations to keep CCW).
            // We’ll reconstruct each face so that its vertices are CCW using the orientation test.

            auto makeCCW = [&](int p, int q, int r)->Face {
                if (ccw(P[p], P[q], P[r])) return Face{p,q,r};
                else return Face{p,r,q};
            };

            Face f1 = makeCCW(x, y, a);
            Face f2 = makeCCW(y, x, b);

            // Remove old edges from incidence
            removeFaceFromIncidence(fA);
            removeFaceFromIncidence(fB);

            // Assign new faces
            F[fA] = f1;
            F[fB] = f2;

            // Re-add edges
            addFaceToIncidence(fA);
            addFaceToIncidence(fB);
        }

        void removeFaceFromIncidence(int f) {
            const auto& T = F[f];
            EdgeKey e1(T[0],T[1]), e2(T[1],T[2]), e3(T[2],T[0]);
            removeEdgeFace(e1, f);
            removeEdgeFace(e2, f);
            removeEdgeFace(e3, f);
        }
        void addFaceToIncidence(int f) {
            const auto& T = F[f];
            addEdgeInc(EdgeKey(T[0],T[1]), f);
            addEdgeInc(EdgeKey(T[1],T[2]), f);
            addEdgeInc(EdgeKey(T[2],T[0]), f);
        }
        void removeEdgeFace(const EdgeKey& e, int f) {
            auto it = edge2faces.find(e);
            if (it == edge2faces.end()) return;
            auto& vec = it->second;
            vec.erase(std::remove(vec.begin(), vec.end(), f), vec.end());
            if (vec.empty()) edge2faces.erase(it);
        }
    };

}
