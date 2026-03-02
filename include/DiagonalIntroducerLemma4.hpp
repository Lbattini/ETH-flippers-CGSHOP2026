#pragma once
#include <vector>
#include <array>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <stdexcept>
#include <cassert>
#include <cmath>
#include <utility>
#include "./TriCommon.hpp"

//
// DiagonalIntroducerLemma4.hpp
//
// Goal: introduce a target diagonal e = (v1, vq) by parallel rounds of flips
// that monotonically reduce the number of crossings with e inside the
// "active polygon" (the union of triangles currently crossing e).
//
// Works with the same basic types as ForwardFlipCanonicalizer.hpp
// (Point, Face, EdgeKey, edge incidence, convexity tests).
//

namespace cgshop {

    using Face = std::array<int,3>;

// ---------- Utility predicates ----------
    static inline double orient2d(const Point& a, const Point& b, const Point& c) {
        long double x1=b.x-a.x, y1=b.y-a.y;
        long double x2=c.x-a.x, y2=c.y-a.y;
        long double v = x1*y2 - y1*x2;
        return (double)v;
    }
    static inline bool ccw(const Point& a, const Point& b, const Point& c, bool flag = 0, double eps=1e-12) {
        if (flag) std::cout << orient2d(a,b,c) << '\n';
        return orient2d(a,b,c) > 0;
    }
    static inline bool segments_strictly_intersect(const Point& a, const Point& b,
                                                   const Point& c, const Point& d) {
        auto sgn = [](double v){ return (v>0)-(v<0); };
        double o1 = orient2d(a,b,c), o2 = orient2d(a,b,d);
        double o3 = orient2d(c,d,a), o4 = orient2d(c,d,b);

        // Proper (non-collinear) intersection
        return (sgn(o1)*sgn(o2) < 0) && (sgn(o3)*sgn(o4) < 0);
    }

    class DiagonalIntroducerLemma4 {
    public:
        struct Stats {
            int rounds = 0;           // # parallel rounds
            int flips_total = 0;      // # single-edge flips
            bool inserted = false;    // did we manage to insert e?
            std::vector<std::vector<EdgeKey>> flips;
        };

        // Construct with a polygon triangulation
        DiagonalIntroducerLemma4(const std::vector<Point>& pts,
                                 std::vector<Face>& faces /* in/out */)
                : P(pts), F(faces)
        {
            buildIncidence();
        }



        // Try to introduce diagonal (v1, vq). Returns stats; mutates F.

        Stats introduce(int v1, int vq) {
            Stats s;

            for (auto p : P)
                std::cout << p << '\n';
            std::cout << v1 << ' ' << vq << '\n';


            if (hasEdge(v1, vq)) { s.inserted = true; return s; }

            auto edgeCrosses = [&](int a,int b)->bool {
                return segments_strictly_intersect(P[a], P[b], P[v1], P[vq]);
            };

            auto count_crossing_edges = [&]()->int {
                int cnt = 0;
                for (const auto& kv : edge2faces) {
                    if (kv.second.size() != 2) continue; // interior only
                    if (edgeCrosses(kv.first.a, kv.first.b)) ++cnt;
                }
                return cnt;
            };

            // safety guard
            int safe_guard = (int)P.size() * 30;

            dbg_round_ = 0;
            if (dbg_out_) {
                dbg("---- Lemma4 introduce start ----");
                dbg("target edge: (" + std::to_string(v1) + "," + std::to_string(vq) + ")");
                dbg("initial crossing edges: " + std::to_string(count_crossing_edges()));
            }

            while (safe_guard-- > 0) {
                ++dbg_round_;

                // stats for rejects
                int reject_boundary = 0, reject_third_bad = 0, reject_still_cross = 0, reject_nonconvex = 0;

                struct Cand { EdgeKey e; int fL=-1, fR=-1; };
                std::vector<Cand> cands; cands.reserve(edge2faces.size());

                for (const auto& kv : edge2faces) {
                    const EdgeKey& e = kv.first;
                    const auto& inc = kv.second;

                    // boundary?
                    if (inc.size() != 2) { if (segments_strictly_intersect(P[e.a], P[e.b], P[v1], P[vq])) ++reject_boundary; continue; }

                    int a = e.a, b = e.b;
                    int fA = inc[0], fB = inc[1];

                    int x = thirdVertexOf(F[fA], a, b);
                    int y = thirdVertexOf(F[fB], a, b);
                    if (x < 0 || y < 0 || x == y) {
                        if (segments_strictly_intersect(P[a], P[b], P[v1], P[vq])) ++reject_third_bad;
                        continue;
                    }

                    bool curCross = edgeCrosses(a,b);
                    if (!curCross) continue; // we only try to remove crossings

                    bool newCross = edgeCrosses(x,y);
                    if (newCross) {
                        ++reject_still_cross;
                        std::cout << "STILL INTERSECTNG " << x << ' ' << y << '\n';
                        continue;
                    }

                    if (!isConvexQuadFlipOK(a,b,x,y)) {
                        std::cout << "NONV CONVEX ";
                        ++reject_nonconvex; continue;
                    }

                    cands.push_back({e, fA, fB});
                }

                // log pre-choice stats
                if (dbg_out_) {
                    dbg("[round " + std::to_string(dbg_round_) + "]");
                    dbg("  crossing edges          : " + std::to_string(count_crossing_edges()));
                    dbg("  candidate flips         : " + std::to_string(cands.size()));
                    dbg("  rejects (boundary)      : " + std::to_string(reject_boundary));
                    dbg("  rejects (third bad)     : " + std::to_string(reject_third_bad));
                    dbg("  rejects (still crossing): " + std::to_string(reject_still_cross));
                    dbg("  rejects (non-convex)    : " + std::to_string(reject_nonconvex));
                }

                if (cands.empty()) {
                    int va, vb;
                    // Try direct insertion (convex quad with (v1,vq) as other diagonal)
                    if (tryInsertDiagonalDirect(v1, vq, va, vb)) {
                        s.inserted = true;
                        s.flips.push_back({{va, vb}});
                        s.flips_total += 1;
                        s.rounds += 1;
                        if (dbg_out_) dbg("  direct insertion succeeded.");
                        return s;
                    }
                    // If no interior crossings remain, we *should* be able to insert; if not, bail with info.
                    int left = count_crossing_edges();
                    if (dbg_out_) {
                        dbg("  NO CANDIDATES; crossing edges remaining: " + std::to_string(left));
                        dbg("  hasEdge(v1,vq) = " + std::to_string(hasEdge(v1,vq)));
                    }
                    s.inserted = hasEdge(v1, vq);
                    return s;
                }

                // choose maximal independent set (no shared face)
                std::unordered_set<int> used;
                std::vector<Cand> chosen; chosen.reserve(cands.size());
                for (const auto& c : cands) {
                    if (used.count(c.fL) || used.count(c.fR)) continue;
                    chosen.push_back(c);
                    used.insert(c.fL);
                    used.insert(c.fR);
                }

                if (dbg_out_) {
                    dbg("  chosen (independent)    : " + std::to_string(chosen.size()));
                }

                std::vector<EdgeKey> tmp;
                for (const auto& c : chosen) {
                    doFlip(c.e, c.fL, c.fR);
                    tmp.push_back(c.e);
                }
                s.flips.push_back(tmp);
                s.rounds += 1;
                s.flips_total += (int)chosen.size();

                if (hasEdge(v1, vq)) { s.inserted = true; if (dbg_out_) dbg("  target edge present after flips."); return s; }
                int va, vb;
                if (tryInsertDiagonalDirect(v1, vq, va, vb)) {
                    s.inserted = true;
                    s.flips.push_back({{va, vb}});
                    s.flips_total += 1;
                    s.rounds += 1;
                    if (dbg_out_) dbg("  target edge inserted directly after flips.");
                    return s; }
            }

            if (dbg_out_) dbg("  hit safety guard; aborting.");
            s.inserted = hasEdge(v1, vq);
            return s;
        }


        void enable_debug(std::ostream& os) { dbg_out_ = &os; }



    private:
        const std::vector<Point>& P;
        std::vector<Face>& F;

        std::unordered_map<EdgeKey, std::vector<int>, EdgeKeyHash> edge2faces;

        // --- debug members ---
        std::ostream* dbg_out_ = nullptr;
        int dbg_round_ = 0;

        void dbg(const std::string& s) const {
            if (dbg_out_) (*dbg_out_) << s << '\n';
        }


        bool anyEdgeCrossesTarget(int v1, int vq) const {
            for (const auto& kv : edge2faces) {
                if (kv.second.size() != 2) continue; // interior only
                int a = kv.first.a, b = kv.first.b;
                if (segments_strictly_intersect(P[a], P[b], P[v1], P[vq])) return true;
            }
            return false;
        }


        // ---------- incidence ----------
        void buildIncidence() {
            edge2faces.clear();
            edge2faces.reserve(F.size()*3);
            for (int i=0;i<(int)F.size();++i) {
                const auto& t = F[i];
                addEdgeInc( EdgeKey(t[0],t[1]), i );
                addEdgeInc( EdgeKey(t[1],t[2]), i );
                addEdgeInc( EdgeKey(t[2],t[0]), i );
            }
        }
        void rebuildIncidence(){ buildIncidence(); }

        void addEdgeInc(const EdgeKey& e, int fIdx) {
            auto& v = edge2faces[e];
            if (v.empty() || (v.size()==1 && v[0]!=fIdx)) v.push_back(fIdx);
        }
        void removeFaceFromIncidence(int f) {
            const auto& T = F[f];
            removeEdgeFace( EdgeKey(T[0],T[1]), f );
            removeEdgeFace( EdgeKey(T[1],T[2]), f );
            removeEdgeFace( EdgeKey(T[2],T[0]), f );
        }
        void addFaceToIncidence(int f) {
            const auto& T = F[f];
            addEdgeInc( EdgeKey(T[0],T[1]), f );
            addEdgeInc( EdgeKey(T[1],T[2]), f );
            addEdgeInc( EdgeKey(T[2],T[0]), f );
        }
        void removeEdgeFace(const EdgeKey& e, int f) {
            auto it = edge2faces.find(e);
            if (it == edge2faces.end()) return;
            auto& vec = it->second;
            vec.erase(std::remove(vec.begin(), vec.end(), f), vec.end());
            if (vec.empty()) edge2faces.erase(it);
        }

        static int thirdVertexOf(const Face& T, int x, int y) {
            int c = -1, cnt=0;
            for (int k=0;k<3;++k) if (T[k]!=x && T[k]!=y) { c=T[k]; ++cnt; }
            return (cnt==1)? c : -1;
        }

        bool quadConvex(int q0, int q1, int q2, int q3, bool flag = 0) const {
            const Point& A = P[q0];
            const Point& B = P[q1];
            const Point& C = P[q2];
            const Point& D = P[q3];

            if (flag)
                std::cout << A << ' ' << B << ' ' << C << ' ' << D << '\n';
            return ccw(A,B,C, flag) && ccw(B,C,D, flag) && ccw(C,D,A, flag) && ccw(D,A,B, flag);
        }
        bool isConvexQuadFlipOK(int a, int b, int x, int y) const {
            // Triangle faces are (a,b,x) and (a,b,y). New diagonal will be (x,y). Check convexity.
            // Accept if either cyclic order yields convex quad.

            if (a == 3 && b == 9 && x == 4 && y == 10) {

                std::cout << "ALI " << quadConvex(x, b, y, a, 1) << " GOTOVO\n";
                exit(0);
            }
            return quadConvex(x,a,y,b) || quadConvex(x,b,y,a);
        }

        void doFlip(const EdgeKey& e, int fA, int fB) {
            int a=e.a, b=e.b;
            int x = thirdVertexOf(F[fA], a, b);
            int y = thirdVertexOf(F[fB], a, b);
            if (x < 0 || y < 0 || x == y) return;

            std::cout << "FLIPPING " << P[a].lbl << ' ' << P[b].lbl << '\n';

            auto makeCCW = [&](int p,int q,int r)->Face {
                if (ccw(P[p],P[q],P[r])) return Face{p,q,r};
                return Face{p,r,q};
            };

            // remove old
            removeFaceFromIncidence(fA);
            removeFaceFromIncidence(fB);

            // new faces share diagonal (x,y)
            F[fA] = makeCCW(x, y, a);
            F[fB] = makeCCW(y, x, b);

           // std::cout << "FLIPPED " << a << ' ' << b << '\n';
           // std::cout << "NEW FACES " <<x << ' ' << y << ' ' << a << ' ' << ", " << y << ' ' << x << ' ' << b << '\n';

            // add new
            addFaceToIncidence(fA);
            addFaceToIncidence(fB);
        }

        bool hasEdge(int u, int v) const {
            EdgeKey e(u,v);
            auto it = edge2faces.find(e);
            return it != edge2faces.end();
        }

        // Try to directly insert (v1,vq) if there is a convex quad around it:
        // meaning: there exist two adjacent triangles whose union forms a convex quadrilateral
        // with (v1,vq) as the other diagonal.
        bool tryInsertDiagonalDirect(int v1, int vq, int& va, int& vb) {
            // Look for a pair of adjacent faces sharing an edge (x,y) such that replacing (x,y) by (v1,vq) is valid.
            // We require v1 and vq to be the third vertices of the two adjacent faces across (x,y).
            for (const auto& kv : edge2faces) {
                if (kv.second.size() != 2) continue;
                int fA = kv.second[0], fB = kv.second[1];
                const auto& T1 = F[fA]; const auto& T2 = F[fB];
                int a = kv.first.a, b = kv.first.b;
                int x = thirdVertexOf(T1, a, b);
                int y = thirdVertexOf(T2, a, b);
                if (x<0||y<0) continue;
                // Want {x,y} = {v1,vq}
                if ((x==v1 && y==vq) || (x==vq && y==v1)) {
                    if (isConvexQuadFlipOK(a,b,v1,vq)) {
                        doFlip( EdgeKey(a,b), fA, fB );
                        va = a;
                        vb = b;
                        return true;
                    }
                }
            }
            return false;
        }

        struct ActiveInfo {
            // faces & interior edges whose *open* geometry intersects segment (v1,vq)
            std::vector<int> faces;
            std::vector<EdgeKey> interiorEdges;
        };

        ActiveInfo makeActiveInfo(int v1, int vq) const {
            ActiveInfo A;

            // Signed distance (orientation) of a point to the line (v1,vq)
            const Point& S = P[v1];
            const Point& T = P[vq];
            auto side = [&](int vid)->double {
                return orient2d(S, T, P[vid]); // >0 left, <0 right of (v1->vq)
            };

            // 1) Mark faces as "active" if (a) any edge strictly intersects (v1,vq)
            //    OR (b) the triangle has vertices on both sides of the line (v1,vq)
            std::vector<char> faceActive(F.size(), 0);
            for (int f = 0; f < (int)F.size(); ++f) {
                const auto& tri = F[f];
                int a = tri[0], b = tri[1], c = tri[2];

                bool edgeCross =
                        segments_strictly_intersect(P[a], P[b], S, T) ||
                        segments_strictly_intersect(P[b], P[c], S, T) ||
                        segments_strictly_intersect(P[c], P[a], S, T);

                if (edgeCross) {
                    faceActive[f] = 1;
                    continue;
                }

                // side-based interior crossing test
                double sa = side(a), sb = side(b), sc = side(c);
                // consider a small epsilon to ignore very-near-line noise
                const double eps = 0.0;
                bool pos = (sa > eps) || (sb > eps) || (sc > eps);
                bool neg = (sa < -eps) || (sb < -eps) || (sc < -eps);

                if (pos && neg) {
                    faceActive[f] = 1;
                }
            }

            // Collect active faces list
            for (int f = 0; f < (int)F.size(); ++f)
                if (faceActive[f]) A.faces.push_back(f);

            // 2) Interior edges whose BOTH incident faces are active
            //    (so flipping them can change the active polygon)
            std::unordered_set<long long> seen;
            auto pack = [](int a, int b)->long long {
                if (a > b) std::swap(a, b);
                return ( (long long)a << 32 ) | (unsigned long long)b;
            };

            for (int f : A.faces) {
                const auto& tri = F[f];
                const int e[3][2] = { {tri[0], tri[1]}, {tri[1], tri[2]}, {tri[2], tri[0]} };
                for (int k = 0; k < 3; ++k) {
                    EdgeKey ek(e[k][0], e[k][1]);
                    auto it = edge2faces.find(ek);
                    if (it == edge2faces.end()) continue;
                    // interior edge must have exactly two incident faces, both active
                    if (it->second.size() == 2) {
                        int fA = it->second[0], fB = it->second[1];
                        if (faceActive[fA] && faceActive[fB]) {
                            long long key = pack(ek.a, ek.b);
                            if (!seen.count(key)) {
                                A.interiorEdges.push_back(ek);
                                seen.insert(key);
                            }
                        }
                    }
                }
            }

            return A;
        }


        // Is (v1,vq) currently a boundary diagonal of the active polygon (i.e., insertable by a single flip)?
        // We approximate this by checking if there exist *two* adjacent faces that currently straddle (v1,vq):
        // one on each side, such that flipping their shared edge to (v1,vq) is convex-valid.
        bool isInsertableNow(int v1, int vq) {
            // Fast path: if edge already exists.
            if (hasEdge(v1, vq)) return true;

            // Look for a convex quad where (v1,vq) is the missing diagonal (as in tryInsertDiagonalDirect).
            for (const auto& kv : edge2faces) {
                if (kv.second.size()!=2) continue;
                int fA = kv.second[0], fB = kv.second[1];
                const auto& T1 = F[fA]; const auto& T2 = F[fB];
                int a = kv.first.a, b = kv.first.b;
                int x = thirdVertexOf(T1, a, b);
                int y = thirdVertexOf(T2, a, b);
                if (x<0||y<0) continue;
                if ((x==v1 && y==vq) || (x==vq && y==v1)) {
                    if (isConvexQuadFlipOK(a,b,v1,vq)) return true;
                }
            }
            return false;
        }
    };

} // namespace cgshop
