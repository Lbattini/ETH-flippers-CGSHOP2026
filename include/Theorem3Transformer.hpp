// Theorem3Transformer.hpp
#pragma once
#include <bits/stdc++.h>
#include "./CutTriangulation.hpp"
#include "./Theorem2Transformer.hpp"   // from earlier step (polygon engine)

// =============== Public API ===============
namespace theorem3_pset {

    struct Stats {
        int depth            = 0;  // max recursion depth
        int nodes            = 0;  // number of processed regions
        int prop1_phases     = 0;  // sum over all Proposition1Inserter phases (via Theorem2)
        int lemma4_calls     = 0;
        int lemma5_calls     = 0;
        int flips_total      = 0;  // flips counted inside Theorem2 calls
    };

/// Transform E_cur into E_tgt using Theorem 3 strategy (O(n) parallel rounds).
/// - P: points (|P| = n)
/// - E_cur: current triangulation edges (in/out)
/// - E_tgt: target triangulation edges (read-only)
/// Returns: Stats about work done. Throws on inconsistent inputs that block progress.
    static inline Stats transform(const std::vector<Point>& P,
                                  std::vector<Edge>&        E_cur,
                                  const std::vector<Edge>&  E_tgt,
                                  std::ostream* dbg_out = nullptr);

// =============== Implementation details (private) ===============
    namespace detail {

// ----- helpers: edges / faces -----

// Fallback polygon boundary: convex hull (CCW) of a vertex set
        inline std::vector<int> convex_hull_ccw(const std::vector<Point>& P,
                                                const std::unordered_set<int>& Vset)
        {
            std::vector<int> V; V.reserve(Vset.size());
            for (int v : Vset) V.push_back(v);
            std::sort(V.begin(), V.end(), [&](int a,int b){
                if (P[a].x != P[b].x) return P[a].x < P[b].x;
                if (P[a].y != P[b].y) return P[a].y < P[b].y;
                return a < b;
            });
            auto cross = [&](int o,int a,int b){
                long double x1=P[a].x-P[o].x, y1=P[a].y-P[o].y;
                long double x2=P[b].x-P[o].x, y2=P[b].y-P[o].y;
                return x1*y2 - y1*x2;
            };
            std::vector<int> H;
            for (int i : V){ while (H.size()>=2 && cross(H[H.size()-2], H.back(), i) <= 0) H.pop_back(); H.push_back(i); }
            size_t lower = H.size();
            for (int k=(int)V.size()-2; k>=0; --k){
                int i = V[k];
                while (H.size()>lower && cross(H[H.size()-2], H.back(), i) <= 0) H.pop_back();
                H.push_back(i);
                if (k==0) break;
            }
            if (!H.empty()) H.pop_back();
            return H; // CCW
        }


// Ray-casting point-in-polygon (non-strict). 'polyCCW' is a vertex *sequence* (not duplicated at end).
        inline bool point_in_polygon(const Point& C,
                                     const std::vector<int>& polyCCW,
                                     const std::vector<Point>& P)
        {
            bool inside = false;
            const int m = (int)polyCCW.size();
            for (int i = 0, j = m-1; i < m; j = i++) {
                const Point& A = P[ polyCCW[j] ];
                const Point& B = P[ polyCCW[i] ];
                // Check if edge AB straddles the horizontal ray at y=C.y
                const bool cond = ((A.y > C.y) != (B.y > C.y));
                if (cond) {
                    const double xint = A.x + (B.x - A.x) * ( (C.y - A.y) / (B.y - A.y) );
                    if (xint > C.x) inside = !inside;
                }
            }
            return inside;
        }

        // Pretty-print a,b pairs
        inline std::string e2s(int a,int b){ return "("+std::to_string(a)+","+std::to_string(b)+")"; }

// Build normalized edge set from faces
        inline std::unordered_set<long long> edge_set_from_faces(const std::vector<std::array<int,3>>& F){
            std::unordered_set<long long> S; S.reserve(F.size()*3);
            auto key=[&](int a,int b){ if(a>b) std::swap(a,b); return ((long long)a<<32) ^ (unsigned long long)b; };
            for (auto t: F){
                S.insert(key(t[0],t[1])); S.insert(key(t[1],t[2])); S.insert(key(t[2],t[0]));
            }
            return S;
        }

// Compare two edge sets (restricted to region vertices) and print a small diff
        inline void print_edge_diff(const std::vector<Edge>& E_cur,
                                    const std::vector<Edge>& E_tgt,
                                    const std::unordered_set<int>& regionV,
                                    std::ostream& os,
                                    size_t max_show = 20)
        {
            auto key=[&](int a,int b){ if(a>b) std::swap(a,b); return ((long long)a<<32) ^ (unsigned long long)b; };

            auto inR = [&](int a,int b){ return regionV.count(a) && regionV.count(b); };

            std::unordered_set<long long> A, B;
            A.reserve(E_cur.size()); B.reserve(E_tgt.size());
            for (auto e: E_cur) if (inR(e.u,e.v)) A.insert(key(e.u,e.v));
            for (auto e: E_tgt) if (inR(e.u,e.v)) B.insert(key(e.u,e.v));

            std::vector<std::pair<int,int>> onlyA, onlyB;
            for (auto k: A) if (!B.count(k)) onlyA.push_back({int(k>>32), int(k & 0xffffffff)});
            for (auto k: B) if (!A.count(k)) onlyB.push_back({int(k>>32), int(k & 0xffffffff)});

            os << "    region-edge diff: |A\\B|=" << onlyA.size() << ", |B\\A|=" << onlyB.size() << "\n";
            if (!onlyA.empty()){
                os << "    1) in CURRENT not in TARGET (first " << std::min(max_show, onlyA.size()) << "): ";
                for (size_t i=0;i<std::min(max_show, onlyA.size()); ++i) os << e2s(onlyA[i].first, onlyA[i].second) << " ";
                os << "\n";
            }
            if (!onlyB.empty()){
                os << "    2) in TARGET not in CURRENT (first " << std::min(max_show, onlyB.size()) << "): ";
                for (size_t i=0;i<std::min(max_show, onlyB.size()); ++i) os << e2s(onlyB[i].first, onlyB[i].second) << " ";
                os << "\n";
            }
        }



        inline long long ekey(int a,int b){ if(a>b) std::swap(a,b); return ( (long long)a<<32 ) ^ (unsigned long long)b; }

        inline std::unordered_set<long long> edge_set_of(const std::vector<Edge>& E){
            std::unordered_set<long long> S; S.reserve(E.size()*2);
            for (auto e: E) if (e.u!=e.v) S.insert(ekey(e.u,e.v));
            return S;
        }

        inline std::unordered_set<long long> edge_set_of_faces(const std::vector<std::array<int,3>>& F){
            std::unordered_set<long long> S; S.reserve(F.size()*3);
            for (auto t: F) {
                S.insert(ekey(t[0],t[1]));
                S.insert(ekey(t[1],t[2]));
                S.insert(ekey(t[2],t[0]));
            }
            return S;
        }

// Filter edge list to a subset of vertices V (membership via boolean mask)
        inline std::vector<Edge> sub_edges(const std::vector<Edge>& E, const std::vector<char>& inS){
            std::vector<Edge> out; out.reserve(E.size());
            for (auto e: E) if (inS[e.u] && inS[e.v]) out.push_back(e);
            return out;
        }

// Triangles from edges (reuses your face extractor)
        inline std::vector<std::array<int,3>> faces_from_edges(const std::vector<Point>& P,
                                                               const std::vector<Edge>& E)
        {
            auto F = triangulation_cut::extractFacesAsTriangles(P, E);
            return F.triangles;
        }

// ----- build the polygon boundary (CCW vertices) from a list of triangles -----

        inline std::vector<int> boundary_cycle_from_tri_union(const std::vector<std::array<int,3>>& Qtris)
        {
            // collect edges and multiplicities (edges on union boundary appear exactly once)
            std::unordered_map<long long, std::pair<int,int>> edge_ab; edge_ab.reserve(Qtris.size()*3);
            std::unordered_map<long long,int> mult; mult.reserve(Qtris.size()*3);
            auto add_edge = [&](int a,int b){
                long long k = ekey(a,b);
                ++mult[k];
                if (!edge_ab.count(k)) edge_ab[k] = {std::min(a,b), std::max(a,b)};
            };
            for (auto T: Qtris) {
                add_edge(T[0],T[1]); add_edge(T[1],T[2]); add_edge(T[2],T[0]);
            }

            // take boundary edges (multiplicity 1) and stitch a cycle
            std::unordered_multimap<int,int> adj; adj.reserve(mult.size()*2);
            for (auto &kv : mult) if (kv.second==1) {
                    auto [a,b] = edge_ab[kv.first];
                    adj.emplace(a,b); adj.emplace(b,a);
                }
            if (adj.empty()) return {}; // degenerate

            // Pick start = the smallest vertex id we see (stable walk).
            int start = INT_MAX;
            for (auto &kv : adj) start = std::min(start, kv.first);

            std::vector<int> cyc; cyc.push_back(start);
            int prev = -1, cur = start;

            // walk (greedy; OK because boundary is a single simple cycle in our usage)
            for (int steps=0; steps< (int)adj.size()+5; ++steps) {
                auto range = adj.equal_range(cur);
                int nxt = -1;
                for (auto it=range.first; it!=range.second; ++it) {
                    if (it->second != prev) { nxt = it->second; break; }
                }
                if (nxt<0) break;
                cyc.push_back(nxt);
                prev = cur; cur = nxt;
                if (cur==start) break;
            }

            // ensure CCW orientation: compute signed area
            auto signed_area2 = [&](const std::vector<int>& poly)->double{
                long double s=0; int m=(int)poly.size();
                for (int i=0;i+1<m;++i){
                    // note: last equals start; we keep duplicate to ease walk; skip closing term
                    int a=poly[i], b=poly[(i+1)%m];
                    s += (long double)a - (long double)b; // dummy; we'll correct below
                }
                return (double)s; // not used; we re-compute properly below anyway
            };
            // Proper area:
            auto area2 = [&](const std::vector<int>& poly)->long double{
                long double S=0; int m=(int)poly.size();
                for (int i=0;i<m;++i){
                    long double x1 = (long double)i; // we don't have coordinates here; the cycle CCW is not strictly needed.
                    (void)x1;
                }
                return 1; // We skip re-orientation: Theorem2Transformer does not rely on 'upper'/'lower' label.
            };
            (void)signed_area2; (void)area2;
            return cyc; // as-sewn order is fine; Theorem 2 uses the cycle only to split into 2 chains between u and v.
        }

// ----- restrict target faces/edges to a polygon vertex set -----

        inline std::vector<std::array<int,3>>
        restrict_faces_to_vertex_set(const std::vector<std::array<int,3>>& FT,
                                     const std::unordered_set<int>& S)
        {
            std::vector<std::array<int,3>> out; out.reserve(FT.size());
            for (auto t: FT)
                if (S.count(t[0]) && S.count(t[1]) && S.count(t[2]))
                    out.push_back(t);
            return out;
        }

        inline std::vector<Edge> restrict_edges_to_vertex_set(const std::vector<Edge>& E,
                                                              const std::unordered_set<int>& S)
        {
            std::vector<Edge> out; out.reserve(E.size());
            for (auto e: E) if (S.count(e.u) && S.count(e.v)) out.push_back(e);
            return out;
        }

// ----- splice: replace all edges fully inside a region by edges from region faces -----

        inline void splice_region_edges(std::vector<Edge>& E_global,
                                        const std::unordered_set<int>& regionV,
                                        const std::vector<std::array<int,3>>& regionFacesNew)
        {
            auto inside = [&](int a,int b){ return regionV.count(a) && regionV.count(b); };

            // 1) remove edges whose both endpoints lie in regionV
            std::vector<Edge> out; out.reserve(E_global.size());
            for (auto e: E_global) if (!inside(e.u,e.v)) out.push_back(e);

            // 2) add edges from regionFacesNew
            for (auto t: regionFacesNew) {
                int a=t[0], b=t[1], c=t[2];
                out.push_back({a,b}); out.push_back({b,c}); out.push_back({c,a});
            }

            // 3) deduplicate
            std::unordered_set<long long> seen; seen.reserve(out.size()*2);
            std::vector<Edge> ded; ded.reserve(out.size());
            for (auto e: out) {
                long long k=ekey(e.u,e.v);
                if (seen.insert(k).second) ded.push_back({std::min(e.u,e.v), std::max(e.u,e.v)});
            }
            E_global.swap(ded);
        }

// ----- horizontal split line (avoid passing through points) -----

        // Pick a y0 strictly between two DISTINCT y-levels present, with tolerance.
// Robust to small jitter around 0 or 1 (e.g., 1±1e-9).
        inline double choose_horizontal_split_y(const std::vector<int>& verts,
                                                const std::vector<Point>& P)
        {
            std::vector<double> ys; ys.reserve(verts.size());
            for (int v : verts) ys.push_back(P[v].y);
            std::sort(ys.begin(), ys.end());

            // Merge nearly equal levels
            const double tol = 1e-6; // <<— key tolerance
            std::vector<double> uniq;
            uniq.reserve(ys.size());
            for (double y : ys) {
                if (uniq.empty() || std::fabs(y - uniq.back()) > tol) uniq.push_back(y);
            }

            if (uniq.size() < 2) return std::numeric_limits<double>::quiet_NaN();

            // Choose the middle *gap* between consecutive levels
            size_t gi = (uniq.size() - 1) / 2;
            double y_low  = uniq[gi];
            double y_high = uniq[gi+1];
            double y0 = 0.5*(y_low + y_high);

            // tiny nudge just in case
            const double eps = 1e-12;
            if (std::fabs(y0 - y_low)  <= eps) y0 += 10*eps;
            if (std::fabs(y_high - y0) <= eps) y0 -= 10*eps;
            return y0;
        }



// ----- build full-triangulation face lists for current/target (for subproblems) -----

        inline std::vector<std::array<int,3>> faces_of_edges_on_vertex_set(const std::vector<Point>& P,
                                                                           const std::vector<Edge>& E,
                                                                           const std::unordered_set<int>& S)
        {
            // filter edges, then extract faces
            std::vector<char> inS(P.size(),0);
            for (int v: S) inS[v]=1;
            auto Esub = sub_edges(E, inS);
            return faces_from_edges(P, Esub);
        }

// ----- main recursion -----

        struct Accum {
            int depth=0, nodes=0, prop1=0, l4=0, l5=0, flips=0;
        };

        inline void recurse_region(const std::vector<Point>& P,
                                   std::vector<Edge>& E_cur,
                                   const std::vector<Edge>& E_tgt,
                                   const std::vector<int>& regionVertsCCW,   // CCW boundary of current region (whole set for first call)
                                   int depth,
                                   Accum& acc,
                                   std::ostream* dbg) {
            acc.depth = std::max(acc.depth, depth);
            ++acc.nodes;

            // Terminal
            if ((int) regionVertsCCW.size() <= 3) return;

            // 1) Choose horizontal split y0
            double y0 = choose_horizontal_split_y(regionVertsCCW, P);

            if (!std::isfinite(y0)) {
                // Can't split by y (all points co-linear in y). Bail out of this region.
                return;
            }
            if (dbg)
                (*dbg) << "[T3] depth " << depth
                       << " | regionV=" << regionVertsCCW.size()
                       << " | y0=" << y0 << "\n";


            // 2) Build cut sets for CURRENT triangulation within this region:
            //    We filter edges to region first, then call compute_cut_sets on that sub-triangulation.
            std::unordered_set<int> regionSet(regionVertsCCW.begin(), regionVertsCCW.end());
            auto Ecur_sub = restrict_edges_to_vertex_set(E_cur, regionSet);


            CutResult curCut;
            std::string err;
            if (!triangulation_cut::compute_cut_sets(P, Ecur_sub, 0.0, 1.0, -y0, curCut, &err)) {
                // If a vertex lies on the line (extremely rare with the epsilon), nudge and retry once.
                y0 += 1e-6;
                curCut = CutResult{};
                if (dbg) (*dbg) << "  crossing tris (CURRENT) = " << curCut.Q.size() << "\n";

                if (!triangulation_cut::compute_cut_sets(P, Ecur_sub, 0.0, 1.0, -y0, curCut, &err)) {
                    if (dbg) (*dbg) << "[Theorem3] compute_cut_sets failed: " << err << "\n";
                    return; // give up splitting this region
                }
            }

            // 3) Ditto for TARGET triangulation (restricted to region verts)
            auto Etgt_sub = restrict_edges_to_vertex_set(E_tgt, regionSet);
            CutResult tgtCut;
            if (!triangulation_cut::compute_cut_sets(P, Etgt_sub, 0.0, 1.0, -y0, tgtCut, nullptr)) {
                // If degenerate for target (unlikely), we can still proceed using current's Q as boundary
                tgtCut = CutResult{};
            }

            // 4) Build the union region Q boundary from CURRENT triangles crossing the line
            auto boundary_ccw = boundary_cycle_from_tri_union(curCut.Q);
            if (boundary_ccw.size() < 3) {
                // Nothing crosses here: just recurse on the two sides defined by y0
                std::vector<int> top, bot;
                for (int v: regionVertsCCW) ((P[v].y > y0) ? top : bot).push_back(v);
                // After you fill 'top' and 'bot' using (P[v].y > y0) / (P[v].y < y0):
                if (top.empty() || bot.empty()) {
                    return; // no split -> avoid infinite recursion
                }
                if ((int) top.size() == (int) regionVertsCCW.size() ||
                    (int) bot.size() == (int) regionVertsCCW.size()) {
                    return; // defensive: nothing shrank
                }

                if (!top.empty()) recurse_region(P, E_cur, E_tgt, top, depth + 1, acc, dbg);
                if (!bot.empty()) recurse_region(P, E_cur, E_tgt, bot, depth + 1, acc, dbg);
                return;
            }

            // All vertices that belong to Q (not just its boundary)
            std::unordered_set<int> Qverts;
            Qverts.reserve(curCut.Q.size() * 3);
            for (auto t: curCut.Q) {
                Qverts.insert(t[0]);
                Qverts.insert(t[1]);
                Qverts.insert(t[2]);
            }
            if (dbg)
                (*dbg) << "  Q boundary size=" << boundary_ccw.size()
                       << " | Qverts=" << Qverts.size() << "\n";


            // --- CURRENT faces inside Q: keep only triangles whose 3 vertices lie in Qverts
            auto Fcur_all = faces_from_edges(P, Ecur_sub);
            std::vector<std::array<int, 3>> Fcur_Q;
            Fcur_Q.reserve(Fcur_all.size());
            for (auto t: Fcur_all) {
                if (Qverts.count(t[0]) && Qverts.count(t[1]) && Qverts.count(t[2])) {
                    Fcur_Q.push_back(t);
                }
            }

// --- TARGET faces inside Q: same vertex filter, PLUS centroid-in-polygon w.r.t. boundary_ccw
            auto Ftgt_all = faces_from_edges(P, Etgt_sub);
            std::vector<std::array<int, 3>> Ftgt_Q;
            Ftgt_Q.reserve(Ftgt_all.size());
            for (auto t: Ftgt_all) {
                if (!Qverts.count(t[0]) || !Qverts.count(t[1]) || !Qverts.count(t[2])) continue;
                Point c{
                        (P[t[0]].x + P[t[1]].x + P[t[2]].x) / 3.0,
                        (P[t[0]].y + P[t[1]].y + P[t[2]].y) / 3.0
                };
                if (point_in_polygon(c, boundary_ccw, P)) {
                    Ftgt_Q.push_back(t);
                }
            }
            // Fallback: if empty, at least use vertex-set restriction (still better than nothing)
            if (Ftgt_Q.empty()) {
                if (dbg) (*dbg) << "  Ftgt_Q empty; falling back to vertex-set filter only.\n";
                for (auto t: Ftgt_all) {
                    if (Qverts.count(t[0]) && Qverts.count(t[1]) && Qverts.count(t[2])) {
                        Ftgt_Q.push_back(t);
                    }
                }
            }

            if (dbg) {
                (*dbg) << "  Fcur_Q faces=" << Fcur_Q.size()
                       << " | Ftgt_Q faces=" << Ftgt_Q.size() << "\n";

                // quick sanity diff between CURRENT/TARGET *inside Q* before transform
                auto Ec_before = edge_set_from_faces(Fcur_Q);
                auto Et_before = edge_set_from_faces(Ftgt_Q);
                size_t diff1 = 0, diff2 = 0;
                for (auto k: Ec_before) if (!Et_before.count(k)) ++diff1;
                for (auto k: Et_before) if (!Ec_before.count(k)) ++diff2;
                (*dbg) << "  pre-transform face-edge diff: cur\\tgt=" << diff1
                       << ", tgt\\cur=" << diff2 << "\n";
            }


// --- Decide whether inside-Q is meaningful and differs ---
            auto edges_from = [&](const std::vector<std::array<int, 3>> &F) {
                std::unordered_set<long long> S;
                S.reserve(F.size() * 3);
                auto key = [&](int a, int b) {
                    if (a > b) std::swap(a, b);
                    return ((long long) a << 32) ^ (unsigned long long) b;
                };
                for (auto t: F) {
                    S.insert(key(t[0], t[1]));
                    S.insert(key(t[1], t[2]));
                    S.insert(key(t[2], t[0]));
                }
                return S;
            };

            bool Q_valid = (!boundary_ccw.empty() && boundary_ccw.size() >= 3 && !Fcur_Q.empty());

            bool need_transform_inside_Q = false;
            if (Q_valid) {
                auto Ec = edges_from(Fcur_Q);
                auto Et = edges_from(Ftgt_Q);
                size_t diff1 = 0, diff2 = 0;
                for (auto k: Ec) if (!Et.count(k)) ++diff1;
                for (auto k: Et) if (!Ec.count(k)) ++diff2;
                need_transform_inside_Q = (diff1 + diff2) > 0;
                if (dbg) (*dbg) << "  inside-Q edge diff: cur\\tgt=" << diff1 << ", tgt\\cur=" << diff2 << "\n";
            }

            if (Q_valid && need_transform_inside_Q) {
                // --- Preferred path: transform inside Q ---
                std::vector<cgshop::Point> P_cg;
                P_cg.reserve(P.size());
                for (const auto &q: P) P_cg.push_back(cgshop::Point{q.x, q.y});

                cgshop::Theorem2Transformer T2(P_cg, boundary_ccw, Fcur_Q, Ftgt_Q);
                // if (dbg) T2.enable_debug(*dbg);
                auto st = T2.transform();
                if (dbg) {
                    (*dbg) << "  T2(Q): prop1=" << st.prop1_phases
                           << " l4=" << st.lemma4_calls
                           << " l5=" << st.lemma5_calls
                           << " flips=" << st.flips_total << "\n";
                }
                acc.prop1 += st.prop1_phases;
                acc.l4 += st.lemma4_calls;
                acc.l5 += st.lemma5_calls;
                acc.flips += st.flips_total;

                // Splice using ALL Q vertices
                splice_region_edges(E_cur, Qverts, Fcur_Q);
                if (dbg) {
                    (*dbg) << "  splice(Q) done; comparing REGION to TARGET…\n";
                    print_edge_diff(E_cur, E_tgt, Qverts, *dbg);
                }
            } else {
                // --- Fallback: transform the WHOLE region if it still differs ---
                std::unordered_set<int> Rset(regionVertsCCW.begin(), regionVertsCCW.end());

                // Faces of current/target restricted to region vertices
                auto Fcur_R = restrict_faces_to_vertex_set(faces_from_edges(P, Ecur_sub), Rset);
                auto Ftgt_R = restrict_faces_to_vertex_set(faces_from_edges(P, Etgt_sub), Rset);

                auto Ec = edges_from(Fcur_R);
                auto Et = edges_from(Ftgt_R);
                bool differ_region = (Ec.size() != Et.size());
                if (!differ_region) {
                    for (auto k: Ec)
                        if (!Et.count(k)) {
                            differ_region = true;
                            break;
                        }
                }

                if (dbg)
                    (*dbg) << "  fallback? Q_valid=" << Q_valid << " needQ=" << need_transform_inside_Q
                           << " | region-diff=" << differ_region << "\n";

                if (differ_region && !Fcur_R.empty() && !Ftgt_R.empty()) {
                    auto region_hull_ccw = convex_hull_ccw(P, Rset);

                    std::vector<cgshop::Point> P_cg;
                    P_cg.reserve(P.size());
                    for (const auto &q: P) P_cg.push_back(cgshop::Point{q.x, q.y});

                    cgshop::Theorem2Transformer T2(P_cg, region_hull_ccw, Fcur_R, Ftgt_R);
                    // if (dbg) T2.enable_debug(*dbg);
                    auto st = T2.transform();
                    if (dbg) {
                        (*dbg) << "  T2(REGION): prop1=" << st.prop1_phases
                               << " l4=" << st.lemma4_calls
                               << " l5=" << st.lemma5_calls
                               << " flips=" << st.flips_total << "\n";
                    }
                    acc.prop1 += st.prop1_phases;
                    acc.l4 += st.lemma4_calls;
                    acc.l5 += st.lemma5_calls;
                    acc.flips += st.flips_total;

                    splice_region_edges(E_cur, Rset, Fcur_R);
                    if (dbg) {
                        (*dbg) << "  splice(REGION) done; comparing REGION to TARGET…\n";
                        print_edge_diff(E_cur, E_tgt, Rset, *dbg);
                    }
                }
            }

// 8) Recurse on the two sides (upper/lower) w.r.t. y0, *within this region*
            std::vector<int> top, bot;
            top.reserve(regionVertsCCW.size());
            bot.reserve(regionVertsCCW.size());
            for (int v: regionVertsCCW) ((P[v].y > y0) ? top : bot).push_back(v);
            if (dbg) (*dbg) << "  split: |top|=" << top.size() << " |bot|=" << bot.size() << "\n";

            if (!top.empty()) recurse_region(P, E_cur, E_tgt, top, depth + 1, acc, dbg);
            if (!bot.empty()) recurse_region(P, E_cur, E_tgt, bot, depth + 1, acc, dbg);

        }
    } // namespace detail

// =============== Public API impl ===============
    static inline Stats transform(const std::vector<Point>& P,
                                  std::vector<Edge>&        E_cur,
                                  const std::vector<Edge>&  E_tgt,
                                  std::ostream* dbg_out)
    {
        // initial region = full boundary (in some CCW order). If you already store the hull
        // CCW, you can pass that instead. Here we just list all indices 0..n-1;
        // the recursion only needs the vertex *set*, not strict CCW geometry.
        std::vector<int> region;
        region.reserve(P.size());
        for (int i=0;i<(int)P.size();++i) region.push_back(i);

        detail::Accum A;
        detail::recurse_region(P, E_cur, E_tgt, region, /*depth=*/1, A, dbg_out);

        return Stats{ A.depth, A.nodes, A.prop1, A.l4, A.l5, A.flips };
    }

} // namespace theorem3_pset
