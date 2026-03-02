#pragma once
#include <vector>
#include <array>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <stdexcept>
#include <limits>
#include <ostream>
#include <cmath>

#include "./TriCommon.hpp"
#include "./ConvexChainTangent.hpp"
#include "./DiagonalIntroducerLemma4.hpp"
#include "./Lemma5Transformer.hpp"

namespace cgshop {

    class Proposition1Inserter {
    public:
        struct Stats {
            int phases = 0;            // number of outer phases
            int lemma4_calls = 0;
            int lemma5_calls = 0;
            int flips_total = 0;       // total flips consumed across all subcalls
            bool inserted = false;     // did we insert the target (u,v)?
        };

        Proposition1Inserter(const std::vector<Point>& pts,
                             const std::vector<int>& boundary_ccw,
                             std::vector<Face>& faces /* in/out */)
                : P(pts), B(boundary_ccw), F(faces)
        {
            if (B.empty()) throw std::runtime_error("Proposition1Inserter: boundary must be non-empty.");
            buildIncidence();
        }

        void enable_debug(std::ostream& os) { dbg_out_ = &os; }

        // Insert target diagonal (u,v) using Proposition-1 style phases.
        Stats insert(int u, int v) {
            Stats S;
            if (u == v) return S;
            if (hasEdge(u,v)) { S.inserted = true; return S; }

            // Cheap attempt: Lemma 4 directly
            {
                DiagonalIntroducerLemma4 L4(P, F);
                if (dbg_out_) L4.enable_debug(*dbg_out_);
                auto st = L4.introduce(u,v);
                S.lemma4_calls += 1;
                S.flips_total  += st.flips_total;
                rebuildIncidence();
                if (hasEdge(u,v)) { S.inserted = true; return S; }
            }

            // Phase loop: reduce active polygon until (u,v) is insertable
            const int MAX_PHASES = std::max<size_t>(8, B.size()); // generous guard
            for (int phase = 1; phase <= MAX_PHASES; ++phase) {
                ++S.phases;
                if (dbg_out_) dbg("== Proposition1 Phase " + std::to_string(phase) + " ==");

                // Build active polygon Pe and its chains from u to v
                ActivePolygon Pe = build_active_polygon(u, v);
                if (Pe.crossing_faces.empty()) {
                    // Nothing crosses; try once more to insert directly
                    if (try_direct_uv(u,v, S)) return S;
                    break;
                }

                // Split Pe boundary at u and v into two CCW chains (u->v along each side)
                Chains CV = chains_between_u_v(Pe, u, v);
                if (CV.upper.empty() || CV.lower.empty()) {
                    // Degenerate Pe — try direct insertion and bail if not possible
                    if (try_direct_uv(u,v, S)) return S;
                    break;
                }

                // Identify convex indices on the lower chain (in Pe boundary order)
                std::vector<int> lower_conv_idx = convex_vertices_on_chain(Pe, CV.lower);
                // Sample every second convex (per the proposition)
                std::vector<int> sample_lower;
                for (size_t i = 0; i < lower_conv_idx.size(); i += 2) sample_lower.push_back(lower_conv_idx[i]);

                if (dbg_out_) {
                    dbg("  Pe: |faces|=" + std::to_string(Pe.crossing_faces.size())
                        + ", upper=" + std::to_string(CV.upper.size())
                        + ", lower=" + std::to_string(CV.lower.size())
                        + ", sampled lower convex = " + std::to_string(sample_lower.size()));
                }

                // For each sampled y (on lower), compute tangent to upper and try to introduce it
                int successes = 0;
                for (int id_in_chain : sample_lower) {
                    int y = CV.lower[id_in_chain];
                    // Compute **right** tangent from y to upper chain (CCW order)
                    int j = right_tangent_to_convex_chain(P, y, CV.upper);
                    if (j < 0) continue;
                    int a = y, b = CV.upper[j];
                    if (hasEdge(a,b)) continue;

                    // Try Lemma 4 on (a,b)
                    {
                        DiagonalIntroducerLemma4 L4(P, F);
                        if (dbg_out_) L4.enable_debug(*dbg_out_);
                        auto st = L4.introduce(a,b);
                        S.lemma4_calls += 1;
                        S.flips_total  += st.flips_total;
                        rebuildIncidence();
                        if (st.inserted) { ++successes; continue; }
                    }

                    // If that failed (rare), nudge using Lemma 5 transformer between a-y and a-b / b-y
                    // Build an auxiliary UV target across the “L4-stuck” pocket
                    {
                        Lemma5Transformer X(P, B, F);
                        if (dbg_out_) X.enable_debug(*dbg_out_);
                        auto st5 = X.transform_from_xy_to_uv(a, b, y, u /*any convex corner works; use u*/);
                        S.lemma5_calls += 1;
                        S.flips_total  += st5.flips;
                        rebuildIncidence();
                        if (hasEdge(a,b)) { ++successes; continue; }
                    }
                }

                // After tangent insertions, try (u,v) again (Lemma 4)
                if (try_direct_uv(u, v, S)) return S;

                // If no progress, try sampling the **upper** chain symmetrically
                if (successes == 0) {
                    std::vector<int> upper_conv_idx = convex_vertices_on_chain(Pe, CV.upper);
                    std::vector<int> sample_upper;
                    for (size_t i = 0; i < upper_conv_idx.size(); i += 2) sample_upper.push_back(upper_conv_idx[i]);

                    for (int id_in_chain : sample_upper) {
                        int y = CV.upper[id_in_chain];
                        int j = right_tangent_to_convex_chain(P, y, CV.lower);
                        if (j < 0) continue;
                        int a = y, b = CV.lower[j];
                        if (hasEdge(a,b)) continue;

                        // Lemma 4
                        {
                            DiagonalIntroducerLemma4 L4(P, F);
                            if (dbg_out_) L4.enable_debug(*dbg_out_);
                            auto st = L4.introduce(a,b);
                            S.lemma4_calls += 1;
                            S.flips_total  += st.flips_total;
                            rebuildIncidence();
                            if (st.inserted) ++successes;
                        }
                        if (hasEdge(u,v)) { S.inserted = true; return S; }
                    }

                    // Try (u,v) once more
                    if (try_direct_uv(u, v, S)) return S;
                }

                // If still nothing changed, we bail to avoid infinite loop
                if (successes == 0) {
                    if (dbg_out_) dbg("  No successful tangents this phase; stopping.");
                    break;
                }
            }

            S.inserted = hasEdge(u,v);
            return S;
        }

    private:
        const std::vector<Point>& P;
        const std::vector<int>& B; // boundary in CCW order
        std::vector<Face>& F;

        std::unordered_map<EdgeKey, std::vector<int>, EdgeKeyHash> edge2faces;

        std::ostream* dbg_out_ = nullptr;
        void dbg(const std::string& s) const { if (dbg_out_) (*dbg_out_) << s << '\n'; }

        // -------- incidence --------
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
        bool hasEdge(int a, int b) const {
            return edge2faces.find(EdgeKey(a,b)) != edge2faces.end();
        }

        // Try a direct Lemma-4 insertion of (u,v)
        bool try_direct_uv(int u, int v, Stats& S) {
            if (hasEdge(u,v)) { S.inserted = true; return true; }
            DiagonalIntroducerLemma4 L4(P, F);
            if (dbg_out_) L4.enable_debug(*dbg_out_);
            auto st = L4.introduce(u,v);
            S.lemma4_calls += 1;
            S.flips_total  += st.flips_total;
            rebuildIncidence();
            S.inserted = hasEdge(u,v);
            return S.inserted;
        }

        // ---------- active polygon ----------
        struct ActivePolygon {
            std::vector<int> crossing_faces;       // indices of faces whose *edge* crosses (u,v) strictly
            std::vector<int> boundary_cycle;       // boundary vertices of the active region (CCW), u..v..u
            // We’ll build the chains later when u and v are known
        };

        // strict segment intersection
        static bool seg_x_seg(const Point& A, const Point& B, const Point& C, const Point& D) {
            auto sgn = [](double v){ return (v>0)-(v<0); };
            auto orient = [](const Point& a, const Point& b, const Point& c){
                long double x1 = (long double)b.x-a.x, y1 = (long double)b.y-a.y;
                long double x2 = (long double)c.x-a.x, y2 = (long double)c.y-a.y;
                long double v = x1*y2 - y1*x2;
                return (double)v;
            };
            double o1 = orient(A,B,C), o2 = orient(A,B,D);
            double o3 = orient(C,D,A), o4 = orient(C,D,B);
            return (sgn(o1)*sgn(o2) < 0) && (sgn(o3)*sgn(o4) < 0);
        }

        ActivePolygon build_active_polygon(int u, int v) const {
            ActivePolygon R;
            // mark faces that cross (u,v)
            for (int f=0; f<(int)F.size(); ++f) {
                const auto& T = F[f];
                int a=T[0], b=T[1], c=T[2];
                if (seg_x_seg(P[a],P[b], P[u],P[v]) ||
                    seg_x_seg(P[b],P[c], P[u],P[v]) ||
                    seg_x_seg(P[c],P[a], P[u],P[v]))
                    R.crossing_faces.push_back(f);
            }
            // A simple boundary cycle approximation: walk all edges of crossing faces that are boundary
            // of the union (appear exactly once among crossing faces). We then sort them into a CCW cycle.
            std::unordered_map< EdgeKey, int, EdgeKeyHash > multiplicity;
            for (int f : R.crossing_faces) {
                const auto& T = F[f];
                EdgeKey e1(T[0],T[1]), e2(T[1],T[2]), e3(T[2],T[0]);
                multiplicity[e1]++; multiplicity[e2]++; multiplicity[e3]++;
            }
            // collect boundary edges of active region (only 1 occurrence within crossing set)
            std::vector<std::pair<int,int>> boundary_edges;
            boundary_edges.reserve(multiplicity.size());
            for (const auto& kv : multiplicity) {
                if (kv.second == 1) {
                    boundary_edges.emplace_back(kv.first.a, kv.first.b);
                }
            }
            if (boundary_edges.empty()) return R;

            // Stitch boundary edges into a CCW cycle of vertices (naive O(E^2) is fine here)
            // Build adjacency
            std::unordered_multimap<int,int> adj;
            for (auto [a,b] : boundary_edges) {
                adj.emplace(a,b);
                adj.emplace(b,a);
            }
            // pick a start (prefer u if present)
            int start = boundary_edges.front().first;
            for (auto [a,b] : boundary_edges) { if (a==u || b==u) { start = u; break; } }
            // walk
            std::vector<int> cyc; cyc.push_back(start);
            int prev = -1, cur = start;
            for (size_t steps=0; steps < boundary_edges.size()*2+5; ++steps) {
                // pick a neighbor that's not 'prev' and that keeps us CCW-ish; we settle for any
                auto range = adj.equal_range(cur);
                int nxt = -1;
                for (auto it = range.first; it != range.second; ++it) {
                    if (it->second != prev) { nxt = it->second; break; }
                }
                if (nxt < 0) break;
                cyc.push_back(nxt);
                if (nxt == start) break;
                prev = cur; cur = nxt;
            }
            R.boundary_cycle = std::move(cyc);
            return R;
        }

        struct Chains {
            std::vector<int> upper; // CCW from u to v along one side
            std::vector<int> lower; // CCW from u to v along the other side
        };

        Chains chains_between_u_v(const ActivePolygon& Pe, int u, int v) const {
            Chains C;
            if (Pe.boundary_cycle.size() < 3) return C;
            // find all occurrences of u and v on the cycle
            auto idxs = [&](int w){
                std::vector<int> res;
                for (int i=0;i<(int)Pe.boundary_cycle.size();++i)
                    if (Pe.boundary_cycle[i]==w) res.push_back(i);
                return res;
            };
            auto iu = idxs(u), iv = idxs(v);
            if (iu.empty() || iv.empty()) return C;

            // choose nearest pair along cycle
            int a = iu.front(), b = iv.front();
            // walk a->b CCW
            std::vector<int> path1;
            int i=a; path1.push_back(Pe.boundary_cycle[i]);
            while (i != b) { i = (i+1) % Pe.boundary_cycle.size(); path1.push_back(Pe.boundary_cycle[i]); if (path1.size()>Pe.boundary_cycle.size()+2) break; }
            // walk b->a CCW (which is other side from u to v)
            std::vector<int> path2;
            int j=b; path2.push_back(Pe.boundary_cycle[j]);
            while (j != a) { j = (j+1) % Pe.boundary_cycle.size(); path2.push_back(Pe.boundary_cycle[j]); if (path2.size()>Pe.boundary_cycle.size()+2) break; }

            // Decide which one is "upper". We don't need true geometric upper/lower — just two disjoint chains from u to v.
            C.upper = std::move(path1);
            C.lower = std::move(path2);

            // ensure both start at u and end at v (strip repeated endpoints)
            auto normalize = [&](std::vector<int>& ch){
                if (!ch.empty() && ch.front()!=u) return;
                if (!ch.empty() && ch.back()!=v)  return;
                // Unique consecutive
                std::vector<int> out; out.reserve(ch.size());
                int last = -1;
                for (int w: ch) {
                    if (out.empty() || w!=last) out.push_back(w);
                    last = w;
                }
                ch.swap(out);
            };
            normalize(C.upper);
            normalize(C.lower);
            return C;
        }

        // Convex vertices on a chain (by signed area test on consecutive triples in that chain order)
        std::vector<int> convex_vertices_on_chain(const ActivePolygon& Pe, const std::vector<int>& chain) const {
            std::vector<int> idx;
            if (chain.size() < 3) return idx;
            for (size_t k=1; k+1<chain.size(); ++k) {
                int a = chain[k-1], b = chain[k], c = chain[k+1];
                long double x1 = (long double)P[b].x - P[a].x, y1 = (long double)P[b].y - P[a].y;
                long double x2 = (long double)P[c].x - P[b].x, y2 = (long double)P[c].y - P[b].y;
                long double v = x1*y2 - y1*x2; // orient(a,b,c)
                if (v > -1e-12) idx.push_back((int)k); // weakly convex
            }
            return idx;
        }
    };

} // namespace cgshop
