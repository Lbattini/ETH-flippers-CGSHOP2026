#pragma once
#include <vector>
#include <array>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <stdexcept>
#include <ostream>
#include <limits>

#include "./TriCommon.hpp"
#include "./Proposition1Inserter.hpp"   // uses Lemma 4 + Lemma 5 and tangents

namespace cgshop {

    class Theorem2Transformer {
    public:
        struct Stats {
            int recursion_nodes = 0;     // number of recursive regions processed
            int prop1_phases    = 0;     // sum over all Proposition1Inserter phases
            int lemma4_calls    = 0;     // sum over all Lemma 4 calls
            int lemma5_calls    = 0;     // sum over all Lemma 5 calls
            int flips_total     = 0;     // total flips spent inside prop.1 calls
            bool finished       = false; // true if all target diagonals present
        };

        Theorem2Transformer(const std::vector<Point>& pts,
                            const std::vector<int>&   boundary_ccw,   // polygon boundary, CCW
                            std::vector<Face>&        faces,          // current triangulation (in/out)
                            const std::vector<Face>&  target_faces)   // target triangulation (read-only)
                : P(pts), B(boundary_ccw), F(faces), FT(target_faces)
        {
            if (B.empty()) throw std::runtime_error("Theorem2Transformer: boundary must not be empty.");
            build_boundary_index();
            collect_target_diagonals_global();
        }

        void enable_debug(std::ostream& os) { dbg_out_ = &os; }

        Stats transform() {
            Stats S;
            // recurse on the whole polygon
            std::vector<int> subB = B; // CCW order
            recurse_region(subB, S);
            // Done if we inserted all target diagonals
            S.finished = all_target_diagonals_present();
            return S;
        }

    private:
        const std::vector<Point>& P;
        const std::vector<int>&   B;     // CCW boundary (global)
        std::vector<Face>&        F;     // current triangulation (global, in/out)
        const std::vector<Face>&  FT;    // target triangulation (global, read-only)

        std::unordered_map<int,int> bIndex;  // vertex id -> index in global boundary
        std::unordered_set<EdgeKey, EdgeKeyHash> boundary_edges_global; // for convenience

        // global target diagonal set (unordered, a<b normalized)
        std::unordered_set<EdgeKey, EdgeKeyHash> target_diagonals;

        // debug
        std::ostream* dbg_out_ = nullptr;
        void dbg(const std::string& s) const { if (dbg_out_) (*dbg_out_) << s << '\n'; }

        // ---------- boundary helpers ----------
        void build_boundary_index() {
            bIndex.clear();
            for (int i=0;i<(int)B.size();++i) bIndex[B[i]] = i;

            boundary_edges_global.clear();
            for (int i=0;i<(int)B.size();++i) {
                int a = B[i], b = B[(i+1)%B.size()];
                boundary_edges_global.insert(EdgeKey(a,b));
            }
        }

        bool is_boundary_edge_global(int a, int b) const {
            return boundary_edges_global.count(EdgeKey(a,b)) > 0;
        }

        // Collect all non-boundary edges from FT (target) as target diagonals
        void collect_target_diagonals_global() {
            target_diagonals.clear();
            for (const auto& T : FT) {
                int x=T[0], y=T[1], z=T[2];
                EdgeKey e1(x,y), e2(y,z), e3(z,x);
                if (!is_boundary_edge_global(e1.a, e1.b)) target_diagonals.insert(e1);
                if (!is_boundary_edge_global(e2.a, e2.b)) target_diagonals.insert(e2);
                if (!is_boundary_edge_global(e3.a, e3.b)) target_diagonals.insert(e3);
            }
        }

        // Are all target diagonals present in current triangulation?
        bool all_target_diagonals_present() const {
            // Build current edge set quickly (only to check membership)
            // We only care that each target diagonal exists; boundary edges ignored.
            std::unordered_set<EdgeKey, EdgeKeyHash> cur;
            for (const auto& t : F) {
                cur.insert(EdgeKey(t[0],t[1]));
                cur.insert(EdgeKey(t[1],t[2]));
                cur.insert(EdgeKey(t[2],t[0]));
            }
            for (const auto& e : target_diagonals) {
                if (cur.count(e)==0) return false;
            }
            return true;
        }

        // ---------- subproblem helpers ----------
        // Return indices of a and b along subB (a,b must be in subB)
        static std::pair<int,int> idx_in_subB(const std::vector<int>& subB, int a, int b) {
            int ia=-1, ib=-1;
            for (int i=0;i<(int)subB.size();++i) {
                if (subB[i]==a) ia=i;
                if (subB[i]==b) ib=i;
            }
            return {ia,ib};
        }

        // Return the two CCW chains from a to b along subB (both inclusive), as (chain1, chain2).
        static std::pair<std::vector<int>, std::vector<int>>
        split_chains_on_subB(const std::vector<int>& subB, int a, int b) {
            auto [ia,ib] = idx_in_subB(subB,a,b);
            if (ia<0 || ib<0) return {{},{}};
            int n = (int)subB.size();

            std::vector<int> c1; c1.push_back(subB[ia]);
            int i = ia;
            while (i != ib) { i = (i+1)%n; c1.push_back(subB[i]); if (c1.size()>n+2) break; }

            std::vector<int> c2; c2.push_back(subB[ib]);
            int j = ib;
            while (j != ia) { j = (j+1)%n; c2.push_back(subB[j]); if (c2.size()>n+2) break; }
            // normalize to start at a and end at b
            std::reverse(c2.begin(), c2.end()); // now c2 also goes a..b
            return {c1,c2};
        }

        // Among all target diagonals fully inside subB, pick a "long" one:
        // minimize max(len(chain1), len(chain2)) along subB.
        EdgeKey pick_long_diagonal_in_subB(const std::vector<int>& subB) const {
            EdgeKey best(-1,-1);
            int best_balance = std::numeric_limits<int>::max();

            // Build a quick membership for subB
            std::unordered_set<int> inSub(subB.begin(), subB.end());

            for (const auto& e : target_diagonals) {
                int a=e.a, b=e.b;
                if (!inSub.count(a) || !inSub.count(b)) continue;

                auto [ia, ib] = idx_in_subB(subB, a, b);
                if (ia<0 || ib<0) continue;

                int n = (int)subB.size();
                int len1 = (ib - ia + n) % n + 1; // a..b inclusive along one side
                int len2 = n - len1 + 2;          // the other side also inclusive
                int bal  = std::max(len1, len2);
                if (bal < best_balance) {
                    best_balance = bal;
                    best = e;
                }
            }
            return best; // possibly (-1,-1) if none found
        }

        // Restrict target diagonal set to a subboundary (both endpoints in subB)
        std::unordered_set<EdgeKey, EdgeKeyHash>
        target_diagonals_in_subB(const std::vector<int>& subB) const {
            std::unordered_set<int> inSub(subB.begin(), subB.end());
            std::unordered_set<EdgeKey, EdgeKeyHash> S;
            for (const auto& e : target_diagonals) {
                if (inSub.count(e.a) && inSub.count(e.b)) S.insert(e);
            }
            return S;
        }

        // ---------- recursion ----------
        void recurse_region(const std::vector<int>& subB, Stats& S) {
            ++S.recursion_nodes;

            // Compute set of target diagonals local to subB
            auto TD = target_diagonals_in_subB(subB);
            if (TD.empty()) return; // nothing to do in this region

            // If all local target diagonals are already present, done.
            if (subregion_done(subB, TD)) return;

            // Pick a balanced diagonal in this region
            EdgeKey e = pick_long_diagonal_in_subB(subB);
            if (e.a < 0) {
                // Fallback: pick any remaining target diagonal in subB
                e = *TD.begin();
            }

            const int u = e.a, v = e.b;

            if (dbg_out_) {
                dbg("Theorem2: inserting diagonal ("+std::to_string(u)+","+std::to_string(v)+") "
                                                                                             "in region of size " + std::to_string(subB.size()));
            }

            // Use Proposition1Inserter on THIS region:
            // We call it with the full face set F (global) but with sub-boundary 'subB' so
            // the controller computes chains and tangents relative to the region.
            {
                Proposition1Inserter ins(P, subB, F);
                if (dbg_out_) ins.enable_debug(*dbg_out_);
                auto st = ins.insert(u, v);
                S.prop1_phases += st.phases;
                S.lemma4_calls += st.lemma4_calls;
                S.lemma5_calls += st.lemma5_calls;
                S.flips_total  += st.flips_total;
            }

            // If insertion somehow failed, stop descending this branch
            // (should not happen on well-formed instances).
            // Try to continue anyway on the assumption some progress was made.
            // But if (u,v) is still absent, we cannot split the region reliably.
            if (!edge_present(u,v)) return;

            // Split subB into two subregions along (u,v)
            auto [c1, c2] = split_chains_on_subB(subB, u, v);

            // Build sub-boundaries that include the diagonal (u,v) as a boundary edge
            auto make_sub = [&](const std::vector<int>& chain)->std::vector<int> {
                // chain is u..v inclusive; add back v..u as a single edge
                std::vector<int> out = chain;
                // ensure uniqueness of consecutive vertices
                if (!out.empty()) {
                    std::vector<int> uniq; uniq.reserve(out.size());
                    int last=-1;
                    for (int w : out) { if (uniq.empty() || w!=last) uniq.push_back(w); last=w; }
                    out.swap(uniq);
                }
                return out;
            };
            std::vector<int> sub1 = make_sub(c1);
            std::vector<int> sub2 = make_sub(c2);

            // Recurse into both sides
            recurse_region(sub1, S);
            recurse_region(sub2, S);
        }

        // Check if all target diagonals inside subB are present in F
        bool subregion_done(const std::vector<int>& subB,
                            const std::unordered_set<EdgeKey, EdgeKeyHash>& TD) const
        {
            // Build current edges present restricted to subB
            std::unordered_set<EdgeKey, EdgeKeyHash> cur;
            for (const auto& t : F) {
                EdgeKey e1(t[0],t[1]), e2(t[1],t[2]), e3(t[2],t[0]);
                cur.insert(e1); cur.insert(e2); cur.insert(e3);
            }
            for (const auto& e : TD) {
                if (cur.count(e)==0) return false;
            }
            return true;
        }

        // Is edge (a,b) currently present in F?
        bool edge_present(int a, int b) const {
            EdgeKey e(a,b);
            for (const auto& t : F) {
                if (EdgeKey(t[0],t[1])==e) return true;
                if (EdgeKey(t[1],t[2])==e) return true;
                if (EdgeKey(t[2],t[0])==e) return true;
            }
            return false;
        }
    };

} // namespace cgshop
