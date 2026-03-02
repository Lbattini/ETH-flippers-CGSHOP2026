#pragma once
#include <vector>
#include <limits>
#include <cmath>
#include <cassert>
#include "./TriCommon.hpp"   // cgshop::Point, Face, EdgeKey, ...



namespace cgshop {

    // (near the top, inside namespace cgshop)
    namespace detail_tangent {
        inline double orient2d_cc(const Point& a, const Point& b, const Point& c) {
            long double x1 = (long double)b.x - a.x, y1 = (long double)b.y - a.y;
            long double x2 = (long double)c.x - a.x, y2 = (long double)c.y - a.y;
            long double v = x1*y2 - y1*x2;
            return (double)v;
        }
    }


    /*
    Tangent to a *convex* chain (in CCW order) from an external point S.

    Definitions:
      - RIGHT tangent at chain[j]: both neighbors (if they exist) lie on or to the RIGHT
        of the line S->chain[j], i.e. orient(S, chain[j], neighbor) <= +eps.
      - LEFT tangent  at chain[j]: both neighbors lie on or to the LEFT
        of the line S->chain[j], i.e. orient(S, chain[j], neighbor) >= -eps.

    Notes:
      - The chain is OPEN (not wrapped). Endpoints have only one neighbor — we check the one they have.
      - The chain must be convex in CCW order. (Upper or lower hull segment of your boundary fits.)
      - Complexity: O(|chain|). For your N ~ 1e4 this is fine and very robust.

    Returns: index in 'chain' (0..chain.size()-1) of a valid tangent vertex, or -1 if none.
*/
    inline int tangent_to_convex_chain(const std::vector<Point>& P,
                                       int s,
                                       const std::vector<int>& chain_ccw,
                                       bool right_tangent,
                                       double eps = 1e-12)
    {
        const int m = (int)chain_ccw.size();
        if (m == 0) return -1;
        if (m == 1) return 0;

        auto ok_at = [&](int j)->bool {
            const int v = chain_ccw[j];
            // check previous neighbor (if any)
            if (j-1 >= 0) {
                const int pv = chain_ccw[j-1];
                double o = detail_tangent::orient2d_cc(P[s], P[v], P[pv]);
                if (right_tangent) { if (o > +eps) return false; } // pv is to the LEFT -> violates right-tangent
                else               { if (o < -eps) return false; } // pv is to the RIGHT -> violates left-tangent
            }
            // check next neighbor (if any)
            if (j+1 < m) {
                const int nv = chain_ccw[j+1];
                double o = detail_tangent::orient2d_cc(P[s], P[v], P[nv]);
                if (right_tangent) { if (o > +eps) return false; }
                else               { if (o < -eps) return false; }
            }
            return true;
        };

        // First pass: collect all valid tangent candidates
        int best = -1;
        double best_score = -std::numeric_limits<double>::infinity();

        for (int j = 0; j < m; ++j) {
            if (!ok_at(j)) continue;

            // Prefer the "most supporting" one: for right tangent minimize max(orient to neighbors),
            // for left tangent maximize min(orient to neighbors). Here we use a simple proxy score:
            // project along +x direction for stability on grid-like inputs.
            const int v = chain_ccw[j];
            double score = (double)P[v].x - (double)P[s].x; // rightward progress
            if (!right_tangent) score = -score;             // flip for left tangent preference

            if (score > best_score) {
                best_score = score;
                best = j;
            }
        }

        // If none satisfied strict checks (rare due to eps), fall back to endpoints as safe supports
        if (best == -1) {
            // Endpoints are always valid “outer” supports for an open convex chain viewed from outside.
            // Pick the one that matches the requested side (right/left) by comparing orientation to inside neighbor.
            auto endpoint_ok = [&](int j)->bool { return ok_at(j); };
            if (endpoint_ok(0)) best = 0;
            if (endpoint_ok(m-1)) {
                if (best == -1) best = m-1;
                else {
                    // choose by side preference
                    double o0 = detail_tangent::orient2d_cc(P[s], P[chain_ccw[best]], P[chain_ccw[best + (best+1<m ? 1 : -1)]]);
                    double o1 = detail_tangent::orient2d_cc(P[s], P[chain_ccw[m-1]],    P[chain_ccw[m-2]]);
                    // prefer the one with "more" support (more negative for right tangent, more positive for left)
                    if (right_tangent) { if (o1 < o0) best = m-1; }
                    else               { if (o1 > o0) best = m-1; }
                }
            }
        }

        return best;
    }

// Convenience wrappers:
    inline int right_tangent_to_convex_chain(const std::vector<Point>& P,
                                             int s,
                                             const std::vector<int>& chain_ccw,
                                             double eps = 1e-12)
    {
        return tangent_to_convex_chain(P, s, chain_ccw, /*right_tangent=*/true, eps);
    }

    inline int left_tangent_to_convex_chain(const std::vector<Point>& P,
                                            int s,
                                            const std::vector<int>& chain_ccw,
                                            double eps = 1e-12)
    {
        return tangent_to_convex_chain(P, s, chain_ccw, /*right_tangent=*/false, eps);
    }

} // namespace cgshop
