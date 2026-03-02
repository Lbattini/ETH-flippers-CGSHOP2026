#pragma once
#include <bits/stdc++.h>

struct Point { double x{}, y{}; };
struct Edge  { int u, v; };

struct CutResult {
    std::vector<int> Ppos, Pneg;              // indices of points on positive/negative side
    std::vector<std::array<int,3>> Q;         // triangles intersecting the line
};

namespace triangulation_cut
{
    static inline int sgn(double val, double eps=1e-12) {
        if (val >  eps) return  1;
        if (val < -eps) return -1;
        return 0;
    }

    // ---------------- Half-edge face extraction ----------------
    struct HalfEdge { int from, to; };
    struct FacesOut { std::vector<std::array<int,3>> triangles; };

    static FacesOut extractFacesAsTriangles(
            const std::vector<Point>& P,
            const std::vector<Edge>& undirectedEdges
    ){
        const int n = (int)P.size();

        std::vector<HalfEdge> he; he.reserve(2*undirectedEdges.size());
        std::unordered_map<long long,int> idOf; idOf.reserve(2*undirectedEdges.size()*2);
        auto key = [](int a,int b)->long long { return ( (long long)a<<32 ) ^ (unsigned long long)b; };

        auto add_half = [&](int u,int v){
            int id = (int)he.size();
            he.push_back({u,v});
            idOf[key(u,v)] = id;
        };

        std::vector<std::vector<int>> out(n);
        for (auto e: undirectedEdges) {
            if (e.u==e.v) continue;
            add_half(e.u, e.v);
            add_half(e.v, e.u);
        }
        for (int id=0; id<(int)he.size(); ++id) out[he[id].from].push_back(id);

        auto angleAt = [&](int from, int to){
            return atan2(P[to].y - P[from].y, P[to].x - P[from].x);
        };
        for (int u=0; u<n; ++u) {
            auto &vec = out[u];
            sort(vec.begin(), vec.end(), [&](int a, int b){
                double aa = angleAt(u, he[a].to);
                double bb = angleAt(u, he[b].to);
                if (aa == bb) return he[a].to < he[b].to;
                return aa < bb; // CCW
            });
        }

        // neighbor->position lookup
        std::vector<std::unordered_map<int,int>> pos(n);
        for (int u=0; u<n; ++u) {
            auto &vec = out[u];
            auto &mp  = pos[u];
            mp.reserve(vec.size()*2+1);
            for (int i=0;i<(int)vec.size();++i) mp[he[vec[i]].to] = i;
        }

        std::vector<char> used(he.size(), 0);
        std::vector<std::vector<int>> faces;

        auto next_left = [&](int hid)->int {
            int u = he[hid].from;
            int v = he[hid].to;
            auto it = pos[v].find(u);
            if (it == pos[v].end()) return -1;
            int idx = it->second;
            const auto &adj = out[v];
            int prev_idx = (idx - 1 + (int)adj.size()) % (int)adj.size();
            return adj[prev_idx];
        };

        for (int h=0; h<(int)he.size(); ++h) {
            if (used[h]) continue;
            int start = h, cur = h;
            std::vector<int> poly; poly.reserve(8);
            while (true) {
                used[cur] = 1;
                int u = he[cur].from;
                poly.push_back(u);
                int nxt = next_left(cur);
                if (nxt<0) { poly.clear(); break; }
                cur = nxt;
                if (cur == start) break;
            }
            if (poly.size() >= 3) faces.push_back(poly);
        }

        auto signedArea2 = [&](const std::vector<int>& poly)->double {
            long double s = 0;
            int m = (int)poly.size();
            for (int i=0;i<m;++i) {
                const Point& A = P[ poly[i] ];
                const Point& B = P[ poly[(i+1)%m] ];
                s += (long double)A.x*B.y - (long double)A.y*B.x;
            }
            return (double)s;
        };

        // find outer face
        int outer_id = -1; double outer_area2 = 0;
        for (int i=0;i<(int)faces.size();++i) {
            double a2 = signedArea2(faces[i]);
            if (i==0 || a2 < outer_area2) { outer_area2 = a2; outer_id = i; }
        }

        FacesOut outF;
        for (int i=0;i<(int)faces.size();++i) {
            if (i == outer_id) continue;
            if (faces[i].size() == 3) {
                std::array<int,3> tri = { faces[i][0], faces[i][1], faces[i][2] };
                int a=tri[0], b=tri[1], c=tri[2];
                if (b<a && b<c) tri = {b,c,a};
                else if (c<a && c<b) tri = {c,a,b};
                const Point &A=P[tri[0]], &B=P[tri[1]], &C=P[tri[2]];
                double cross = (B.x-A.x)*(C.y-A.y) - (B.y-A.y)*(C.x-A.x);
                if (cross < 0) std::swap(tri[1], tri[2]);
                outF.triangles.push_back(tri);
            }
        }
        sort(outF.triangles.begin(), outF.triangles.end());
        outF.triangles.erase(unique(outF.triangles.begin(), outF.triangles.end()), outF.triangles.end());
        return outF;
    }

    // ---------------- Main API: returns false if any vertex lies on the line ----------------
    static bool compute_cut_sets(
            const std::vector<Point>& P,
            const std::vector<Edge>& E,
            double a, double b, double c,
            CutResult& out,
            std::string* error_msg = nullptr,
            double eps = 1e-12
    ){
        const int n = (int)P.size();
        std::vector<int> sign(n,0);

        // 1) classify points; fail-fast if any lies exactly on the line
        for (int i=0;i<n;++i) {
            double v = a*P[i].x + b*P[i].y + c;
            int s = sgn(v, eps);
            if (s == 0) {
                if (error_msg)
                    *error_msg = "Degenerate cut: vertex " + std::to_string(i) + " lies on the line.";
                return false;
            }
            sign[i] = s;
            if (s > 0) out.Ppos.push_back(i);
            else       out.Pneg.push_back(i);
        }

        // 2) extract triangular faces
        FacesOut F = extractFacesAsTriangles(P, E);

        // 3) collect intersecting triangles
        auto tri_intersects = [&](const std::array<int,3>& T)->bool{
            int s0 = sign[T[0]], s1 = sign[T[1]], s2 = sign[T[2]];
            return ((s0>0||s1>0||s2>0) && (s0<0||s1<0||s2<0));
        };

        for (auto &tri : F.triangles)
            if (tri_intersects(tri))
                out.Q.push_back(tri);

        return true;
    }
} // namespace triangulation_cut
