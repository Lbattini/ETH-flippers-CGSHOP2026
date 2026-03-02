#include <bits/stdc++.h>

using namespace std;

struct Point{
    double x;
    double y;
    int ind;
};


struct Quadrilateral{
    Point a, b, c, d;
};

struct Edge{
    int u;
    int v;
};

struct Triangulation{
    int n;
    vector<Point> pts;
    vector<Edge> edges;
};

struct combQuadrilateral{
    int a;
    int b;
    int c;
    int d;
    int ind;

    combQuadrilateral(int a1, int b1, int c1, int d1, int ind) : ind(ind) {
        vector<int> t = {a1,b1,c1,d1};
        sort(t.begin(), t.end());

        a = t[0]; b = t[1]; c = t[2]; d = t[3];
    }

    bool operator<(const combQuadrilateral& other) const{

        if (a != other.a)
            return a < other.a;
        if (b != other.b)
            return b < other.b;
        if (c != other.c)
            return c < other.c;

        return d < other.d;
    }
};


bool operator!=(const Quadrilateral& a, const combQuadrilateral& b){
    if (a.a.ind != b.a) return 1;
    if (a.b.ind != b.b) return 1;
    if (a.c.ind != b.c) return 1;
    return (a.d.ind != b.d);
}

double orient(const Point& a, const Point& b, const Point& c) {
    return (b.x - a.x) * (c.y - a.y) -
           (b.y - a.y) * (c.x - a.x);
}

double cross(const Point& a, const Point& b, const Point& c) {
    return orient(a, b, c);
}

bool isConvexQuadrilateral(std::vector<Point> p) {
    if (p.size() != 4) return false;

    // ---- Compute convex hull (Graham scan for 4 points) ----
    // Find lowest-leftmost point
    int start = 0;
    for (int i = 1; i < 4; i++) {
        if (p[i].y < p[start].y ||
            (p[i].y == p[start].y && p[i].x < p[start].x))
            start = i;
    }
    std::swap(p[0], p[start]);
    Point base = p[0];

    // Sort by polar angle around base
    std::sort(p.begin() + 1, p.end(), [&](const Point& a, const Point& b) {
        double o = orient(base, a, b);
        if (o == 0) {
            // closer one first
            double da = (a.x - base.x)*(a.x-base.x) + (a.y-base.y)*(a.y-base.y);
            double db = (b.x - base.x)*(b.x-base.x) + (b.y-base.y)*(b.y-base.y);
            return da < db;
        }
        return o > 0; // CCW first
    });

    // Build hull
    std::vector<Point> h;
    for (auto &pt : p) {
        while (h.size() >= 2 && orient(h[h.size()-2], h.back(), pt) <= 0)
            h.pop_back();
        h.push_back(pt);
    }

    // Hull must contain all 4 points to form a convex quadrilateral
    if (h.size() != 4) return false;

    // ---- Check all turns have same orientation sign ----
    bool all_pos = true, all_neg = true;

    for (int i = 0; i < 4; i++) {
        double o = orient(h[i], h[(i+1)%4], h[(i+2)%4]);
        if (o <= 0) all_pos = false;
        if (o >= 0) all_neg = false;
    }

    return all_pos || all_neg;
}

double dist2(const Point& a, const Point& b) {
    return (a.x - b.x)*(a.x - b.x) +
           (a.y - b.y)*(a.y - b.y);
}

// -----------------------------------------------------
// 1) Sort 4 points into cyclic convex order
// -----------------------------------------------------
vector<Point> convexOrder4(vector<Point> p) {
    // Step 1: Pick lowest-leftmost point
    int start = 0;
    for (int i = 1; i < 4; i++)
        if (p[i].y < p[start].y ||
            (p[i].y == p[start].y && p[i].x < p[start].x))
            start = i;

    swap(p[0], p[start]);
    Point base = p[0];

    // Step 2: Sort by polar angle around base
    sort(p.begin() + 1, p.end(),
         [&](const Point& a, const Point& b) {
             double o = orient(base, a, b);
             if (o == 0) {
                 return dist2(base, a) < dist2(base, b);  // closer first
             }
             return o > 0; // CCW first
         });

    // Step 3: Graham-style hull (should return size 4 for convex quad)
    vector<Point> h;
    for (auto& pt : p) {
        while (h.size() >= 2 && orient(h[h.size()-2], h.back(), pt) <= 0)
            h.pop_back();
        h.push_back(pt);
    }

    // h now contains vertices in CCW order
    return h;  // guaranteed convex order (since input quad is convex)
}

// -----------------------------------------------------
// 2) Point inside-or-on-boundary of convex quad (in CCW order)
// -----------------------------------------------------
bool pointInConvexQuadOrdered(const Point& p, const vector<Point>& q) {
    // Determine polygon orientation
    double o0 = orient(q[0], q[1], q[2]);
    bool isCCW = o0 > 0;

    for (int i = 0; i < 4; i++) {
        const Point& a = q[i];
        const Point& b = q[(i + 1) % 4];

        double o = orient(a, b, p);

        if (o == 0) {
            // Check if p lies on segment ab
            double dot1 = (p.x - a.x) * (b.x - a.x) + (p.y - a.y) * (b.y - a.y);
            double dot2 = (p.x - b.x) * (a.x - b.x) + (p.y - b.y) * (a.y - b.y);
            if (dot1 >= 0 && dot2 >= 0)
                continue;  // on edge → inside
            else
                return false;
        }

        // Strict inside test
        if (isCCW && o < 0) return false;
        if (!isCCW && o > 0) return false;
    }

    return true;
}

// -----------------------------------------------------
// 3) Public function: does NOT require ordered quad vertices
// -----------------------------------------------------
bool pointInConvexQuadUnordered(const Point& p, vector<Point> quad) {
    // Reorder quad into CCW convex order
    vector<Point> cq = convexOrder4(quad);

    // Now perform correct inside/boundary test
    return pointInConvexQuadOrdered(p, cq);
}



vector<Quadrilateral> findAllQuadrilaterals(Triangulation T){
    vector<Quadrilateral> q;
    for (int i = 0; i < T.n; ++i){
    for (int i1 = i+1; i1 < T.n; ++i1){
    for (int i2 = i1+1; i2 < T.n; ++i2){
    for (int i3 = i2 + 1; i3 < T.n; ++i3){
        if (isConvexQuadrilateral({T.pts[i], T.pts[i1], T.pts[i2], T.pts[i3]})){
            bool ok = 1;
            for (int j = 0; j < T.n; ++j){
                if (j == i || j == i1 || j == i2 || j == i3) continue;
                if (pointInConvexQuadUnordered(T.pts[j], {T.pts[i], T.pts[i1], T.pts[i2], T.pts[i3]})){
                    ok = 0;
                    break;
                }
            }

            if (ok)
                q.push_back({T.pts[i], T.pts[i1], T.pts[i2], T.pts[i3]});
        }
    }}}}

    return q;
}

struct Triangle{
    Point a, b, c;
};

bool pointInOrOnTriangle(const Point& p,
                         const Point& a,
                         const Point& b,
                         const Point& c)
{
    double o1 = orient(a, b, p);
    double o2 = orient(b, c, p);
    double o3 = orient(c, a, p);

    bool has_pos = (o1 > 0) || (o2 > 0) || (o3 > 0);
    bool has_neg = (o1 < 0) || (o2 < 0) || (o3 < 0);

    // Inside or on boundary iff all orientations have same sign or zero
    return !(has_pos && has_neg);
}
bool sameSideOfLine(const Point& a, const Point& b,
                    const Point& p, const Point& q)
{
    double o1 = orient(a, b, p);
    double o2 = orient(a, b, q);

    // Same side or on the line
    return (o1 >= 0 && o2 >= 0) || (o1 <= 0 && o2 <= 0);
}


vector<Quadrilateral> findAllQuadrilateralsFast(Triangulation T){
    vector<Quadrilateral> q;

    vector<vector<vector<int>>> t(T.n, vector<vector<int>>(T.n));
    for (int i = 0; i < T.n; ++i){
        for (int i1 = i+1; i1 < T.n; ++i1){
            for (int i2 = i1+1; i2 < T.n; ++i2){
                bool empty = 1;
                for (int i3 = 0; i3 < T.n; ++i3){
                    if (i3 == i2 || i3 == i1 || i3 == i) continue;

                    if (pointInOrOnTriangle(T.pts[i3], T.pts[i], T.pts[i1], T.pts[i2])){
                        empty = 0;
                        break;
                    }


                }

                if (!empty) continue;

                t[i][i1].push_back(i2);
                t[i][i2].push_back(i1);
                t[i1][i2].push_back(i);
            }
        }
    }

    for (int i = 0; i < T.n; ++i){
        for (int j = i+1; j < T.n; ++j){
            for (int l = 0; l < t[i][j].size(); ++l){
                for (int r = l+1; r < t[i][j].size(); ++r){
                    if (make_pair(i,j) > make_pair(t[i][j][l], t[i][j][r])) continue;
                    if (sameSideOfLine(T.pts[i], T.pts[j], T.pts[t[i][j][l]], T.pts[t[i][j][r]]))
                        continue;

                    if (isConvexQuadrilateral({T.pts[i], T.pts[j], T.pts[t[i][j][l]], T.pts[t[i][j][r]]})){
                        vector<int> tmp = {i, j, t[i][j][l], t[i][j][r]};
                        sort(tmp.begin(),tmp.end());
                        q.push_back({T.pts[tmp[0]], T.pts[tmp[1]], T.pts[tmp[2]],T.pts[tmp[3]]});
                    }
                }
            }
        }
    }

    return q;
}
