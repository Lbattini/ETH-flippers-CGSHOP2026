#pragma once
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include "../include/CutTriangulation.hpp"

namespace triangulation_cut_tests {

    using namespace std;
    using namespace triangulation_cut;

    static void require(bool cond, const string& msg) {
        if (!cond) throw runtime_error("ASSERT FAIL: " + msg);
    }

    static void print_vec(const vector<int>& v) {
        cout << "{";
        for (size_t i=0;i<v.size();++i){ if(i) cout<<","; cout<<v[i]; }
        cout << "}";
    }

// ---------------------- INDIVIDUAL TESTS ----------------------

    static void run_simple_horizontal_cut() {
        cout << "\n[TEST] Simple square, cut by y = 0.5\n";
        vector<Point> pts = { {0,0},{1,0},{1,1},{0,1} };
        vector<Edge>  E   = { {0,1},{1,2},{2,3},{3,0},{0,2} };
        double a=0,b=1,c=-0.5;

        CutResult res; string err;
        bool ok = compute_cut_sets(pts, E, a,b,c, res, &err);
        require(ok, "unexpected failure: " + err);

        sort(res.Ppos.begin(), res.Ppos.end());
        sort(res.Pneg.begin(), res.Pneg.end());
        require(res.Ppos==vector<int>({2,3}), "P+ expected {2,3}");
        require(res.Pneg==vector<int>({0,1}), "P- expected {0,1}");
        require(res.Q.size()==2, "Q expected 2 triangles");

        cout << "  P+ = "; print_vec(res.Ppos); cout << "\n";
        cout << "  P- = "; print_vec(res.Pneg); cout << "\n";
        cout << "  Q count = " << res.Q.size() << "\n";
    }

    static void run_diagonal_cut() {
        cout << "\n[TEST] Tilted diagonal: y = 0.2x + 0.4\n";
        vector<Point> pts = { {0,0},{1,0},{1,1},{0,1} };
        vector<Edge>  E   = { {0,1},{1,2},{2,3},{3,0},{0,2} };

        // Line: y = 0.2x + 0.4  ->  (-0.2)x + 1*y - 0.4 = 0
        double a = -0.2, b = 1.0, c = -0.4;

        CutResult res; string err;
        bool ok = compute_cut_sets(pts, E, a,b,c, res, &err);
        require(ok, "unexpected failure: " + err);

        // Expect P+ = {2,3}, P- = {0,1}, and it intersects both triangles
        sort(res.Ppos.begin(), res.Ppos.end());
        sort(res.Pneg.begin(), res.Pneg.end());
        require(res.Ppos==vector<int>({2,3}), "P+ expected {2,3}");
        require(res.Pneg==vector<int>({0,1}), "P- expected {0,1}");
        require(res.Q.size()==2, "Q expected 2 triangles");

        cout << "  P+ = "; print_vec(res.Ppos); cout << "\n";
        cout << "  P- = "; print_vec(res.Pneg); cout << "\n";
        cout << "  Q count = " << res.Q.size() << "\n";
    }



    static void run_no_intersection() {
        cout << "\n[TEST] Line above everything: y = 2\n";
        vector<Point> pts = { {0,0},{1,0},{1,1},{0,1} };
        vector<Edge>  E   = { {0,1},{1,2},{2,3},{3,0},{0,2} };
        double a=0,b=1,c=-2; // y = 2

        CutResult res; string err;
        bool ok = compute_cut_sets(pts, E, a,b,c, res, &err);
        require(ok, "unexpected failure: " + err);

        sort(res.Ppos.begin(), res.Ppos.end());
        sort(res.Pneg.begin(), res.Pneg.end());
        require(res.Ppos.empty(), "P+ expected empty");
        require(res.Pneg==vector<int>({0,1,2,3}), "P- expected {0,1,2,3}");
        require(res.Q.empty(), "Q expected 0 triangles");

        cout << "  P+ = "; print_vec(res.Ppos); cout << "\n";
        cout << "  P- = "; print_vec(res.Pneg); cout << "\n";
        cout << "  Q count = " << res.Q.size() << "\n";
    }

    static void run_degenerate_vertex_on_line() {
        cout << "\n[TEST] Degenerate: vertex lies on line y=0\n";
        vector<Point> pts = { {0,0},{1,0},{1,1},{0,1} };
        vector<Edge>  E   = { {0,1},{1,2},{2,3},{3,0},{0,2} };
        double a=0,b=1,c=0;

        CutResult res; string err;
        bool ok = compute_cut_sets(pts, E, a,b,c, res, &err);
        require(!ok, "expected failure when vertex lies on line");
        cout << "  Expected failure: " << err << "\n";
    }

// ---------------------- MASTER RUNNER ----------------------

    inline void run_all_cut_triangulation_tests() {
        cout << "\n===== Running Triangulation Cut Tests =====\n";
        run_simple_horizontal_cut();
        run_diagonal_cut();
        run_no_intersection();
        run_degenerate_vertex_on_line();
        cout << "\n✅ Triangulation Cut tests passed.\n";
    }

} // namespace triangulation_cut_tests
