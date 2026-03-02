#pragma once
#include <bits/stdc++.h>
#include "Lemma4ManualHarness.hpp"

namespace lemma4_manual_tests {

    inline void run_all_lemma4_manual_tests() {
        using namespace lemma4_manual_harness;

        // --- your manual instance ---
        std::vector<Point> P = {
                {0,0},{3,3},{7,4},{11,3},{13,0},{0,8},
                {4,9},{7,15},{9,9},{13,8}
        };
        std::vector<Edge> E = {
                {0,1},{1,2},{2,3},{3,4},{4,9},{8,9},{7,8},{6,7},{5,6},{0,5},
                {5,1},{1,6},{2,6},{2,7},{2,8},{3,8},{3,9}
        };
        int v1 = 5, vq = 9;

        std::cout << "===== Lemma 4 — Manual Instance Test =====\n";

        L4RunResult res;
        // require_full_triangulation = false (as requested)
        bool ok = run_lemma4_on(P, E, v1, vq, res, &std::cout, /*require_full_triangulation=*/false);

        if (!ok) {
            std::cerr << "❌ Lemma4 failed: " << res.error << "\n";
            throw std::runtime_error(res.error);
        }

        std::cout << "✅ Lemma4 OK: flips=" << res.flips
                  << ", rounds=" << res.rounds
                  << ", rejects=" << res.rejects << "\n";

        // Final triangulation edges (derived from faces)
        print_faces_as_edges(res.faces_out, std::cout);

        // Optional flip log (if your Lemma4 exposes export_flip_log)
        if (!res.flip_seq.empty()) {
            std::cout << "\nFLIPS\n";
            for (const auto& round : res.flip_seq) {
                for (const auto& f : round) std::cout << f.a << ' ' << f.b << '\n';
                std::cout << "------------------\n";
            }
        } else {
            std::cout << "\n(no flip sequence exported — implement export_flip_log(...) to enable)\n";
        }


        //OBSERVATION 2 TEST


        // --- your manual instance ---
        P = {
                {0,0},{3,1},{6,4},{7,8},{12,2},{8,-6},{6,-3},{3, -1}
        };
        E = {
                {0,1},{1,2},{2,3},{3,4},{4,5}, {5, 6}, {6, 7}, {7,0},
                {1,7},{2,7},{2,6},{3,6},{4,6}
        };
        v1 = 4;
        vq = 0;

        std::cout << "===== Lemma 4 — Manual Instance Test - OBSERVATION 2=====\n";

        L4RunResult res1;
        // require_full_triangulation = false (as requested)
        ok = run_lemma4_on(P, E, v1, vq, res1, &std::cout, /*require_full_triangulation=*/false);

        if (!ok) {
            std::cerr << "❌ Lemma4 failed: " << res1.error << "\n";
            throw std::runtime_error(res1.error);
        }

        std::cout << "✅ Lemma4 OK: flips=" << res1.flips
                  << ", rounds=" << res1.rounds
                  << ", rejects=" << res1.rejects << "\n";

        // Final triangulation edges (derived from faces)
        print_faces_as_edges(res1.faces_out, std::cout);

        // Optional flip log (if your Lemma4 exposes export_flip_log)
        if (!res1.flip_seq.empty()) {
            std::cout << "\nFLIPS\n";
            for (const auto& round : res1.flip_seq) {
                for (const auto& f : round) std::cout << f.a << ' ' << f.b << '\n';
                std::cout << "------------------\n";
            }
        } else {
            std::cout << "\n(no flip sequence exported — implement export_flip_log(...) to enable)\n";
        }
    }

} // namespace lemma4_manual_tests
