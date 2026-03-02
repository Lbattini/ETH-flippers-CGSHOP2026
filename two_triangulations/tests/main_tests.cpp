#include <iostream>
#include "./CutTriangulationTests.hpp"
#include "./ForwardFlipCanonicalizerTests.hpp"
#include "./DiagonalIntroducerLemma4Tests.hpp"
#include "./Lemma5Tests.hpp"
#include "./Proposition1Tests.hpp"
#include "./Theorem2Tests.hpp"
#include "./Observation3Tests.hpp"
#include "./Theorem3Tests.hpp"
#include "./Lemma4ManualTest.hpp"
using namespace std;

int main() {
    cout << "========== PROJECT TEST RUNNER ==========\n";

    try {
       // triangulation_cut_tests::run_all_cut_triangulation_tests();
       // forward_flip_tests::run_all_forward_flip_tests();

   //     lemma4_manual_tests::run_all_lemma4_manual_tests();
        lemma5_tests::run_all_lemma5_tests();
        return 0;
        //proposition1_tests::run_all_proposition1_tests();
        //theorem2_tests::run_all_theorem2_tests();
        //observation3_tests::run_all_observation3_tests();
        //theorem3_tests::run_all_theorem3_tests();
    }
    catch (const exception& ex) {
        cerr << "\n❌ Tests failed: " << ex.what() << "\n";
        return 1;
    }

    cout << "\n✅ All tests finished successfully.\n";
    return 0;
}
