#include <bits/stdc++.h>
#include "all_pairs.hpp"

using namespace std;

void empty_folder(const filesystem::path& folder)
{
    for (const auto& entry : filesystem::directory_iterator(folder)) {
        filesystem::remove_all(entry.path());
    }
}

int file_cnt = 0;
string nxt(){
    ++file_cnt;

    string path = "C:\\libs\\cryptominisat5-win\\formulations\\";

    if (file_cnt >= 10)
        path += to_string(file_cnt);
    else
        path += "0" + to_string(file_cnt);

    path += ".cnf";

    return path;
}

string edge_nxt(){

    string path = "C:\\libs\\cryptominisat5-win\\edges\\";

    if (file_cnt >= 10)
        path += to_string(file_cnt);
    else
        path += "0" + to_string(file_cnt);

    path += ".cnf";

    return path;
}


int main(int argc, char *argv[]) {

    assert(argc >= 2);
    omp_set_num_threads(8);

    string input_path = argv[1];
  //string input_path = "random_instance_840_15_3";
    ofstream validate("C:/libs/cryptominisat5-win/validate_path.txt");
    validate << input_path;
    validate.close();

    ifstream in("D:/ETH/CG_SHOP/benchmark_instances_rev1/benchmark_instances/" + input_path + ".in");
    if (!in.is_open()) {
        cerr << "ERROR: Could not open input file for: " << input_path << endl;
        return 1;
    }

    int num_triangulations = 2;
    int n, m;
    in >> n >> m >> num_triangulations;
    //num_triangulations = 4;
    vector<Point> p(n);
    for (int i = 0; i < n; ++i) {
        in >> p[i].x >> p[i].y;
        p[i].ind = i;
    }

    empty_folder("C:\\libs\\cryptominisat5-win\\formulations");

    ifstream ind("D:\\ETH\\CG_SHOP\\distances.txt");

    vector<int> distances;
    for (int i = 0; i < num_triangulations; ++i){
        int tmp;
        ind >> tmp;
        distances.push_back(tmp);
    }

    ind.close();
    auto start = chrono::high_resolution_clock::now();
    int tt = 0;
    while(num_triangulations--) {
        vector<Edge> eA;
        ++tt;

        for (int i = 0; i < m; ++i) {
            int u, v;
            in >> u >> v;

            eA.push_back({u, v});
        }

        Triangulation A;
        A.n  = n;
        A.pts = p;
        A.edges = eA;

        string path = nxt();
        string edge_path = edge_nxt();
        formulate(A, distances[file_cnt-1], path, edge_path, input_path);
    }

    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);


    ofstream out_main("C:\\libs\\cryptominisat5-win\\formulations\\00.cnf");
    cout << "Number of clauses " << count_clauses << '\n';
    cout << "Number of variables " << cur-1 << '\n';
    cout << "Time: " << duration.count() << endl;


    out_main << "p cnf " << cur.load()-1 << ' ' << count_clauses << '\n';

    return 0;
}
