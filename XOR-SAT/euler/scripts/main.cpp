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

    string path = "../formulations/";

    if (file_cnt >= 10)
        path += to_string(file_cnt);
    else
        path += "0" + to_string(file_cnt);

    path += ".cnf";

    return path;
}

string edge_nxt(){

    string path = "../edges/";

    if (file_cnt >= 10)
        path += to_string(file_cnt);
    else
        path += "0" + to_string(file_cnt);

    path += ".cnf";

    return path;
}


int main(int argc, char *argv[]) {

    assert(argc >= 2);
    omp_set_num_threads(48);

    string input_path = argv[1];
    //string id_path = argv[2];
    //string validate_path = argv[3];
    //string input_path = "random_instance_840_15_3";
    ofstream validate("../validate/validate_path.txt");
    validate << input_path << '\n';
    validate.close();

    ifstream in("../data/" + input_path + ".in");
    if (!in.is_open()) {
        cerr << "ERROR: Could not open input file for: " << input_path << endl;
        return 1;
    }
    double maxX = INT32_MIN, maxY = INT32_MIN;

    double minX = INT32_MAX, minY = INT32_MIN;

    int num_triangulations = 2;
    int n, m;
    in >> n >> m >> num_triangulations;
    //num_triangulations = 2;
    vector<Point> p(n);
    for (int i = 0; i < n; ++i) {
        in >> p[i].x >> p[i].y;
        p[i].ind = i;
        minX = min(minX, p[i].x);
        minY = min(minY, p[i].y);
        maxX = max(maxX, p[i].x);
        maxY = max(maxY, p[i].y);
    }

    empty_folder("../formulations/");
    empty_folder("../edges/");

    ifstream ind("../config/distances.txt");

    ifstream ins("../config/selection.txt");
    unordered_set<int> selection;
    int lensel;
    ins >> lensel;
    cout << "AJME " << lensel << '\n';

    if (!lensel){
	for (int i = 0; i < num_triangulations; ++i){
		selection.insert(i);
	}
    } else {
	for (int i = 0; i < lensel; ++i){
		int x;
		ins >> x;
		selection.insert(x);
	}
    }

    vector<int> distances;
    for (int i = 0; i < num_triangulations; ++i){
        int tmp;
        ind >> tmp;
        distances.push_back(tmp);
    }

    ind.close();

    auto start = chrono::high_resolution_clock::now();
    int tt = 0;
    ofstream izlaz("/izlaz.txt");
    while(num_triangulations--) {
        vector<Edge> eA;
	++tt;

        for (int i = 0; i < m; ++i) {
            int u, v;
            in >> u >> v;

            eA.push_back({u, v});
        }

	if (!selection.count(tt-1)) continue;

        Triangulation A;
        A.n  = n;
        A.pts = p;
        A.edges = eA;
        A.maxX = maxX;
        A.maxY = maxY;
        A.minX = minX;
        A.minY = minY;

        string path = nxt();
        string edge_path = edge_nxt();

        formulate(A, distances[tt-1], path, edge_path, input_path);
	izlaz << tt-1 << ' ' << distances[tt-1] << '\n';
    }

    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);


    ofstream out_main("../formulations/00.cnf");
    izlaz << "Number of clauses " << count_clauses << '\n';
    izlaz << "Number of variables " << cur-1 << '\n';
    izlaz << "Time: " << duration.count() << endl;
    izlaz.close();


    out_main << "p cnf " << cur.load()-1 << ' ' << count_clauses << '\n';

    return 0;
}
