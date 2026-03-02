//#include <bits/stdc++.h>
#include "Utils.hpp"
#include <atomic>
#include <omp.h>

using namespace std;

atomic<int> cur{1};
ofstream out;
int count_clauses = 0;
vector<vector<int>> c;
vector<vector<vector<Edge>>> qg;

void single_clause(int a){
    ++count_clauses;
    out << a << " 0\n";
}

void clause(int a, int b){
    ++count_clauses;
    out << a << ' ' << b << " 0\n";
}

void implication(int a, int b){
    clause(-a, b);
}

void equiv(int a, int b){
    implication(a, b);
    implication(b, a);
}

void xorclause(vector<int> a, int res){
    int n = a.size();

    ++count_clauses;
    out << "x";
    for (int i = 0; i < n-1; ++i)
        out << a[i] << ' ';

    if (!res)
        out << -a[n-1] << " 0\n";
    else
        out << a[n-1] << " 0\n";
}

vector<Quadrilateral> qs;
//ofstream centar_out;
bool debug = 0;
void formulate(Triangulation A, int k, string path, string edge_path, string input_path){

    int second = 0;
    if (!qs.empty())
        second++;
    else if (filesystem::exists("../quads/" + input_path + ".txt")){
        int lenqs;
        ifstream quads("../quads/" + input_path + ".txt");

        quads >> lenqs;
        qs = vector<Quadrilateral>(lenqs);
        for (int i = 0; i < lenqs; ++i){
            int a1, b1, c1, d1;
            quads >> a1 >> b1 >> c1 >> d1;

            qs[i] = {A.pts[a1], A.pts[b1], A.pts[c1], A.pts[d1]};
        }
        quads.close();
        qg = intersecting_diagonals(qs, A.n);
    } else{
        qs = findAllQuadrilateralsSuperFast(A);
        ofstream quads("../quads/" + input_path + ".txt");

        quads << (int)qs.size() << '\n';
        for (auto q : qs){
            quads << q.a.ind << ' ' << q.b.ind << ' ' << q.c.ind << ' ' << q.d.ind << '\n';
        }

        quads.close();
        qg = intersecting_diagonals(qs, A.n);
    }

    auto might = computeMightBePresentOne(A.n, qg, A, k);

    cout << "DONE PREPROCESS\n";
    int nqs = qs.size();

    vector<int> qmight(nqs, 0);
    for (int i = 0; i < nqs; ++i){
        auto h = convexOrder4({qs[i].a, qs[i].b, qs[i].c, qs[i].d});

        qmight[i] = max(qmight[i], might[h[0].ind][h[1].ind]);
        qmight[i] = max(qmight[i], might[h[1].ind][h[2].ind]);
        qmight[i] = max(qmight[i], might[h[2].ind][h[3].ind]);
        qmight[i] = max(qmight[i], might[h[3].ind][h[0].ind]);
    }

    vector<vector<int>> f(nqs, vector<int>(k+1, 0));
    for (int i = 1; i <= k; ++i){
        for (int j = 0; j < nqs; ++j){

            if (qmight[j] >= i) continue;

            f[j][i] = cur++;

            if (debug)
            cout << "FLIP (" << qs[j].a.ind << ' ' << qs[j].b.ind << ' ' <<
                qs[j].c.ind << ' ' << qs[j].d.ind << ") " << cur-1 << '\n';
        }
    }

    if (!second) {
        //centar_out.open("C:\\libs\\cryptominisat5-win\\centar.txt");
        c = vector<vector<int>>(A.n, vector<int>(A.n));
    }
    vector<vector<vector<int>>> e(A.n, vector<vector<int>>(A.n, vector<int>(k+1, 0)));
    for (int i = 0; i < A.n; ++i){
        for (int j = i+1; j < A.n; ++j){
            for (int l = 0; l <= k; ++l){
                if (might[i][j] > l)
                    continue;
                e[j][i][l] = cur;
                if (debug)
                cout << "EDGE " << i << ' ' << j << " AT MOMENT " << l << ": " << cur << '\n';
                e[i][j][l] = cur++;
            }

            if (!second){
                //centar_out << "CENTAR: " << i << ' ' << j << " " << cur << '\n';
                c[i][j] = cur;
                c[j][i] = cur++;
            }
        }
    }

    out.open(path);

    // geometric precondition
    for (int i = 0; i < nqs; ++i){
        auto ps = convexOrder4({qs[i].a, qs[i].b, qs[i].c, qs[i].d});

        int v1 = ps[0].ind;
        int v2 = ps[1].ind;
        int v3 = ps[2].ind;
        int v4 = ps[3].ind;
      //  cout << "Q " << v1 << ' ' << v2 << ' ' << v3 << ' ' << v4 << '\n';
        for (int j = 1; j <= k; ++j) {
            if (!f[i][j]) continue;

            if (!e[v1][v2][j-1] || !e[v2][v3][j-1] || !e[v3][v4][j-1] || !e[v4][v1][j-1]){
                single_clause(-f[i][j]);
                continue;
            }
            implication(f[i][j],e[v1][v2][j - 1] );
            implication(f[i][j],e[v2][v3][j - 1] );
            implication(f[i][j],e[v3][v4][j - 1] );
            implication(f[i][j],e[v4][v1][j - 1] );
        }
    }


    // target
    for (int i = 0; i < A.n; ++i){
        for (int j = i + 1; j < A.n; ++j){
            if (e[i][j][k]) {
                equiv(e[i][j][k], c[i][j]);
            } else {
                single_clause(-c[i][j]);
            }
        }
    }

    // Mechanics of State Evolution

    vector<vector<vector<int>>> c_edge(A.n, vector<vector<int>>(A.n));

    for (int i = 0; i < nqs; ++i){
        auto ps = convexOrder4({qs[i].a, qs[i].b, qs[i].c, qs[i].d});

        int v1 = ps[0].ind;
        int v2 = ps[1].ind;
        int v3 = ps[2].ind;
        int v4 = ps[3].ind;

        c_edge[min(v1, v3)][max(v1,v3)].push_back(i);
        c_edge[min(v2, v4)][max(v2,v4)].push_back(i);
    }

    vector<Edge> hull;
    for (int u = 0; u < A.n; ++u){
        for (int v = u+1; v < A.n; ++v){
            if (c_edge[u][v].empty()){
                hull.push_back({u, v});
                continue;
            }

            for (int t = 1; t <= k; ++t){
                vector<int> a = {};

                if (e[u][v][t-1])
                    a.push_back(e[u][v][t-1]);
                if (e[u][v][t])
                    a.push_back(e[u][v][t]);

                for (int q : c_edge[u][v]){
                    if (!f[q][t]) continue;
                    a.push_back(f[q][t]);
                }

                if ((int)a.size() == 1){
                    single_clause(-a[0]);
                    continue;
                }

                if (!a.empty()) {
                    xorclause(a, 0);
                }


            }
        }
    }

    // Parallel Independence

    set<combQuadrilateral> s;
    for (int i = 0; i < nqs; ++i){
        s.insert({qs[i].a.ind, qs[i].b.ind, qs[i].c.ind, qs[i].d.ind, i});
    }

    //cout << "START PAIRWISE\n";
    vector<string> thread_outputs(omp_get_max_threads());
    vector<int> thread_clause_counts(omp_get_max_threads(), 0);
#pragma omp parallel for schedule(dynamic)
    for (int i1 = 0; i1 < A.n; ++i1){
        int tid = omp_get_thread_num();
        string& out_buf = thread_outputs[tid];
        int& clause_cnt = thread_clause_counts[tid];
        for (int i2 = i1 + 1; i2 < A.n; ++i2){
            for (int i3 = i2 + 1; i3 < A.n; ++i3){
                vector<int> conf_qs;
                for (int i4 = 0; i4 < A.n; ++i4){
                    if (i4 == i1 || i4 == i2 || i4 == i3) continue;
                    auto it = s.lower_bound({i1,i2,i3,i4,-1});
                    if (it == s.end()) continue;

                    combQuadrilateral cand = *it;

                    int ind = cand.ind;
                    Quadrilateral qq = qs[ind];

                    if (qs[ind] != combQuadrilateral{i1,i2,i3,i4,-1}) continue;

                    if (!s.count({i1,i2,i3,i4,ind})) continue;

                    conf_qs.push_back(ind);
                }

                for (int t = 1; t <= k; ++t) {

                    vector<int> upd_conf_qs;
                    for (auto q : conf_qs)
                        if (f[q][t])
                            upd_conf_qs.push_back(q);

                    int len = upd_conf_qs.size();
                    if (len <= 3) {
                        for (int i = 0; i < len; ++i) {
                            for (int j = i + 1; j < len; ++j) {
                                //clause(-f[conf_qs[i]][t], -f[conf_qs[j]][t]);
                                out_buf += to_string(-f[upd_conf_qs[i]][t]) + "  " + std::to_string(-f[upd_conf_qs[j]][t]) + " 0\n";
                                clause_cnt++;
                            }
                        }
                    } else {
                        int base = cur.fetch_add(len);
                        vector<int> state(len);
                        for (int i = 0; i < len; ++i) {
                            state[i] = base + i;
                            //clause(-f[conf_qs[i]][t], state[i]);
                            out_buf += to_string(-f[upd_conf_qs[i]][t]) + "  " + to_string(state[i]) + " 0\n";
                            clause_cnt++;

                            if (i > 0) {
                                //clause(-f[conf_qs[i]][t], -state[i - 1]);
                                //clause(-state[i - 1], state[i]);
                                out_buf += to_string(-f[upd_conf_qs[i]][t]) + "  " + to_string(-state[i-1]) + " 0\n";
                                out_buf += to_string(state[i]) + "  " + to_string(-state[i-1]) + " 0\n";

                                clause_cnt += 2;
                            }
                        }
                    }
                }
            }
        }
    }

    for (int x : thread_clause_counts)
        count_clauses += x;


    //cout << "DONE PAIRWISE\n";

    // Set of Existing Edges
    vector<vector<int>> vis(A.n, vector<int>(A.n, 0));
    for (auto ed : A.edges){
        vis[ed.u][ed.v] = vis[ed.v][ed.u] = 1;
        single_clause(e[ed.u][ed.v][0]);
    }

    // should be redundant
    for (int i = 0; i < A.n; ++i){
        for (int j = i + 1; j < A.n; ++j){
            if (!vis[i][j] && e[i][j][0]){
                single_clause(-e[i][j][0]);
            }
        }
    }

    for (auto ee : hull){
        int z = -1;
        if (vis[ee.u][ee.v]){
            z = 1;
        }

        for (int t = 1; t <= k; ++t)
            if (e[ee.u][ee.v][t])
                single_clause(e[ee.u][ee.v][t] * z);
    }

    for (auto& s : thread_outputs)
        out << s;

    out.close();

    out.open(edge_path);

    out << k+1 << '\n';
    for (int i = 0; i < A.n; ++i){
        for (int j = i + 1; j < A.n; ++j){
            out << "EDGE " << i << ' ' << j << ' ';
            for (int t = 0; t <= k; ++t){
                out << e[i][j][t] << ' ';
            }
            out << '\n';
        }
    }

    out.close();

    edge_path += "m";
    cout << edge_path << '\n';
    out.open(edge_path);

    for (int i = 0; i < A.n; ++i){
        for (int j = i + 1; j < A.n; ++j){
            out << "EDGE " << i << ' ' << j << ' ' << might[i][j] << '\n';
            out << '\n';
        }
    }
    out.close();
}