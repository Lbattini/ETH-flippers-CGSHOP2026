// main.cpp
#include <bits/stdc++.h>
#include "./include/CutTriangulation.hpp"
#include "./include/Theorem2Transformer.hpp"
#include "./include/Theorem3Transformer.hpp"

using std::string;
using std::vector;

// ============ sitni utili ============
static inline long long ekey(int a,int b){ if(a>b) std::swap(a,b); return ((long long)a<<32) ^ (unsigned long long)b; }
static inline std::unordered_set<long long> edge_set_of(const vector<Edge>& E){
    std::unordered_set<long long> S; S.reserve(E.size()*2);
    for (auto e: E) if (e.u!=e.v) S.insert(ekey(e.u,e.v));
    return S;
}
static void print_edges(const vector<Edge>& E, const char* title){
    vector<Edge> sorted=E;
    std::sort(sorted.begin(), sorted.end(), [](const Edge& a,const Edge& b){ return (a.u==b.u)? a.v<b.v : a.u<b.u; });
    std::cout << title << " (" << sorted.size() << "): ";
    for (auto e: sorted) std::cout << "("<<e.u<<","<<e.v<<") ";
    std::cout << "\n";
}
static void print_edge_diff_sets(const vector<Edge>& A, const vector<Edge>& B, const char* label){
    auto SA=edge_set_of(A), SB=edge_set_of(B);
    vector<std::pair<int,int>> onlyA, onlyB;
    for (auto k: SA) if (!SB.count(k)) onlyA.push_back({int(k>>32), int(k&0xffffffff)});
    for (auto k: SB) if (!SA.count(k)) onlyB.push_back({int(k>>32), int(k&0xffffffff)});
    std::sort(onlyA.begin(), onlyA.end()); std::sort(onlyB.begin(), onlyB.end());
    std::cout << label << "  A\\B=" << onlyA.size() << ", B\\A=" << onlyB.size() << "\n";
    if (!onlyA.empty()) { std::cout << "   in A not in B: "; for (auto [u,v]: onlyA) std::cout<<"("<<u<<","<<v<<") "; std::cout<<"\n"; }
    if (!onlyB.empty()) { std::cout << "   in B not in A: "; for (auto [u,v]: onlyB) std::cout<<"("<<u<<","<<v<<") "; std::cout<<"\n"; }
}

// ============ čitanje fajla / JSON parser za očekivani format ============
static string read_file(const string& path){
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("Cannot open file: " + path);
    std::ostringstream ss; ss << in.rdbuf();
    return ss.str();
}
static void skip_ws(const string& s, size_t& i){ while (i<s.size() && isspace((unsigned char)s[i])) ++i; }

static vector<double> parse_number_array(const string& s, size_t& i){
    vector<double> out; skip_ws(s,i); if (i>=s.size() || s[i]!='[') throw std::runtime_error("Expected '[' num array");
    ++i; skip_ws(s,i);
    while (i<s.size() && s[i]!=']'){
        size_t j=i;
        while (j<s.size() && (isdigit((unsigned char)s[j])||s[j]=='-'||s[j]=='+'||s[j]=='.'||s[j]=='e'||s[j]=='E')) ++j;
        if (j==i) throw std::runtime_error("Bad number");
        out.push_back(std::stod(s.substr(i,j-i)));
        i=j; skip_ws(s,i);
        if (s[i]==','){++i; skip_ws(s,i);} else if (s[i]==']') break; else throw std::runtime_error("num array sep");
    }
    if (i>=s.size() || s[i]!=']') throw std::runtime_error("Unclosed num array");
    ++i; return out;
}

static vector<Edge> parse_edge_array(const string& s, size_t& i){
    vector<Edge> E; skip_ws(s,i); if (i>=s.size() || s[i]!='[') throw std::runtime_error("Expected '[' edge array");
    ++i; skip_ws(s,i);
    while (i<s.size() && s[i]!=']'){
        if (s[i]!='[') throw std::runtime_error("Expected '[' pair");
        ++i; skip_ws(s,i);
        size_t j=i; while (j<s.size() && (isdigit((unsigned char)s[j])||s[j]=='-')) ++j; if (j==i) throw std::runtime_error("Bad u");
        int u=std::stoi(s.substr(i,j-i)); i=j; skip_ws(s,i); if (s[i]!=',') throw std::runtime_error("Missing ,"); ++i; skip_ws(s,i);
        j=i; while (j<s.size() && (isdigit((unsigned char)s[j])||s[j]=='-')) ++j; if (j==i) throw std::runtime_error("Bad v");
        int v=std::stoi(s.substr(i,j-i)); i=j; skip_ws(s,i); if (s[i]!=']') throw std::runtime_error("Missing ]"); ++i; skip_ws(s,i);
        if (u!=v) E.push_back({std::min(u,v), std::max(u,v)});
        if (s[i]==','){++i; skip_ws(s,i);} else if (s[i]==']') break; else throw std::runtime_error("edge array sep");
    }
    if (i>=s.size() || s[i]!=']') throw std::runtime_error("Unclosed edge array");
    ++i;
    auto S=edge_set_of(E); vector<Edge> ded; ded.reserve(S.size());
    for (auto k:S) ded.push_back({int(k>>32), int(k&0xffffffff)});
    return ded;
}

static void parse_input_json(const string& text, vector<Point>& P, vector<Edge>& E1, vector<Edge>& E2){
    size_t i=0; auto find_key=[&](const string& key){ size_t pos=text.find("\""+key+"\"",i); if(pos==string::npos) throw std::runtime_error("Key not found: "+key); return pos+key.size()+2; };
    i=find_key("points_x"); while (i<text.size() && text[i]!='[') ++i; auto xs=parse_number_array(text,i);
    i=find_key("points_y"); while (i<text.size() && text[i]!='[') ++i; auto ys=parse_number_array(text,i);
    if (xs.size()!=ys.size()) throw std::runtime_error("points_x/points_y length mismatch");
    P.resize(xs.size()); for (size_t k=0;k<xs.size();++k){ P[k].x=xs[k]; P[k].y=ys[k]; }

    i=find_key("triangulations"); while (i<text.size() && text[i]!='[') ++i; ++i; // ulazimo u listu triangulacija
    while (i<text.size() && text[i]!='[') ++i; E1=parse_edge_array(text,i);
    while (i<text.size() && text[i]!='[') { if (text[i]==']') break; ++i; }
    if (i<text.size() && text[i]=='[') E2=parse_edge_array(text,i); else throw std::runtime_error("Second triangulation not found");
}

// ============ kanonska meta triangulacija (zavisna samo od tačaka) ============
static double choose_y0(const vector<int>& V, const vector<Point>& P){
    vector<double> ys; ys.reserve(V.size()); for (int v: V) ys.push_back(P[v].y);
    std::sort(ys.begin(), ys.end());
    const double tol=1e-6;
    vector<double> uniq; uniq.reserve(ys.size());
    for (double y: ys) if (uniq.empty() || std::fabs(y-uniq.back())>tol) uniq.push_back(y);
    if (uniq.size()<2) return std::numeric_limits<double>::quiet_NaN();
    size_t gi=(uniq.size()-1)/2;
    return 0.5*(uniq[gi]+uniq[gi+1]);
}
static vector<int> convex_hull_indices(const vector<Point>& P, const vector<int>& V){
    vector<int> I=V;
    std::sort(I.begin(), I.end(), [&](int a,int b){
        if (P[a].x!=P[b].x) return P[a].x<P[b].x;
        if (P[a].y!=P[b].y) return P[a].y<P[b].y;
        return a<b;
    });
    auto cross=[&](int o,int a,int b){ long double x1=P[a].x-P[o].x, y1=P[a].y-P[o].y; long double x2=P[b].x-P[o].x, y2=P[b].y-P[o].y; return x1*y2 - y1*x2; };
    vector<int> H;
    for (int idx: I){ while (H.size()>=2 && cross(H[H.size()-2], H.back(), idx) <= 0) H.pop_back(); H.push_back(idx); }
    size_t lower_sz=H.size();
    for (int k=(int)I.size()-2; k>=0; --k){
        int idx=I[k];
        while (H.size()>lower_sz && cross(H[H.size()-2], H.back(), idx) <= 0) H.pop_back();
        H.push_back(idx);
        if (k==0) break;
    }
    if (!H.empty()) H.pop_back();
    return H; // CCW
}
static void strip_triangulate(const vector<int>& chain_low_x,
                              const vector<int>& chain_up_x,
                              vector<Edge>& outE)
{
    auto add=[&](int a,int b){ if (a!=b) outE.push_back({std::min(a,b),std::max(a,b)}); };
    int p=(int)chain_low_x.size(), m=(int)chain_up_x.size(), k=std::min(p,m);
    for (int i=1;i<p;++i) add(chain_low_x[i-1], chain_low_x[i]);
    for (int j=1;j<m;++j) add(chain_up_x[j-1],  chain_up_x[j]);
    if (p && m){ add(chain_low_x.front(), chain_up_x.front()); add(chain_low_x.back(), chain_up_x.back()); }
    for (int i=1;i<k;++i){
        int L_im1=chain_low_x[i-1], L_i=chain_low_x[i], U_im1=chain_up_x[i-1], U_i=chain_up_x[i];
        // fiksna orijentacija: (L_{i-1}, U_i)
        add(L_im1, U_i); add(U_i, U_im1); add(U_im1, L_im1);
        add(L_im1, L_i); add(L_i, U_i);   add(U_i, L_im1);
    }
    if (m>p){ int base=chain_low_x.back(); for (int j=p;j<m;++j){ add(base, chain_up_x[j-1]); add(chain_up_x[j-1], chain_up_x[j]); add(chain_up_x[j], base); } }
    if (p>m){ int base=chain_up_x.back();  for (int i=m;i<p;++i){ add(chain_low_x[i-1], chain_low_x[i]); add(chain_low_x[i], base); add(base, chain_low_x[i-1]); } }
}
static void build_canonical_target_rec(const vector<Point>& P, const vector<int>& V, vector<Edge>& Eout)
{
    if ((int)V.size()<=3) return;
    double y0 = choose_y0(V, P);
    if (!std::isfinite(y0)) return;

    vector<int> Vneg, Vpos; Vneg.reserve(V.size()); Vpos.reserve(V.size());
    for (int v: V){ ((P[v].y < y0) ? Vneg : Vpos).push_back(v); }
    if (Vneg.empty() || Vpos.empty() || Vneg.size()==V.size() || Vpos.size()==V.size()) return;

    auto Hpos = convex_hull_indices(P, Vpos);
    auto Hneg = convex_hull_indices(P, Vneg);
    auto xsort=[&](vector<int> h){ std::sort(h.begin(), h.end(),
                                             [&](int a,int b){ if (P[a].x!=P[b].x) return P[a].x<P[b].x; if (P[a].y!=P[b].y) return P[a].y<P[b].y; return a<b; }); return h; };
    vector<int> chain_up_x  = xsort(Hpos);
    vector<int> chain_low_x = xsort(Hneg);

    strip_triangulate(chain_low_x, chain_up_x, Eout);
    build_canonical_target_rec(P, Vpos, Eout);
    build_canonical_target_rec(P, Vneg, Eout);
}
static vector<Edge> build_canonical_target(const vector<Point>& P){
    vector<int> V(P.size()); std::iota(V.begin(), V.end(), 0);
    vector<Edge> E; E.reserve(P.size()*3);
    build_canonical_target_rec(P, V, E);
    auto S=edge_set_of(E); vector<Edge> ded; ded.reserve(S.size());
    for (auto k:S) ded.push_back({int(k>>32), int(k&0xffffffff)});
    return ded;
}

// ============ MAIN ============
int main(int argc, char** argv){
    try{
        const string path = (argc>=2 ? argv[1] : "../data/example_ps_10_nt2_pdf_random.json");
        std::cout << "========== PROJECT RUNNER (Theorem 3 — print outputs) ==========\n";
        std::cout << "Reading: " << path << "\n";

        vector<Point> P; vector<Edge> E1, E2;
        {
            string text = read_file(path);
            parse_input_json(text, P, E1, E2);
        }
        std::cout << "Loaded " << P.size() << " points; T1 edges: " << E1.size() << ", T2 edges: " << E2.size() << "\n";

        // Kanonska meta (zavisi samo od P)
        vector<Edge> Ecanon = build_canonical_target(P);
        if (Ecanon.empty()) {
            std::cerr << "Warning: canonical target is empty (degenerate split). We'll still print T1→T2.\n";
        }

        // Pokreni Theorem 3 ka kanonskoj za obe triangulacije
        vector<Edge> T3_1 = E1;
        vector<Edge> T3_2 = E2;
        auto stats1 = theorem3_pset::transform(P, T3_1, Ecanon, /*dbg_out=*/nullptr);
        auto stats2 = theorem3_pset::transform(P, T3_2, Ecanon, /*dbg_out=*/nullptr);

        // Odštampaj rezultate (edge liste) i truglove koje induciraju
        std::cout << "\n===== OUTPUTS FROM THEOREM 3 =====\n";
        print_edges(T3_1, "T1 after Theorem 3");
        auto F1 = triangulation_cut::extractFacesAsTriangles(P, T3_1).triangles;
        std::cout << "Triangles (T1 after T3) [" << F1.size() << "]: ";
        for (auto t: F1) std::cout << "["<<t[0]<<","<<t[1]<<","<<t[2]<<"] ";
        std::cout << "\n";

        print_edges(T3_2, "T2 after Theorem 3");
        auto F2 = triangulation_cut::extractFacesAsTriangles(P, T3_2).triangles;
        std::cout << "Triangles (T2 after T3) [" << F2.size() << "]: ";
        for (auto t: F2) std::cout << "["<<t[0]<<","<<t[1]<<","<<t[2]<<"] ";
        std::cout << "\n";

        // Uporedi (ako želiš)
        auto S1=edge_set_of(T3_1), S2=edge_set_of(T3_2);
        std::cout << "\nEqual after Theorem 3? " << (S1==S2 ? "YES" : "NO") << "\n";
        if (S1!=S2) print_edge_diff_sets(T3_1, T3_2, "Diff(T3_1, T3_2):");

        // Info o radu
        std::cout << "\nStats T1→canon: depth="<<stats1.depth<<", nodes="<<stats1.nodes
                  <<", l4="<<stats1.lemma4_calls<<", l5="<<stats1.lemma5_calls
                  <<", flips_total="<<stats1.flips_total<<"\n";
        std::cout << "Stats T2→canon: depth="<<stats2.depth<<", nodes="<<stats2.nodes
                  <<", l4="<<stats2.lemma4_calls<<", l5="<<stats2.lemma5_calls
                  <<", flips_total="<<stats2.flips_total<<"\n";

        // (sanity) direktno T1→T2, čisto da vidiš da li driver išta radi kada target != current
        vector<Edge> T12 = E1;
        auto stats12 = theorem3_pset::transform(P, T12, E2, /*dbg_out=*/nullptr);
        auto Sd=edge_set_of(T12), S2t=edge_set_of(E2);
        std::cout << "\nSanity T1→T2 equal? " << (Sd==S2t ? "YES" : "NO") << "\n";
        if (Sd!=S2t) print_edge_diff_sets(T12, E2, "Diff(T1→T2 result, T2):");
        std::cout << "Stats T1→T2: depth="<<stats12.depth<<", nodes="<<stats12.nodes
                  <<", l4="<<stats12.lemma4_calls<<", l5="<<stats12.lemma5_calls
                  <<", flips_total="<<stats12.flips_total<<"\n";

        return 0;
    }
    catch (const std::exception& ex){
        std::cerr << "\n❌ Error: " << ex.what() << "\n";
        return 1;
    }
}
