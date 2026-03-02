#ifndef GENERATECNF_HPP
#define GENERATECNF_HPP

#include<flippingDS.hpp>

typedef std::vector<int> vi;
typedef std::vector<vi> vvi;
typedef std::vector<bool> vb;
typedef std::vector<vb> vvb;
typedef std::vector<vvb> vvvb;

void genCNF(const vP& points, const vvi& qg, const vvvb& emptyTriangle, vvi& crossingEdgesVec, std::vector<facesDS>& tVec, int n, vi& distCenter, const std::string& outFileName, const std::string& idxFileName);

void genCNF(const vP& points, const vvi& qg, vvi& crossingEdgesVec, std::vector<facesDS>& tVec, int n, int distCenter, const std::string& outFileName);

void cnfToFlipList(const std::string& inFileName, const std::string& idxFileName, std::vector<std::vector<vi2>>& flipsPerRound);

void cnfToFlipList(const std::string& inFileName, size_t n, size_t m, int distCenter, std::vector<std::vector<vi2>>& flipsPerRound);
#endif //GENERATECNF_HPP
