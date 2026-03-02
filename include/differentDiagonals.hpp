#ifndef DIFFERENTDIAGONALS_HPP
#define DIFFERENTDIAGONALS_HPP

#endif //DIFFERENTDIAGONALS_HPP

#include <vector>

#include "../include/inOut.hpp"
#include "../include/flippingDS.hpp"

typedef std::vector<int> vi;
typedef std::unordered_map<int,int> mii;

//try to change t into t1 by flipping diagonals in t which are not present in t1
//this is in general not enough
int flipDifferentDiagonals(flippingDS& t, flippingDS& t1, size_t n, const vP& points, std::vector<vi2>& flipsPerRound);

//O(n) memory
int flipToReduceIntesections(flippingDS& t, flippingDS& t1, size_t n, const vP& points, std::vector<vi2>& flipsPerRound, mii& intersectionsWithTarget);

int flipToReduceIntesections(facesDS& t, facesDS& t1, size_t n, const vP& points, std::vector<vi2>& flipsPerRound, mii& intersectionsWithTarget);

//For at most cntTtoT1 perform flips to transform t into t1; then use revDir to transform t1 into the (partially transformed) t
int flipToReduceIntesections(facesDS& t, facesDS& t1, size_t n, const vP& points, std::vector<vi2>& flipsPerRound, mii& intersectionsWithTarget, int cntTtoT1);

int flipToReduceIntesectionsBiDir(facesDS& t, facesDS& t1, size_t n, const vP& points, std::vector<vi2>& flipsPerRound);

int flipToReduceIntesections(facesDS& t, facesDS& t1, size_t n, const vP& points, std::vector<vi2>& flipsPerRound, int dir);

int flipToReduceIntesectionsRevDir(facesDS& t, facesDS& t1, size_t n, const vP& points, std::vector<vi2>& flipsPerRound, mii& intersectionsWithTarget);

//O(n) memory, slower but less memory TOO SLOW and you don't save that much memory
// int flipToReduceIntesections(flippingDS& t, flippingDS& t1, size_t n, const vP& points, std::vector<vi2>& flipsPerRound);

//faster but requires O(n^2) memory
int flipToReduceIntesections(flippingDS& t, flippingDS& t1, size_t n, const vP& points, std::vector<vi2>& flipsPerRound, vi& intersectionsWithTarget);


std::vector<simpleDS> flipToReduceIntesectionsWithHistory(flippingDS& t, flippingDS& t1, size_t n, const vP& points, std::vector<vi2>& flipsPerRound, mii& intersectionsWithTarget);

std::vector<facesDS> flipToReduceIntesectionsWithHistory(facesDS& t, facesDS& t1, size_t n, const vP& points, std::vector<vi2>& flipsPerRound, mii& intersectionsWithTarget);

int improveCenter(facesDS& center, std::vector<facesDS>& tVec, size_t n, const vP& points, std::vector<vi2>& flipsPerRound);

//flip all triangulations at once
//at each step, do the following in this order:
//1) if possible, flip center to reduce its total number of intersections
//2) for all other triangulations, do all flips that reduce intersections with the center from 1 to 0, ie introduce a diagonal present in the center
//3) for all other triangulations, do all flips that reduce intersections with the center
int omniDirFlips(int centerIdx, std::vector<facesDS>& tVec, size_t n, const vP& points, std::vector<std::vector<vi2>>& flipsPerRound);