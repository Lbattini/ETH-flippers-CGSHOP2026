#ifndef QUADRILATERALGRAPH_HPP
#define QUADRILATERALGRAPH_HPP

#include<flippingDS.hpp>

typedef std::vector<int> vi;
typedef std::vector<vi> vvi;
typedef std::vector<bool> vb;
typedef std::vector<vb> vvb;
typedef std::vector<vvb> vvvb;

//for each edge, store all the edges with a larger index that it intersects
vvi crossingEdges(vP points);

//build the quadrilateral graph corresponding to a set of points
vvi buildQG(vP points);

//build the quadrilateral graph corresponding to a set of points PLUS empty non-convex quadrilaterals
//(stored with negative sign)
//ALSO compute and pass by reference the empty triangles
vvi buildQGPlusNonConvex(vP points, const std::string& fileName, vvvb& emptyTriangle);

//compute the maximum among the shortest path lengths between a diagonal of t0 and any diagonal of t1
//inside the quadrilateral graph
//todo this can be improved to min max perfect matching (a perfect matching that minimizes the maximum weight)
/* todo also, if you need to call this multiple times with different t0 and t1 in the same point set,
  * it would be worth it to precompute all shortest distances (APSP)*/
int longestSP( flippingDS& t0,  flippingDS& t1, vvi& qg, int nPoints);
#endif //QUADRILATERALGRAPH_HPP
