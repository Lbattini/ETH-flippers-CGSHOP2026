//The functions here were written before the examples were released, in part just to get more used to CGAL, in part
//to be able to generate more examples if needed

#ifndef GENRANDOMTRIANGULATIONS_HPP
#define GENRANDOMTRIANGULATIONS_HPP

#endif //GENRANDOMTRIANGULATIONS_HPP

#include <vector>
#include "../include/inOut.hpp"

//CGAL includes
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Triangulation_2.h>

#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/Triangulation_data_structure_2.h>

#include <CGAL/Constrained_triangulation_2.h>

// #include <CGAL/Exact_predicates_exact_constructions_kernel.h>
// typedef CGAL::Exact_predicates_exact_constructions_kernel EK;

/*
 * TODO For now, I'll try using EPIK, it should suffice for what I have to do. Check that input coordinates are not too
 * big (they shouldn't be)
 * If needed, I'll find and replace
 */
// From algolab examples: Epic kernel is enough, no constructions needed, provided the squared distance
// fits into a double (!)
typedef CGAL::Exact_predicates_inexact_constructions_kernel IK;

// we want to store an index with each vertex
typedef std::size_t                                            Index;
typedef CGAL::Triangulation_vertex_base_with_info_2<Index,IK>   Vb;
typedef CGAL::Triangulation_face_base_2<IK>                     Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>            Tds;

// As edges are not explicitly represented in the triangulation, we extract them
// from the triangulation to be able to sort and process them. We store the
// indices of the two endpoints, first the smaller, second the larger, and third
// the squared length of the edge. The i-th entry, for i=0,... of a tuple t can
// be accessed using std::get<i>(t).
//TODO is the squared length of the edge really necessary for my purposes? It shouldn't, but removing should be easy

/*edge index for edge (i,j), with i<j, is i*n+j.
 *this is not an index among the edges of the current triangulation, but rather an index
 *among all possible edges for the point set
 */
typedef std::tuple<Index,Index,IK::FT> Edge;
typedef std::vector<Edge> EdgeV;

typedef IK::Point_2 P;
typedef std::pair<IK::Point_2,Index> IPoint;

typedef CGAL::Triangulation_2<IK,Tds> Triangulation;

/*TODO (note to myself) CGAL::No_constraint_intersection_requiring_constructions_tag should suffice
 *since the triangulation should be valid. If intersections exception came out anyway, consider switching to Itag
 */
typedef CGAL::Exact_predicates_tag                               Itag;
typedef CGAL::No_constraint_intersection_requiring_constructions_tag NoIntersectionTag;
typedef CGAL::Constrained_triangulation_2<IK,CGAL::Default,NoIntersectionTag> ConstrainedTriangulation;

//For initial testing, something similar will probably be inside the utils we'll be given
std::vector<IPoint> genRandomPoints(int n, double minCoord, double maxCoord);

//genarate random points with integer coordinates (for debugging/testing)
std::vector<IPoint> genRandomPoints(int n, int minCoord, int maxCoord);

//TODO switch to ConstrainedTriangulation version (for compatibility with the rest of the code)
//This is also for initial testing. Scans all the edges and
//if the edge can be flipped, sample uar in [0,1], if >=(1-flipProbability) flips the edge
//This is repeated for nIterations
void randomFlips(Triangulation& t, int nIterations, double flipProbability);
