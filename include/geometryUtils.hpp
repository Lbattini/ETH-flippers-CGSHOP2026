#ifndef GEOMETRYUTILS_HPP
#define GEOMETRYUTILS_HPP

#endif //GEOMETRYUTILS_HPP

#include <vector>

//CGAL includes
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Triangulation_2.h>

#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/Triangulation_data_structure_2.h>

#include <CGAL/Constrained_triangulation_2.h>

// #include <CGAL/Kernel/global_functions.h>
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

// TODO This version is only used for my own random input generation (which isn't compatible with all the rest, so remove this)

bool canBeFlipped(const Triangulation& t, Triangulation::Finite_edges_iterator e);

bool canBeFlipped(const ConstrainedTriangulation& t, ConstrainedTriangulation::Finite_edges_iterator e) ;

//checks if the diagonal of the quadrilateral can be flipped
//the points of the quadrilateral can be passed either in clockwise or counterclockwise order
bool canBeFlipped(const std::array<P,4>& quadrilateral);


bool moreThanKInside(const ConstrainedTriangulation& t, ConstrainedTriangulation::Face_handle& face, int k, bool print=false);

void printEdges(const Triangulation& t);

void printEdges(const ConstrainedTriangulation& t);

//check if the edges of two triangulations are all equal
bool sameEdges(const ConstrainedTriangulation& t0, const ConstrainedTriangulation& t1);