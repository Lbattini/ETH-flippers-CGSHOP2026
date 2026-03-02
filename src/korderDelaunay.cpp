/*
Useful functions:
- NOT THIS
side_of_oriented_circle, returns ON_POSITIVE_SIDE (bounded side,
the circle is assumed to be counterclockwise oriented), ON_NEGATIVE_SIDE or,
ON_ORIENTED_BOUNDARY
BUT INSTEAD side_of_bounded_circle Precondition p, q and r are not collinear.
For the first implementation, I plan to test all vertices to count how many are inside
(at least for k>0). Later on, I can just explore layer by layer starting from
the current vertex (if in a certain layer, ie bfs distance in the graph that has
faces as nodes, I find 0 vertices there is no need to go further).

- flip, you need to check that The faces f and f->neighbor(i) are
finite faces and their union forms a convex quadrilateral by using is_ccw_strongly_convex_2().

To generate a random triangulation (no guarantee to be uar), start from scan
triangulation; then do a certain number of random edge flips.

CGAL constrained triangulation+map from point to index should contain all the information I need.
But for efficiency, it might be better to maintain a separate data structure to store which edges can be flipped.
If that's the case, refer to the A* repo for a possible way of doing that.

*/

//My local installation uses CGAL 6.0.1, so that's what this code is tested with

//WARNING: behavior of Lawson flip and its variations might change from run to run, since the CGAL internal ordering of edges
//(and therefore the order in which the flips are performed) might change. This is particularly relevant for greedy parallel flips

#include <iostream>
#include <algorithm>
#include "../include/korderDelaunay.hpp"
#include "../include/geometryUtils.hpp"

//Simple function for the scan triangulation that sorts the points (lexicographically) and then inserts them
Triangulation scanTriangulation(std::vector<IPoint> &points) {
  Triangulation triangulation;
  std::sort(points.begin(),points.end());
  triangulation.insert(points.begin(),points.end());
  return triangulation;
}

bool isKOrderDelaunay(ConstrainedTriangulation& t, int k, size_t n, const pointToIdx& idxMap, bool print=false) {
  bool flag=true;
  for (ConstrainedTriangulation::Face_iterator f=t.finite_faces_begin();f!=t.finite_faces_end();++f) {
    ConstrainedTriangulation::Face_handle face=f;
    if (!t.is_infinite(f) && moreThanKInside(t, face, k)) {
      if (print) {
        P p0=f->vertex(0)->point(), p1=f->vertex(1)->point(),
        p2=f->vertex(2)->point();
        // int a=idxMap.at(p0), b=idxMap.at(p1);
        // if (a>b) std::swap(a,b);
        std::cout<<"Triangle with points inside: \n";
        std::cout<<p0<<std::endl;
        std::cout<<p1<<std::endl;
        std::cout<<p2<<std::endl;
        std::cout<<"____\n";
      }
      flag=false;
    }
  }
  return flag;
}

// bool shareATriangle(const ConstrainedTriangulation::Finite_edges_iterator& e0, const ConstrainedTriangulation::Finite_edges_iterator& e1) {
//   bool flag=false;
//   if (e1->first==e0->first || e1->first==e0->first->neighbor(e0->second)) {
//
//   }
//   return flag;
// }

void printEdge(const ConstrainedTriangulation::Finite_edges_iterator& e, const pointToIdx& idxMap) {
  ConstrainedTriangulation::Face_handle currentFace=e->first;
  P p0=currentFace->vertex((e->second+1)%3)->point(), p1=currentFace->vertex((e->second+2)%3)->point();
  int a=idxMap.at(p0), b=idxMap.at(p1);
  if (a>b) std::swap(a,b);
  std::cout<<a<<' '<<b<<std::endl;
}

//Lawson flip algorithm with parallel flips
//TODO deal with collinear/cocircular points (at least print a warning)
int lawsonFlipParallel(ConstrainedTriangulation& t, size_t n, const pointToIdx& idxMap) {
  int F=0;
  //stop when no edge has been flipped in the previous round
  bool toFlip=true;

  int curK=0;
  int totalFlips=0;
  while (toFlip && totalFlips<n*n) {
    toFlip=false;
    std::vector<ConstrainedTriangulation::Finite_edges_iterator> delEdges;
    for (auto e = t.finite_edges_begin(); e != t.finite_edges_end(); ++e) {
      if (!t.is_infinite(e->first) && !t.is_infinite(e->first->neighbor(e->second)) && moreThanKInside(t, e->first, curK)) {
        if (canBeFlipped(t,e)) {
          delEdges.push_back(e);
          toFlip=true;
        }
      }
    }

    //2 means don't flip because an edge that shares a triangle with it was flipped
    std::vector<int> flipped(delEdges.size(),0);
    int flipsThisRound=0;
    for (int i=0;i<delEdges.size();i++) {
      auto e=delEdges[i];
      if (flipped[i]==0) {
        //mark the edges that share a triangle with e
        for (int j=i+1;j<delEdges.size();j++) {
          auto e1=delEdges[j];
          if (e1->first==e->first || e1->first==e->first->neighbor(e->second)
            || e1->first->neighbor(e1->second)==e->first || e1->first->neighbor(e1->second)==e->first->neighbor(e->second)) {

            flipped[j]=2;
          }
        }

        //flip e
        // printEdge(e,idxMap);
        t.flip(e->first,e->second);
        totalFlips++;
        flipsThisRound++;
      }
    }

    // if (delEdges.size()>1)
    //   std::cout<<delEdges.size()<<' '<<flipsThisRound<<std::endl;
    // std::cout<<"____\n";

    if (toFlip)
      F++;
  }
  bool validDelaunay=isKOrderDelaunay(t,0,n,idxMap,true);
  // assert(validDelaunay);
  if (!validDelaunay)
    std::cout<<"Lawson failed, check for collinear or circular points\n";
  std::cout<<totalFlips<<' '<<F<<std::endl;
  return F;
}

//Lawson flip algorithm
//TODO deal with collinear/cocircular points (at least print a warning)
int lawsonFlip(ConstrainedTriangulation& t, size_t n, const pointToIdx& idxMap) {
  int F=0;
  //stop when no edge has been flipped in the previous round
  bool toFlip=true;
  int curK=0;
  while (toFlip && F<n*n) {
    toFlip=false;
    for (auto e = t.finite_edges_begin(); e != t.finite_edges_end() && !toFlip; ++e) {
      if (!t.is_infinite(e->first) && !t.is_infinite(e->first->neighbor(e->second)) && moreThanKInside(t, e->first, curK)) {
        if (canBeFlipped(t,e)) {
          t.flip(e->first,e->second);
          F++;
          toFlip=true;
        }
      }
    }
  }
  bool validDelaunay=isKOrderDelaunay(t,0,n,idxMap,true);
  // assert(validDelaunay);
  if (!validDelaunay)
    std::cout<<"Lawson failed, check for collinear or circular points\n";

  return F;
}

//k-order Delaunay using Lawson Flip algorithm with parallel flips
//note this specific version has no theoretical analysis behind it
//Known results are:
//- O(n) parallel flips are enough to turn a triangulation into any other (while O(n^2) regular flips might be necessary)
//- Lawson's flip algorithm with regular flips takes O(n^2)
//- exactly (n-4)/5 edges can be flipped in parallel. But not all of these are needed for a Delaunay triangulation
//- O(nk) edges are possible in a k order triangulation (check hypothesis)
//- k order triangulations are very close for constant k<=7 (check precise statement)
//What I would like to empirically study (with no pretense of being rigorous) is:
//- does it take fewer parallel flips to reach a k order triangulation, for constant 0<k<=7, compared to a standard Delaunay triangulation (k=0)?
//Limits of planned implementation:
//- handling collinear and cocircular points;
//- the set of parallel edge flips at each step might not be of maximal cardinality
//- is choosing a set of Delaunay edge flips of maximal cardinality for each parallel flip optimal?

/*
 *
 */
//TODO add parallel flips
//TODO Lawson might not work for k>0... Use point to idx map to get edge idx, and avoid flipping the same edge twice
int korderDelaunay(ConstrainedTriangulation& t, int k, size_t n, const pointToIdx& idxMap) {
  if (k==0)
    return lawsonFlipParallel(t,n,idxMap);
  int F=0;
  bool toFlip=true;
  //only flipping edges when >k points inside might not be enough for k>0
  //as a workaround, decrease curK when all those edges have been flipped once
  int curK=k;
  bool printStuff=false;
  //flip each edge only once
  std::vector<int> flipped(n*n,0);
  while (toFlip && curK>=0) {
    // toFlip=false;
    bool noneFlipped=true,someHaveToBeFlipped=false;;
    for (auto e = t.finite_edges_begin(); e != t.finite_edges_end() && noneFlipped; ++e) {
      ConstrainedTriangulation::Face_handle currentFace=e->first;
      P p0=currentFace->vertex((e->second+1)%3)->point(), p1=currentFace->vertex((e->second+2)%3)->point();
      int a=idxMap.at(p0), b=idxMap.at(p1);
      if (a>b) std::swap(a,b);
      if (!t.is_infinite(e->first) && !t.is_infinite(e->first->neighbor(e->second)) && moreThanKInside(t, currentFace, 0)) {
        // toFlip=true;
        // someHaveToBeFlipped=true;
        if (/*flipped[a*n+b]<2 &&*/ canBeFlipped(t,e)) {
          t.flip(e->first,e->second);
          if (flipped[a*n+b]>1) {
            std::cout<<"___\n";
            std::cout<<F<<std::endl;
            std::cout<<a<<' '<<b<<std::endl;
            std::cout<<p0.x()<<' '<<p0.y()<<std::endl;
            std::cout<<p1.x()<<' '<<p1.y()<<std::endl;
            printStuff=true;
            // CGAL::draw(t);
          }
          flipped[a*n+b]++;
          F++;

          noneFlipped=false;
        }
      }
    }
    // if (printStuff) {
    //   std::cout<<"__________\n";
    //   CGAL::draw(t);
    // }
    for (auto e = t.finite_edges_begin(); e != t.finite_edges_end() && !someHaveToBeFlipped /*&& noneFlipped*/; ++e) {
      if (!t.is_infinite(e->first) && !t.is_infinite(e->first->neighbor(e->second)) && moreThanKInside(t, e->first, curK)) {
        someHaveToBeFlipped=true;
      }
    }
    // bool someCanBeFlipped=false,someHaveToBeFlipped=false;
    // for (auto e = t.finite_edges_begin(); e != t.finite_edges_end(); ++e) {
    //   ConstrainedTriangulation::Face_handle currentFace=e->first;
    //   P p0=currentFace->vertex((e->second+1)%3)->point(), p1=currentFace->vertex((e->second+2)%3)->point();
    //   int a=idxMap.at(p0), b=idxMap.at(p1);
    //   if (moreThanKInside(t, currentFace, curK)) {
    //     // std::cout<<a*n+b<<'\n';
    //     if (!flipped[a*n+b] && canBeFlipped(t,e)) {
    //       someCanBeFlipped=true;
    //     }
    //
    //     someHaveToBeFlipped=true;
    //   }
    // }
    if (noneFlipped && someHaveToBeFlipped) {
      curK--;
      std::cout<<curK<<std::endl;
    }
    toFlip=someHaveToBeFlipped;
    // std::cout<<F<<'\n';
  }
  // while (toFlip && F<n*n) {
  //   toFlip=false;
  //   for (auto e = t.finite_edges_begin(); e != t.finite_edges_end(); ++e) {
  //     ConstrainedTriangulation::Face_handle currentFace=e->first;
  //     if (canBeFlipped(t,e) && moreThanKInside(t, currentFace, k)) {
  //       t.flip(e->first,e->second);
  //       F++;
  //     }
  //   }
  //   if (k>0) {
  //     for (auto e = t.finite_edges_begin(); e != t.finite_edges_end() && !toFlip; ++e) {
  //       ConstrainedTriangulation::Face_handle currentFace=e->first;
  //       if (canBeFlipped(t,e) && moreThanKInside(t, currentFace, k)) {
  //         toFlip=true;
  //       }
  //     }
  //     if (toFlip) {
  //       for (auto e = t.finite_edges_begin(); e != t.finite_edges_end(); ++e) {
  //         ConstrainedTriangulation::Face_handle currentFace=e->first;
  //         if (canBeFlipped(t,e)) {
  //           t.flip(e->first,e->second);
  //           F++;
  //           if (moreThanKInside(t, currentFace, k))
  //             toFlip=true;
  //           // std::cout<<F<<'\n';
  //         }
  //       }
  //     }
  //   }
  // }
  // if (F==n*n)
  //   std::cout<<"Lawson failed\n";
  if (!isKOrderDelaunay(t,k,n,idxMap))
    std::cout<<"Failed to turn into k order Delaunay"<<std::endl;
  return F;
}
