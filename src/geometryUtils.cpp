#include "../include/geometryUtils.hpp"

#include <CGAL/convexity_check_2.h>

bool canBeFlipped(const Triangulation& t, Triangulation::Finite_edges_iterator e) {
  //edges are represented as std::pair<Face_handle,int> (face_handle,edge_index)
  bool flag=false;
  auto f=e->first;
  int j=e->second;
  auto neighborJ=f->neighbor(j);
  if (!t.is_infinite(f) && !t.is_infinite(neighborJ)) {
    //TODO according to docs is_ccw_strongly_convex takes beyond as a parameter https://doc.cgal.org/latest/Convex_hull_2/group__PkgConvexHull2Convexity.html
    //TODO but github implementation (same version) takes last instead https://github.com/CGAL/cgal/blob/v6.0.1/Convex_hull_2/include/CGAL/convexity_check_2.h
    std::array<P,4> quadrilateral={ f->vertex(j)->point(),f->vertex((j+1)%3)->point(),t.mirror_vertex(f,j)->point(),f->vertex((j+2)%3)->point()};//,P(0,0)};

    if (CGAL::is_ccw_strongly_convex_2(quadrilateral.begin(),quadrilateral.end())) {
      flag=true;
    }
  }
  // std::cout<<flag<<'\n';
  return flag;
}

bool canBeFlipped(const ConstrainedTriangulation& t, ConstrainedTriangulation::Finite_edges_iterator e) {
  //edges are represented as std::pair<Face_handle,int> (face_handle,edge_index)
  bool flag=false;
  auto f=e->first;
  int j=e->second;
  auto neighborJ=f->neighbor(j);
  if (!t.is_infinite(f) && !t.is_infinite(neighborJ)) {
    //TODO according to docs is_ccw_strongly_convex takes beyond as a parameter https://doc.cgal.org/latest/Convex_hull_2/group__PkgConvexHull2Convexity.html
    //TODO but github implementation (same version) takes last instead https://github.com/CGAL/cgal/blob/v6.0.1/Convex_hull_2/include/CGAL/convexity_check_2.h
    std::array<P,4> quadrilateral={ f->vertex(j)->point(),f->vertex((j+1)%3)->point(),t.mirror_vertex(f,j)->point(),f->vertex((j+2)%3)->point()};//,P(0,0)};

    if (!CGAL::collinear(quadrilateral[1],quadrilateral[2],quadrilateral[3]) && CGAL::is_ccw_strongly_convex_2(quadrilateral.begin(),quadrilateral.end())) {
      flag=true;
    }
  }
  return flag;
}

bool canBeFlipped(const std::array<P,4>& quadrilateral) {
  //TODO according to docs is_ccw_strongly_convex takes beyond as a parameter https://doc.cgal.org/latest/Convex_hull_2/group__PkgConvexHull2Convexity.html
  //TODO but github implementation (same version) takes last instead https://github.com/CGAL/cgal/blob/v6.0.1/Convex_hull_2/include/CGAL/convexity_check_2.h
  if (!CGAL::collinear(quadrilateral[1],quadrilateral[2],quadrilateral[3]) &&
    (CGAL::is_ccw_strongly_convex_2(quadrilateral.begin(),quadrilateral.end()) ||
      CGAL::is_cw_strongly_convex_2(quadrilateral.begin(),quadrilateral.end()))
      )
    return true;

  return false;
}

//TODO for k>0, this is a very bad/slow initial implementation, improve it by using a BFS
bool moreThanKInside(const ConstrainedTriangulation& t, ConstrainedTriangulation::Face_handle& face, int k, bool print) {
  int cnt=0;
  P p0=face->vertex(0)->point(), p1= face->vertex(1)->point(), p2=face->vertex(2)->point();
  for (int i=0;i<3 && cnt<=k;i++) {
    if (!t.is_infinite(t.mirror_vertex(face,i)) && CGAL::side_of_bounded_circle(p0,p1,p2,t.mirror_vertex(face,i)->point())==CGAL::ON_BOUNDED_SIDE) {
      cnt++;
      if (print)
        std::cout<<t.mirror_vertex(face,i)->point()<<'\n';
    }
  }
  if (k>0 && cnt<=k) {
    cnt=0;
    for (auto v=t.finite_vertices_begin(); v!=t.finite_vertices_end() && cnt<=k; ++v) {
      if (/*v->point()!=p0 && v->point()!=p1 && v->point()!=p2 &&*/ CGAL::side_of_bounded_circle(p0,p1,p2,v->point())==CGAL::ON_BOUNDED_SIDE)
        cnt++;
    }
  }

  return cnt>k;
}

void printEdges(const Triangulation& t) {
  for (auto e = t.finite_edges_begin(); e != t.finite_edges_end(); ++e) {
    // get the vertices of e
    Triangulation::Vertex_handle v1 = e->first->vertex((e->second + 1) % 3);
    Triangulation::Vertex_handle v2 = e->first->vertex((e->second + 2) % 3);
    std::cout << "e = " << v1->point() << " <-> " << v2->point() << std::endl;
  }
}

void printEdges(const ConstrainedTriangulation& t) {
  for (auto e = t.finite_edges_begin(); e != t.finite_edges_end(); ++e) {
    // get the vertices of e
    ConstrainedTriangulation::Vertex_handle v1 = e->first->vertex((e->second + 1) % 3);
    ConstrainedTriangulation::Vertex_handle v2 = e->first->vertex((e->second + 2) % 3);
    std::cout << "e = " << v1->point() << " <-> " << v2->point() << std::endl;
  }
}

bool sameEdges(const ConstrainedTriangulation& t0, const ConstrainedTriangulation& t1) {
  bool allSame=true;
  for (auto e0 = t0.finite_edges_begin(); e0 != t0.finite_edges_end() && allSame; ++e0) {
    bool equal=false;
    P p0=e0->first->vertex((e0->second + 1) % 3)->point(),
    p1=e0->first->vertex((e0->second + 2) % 3)->point();
    for (auto e1 = t1.finite_edges_begin(); e1 != t1.finite_edges_end() && !equal; ++e1) {
      P p2=e1->first->vertex((e1->second + 1) % 3)->point(),
      p3=e1->first->vertex((e1->second + 2) % 3)->point();
      if ((p0==p2 && p1==p3) || (p0==p3 && p1==p2))
        equal=true;
    }
    allSame=equal;
  }
  return allSame;
}