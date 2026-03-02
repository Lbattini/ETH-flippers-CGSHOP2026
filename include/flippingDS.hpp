#ifndef FLIPPINGDS_HPP
#define FLIPPINGDS_HPP

#include <vector>
#include <array>
#include <unordered_map>

//CGAL includes
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Triangulation_2.h>

#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/Triangulation_data_structure_2.h>

#include <CGAL/Constrained_triangulation_2.h>

#include <CGAL/draw_triangulation_2.h>
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

typedef std::unordered_map<P,int> pointToIdx;
typedef std::vector<P> vP;
typedef std::array<int,2> i2;
typedef std::vector<i2> vi2;
typedef std::vector<ConstrainedTriangulation> vConstrainedTriangulation;



class flippingDS {

public:
  flippingDS()= default;

  virtual ~flippingDS()=default;

  //Return the vector with the current edges in the triangulation. Each edge (u,v) is such that u<v
  //This only needs to be called once at the beginning, then the callee can keep track of the updates locally.
  //To do so, before flipping edge i, store the result of quadrilateral(i), and then put that value into edges[i]
  virtual vi2 getEdges() const=0;

  //flip the edge with index i
  virtual void flip(int i)=0;

  //return the indices of the other two vertices of the quadrilateral with edge i as a diagonal
  //no guarantee on the order, so you need to check both is cw strongly convex and is ccw strongly convex
  //returns {-1,-1} if infinite (ie the edge belongs to the convex hull)
  virtual i2 otherDiagonal(int i)=0;

  virtual std::array<int,4> quadrilateralEdges(int i)=0;

  //check if edges i and j share a triangle
  virtual bool shareTriangle(int i, int j)=0;

  //must be called in the innermost method that needs the auxiliary data structure;
  //for some reason I cannot do some things in the constructor (something related to CGAL pointers changing when passing something, idk)
  virtual void init()=0;

  [[nodiscard]] bool equal(const flippingDS& other) const{
    vi2 edges0=this->getEdges();
    vi2 edges1=other.getEdges();
    bool allEqual=true;
    if (edges0.size()!=edges1.size())
      allEqual=false;
    else {
      for (int i=0;i<edges0.size();i++) {
        std::sort(edges0[i].begin(),edges0[i].end());
        std::sort(edges1[i].begin(),edges1[i].end());
      }
      std::sort(edges0.begin(),edges0.end());
      std::sort(edges1.begin(),edges1.end());
      for (int i=0;i<edges0.size() && allEqual;i++)
        for (int j=0;j<2;j++)
          if (edges0[i][j]!=edges1[i][j]) {
            std::cout<<edges0[i][0]<<' '<<edges0[i][1]<<'\n';
            std::cout<<edges1[i][0]<<' '<<edges1[i][1]<<'\n';
            allEqual=false;
          }
    }
    return allEqual;
  }

  [[nodiscard]] int cntIntersections(const flippingDS& other) const{
    vi2 edges0=this->getEdges();
    vi2 edges1=other.getEdges();
    int count=0;
    if (edges0.size()!=edges1.size())
      count=10000000;
    else {
      for (int i=0;i<edges0.size();i++) {
        std::sort(edges0[i].begin(),edges0[i].end());
        std::sort(edges1[i].begin(),edges1[i].end());
      }
      std::sort(edges0.begin(),edges0.end());
      std::sort(edges1.begin(),edges1.end());
      for (int i=0;i<edges0.size();i++)
        for (int j=0;j<2;j++)
          if (edges0[i][j]!=edges1[i][j]) {
            // std::cout<<edges0[i][0]<<' '<<edges0[i][1]<<'\n';
            // std::cout<<edges1[i][0]<<' '<<edges1[i][1]<<'\n';
            count++;
            break;
          }
    }
    return count;
  }

};

class CGALconstrainedTriangulation : public flippingDS{
  ConstrainedTriangulation t;
  pointToIdx idxMap;
  //store the ConstrainedTriangulation::Edge corresponding to each edge index
  std::vector<ConstrainedTriangulation::Edge> edgesTDS;
  vi2 edges;
  //store the (up to 4) edges that share a triangle with the edge (u,v)
  // std::vector<std::array<int,4>> shareTriangleVec;
  std::map<i2,int> edgeIdxMap;

  //store the idx of each with endpoints (u,v)
  // std::vector<std::vector<int>> edgeIdx;

public:
  CGALconstrainedTriangulation(ConstrainedTriangulation& t, pointToIdx& idxMap) {
    this->t = t;
    this->idxMap = idxMap;
    // this->CGALconstrainedTriangulation::init();
  }

  ~CGALconstrainedTriangulation() override =default;

  void init() override{
    if (!edges.empty()) {
      edges.clear();
      edgesTDS.clear();
      edgeIdxMap.clear();
    }
    for (auto e : t.finite_edges()) {
      P p0=e.first->vertex((e.second+1)%3)->point(), p1=e.first->vertex((e.second+2)%3)->point();
      edges.push_back({std::min(idxMap[p0],idxMap[p1]),std::max(idxMap[p0],idxMap[p1])});
      edgeIdxMap[edges[edges.size()-1]]=edges.size()-1;
      edgesTDS.push_back(e);
    }
  }

  vi2 getEdges() const override {
    if (edges.empty())
      std::cout<<"something went wrong\n";
    // if (edges.empty()) {
    //   initDS();
    // }

    return edges;
  }

  //todo this is currently O(log n) because of the edgeIdxMap it uses. It can be replaced either by a O(n^2) memory matrix or an undordered map

  void flip(int i) override {
    assert(!edges.empty());
    auto neigh2=edgesTDS[i].first->neighbor((edgesTDS[i].second+2)%3);
    int j=t.mirror_index(edgesTDS[i].first,(edgesTDS[i].second+2)%3);
    i2 quadr=otherDiagonal(i);
    if (quadr[0]>quadr[1]) std::swap(quadr[0],quadr[1]);

    std::vector<int> shareTriangleVec;
    for (int l=0;l<2;l++)
      for (int m=0;m<2;m++) {
        i2 e={edges[i][l],quadr[m]};
        if (e[0]>e[1])
          std::swap(e[0],e[1]);

        if (edgeIdxMap.count(e)>0)
          shareTriangleVec.push_back(edgeIdxMap[e]);
      }

    t.flip(edgesTDS[i].first,edgesTDS[i].second);
    // CGAL::draw(t);

    bool found=false;

    for (int k=0;k<3 && !found;k++) {
      int a=idxMap[neigh2->neighbor(j)->vertex((k+1)%3)->point()], b=idxMap[neigh2->neighbor(j)->vertex((k+2)%3)->point()];
      if (a>b) std::swap(a,b);
      if (a==quadr[0] && b==quadr[1]) {
        edgesTDS[i]={neigh2->neighbor(j),k};
        found=true;
      }
    }

    edges[i]=quadr;
    edgeIdxMap[edges[i]]=i;
    int cnt=0;

    for (int k : shareTriangleVec){
      bool found1=false;
      for (int l=0;l<3 && !found1;l++) {
        int a=idxMap[edgesTDS[i].first->vertex((l+1)%3)->point()], b=idxMap[edgesTDS[i].first->vertex((l+2)%3)->point()];
        if (a>b) std::swap(a,b);
        if (a==edges[k][0] && b==edges[k][1]) {
          edgesTDS[k]={edgesTDS[i].first,l};
          found1=true;
          cnt++;
        }
      }
      for (int l=0;l<3 && !found1;l++) {
        int a=idxMap[edgesTDS[i].first->neighbor(edgesTDS[i].second)->vertex((l+1)%3)->point()],
        b=idxMap[edgesTDS[i].first->neighbor(edgesTDS[i].second)->vertex((l+2)%3)->point()];
        if (a>b) std::swap(a,b);
        if (a==edges[k][0] && b==edges[k][1]) {
          edgesTDS[k]={edgesTDS[i].first->neighbor(edgesTDS[i].second),l};
          found1=true;
          cnt++;
        }
      }
    }

    assert(cnt==4);

    //if this fails, it means that the neighboring face also changes for some reason
    assert(found);
  }

  i2 otherDiagonal(int i) override {
    assert(!edges.empty());
    if (t.is_infinite(edgesTDS[i].first))
      return {-1,-1};
    if (t.is_infinite(edgesTDS[i].first) || t.is_infinite(edgesTDS[i].first->neighbor(edgesTDS[i].second)))
      return {-1,-1};
    auto p0=edgesTDS[i].first->vertex(edgesTDS[i].second),
    p1=t.mirror_vertex(edgesTDS[i].first,edgesTDS[i].second);

    if (t.is_infinite(p0) || t.is_infinite(p1))
      return {-1,-1};

    return {idxMap[p0->point()],idxMap[p1->point()]};
  }

  bool shareTriangle(int i, int j) override {
    assert(!edges.empty());
    auto e=edgesTDS[i],e1=edgesTDS[j];
    return e1.first==e.first || e1.first==e.first->neighbor(e.second)
            || e1.first->neighbor(e1.second)==e.first || e1.first->neighbor(e1.second)==e.first->neighbor(e.second);
  }

  std::array<int,4> quadrilateralEdges(int i) override {
    std::cout<<"Not implemented";
    return {-1,-1,-1,-1};
  }

};

class simpleDS: public flippingDS {
  vi2 edges;
  //store for edge i the indices (u,v) of the other diagonal, or {-1,-1} if that diagonal doesn't exist
  vi2 otherDiag;
  //store for each edge i the 4 edges in the quadrilateral that has i as a diagonal (-1 if the edge is in the convex hull)
  //invariant: it's (firstEdgeFirstFace, secondEdgeFirstFace, firstEdgeSecondFace, secondEdgeSecondFace)
  //In other words: the edges of the quadrilateral that share a triangle are kept adjacent
  std::vector<std::array<int,4>> quadrilaterals;


public:

  bool compatible(int j) {
    bool res=true;
    if (quadrilaterals[j][0]!=-1)
      for (int e : quadrilaterals[j])
        for (int ii=0;ii<2;ii++)
          if (edges[e][ii]!=edges[j][0] && edges[e][ii]!=edges[j][1] &&
            edges[e][ii]!=otherDiag[j][0] && edges[e][ii]!=otherDiag[j][1]) {
            std::cout<<edges[j][0]<<' '<<edges[j][1]<<' '<<otherDiag[j][0]<<' '<<otherDiag[j][1]<<'\n';
            std::cout<<e<<' '<<edges[e][0]<<' '<<edges[e][1]<<'\n';
            res=false;
          }
    return res;
  }

  explicit simpleDS(flippingDS& t1) {
    edges=t1.getEdges();
    otherDiag.resize(edges.size());
    quadrilaterals.resize(edges.size());
    for (int i=0;i<edges.size();i++) {
      otherDiag[i]=t1.otherDiagonal(i);
      quadrilaterals[i]=t1.quadrilateralEdges(i);
    }
  }
  simpleDS(ConstrainedTriangulation& t, pointToIdx& idxMap) {
    int cnt=0;
    for (auto f=t.finite_faces_begin();f!=t.finite_faces_end();++f)
      cnt++;
    std::vector<std::array<int,7>> edgesWithInfo(cnt*3);
    int currentEdgeIdx=0;
    //store edges face by face: since this a triangulation, each edge not in the convex hull will be inserted twice
    for (auto f=t.finite_faces_begin();f!=t.finite_faces_end();++f) {
      for (int i=0;i<3;i++) {
        int a=idxMap[f->vertex((i+1)%3)->point()],b=idxMap[f->vertex((i+2)%3)->point()];
        if (a>b) std::swap(a,b);
        std::array<int,7> tmp={-1};
        tmp[0]=a;
        tmp[1]=b;

        auto v0=f->vertex(i),
        v1=t.mirror_vertex(f,i);

        if (t.is_infinite(v0) || t.is_infinite(v1)) {
          tmp[2]=-1;
          tmp[3]=-1;
        }
        else {
          tmp[2]=idxMap[v0->point()];
          tmp[3]=idxMap[v1->point()];
          if (tmp[2]>tmp[3])
            std::swap(tmp[2],tmp[3]);
          // if (tmp[2]==tmp[3])
          //   std::cout<<v0->point()<<' '<<v1->point()<<'\n';
          assert(tmp[2]!=tmp[3]);
        }

        tmp[4]=currentEdgeIdx+(i+1)%3;
        tmp[5]=currentEdgeIdx+(i+2)%3;
        tmp[6]=currentEdgeIdx+i;
        edgesWithInfo[currentEdgeIdx+i]=tmp;
      }
      currentEdgeIdx+=3;
    }

    std::sort(edgesWithInfo.begin(),edgesWithInfo.end());
    std::vector<int> newIdx(edgesWithInfo.size());
    int cntEdges=0;
    for (int i=0;i<edgesWithInfo.size();i++) {
      //check that the edge is stored twice, which means it is not in the boundary of the convex hull
      newIdx[edgesWithInfo[i][6]]=cntEdges;
      if (i+1<edgesWithInfo.size() && edgesWithInfo[i][0]==edgesWithInfo[i+1][0] && edgesWithInfo[i][1]==edgesWithInfo[i+1][1]) {
        i++;
        newIdx[edgesWithInfo[i][6]]=cntEdges;
      }
      cntEdges++;
    }
    edges.resize(cntEdges);
    otherDiag.resize(cntEdges);
    quadrilaterals.resize(cntEdges);
    int j=0;
    for (int i=0;i<edgesWithInfo.size();i++) {
      //check that the edge is stored twice, which means it is not in the boundary of the convex hull
      if (i+1<edgesWithInfo.size() && edgesWithInfo[i][0]==edgesWithInfo[i+1][0] && edgesWithInfo[i][1]==edgesWithInfo[i+1][1]) {
        edges[j]={edgesWithInfo[i][0],edgesWithInfo[i][1]};
        otherDiag[j]={edgesWithInfo[i][2],edgesWithInfo[i][3]};
        assert(otherDiag[j][0]<otherDiag[j][1]);
        quadrilaterals[j]={newIdx[edgesWithInfo[i][4]],newIdx[edgesWithInfo[i][5]],newIdx[edgesWithInfo[i+1][4]],newIdx[edgesWithInfo[i+1][5]]};
        j++;
        i++;
      }
      else {
        edges[j]={edgesWithInfo[i][0],edgesWithInfo[i][1]};
        otherDiag[j]={-1,-1};
        quadrilaterals[j]={-1,-1,-1,-1};
        j++;
      }
      assert(edges[j-1][0]<edges[j-1][1]);

    }
    for (int k=0;k<j;k++) {
      assert(compatible(k));
      assert(edges[k][0]!=edges[k][1]);
    }
  }

  //nothing to do here
  void init() override {
    ;
  }

  [[nodiscard]] vi2 getEdges() const override {
    if (edges.empty())
      std::cout<<"something went wrong\n";

    return edges;
  }

  //14 array values are changed every time this is called
  //On top of the assignments, up to 20 comparisons are performed
  void flip(int i) override {

    //update the other edges in the quadrilateral
    for (int j=0;j<4;j++) {
      int currentEdgeIdx=quadrilaterals[i][j];
      if (otherDiag[currentEdgeIdx][0]==-1)
        continue;

      //fix otherDiag
      //the new vertex on the diagonal belongs to otherDiag[i], but not to edges[currentEdgeIdx]
      int newDiagVertex=(otherDiag[i][0]==edges[currentEdgeIdx][0] || otherDiag[i][0]==edges[currentEdgeIdx][1]) ? otherDiag[i][1] : otherDiag[i][0];
      if (otherDiag[currentEdgeIdx][0]==edges[i][0] || otherDiag[currentEdgeIdx][0]==edges[i][1]) {
      // if (otherDiag[currentEdgeIdx][0]!=newDiagVertex){
        otherDiag[currentEdgeIdx][0]=newDiagVertex;
      }
      else {
        if (!(otherDiag[currentEdgeIdx][1]==edges[i][0] || otherDiag[currentEdgeIdx][1]==edges[i][1])) {
          std::cout<<edges[i][0]<<' '<<edges[i][1]<<'\n';
          std::cout<<otherDiag[i][0]<<' '<<otherDiag[i][1]<<'\n';
          std::cout<<"Adjacent quadrilateral:\n"<<edges[currentEdgeIdx][0]<<' '<<edges[currentEdgeIdx][1]<<'\n';
          std::cout<<otherDiag[currentEdgeIdx][0]<<' '<<otherDiag[currentEdgeIdx][1]<<'\n';
        }
        assert(otherDiag[currentEdgeIdx][1]==edges[i][0] || otherDiag[currentEdgeIdx][1]==edges[i][1]);
        otherDiag[currentEdgeIdx][1]=newDiagVertex;
      }
      assert(otherDiag[currentEdgeIdx][0]!=otherDiag[currentEdgeIdx][1]);

      //fix quadrilaterals
      bool found=false;
      for (int k=0;k<4 && !found;k++) {
        //see invariant: the edge we need to update is the one that is adjacent to e (before flipping)
        //if e is the k-th edge, then the one adjacent to it is at k+1 if k%2==0, k-1 otherwise
        if (quadrilaterals[currentEdgeIdx][k]==i) {
          int vertexToFind=(edges[currentEdgeIdx][0]==edges[i][0] || edges[currentEdgeIdx][0]==edges[i][1]) ? edges[currentEdgeIdx][0] : edges[currentEdgeIdx][1];
          //the new edge in the quadrilateral of currentEdgeIdx belongs to the quadrilateral of i;

          bool found2=false;
          for (int l=0;l<4;l++) {
            if (quadrilaterals[i][l]!=currentEdgeIdx &&
              (edges[quadrilaterals[i][l]][0]==vertexToFind ||
            edges[quadrilaterals[i][l]][1]==vertexToFind)) {
              quadrilaterals[currentEdgeIdx][k+1-2*(k%2)]=quadrilaterals[i][l];
              found2=true;
              break;
            }
          }
          assert(found2);

          /* todo: the new edge should be in a different triangle than currentEdgeIdx, but this doesn't work
           * try to think about why (the version above works but it's slower)
           */

          // if (edges[quadrilaterals[i][(j+2)%4]][0]==vertexToFind ||
          //   edges[quadrilaterals[i][(j+2)%4]][1]==vertexToFind)
          //   quadrilaterals[currentEdgeIdx][k+1-2*(k%2)]=quadrilaterals[i][(j+2)%4];
          // else
          //   quadrilaterals[currentEdgeIdx][k+1-2*(k%2)]=quadrilaterals[i][(j+3)%4];

          found=true;
        }
      }
      // assert(compatible(currentEdgeIdx));
    }

    i2 tmp=edges[i];
    edges[i]=otherDiag[i];
    assert(edges[i][0]!=edges[i][1]);
    otherDiag[i]=tmp;

    //other possibility: set to {-1,-1} to make it impossible to reverse an action
    // otherDiag[i]={-1,-1};

    //maintain the invariant
    if (edges[quadrilaterals[i][0]][0]==edges[quadrilaterals[i][2]][0] ||
      edges[quadrilaterals[i][0]][0]==edges[quadrilaterals[i][2]][1] ||
      edges[quadrilaterals[i][0]][1]==edges[quadrilaterals[i][2]][0] ||
      edges[quadrilaterals[i][0]][1]==edges[quadrilaterals[i][2]][1])
    std::swap(quadrilaterals[i][1],quadrilaterals[i][2]);
    else
      std::swap(quadrilaterals[i][1],quadrilaterals[i][3]);

    // debugCheck=true;
    // for (int e : quadrilaterals[i])
    //   for (int j=0;j<2;j++)
    //     if (edges[e][j]!=edges[i][0] && edges[e][j]!=edges[i][1] &&
    //       edges[e][j]!=otherDiag[i][0] && edges[e][j]!=otherDiag[i][1])
    //         debugCheck=false;
    // assert(debugCheck);
    assert(compatible(i));
    for (int j=0;j<4;j++) {
      int currentEdgeIdx=quadrilaterals[i][j];
      bool comp=compatible(currentEdgeIdx);
      if (!comp) {
        std::cout<<i<<": "<<edges[i][0]<<' '<<edges[i][1]<<'\n';
        for (int q : quadrilaterals[i])
          std::cout<<edges[q][0]<<' '<<edges[q][1]<<'\n';
      }
      assert(comp);
    }
  }

  ~simpleDS() override =default;

  i2 otherDiagonal(int i) override {
    //todo fix this properly inside the flip operation
    if (otherDiag[i][0]>otherDiag[i][1])
      return {otherDiag[i][1],otherDiag[i][0]};
    return otherDiag[i];
  }

  bool shareTriangle(int i, int j) override {
    //i and j share a triangle iff j belongs to the quadrilateral of i
    for (int k=0;k<4;k++)
      if (quadrilaterals[i][k]==j)
        return true;

    return false;
  }

  std::array<int,4> quadrilateralEdges(int i) override {
    return quadrilaterals[i];
  }
};

class facesDS: public flippingDS {
  typedef std::array<int,3> i3;
  typedef std::array<int,4> i4;
  typedef std::vector<i3> vi3;
  typedef std::vector<i4> vi4;

  vi2 edges;
  //store for edge i the indices (u,v) of the other diagonal, or {-1,-1} if that diagonal doesn't exist
  vi2 otherDiag;

  //store for each face the 3 edges it contains in ccw order
  vi3 faces;

  //store the <=2 faces each edge is incident to, and the idx of the edge along the face
  vi4 edgeIncidentFaces;


public:

  bool compatible(int j) {
    bool res=true;
    i4 quadrilateral=quadrilateralEdges(j);
    if (quadrilateral[0]!=-1)
      for (int e : quadrilateral)
        for (int ii=0;ii<2;ii++)
          if (edges[e][ii]!=edges[j][0] && edges[e][ii]!=edges[j][1] &&
            edges[e][ii]!=otherDiag[j][0] && edges[e][ii]!=otherDiag[j][1]) {
            std::cout<<edges[j][0]<<' '<<edges[j][1]<<' '<<otherDiag[j][0]<<' '<<otherDiag[j][1]<<'\n';
            std::cout<<e<<' '<<edges[e][0]<<' '<<edges[e][1]<<'\n';
            res=false;
          }
    return res;
  }

  // this does not compile and if I understand correctly the compiler should just do it on its own, no need
  // to do it myself
  // facesDS(facesDS& t1) {
  //   edges=t1.getEdges();
  //   otherDiag.resize(edges.size());
  //   for (int i=0;i<edges.size();i++) {
  //     otherDiag[i]=t1.otherDiagonal(i);
  //   }
  //   faces=t1.getFaces();
  //   edgeIncidentFaces=t1.getEdgeIncidentFaces();
  // }

  facesDS(ConstrainedTriangulation& t, pointToIdx& idxMap) {
    int cnt=0;
    for (auto f=t.finite_faces_begin();f!=t.finite_faces_end();++f)
      cnt++;
    std::vector<std::array<int,7>> edgesWithInfo(cnt*3);
    int currentEdgeIdx=0;
    int faceIdx=0;
    faces.resize(t.number_of_faces());
    //store edges face by face: since this a triangulation, each edge not in the convex hull will be inserted twice
    for (auto f=t.finite_faces_begin();f!=t.finite_faces_end();++f) {
      for (int i=0;i<3;i++) {
        int a=idxMap[f->vertex((i+1)%3)->point()],b=idxMap[f->vertex((i+2)%3)->point()];
        if (a>b) std::swap(a,b);
        std::array<int,7> tmp={-1};
        tmp[0]=a;
        tmp[1]=b;

        auto v0=f->vertex(i),
        v1=t.mirror_vertex(f,i);

        if (t.is_infinite(v0) || t.is_infinite(v1)) {
          tmp[2]=-1;
          tmp[3]=-1;
        }
        else {
          tmp[2]=idxMap[v0->point()];
          tmp[3]=idxMap[v1->point()];
          if (tmp[2]>tmp[3])
            std::swap(tmp[2],tmp[3]);
          assert(tmp[2]!=tmp[3]);
        }

        tmp[4]=faceIdx;
        tmp[5]=i;
        tmp[6]=currentEdgeIdx+i;
        edgesWithInfo[currentEdgeIdx+i]=tmp;

        faces[faceIdx][i]=currentEdgeIdx+i;
      }
      currentEdgeIdx+=3;
      faceIdx++;
    }

    std::sort(edgesWithInfo.begin(),edgesWithInfo.end());
    std::vector<int> newIdx(edgesWithInfo.size());
    int cntEdges=0;
    for (int i=0;i<edgesWithInfo.size();i++) {
      //check that the edge is stored twice, which means it is not in the boundary of the convex hull
      newIdx[edgesWithInfo[i][6]]=cntEdges;
      if (i+1<edgesWithInfo.size() && edgesWithInfo[i][0]==edgesWithInfo[i+1][0] && edgesWithInfo[i][1]==edgesWithInfo[i+1][1]) {
        i++;
        newIdx[edgesWithInfo[i][6]]=cntEdges;
      }
      cntEdges++;
    }
    edges.resize(cntEdges);
    otherDiag.resize(cntEdges);
    edgeIncidentFaces.resize(cntEdges);
    int j=0;
    for (int i=0;i<edgesWithInfo.size();i++) {
      //check that the edge is stored twice, which means it is not in the boundary of the convex hull
      if (i+1<edgesWithInfo.size() && edgesWithInfo[i][0]==edgesWithInfo[i+1][0] && edgesWithInfo[i][1]==edgesWithInfo[i+1][1]) {
        edges[j]={edgesWithInfo[i][0],edgesWithInfo[i][1]};
        otherDiag[j]={edgesWithInfo[i][2],edgesWithInfo[i][3]};
        assert(otherDiag[j][0]<otherDiag[j][1]);
        edgeIncidentFaces[j]={edgesWithInfo[i][4],edgesWithInfo[i][5],edgesWithInfo[i+1][4],edgesWithInfo[i+1][5]};
        j++;
        i++;
      }
      else {
        edges[j]={edgesWithInfo[i][0],edgesWithInfo[i][1]};
        otherDiag[j]={-1,-1};
        edgeIncidentFaces[j]={-1,-1,-1,-1};
        j++;
      }
      assert(edges[j-1][0]<edges[j-1][1]);

    }
    for (int i=0;i<faces.size();i++)
      for (int ii=0;ii<3;ii++)
        faces[i][ii]=newIdx[faces[i][ii]];

    for (int k=0;k<j;k++) {
      assert(compatible(k));
      assert(edges[k][0]!=edges[k][1]);
    }
  }

  //nothing to do here
  void init() override {
    ;
  }

  [[nodiscard]] vi2 getEdges() const override {
    if (edges.empty())
      std::cout<<"something went wrong\n";

    return edges;
  }

  //todo old opcounts, compute new ones 14 array values are changed every time this is called
  //On top of the assignments, up to 20 comparisons are performed
  void flip(int i) override {

    assert(compatible(i));

    i4 quadrilateralI=quadrilateralEdges(i);

    //update otherDiag of the other edges in the quadrilateral
    for (int j=0;j<4;j++) {
      int currentEdgeIdx=quadrilateralI[j];
      if (otherDiag[currentEdgeIdx][0]==-1)
        continue;

      //fix otherDiag
      //the new vertex on the diagonal belongs to otherDiag[i], but not to edges[currentEdgeIdx]
      int newDiagVertex=(otherDiag[i][0]==edges[currentEdgeIdx][0] || otherDiag[i][0]==edges[currentEdgeIdx][1]) ? otherDiag[i][1] : otherDiag[i][0];
      if (otherDiag[currentEdgeIdx][0]==edges[i][0] || otherDiag[currentEdgeIdx][0]==edges[i][1]) {
        otherDiag[currentEdgeIdx][0]=newDiagVertex;
      }
      else {
        if (!(otherDiag[currentEdgeIdx][1]==edges[i][0] || otherDiag[currentEdgeIdx][1]==edges[i][1])) {
          std::cout<<edges[i][0]<<' '<<edges[i][1]<<'\n';
          std::cout<<otherDiag[i][0]<<' '<<otherDiag[i][1]<<'\n';
          std::cout<<"Adjacent quadrilateral:\n"<<edges[currentEdgeIdx][0]<<' '<<edges[currentEdgeIdx][1]<<'\n';
          std::cout<<otherDiag[currentEdgeIdx][0]<<' '<<otherDiag[currentEdgeIdx][1]<<'\n';
        }
        assert(otherDiag[currentEdgeIdx][1]==edges[i][0] || otherDiag[currentEdgeIdx][1]==edges[i][1]);
        otherDiag[currentEdgeIdx][1]=newDiagVertex;
      }
      assert(otherDiag[currentEdgeIdx][0]!=otherDiag[currentEdgeIdx][1]);
    }

    i2 tmp=edges[i];
    edges[i]=otherDiag[i];
    assert(edges[i][0]!=edges[i][1]);
    otherDiag[i]=tmp;

    //other possibility: set to {-1,-1} to make it impossible to reverse an action
    // otherDiag[i]={-1,-1};

    //update faces
    i4 eIF=edgeIncidentFaces[i];
    i3 tmpF=faces[eIF[0]];
    faces[eIF[0]][(eIF[1]+1)%3]=faces[eIF[0]][(eIF[1]+2)%3];
    faces[eIF[0]][(eIF[1]+2)%3]=faces[eIF[2]][(eIF[3]+1)%3];

    faces[eIF[2]][(eIF[3]+1)%3]=faces[eIF[2]][(eIF[3]+2)%3];
    faces[eIF[2]][(eIF[3]+2)%3]=tmpF[(eIF[1]+1)%3];

    //TODO most "delicate" part, check this
    if (edgeIncidentFaces[faces[eIF[0]][(eIF[1]+1)%3]][0]!=-1) {
      int toUpdate=(edgeIncidentFaces[faces[eIF[0]][(eIF[1]+1)%3]][0]==eIF[0]) ? 1 : 3;
      assert(edgeIncidentFaces[faces[eIF[0]][(eIF[1]+1)%3]][toUpdate-1]==eIF[0]);
      edgeIncidentFaces[faces[eIF[0]][(eIF[1]+1)%3]][toUpdate]=(eIF[1]+1)%3;
    }

    if (edgeIncidentFaces[faces[eIF[0]][(eIF[1]+2)%3]][0]!=-1) {
      int toUpdate=(edgeIncidentFaces[faces[eIF[0]][(eIF[1]+2)%3]][0]==eIF[2]) ? 0 : 2;
      assert(edgeIncidentFaces[faces[eIF[0]][(eIF[1]+2)%3]][toUpdate]==eIF[2]);
      edgeIncidentFaces[faces[eIF[0]][(eIF[1]+2)%3]][toUpdate]=eIF[0];
      edgeIncidentFaces[faces[eIF[0]][(eIF[1]+2)%3]][toUpdate+1]=(eIF[1]+2)%3;
    }

    if (edgeIncidentFaces[faces[eIF[2]][(eIF[3]+1)%3]][0]!=-1) {
      int toUpdate=(edgeIncidentFaces[faces[eIF[2]][(eIF[3]+1)%3]][0]==eIF[2]) ? 1 : 3;
      assert(edgeIncidentFaces[faces[eIF[2]][(eIF[3]+1)%3]][toUpdate-1]==eIF[2]);
      edgeIncidentFaces[faces[eIF[2]][(eIF[3]+1)%3]][toUpdate]=(eIF[3]+1)%3;
    }

    if (edgeIncidentFaces[faces[eIF[2]][(eIF[3]+2)%3]][0]!=-1) {
      int toUpdate=(edgeIncidentFaces[faces[eIF[2]][(eIF[3]+2)%3]][0]==eIF[0]) ? 0 : 2;
      assert(edgeIncidentFaces[faces[eIF[2]][(eIF[3]+2)%3]][toUpdate]==eIF[0]);
      edgeIncidentFaces[faces[eIF[2]][(eIF[3]+2)%3]][toUpdate]=eIF[2];
      edgeIncidentFaces[faces[eIF[2]][(eIF[3]+2)%3]][toUpdate+1]=(eIF[3]+2)%3;
    }


    quadrilateralI=quadrilateralEdges(i);
    assert(compatible(i));

    for (int j=0;j<4;j++) {
      // std::cout<<j<<'\n';
      int currentEdgeIdx=quadrilateralI[j];
      bool comp=compatible(currentEdgeIdx);
      assert(comp);
    }
  }

  ~facesDS() override =default;

  i2 otherDiagonal(int i) override {
    //todo fix this properly inside the flip operation
    if (otherDiag[i][0]>otherDiag[i][1])
      return {otherDiag[i][1],otherDiag[i][0]};
    return otherDiag[i];
  }

  bool shareTriangle(int i, int j) override {
    //i and j share a triangle iff j belongs to the quadrilateral of i
    i4 quadrilateral=quadrilateralEdges(i);
    for (int k=0;k<4;k++)
      if (quadrilateral[k]==j)
        return true;

    return false;
  }

  std::array<int,4> quadrilateralEdges(int i) override {
    if (edgeIncidentFaces[i][0]==-1 || edgeIncidentFaces[i][2]==-1)
      return {-1,-1,-1,-1};
    auto eIF=edgeIncidentFaces[i];
    return {faces[eIF[0]][(eIF[1]+1)%3],faces[eIF[0]][(eIF[1]+2)%3],faces[eIF[2]][(eIF[3]+1)%3],faces[eIF[2]][(eIF[3]+2)%3]};
  }

  i2 incidentFaces(int i) {
    return {edgeIncidentFaces[i][0],edgeIncidentFaces[i][2]};
  }
  // vi3 getFaces() {
  //   return faces;
  // }
  //
  // vi4 getEdgeIncidentFaces() {
  //   return edgeIncidentFaces;
  // }
};




#endif //FLIPPINGDS_HPP
