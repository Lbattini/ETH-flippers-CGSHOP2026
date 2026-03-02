#include "../include/differentDiagonals.hpp"
#include "../include/geometryUtils.hpp"
#include "../include/flippingDS.hpp"
#include <boost/graph/edge_coloring.hpp>
// #include <CGAL/draw_triangulation_2.h>

typedef IK::Segment_2 S;

struct EdgeProperty {
  size_t color;
  int idx;
};
typedef boost::adjacency_list< boost::vecS, boost::vecS, boost::undirectedS, boost::no_property, EdgeProperty,
        boost::no_property >
        Graph;


int flipDifferentDiagonals(flippingDS& t, flippingDS& t1, size_t n, const vP& points, std::vector<vi2>& flipsPerRound) {
  t.init();
  t1.init();
  int F=0;
  //stop when no edge has been flipped in the previous round
  bool toFlip=true;

  int totalFlips=0;
  std::vector<std::vector<int>> adjMatT1(n*n, std::vector<int>(n,0));

  vi2 edgeL=t.getEdges();
  vi2 edgeL1=t1.getEdges();
  // for (i2 e : edgeL)
  //   std::cout<<e[0]<<' '<<e[1]<<std::endl;
  for (i2 e : edgeL1) {
    int a=e[0],b=e[1];
    if (a>b) std::swap(a,b);
    adjMatT1[a][b]=1;
  }
  // for (auto e : t1.finite_edges()) {
  //   ConstrainedTriangulation::Face_handle currentFace=e.first;
  //   P p0=currentFace->vertex((e.second+1)%3)->point(), p1=currentFace->vertex((e.second+2)%3)->point();
  //   int a=idxMap.at(p0), b=idxMap.at(p1);
  //   if (a>b) std::swap(a,b);
  //   adjMatT1[a][b]=1;
  // }

  while (toFlip && totalFlips<(n*n)) {
    toFlip=false;
    std::vector<ConstrainedTriangulation::Finite_edges_iterator> delEdges;
    //2 means don't flip because an edge that shares a triangle with it was flipped
    std::vector<int> flipped(n*n,0);

    // std::vector<int> wait(n*n,1);
    int flipsThisRound=0;
    for (int i=0;i<edgeL.size();i++) {
      int a=edgeL[i][0],b=edgeL[i][1];
      if (a>b) std::swap(a,b);
      // std::cout<<"candidate: "<<a<<' '<<b<<'\n';
      // std::cout<<"idx edge list: "<<i<<'\n';
      i2 otherDiag=t.otherDiagonal(i);
      if (flipped[a*n+b]==0 && otherDiag[0]!=-1 && adjMatT1[a][b]==0) {
        if (adjMatT1[otherDiag[0]][otherDiag[1]]==1 /*|| wait[a*n+b]==0*/){
          // if (adjMatT1[otherDiag[0]][otherDiag[1]]==0)
            // wait[a*n+b]=1;
          std::array<P,4> quad={points[a],points[otherDiag[0]],points[b],points[otherDiag[1]]};

          // std::cout<<"candidate: "<<a<<' '<<b<<'\n';
          // std::cout<<otherDiag[0]<<' '<<otherDiag[1]<<'\n';
          if (canBeFlipped(quad)) {
            int c=otherDiag[0],d=otherDiag[1];
            flipped[std::min(a,c)*n+std::max(a,c)]=2;
            flipped[std::min(b,c)*n+std::max(b,c)]=2;
            flipped[std::min(b,d)*n+std::max(b,d)]=2;
            flipped[std::min(a,d)*n+std::max(a,d)]=2;
            // for (int j=i+1;j<edgeL.size();j++)
            //   if (t.shareTriangle(i,j))
            //     flipped[j]=2;

            edgeL[i]=otherDiag;

            if (!toFlip) {
              toFlip=true;
              flipsPerRound.push_back(vi2());
            }
            t.flip(i);
            // std::cout<<"flipped: "<<a<<' '<<b<<std::endl;

            flipsPerRound[F].push_back({a,b});
            flipsThisRound++;
          }
        }
        // else
        //   wait[a*n+b]--;
      }
    }

    if (toFlip)
      F++;
    totalFlips+=flipsThisRound;
  }

  //todo check that t and t1 are equal
  return F;
}

int countIntersections(vi2& edgeL, vi2& edgeL1, const vP& points) {
  int count=0;

  for (i2 e : edgeL) {
    S s(points[e[0]],points[e[1]]);
    for (i2 e1 : edgeL1) {
      S s1(points[e1[0]],points[e1[1]]);
      if (e[0]!=e1[0] && e[0]!=e1[1] && e[1]!=e1[0] && e[1]!=e1[1] && CGAL::do_intersect(s,s1))
        count++;
    }
  }
  return count;
}

int countIntersections(i2 e, vi2& edgeL1, const vP& points, vi& intersectionsWithTarget) {
  int count=0;
  int n=points.size();
  int edgeUID=std::min(e[0],e[1])*n+std::max(e[0],e[1]);
  if (intersectionsWithTarget[edgeUID]!=-1)
    return intersectionsWithTarget[edgeUID];

  S s(points[e[0]],points[e[1]]);
  for (i2 e1 : edgeL1) {
    if (std::min(points[e1[0]].x(),points[e1[1]].x()) > std::max(points[e[0]].x(),points[e[1]].x()))
      break;
    S s1(points[e1[0]],points[e1[1]]);
    if (e[0]!=e1[0] && e[0]!=e1[1] && e[1]!=e1[0] && e[1]!=e1[1] && CGAL::do_intersect(s,s1))
      count++;
  }

  return intersectionsWithTarget[edgeUID]=count;
}

int countIntersections(i2 e, vi2& edgeL1, const vP& points, mii& intersectionsWithTarget) {
  int count=0;
  int n=points.size();
  int edgeUID=std::min(e[0],e[1])*n+std::max(e[0],e[1]);
  if (intersectionsWithTarget.count(edgeUID)!=0)
    return intersectionsWithTarget[edgeUID];

  S s(points[e[0]],points[e[1]]);
  for (i2 e1 : edgeL1) {
    if (std::min(points[e1[0]].x(),points[e1[1]].x()) > std::max(points[e[0]].x(),points[e[1]].x()))
      break;
    S s1(points[e1[0]],points[e1[1]]);
    if (e[0]!=e1[0] && e[0]!=e1[1] && e[1]!=e1[0] && e[1]!=e1[1] && CGAL::do_intersect(s,s1))
      count++;
  }

  return intersectionsWithTarget[edgeUID]=count;
}

int countIntersections(i2 e, std::vector<vi2>& edgeLVec, const vP& points) {
  int count=0;
  int n=points.size();
  int edgeUID=std::min(e[0],e[1])*n+std::max(e[0],e[1]);

  S s(points[e[0]],points[e[1]]);
  for (const auto& edgeL1 : edgeLVec)
    for (i2 e1 : edgeL1) {
      // if (std::min(points[e1[0]].x(),points[e1[1]].x()) > std::max(points[e[0]].x(),points[e[1]].x()))
      //   break;
      S s1(points[e1[0]],points[e1[1]]);
      if (e[0]!=e1[0] && e[0]!=e1[1] && e[1]!=e1[0] && e[1]!=e1[1] && CGAL::do_intersect(s,s1))
        count++;
    }

  return count;
}

int countIntersections(i2 e, std::vector<vi2>& edgeLVec, const vP& points,  mii& intersectionsWithTarget) {
  int count=0;
  int n=points.size();
  int edgeUID=std::min(e[0],e[1])*n+std::max(e[0],e[1]);
  if (intersectionsWithTarget.count(edgeUID)!=0)
    return intersectionsWithTarget[edgeUID];

  S s(points[e[0]],points[e[1]]);
  for (const auto& edgeL1 : edgeLVec)
    for (i2 e1 : edgeL1) {
      if (std::min(points[e1[0]].x(),points[e1[1]].x()) > std::max(points[e[0]].x(),points[e[1]].x()))
        break;
      S s1(points[e1[0]],points[e1[1]]);
      if (e[0]!=e1[0] && e[0]!=e1[1] && e[1]!=e1[0] && e[1]!=e1[1] && CGAL::do_intersect(s,s1))
        count++;
    }

  return intersectionsWithTarget[edgeUID]=count;
}

int countIntersections(i2 e, vi2& edgeL1, const vP& points) {
  int count=0;
  int n=points.size();
  // int edgeUID=std::min(e[0],e[1])*n+std::max(e[0],e[1]);

  S s(points[e[0]],points[e[1]]);
  for (i2 e1 : edgeL1) {
    // if (std::min(points[e1[0]].x(),points[e1[1]].x()) > std::max(points[e[0]].x(),points[e[1]].x()))
    //   break;
    S s1(points[e1[0]],points[e1[1]]);
    if (e[0]!=e1[0] && e[0]!=e1[1] && e[1]!=e1[0] && e[1]!=e1[1] && CGAL::do_intersect(s,s1))
      count++;
  }

  return count;
}

//return true if e1 reduces the number of intersections to most triangulations
bool compareItsReductionCardinality(i2 e0, i2 e1, std::vector<vi2>& edgeLVec, const vP& points,  std::vector<mii>& intersectionsWithTarget) {
  int cardinalityDifference=0;
  int i=0;
  for (auto& edgeL1 : edgeLVec) {
    int currentDif=countIntersections(e0,edgeL1,points,intersectionsWithTarget[i])-
            countIntersections(e1,edgeL1,points,intersectionsWithTarget[i]);
    i+=1;
    if (currentDif>0)
      cardinalityDifference+=1;
    else if (currentDif < 0)
      cardinalityDifference-=1;
  }
  return cardinalityDifference>0;
}
int flipToReduceIntesections(flippingDS& t, flippingDS& t1, size_t n, const vP& points, std::vector<vi2>& flipsPerRound, vi& intersectionsWithTarget) {
  t.init();
  t1.init();
  int F=0;
  //stop when no edge has been flipped in the previous round
  bool toFlip=true;

  int totalFlips=0;
  std::vector<std::vector<int>> adjMatT1(n*n, std::vector<int>(n,0));

  vi2 edgeL=t.getEdges();
  vi2 edgeL1=t1.getEdges();

  for (i2& e : edgeL1)
    if (e[0]>e[1])
      std::swap(e[0],e[1]);
  std::sort(edgeL1.begin(),edgeL1.end(),[points](i2 edgeA, i2 edgeB) {
    return std::min(points[edgeA[0]].x(),points[edgeA[1]].x())<std::min(points[edgeB[0]].x(),points[edgeB[1]].x());
  });

  for (i2 e : edgeL1) {
    int a=e[0],b=e[1];
    if (a>b) std::swap(a,b);
    adjMatT1[a][b]=1;
  }

  // std::cout<<"Intersections in the same triangulation (should be 0): "<<countIntersections(edgeL1,edgeL1,points)<<'\n';
  // std::cout<<"Intersections between t0 and t1: "<<countIntersections(edgeL,edgeL1,points)<<'\n';

  while (toFlip && totalFlips<(n*n)) {
    toFlip=false;
    //2 means don't flip because an edge that shares a triangle with it was flipped
    std::vector<int> flipped(n*n,0);

    int flipsThisRound=0;

    for (int i=0;i<edgeL.size();i++) {
      int a=edgeL[i][0],b=edgeL[i][1];
      if (a>b) std::swap(a,b);

      i2 otherDiag=t.otherDiagonal(i);
      if (flipped[a*n+b]==0 && otherDiag[0]!=-1 && adjMatT1[a][b]==0) {
        if (adjMatT1[otherDiag[0]][otherDiag[1]]==1){
          // if (countIntersections(otherDiag,edgeL1,points,intersectionsWithTarget)!=0)
          //   std::cout<<"Error\n";

          std::array<P,4> quad={points[a],points[otherDiag[0]],points[b],points[otherDiag[1]]};

          if (canBeFlipped(quad)) {
            int c=otherDiag[0],d=otherDiag[1];
            flipped[a*n+b]=1;
            flipped[std::min(a,c)*n+std::max(a,c)]=2;
            flipped[std::min(b,c)*n+std::max(b,c)]=2;
            flipped[std::min(b,d)*n+std::max(b,d)]=2;
            flipped[std::min(a,d)*n+std::max(a,d)]=2;

            edgeL[i]=otherDiag;

            if (!toFlip) {
              toFlip=true;
              flipsPerRound.push_back(vi2());
            }
            t.flip(i);

            flipsPerRound[F].push_back({a,b});
            flipsThisRound++;
          }
        }
      }
    }

    for (int i=0;i<edgeL.size();i++) {
      int a=edgeL[i][0],b=edgeL[i][1];
      if (a>b) std::swap(a,b);
      // std::cout<<"candidate: "<<a<<' '<<b<<'\n';
      // std::cout<<"idx edge list: "<<i<<'\n';
      i2 otherDiag=t.otherDiagonal(i);

      if (flipped[a*n+b]==0 && otherDiag[0]!=-1
        && countIntersections(otherDiag,edgeL1,points,intersectionsWithTarget)<countIntersections(edgeL[i],edgeL1,points,intersectionsWithTarget)){
        // if (adjMatT1[otherDiag[0]][otherDiag[1]]==0)
          // wait[a*n+b]=1;
        std::array<P,4> quad={points[a],points[otherDiag[0]],points[b],points[otherDiag[1]]};

        // std::cout<<"candidate: "<<a<<' '<<b<<'\n';
        // std::cout<<otherDiag[0]<<' '<<otherDiag[1]<<'\n';
        if (canBeFlipped(quad)) {
          int c=otherDiag[0],d=otherDiag[1];
          flipped[std::min(a,c)*n+std::max(a,c)]=2;
          flipped[std::min(b,c)*n+std::max(b,c)]=2;
          flipped[std::min(b,d)*n+std::max(b,d)]=2;
          flipped[std::min(a,d)*n+std::max(a,d)]=2;
          // for (int j=i+1;j<edgeL.size();j++)
          //   if (t.shareTriangle(i,j))
          //     flipped[j]=2;

          edgeL[i]=otherDiag;

          if (!toFlip) {
            toFlip=true;
            flipsPerRound.push_back(vi2());
          }
          t.flip(i);
          // std::cout<<"flipped: "<<a<<' '<<b<<std::endl;

          flipsPerRound[F].push_back({a,b});
          flipsThisRound++;
        }

      }
    }

    if (toFlip)
      F++;
    totalFlips+=flipsThisRound;
  }

  return F;
}

// void flipEdge(int a, int b, int c, int d, int edgeUid, size_t n, bool& toFlip, int F,  std::set<int>& flipped, std::vector<vi2>& flipsPerRound, int i, flippingDS& t) {
//   // int c=otherDiag[0],d=otherDiag[1];
//   flipped.insert(edgeUid);
//   flipped.insert(std::min(a,c)*n+std::max(a,c));
//   flipped.insert(std::min(b,c)*n+std::max(b,c));
//   flipped.insert(std::min(b,d)*n+std::max(b,d));
//   flipped.insert(std::min(a,d)*n+std::max(a,d));
//
//   edgeL[i]=otherDiag;
//
//   if (!toFlip) {
//     toFlip=true;
//     flipsPerRound.push_back(vi2());
//   }
//   t.flip(i);
//
//   flipsPerRound[F].push_back({a,b});
//   flipsThisRound++;
// }

int flipToReduceIntesections(flippingDS& t, flippingDS& t1, size_t n, const vP& points, std::vector<vi2>& flipsPerRound, mii& intersectionsWithTarget) {
  t.init();
  t1.init();
  int F=0;
  //stop when no edge has been flipped in the previous round
  bool toFlip=true;

  int totalFlips=0;

  vi2 edgeL=t.getEdges();
  vi2 edgeL1=t1.getEdges();
  vi2 edgeL1IncrX(edgeL1.size());

  for (size_t i=0;i<edgeL1.size();i++){
    assert(edgeL[i][0]<edgeL[i][1]);
    assert(t.otherDiagonal(i)[0]<=t.otherDiagonal(i)[1]);
    if (edgeL1[i][0]>edgeL1[i][1])
      std::swap(edgeL1[i][0],edgeL1[i][1]);
    edgeL1IncrX[i]=edgeL1[i];
  }
  std::sort(edgeL1.begin(),edgeL1.end());
  std::sort(edgeL1IncrX.begin(),edgeL1IncrX.end(),[points](i2 edgeA, i2 edgeB) {
    return std::min(points[edgeA[0]].x(),points[edgeA[1]].x())<std::min(points[edgeB[0]].x(),points[edgeB[1]].x());
  });

  // std::cout<<"Intersections in the same triangulation (should be 0): "<<countIntersections(edgeL1,edgeL1,points)<<'\n';
  // std::cout<<"Intersections between t0 and t1: "<<countIntersections(edgeL,edgeL1,points)<<'\n';

  while (toFlip && totalFlips<(n*n)) {
    toFlip=false;
    std::set<int> flipped;

    int flipsThisRound=0;


    Graph G0(n);
    int cntFlippable=0;
    //todo always flip isolated edges
    for (int i=0;i<edgeL.size();i++) {
      int a=edgeL[i][0],b=edgeL[i][1];

      i2 otherDiag=t.otherDiagonal(i);
      if (otherDiag[0]!=-1
        && std::binary_search(edgeL1.begin(),edgeL1.end(),otherDiag)){

        std::array<P,4> quad={points[a],points[otherDiag[0]],points[b],points[otherDiag[1]]};
        if (canBeFlipped(quad)){
          boost::add_edge(a,b,G0);
          G0[boost::edge(a,b,G0).first].idx=i;
          cntFlippable++;
        }

        }
    }

    if (cntFlippable>0) {
      size_t colors = boost::edge_coloring(G0, get(&EdgeProperty::color, G0));

      vi colorFreq(colors,0);

      boost::graph_traits<Graph>::edge_iterator ei, ei_end;
      for (boost::tie(ei, ei_end) = boost::edges(G0); ei != ei_end; ++ei)
        colorFreq[G0[boost::edge(boost::source(*ei, G0),boost::target(*ei,G0),G0).first].color]++;

      int maxFreq=colorFreq[0],bestColor=0;
      for (int i=1;i<colors;i++)
        if (colorFreq[i]>maxFreq) {
          maxFreq=colorFreq[i];
          bestColor=i;
        }
      for (boost::tie(ei, ei_end) = boost::edges(G0); ei != ei_end; ++ei) {
        int a=boost::source(*ei, G0), b=boost::target(*ei,G0);
        if (G0[boost::edge(a,b,G0).first].color==bestColor || (boost::out_degree(a,G0)==1 && boost::out_degree(b,G0)==1)) {
          int i=G0[boost::edge(a,b,G0).first].idx;
          i2 otherDiag=t.otherDiagonal(i);
          int c=otherDiag[0],d=otherDiag[1];
          if (a>b)
            std::swap(a,b);

          flipped.insert(a*n+b);
          flipped.insert(std::min(a,c)*n+std::max(a,c));
          flipped.insert(std::min(b,c)*n+std::max(b,c));
          flipped.insert(std::min(b,d)*n+std::max(b,d));
          flipped.insert(std::min(a,d)*n+std::max(a,d));

          edgeL[i]=otherDiag;
          assert(c<d);

          if (!toFlip) {
            toFlip=true;
            flipsPerRound.push_back(vi2());
          }
          t.flip(i);

          flipsPerRound[F].push_back({a,b});
          flipsThisRound++;
        }
      }
    }

    Graph G1(n);
    cntFlippable=0;
    for (int i=0;i<edgeL.size();i++) {
      int a=edgeL[i][0],b=edgeL[i][1];
      i2 otherDiag=t.otherDiagonal(i);
      if (flipped.count(a*n+b)==0 && otherDiag[0]!=-1){
        std::array<P,4> quad={points[a],points[otherDiag[0]],points[b],points[otherDiag[1]]};
        if (canBeFlipped(quad)
        && countIntersections(otherDiag,edgeL1IncrX,points,intersectionsWithTarget)<countIntersections(edgeL[i],edgeL1IncrX,points,intersectionsWithTarget)) {
          boost::add_edge(a,b,G1);
          G1[boost::edge(a,b,G1).first].idx=i;
          cntFlippable++;
        }

        }
    }

    if (cntFlippable>0) {
      size_t colors = boost::edge_coloring(G1, get(&EdgeProperty::color, G1));

      vi colorFreq(colors,0);

      boost::graph_traits<Graph>::edge_iterator ei, ei_end;
      for (boost::tie(ei, ei_end) = boost::edges(G1); ei != ei_end; ++ei)
        colorFreq[G1[boost::edge(boost::source(*ei, G1),boost::target(*ei,G1),G1).first].color]++;

      int maxFreq=colorFreq[0],bestColor=0;
      for (int i=1;i<colors;i++)
        if (colorFreq[i]>maxFreq) {
          maxFreq=colorFreq[i];
          bestColor=i;
        }
      // std::cout<<"Flippable "<<cntFlippable<<", number of colors "<<colors<<", best color "<<bestColor<<" frequency "<<maxFreq<<'\n';
      for (boost::tie(ei, ei_end) = boost::edges(G1); ei != ei_end; ++ei) {
        int a=boost::source(*ei, G1), b=boost::target(*ei,G1);
        if (G1[boost::edge(a,b,G1).first].color==bestColor || (boost::out_degree(a,G1)==1 && boost::out_degree(b,G1)==1)) {
          int i=G1[boost::edge(a,b,G1).first].idx;
          i2 otherDiag=t.otherDiagonal(i);
          int c=otherDiag[0],d=otherDiag[1];

          edgeL[i]=otherDiag;
          assert(c<d);

          if (!toFlip) {
            toFlip=true;
            flipsPerRound.push_back(vi2());
          }
          t.flip(i);

          flipsPerRound[F].push_back({a,b});
          flipsThisRound++;
        }
      }
    }


    // for (int i=0;i<edgeL.size();i++) {
    //   int a=edgeL[i][0],b=edgeL[i][1];
    //   // if (a>b)
    //   //   std::cout<<a<<' '<<b<<'\n';
    //   // assert(a<b);
    //   int edgeUid=a*n+b;
    //
    //   i2 otherDiag=t.otherDiagonal(i);
    //   if (flipped.count(edgeUid)==0 && otherDiag[0]!=-1 //&& !std::binary_search(edgeL1.begin(),edgeL1.end(),edgeL[i])
    //     && std::binary_search(edgeL1.begin(),edgeL1.end(),otherDiag)){
    //       std::array<P,4> quad={points[a],points[otherDiag[0]],points[b],points[otherDiag[1]]};
    //
    //       if (canBeFlipped(quad)) {
    //         int c=otherDiag[0],d=otherDiag[1];
    //         flipped.insert(edgeUid);
    //         flipped.insert(std::min(a,c)*n+std::max(a,c));
    //         flipped.insert(std::min(b,c)*n+std::max(b,c));
    //         flipped.insert(std::min(b,d)*n+std::max(b,d));
    //         flipped.insert(std::min(a,d)*n+std::max(a,d));
    //
    //         edgeL[i]=otherDiag;
    //         assert(c<d);
    //
    //         if (!toFlip) {
    //           toFlip=true;
    //           flipsPerRound.push_back(vi2());
    //         }
    //         t.flip(i);
    //
    //         flipsPerRound[F].push_back({a,b});
    //         flipsThisRound++;
    //       }
    //
    //   }
    // }

    // for (int i=0;i<edgeL.size();i++) {
    //   int a=edgeL[i][0],b=edgeL[i][1];
    //   // assert(a<b);
    //   int edgeUid=a*n+b;
    //
    //   i2 otherDiag=t.otherDiagonal(i);
    //
    //   if (flipped.count(edgeUid)==0 && otherDiag[0]!=-1
    //     // && !std::binary_search(edgeL1.begin(),edgeL1.end(),otherDiag)
    //     // && countIntersections(otherDiag,edgeL1IncrX,points,intersectionsWithTarget)<countIntersections(edgeL[i],edgeL1IncrX,points,intersectionsWithTarget)
    //     ){
    //     std::array<P,4> quad={points[a],points[otherDiag[0]],points[b],points[otherDiag[1]]};
    //
    //     if (canBeFlipped(quad)
    //       && countIntersections(otherDiag,edgeL1IncrX,points,intersectionsWithTarget)<countIntersections(edgeL[i],edgeL1IncrX,points,intersectionsWithTarget)) {
    //       int c=otherDiag[0],d=otherDiag[1];
    //       flipped.insert(edgeUid);
    //       flipped.insert(std::min(a,c)*n+std::max(a,c));
    //       flipped.insert(std::min(b,c)*n+std::max(b,c));
    //       flipped.insert(std::min(b,d)*n+std::max(b,d));
    //       flipped.insert(std::min(a,d)*n+std::max(a,d));
    //
    //       edgeL[i]=otherDiag;
    //
    //       if (!toFlip) {
    //         toFlip=true;
    //         flipsPerRound.push_back(vi2());
    //       }
    //       t.flip(i);
    //
    //       flipsPerRound[F].push_back({a,b});
    //       flipsThisRound++;
    //     }
    //
    //   }
    // }

    // std::cout<<flipsThisRound<<'\n';

    if (toFlip)
      F++;
    totalFlips+=flipsThisRound;
  }

  return F;
}

int parallelFlip(int j, facesDS& tJ, std::array<vi2,2>& edgeL, size_t n, const vP& points, vi2& flips, vi& flipsIdx) {
  vi2 edgeLSorted=edgeL[(j+1)%2];

  vi2 edgeLIncrX=edgeL[(j+1)%2];

  std::sort(edgeLSorted.begin(),edgeLSorted.end());
  std::sort(edgeLIncrX.begin(),edgeLIncrX.end(),[points](i2 edgeA, i2 edgeB) {
    return std::min(points[edgeA[0]].x(),points[edgeA[1]].x())<std::min(points[edgeB[0]].x(),points[edgeB[1]].x());
  });



  int flipsThisRound=0;

  bool toFlip=false;
  std::set<int> flipped;

  Graph G0(3*n);
  int cntFlippable=0;
  //todo always flip isolated edges
  for (int i=0;i<edgeL[j].size();i++) {
    int a=edgeL[j][i][0],b=edgeL[j][i][1];

    i2 otherDiag;

    otherDiag=tJ.otherDiagonal(i);


    if (otherDiag[0]!=-1
      && std::binary_search(edgeLSorted.begin(),edgeLSorted.end(),otherDiag)){

      std::array<P,4> quad={points[a],points[otherDiag[0]],points[b],points[otherDiag[1]]};
      if (canBeFlipped(quad)){
        i2 uv=tJ.incidentFaces(i);

        boost::add_edge(uv[0],uv[1],G0);
        G0[boost::edge(uv[0],uv[1],G0).first].idx=i;
        cntFlippable++;
      }

      }
  }

  if (cntFlippable>0) {
    size_t colors = boost::edge_coloring(G0, get(&EdgeProperty::color, G0));

    vi colorFreq(colors,0);

    boost::graph_traits<Graph>::edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = boost::edges(G0); ei != ei_end; ++ei)
      colorFreq[G0[boost::edge(boost::source(*ei, G0),boost::target(*ei,G0),G0).first].color]++;

    int maxFreq=colorFreq[0],bestColor=0;
    for (int i=1;i<colors;i++)
      if (colorFreq[i]>maxFreq) {
        maxFreq=colorFreq[i];
        bestColor=i;
      }
    for (boost::tie(ei, ei_end) = boost::edges(G0); ei != ei_end; ++ei) {
      int u=boost::source(*ei, G0), v=boost::target(*ei,G0);
      if (G0[boost::edge(u,v,G0).first].color==bestColor /*|| (boost::out_degree(u,G0)<=1 && boost::out_degree(v,G0)<=1)*/) {
        int i=G0[boost::edge(u,v,G0).first].idx;
        int a=edgeL[j][i][0],b=edgeL[j][i][1];
        i2 otherDiag;

        otherDiag=tJ.otherDiagonal(i);
        int c=otherDiag[0],d=otherDiag[1];
        if (a>b)
          std::swap(a,b);

        flipped.insert(a*n+b);
        flipped.insert(std::min(a,c)*n+std::max(a,c));
        flipped.insert(std::min(b,c)*n+std::max(b,c));
        flipped.insert(std::min(b,d)*n+std::max(b,d));
        flipped.insert(std::min(a,d)*n+std::max(a,d));

        edgeL[j][i]=otherDiag;
        assert(c<d);

        if (!toFlip) {
          toFlip=true;
        }

        tJ.flip(i);

        if (j==0)
          flips.push_back({a,b});
        else
          flips.push_back({c,d});

        flipsIdx.push_back(i);
        flipsThisRound++;
      }
    }
  }

  vi2 sortedEdges;
  for (int i=0;i<edgeL[j].size();i++) {
    int a=edgeL[j][i][0],b=edgeL[j][i][1];
    if (a>b)
      std::swap(a,b);
    int edgeUid=a*n+b;

    i2 otherDiag;

    otherDiag=tJ.otherDiagonal(i);


    if (flipped.count(edgeUid)==0 && otherDiag[0]!=-1) {
      std::array<P,4> quad={points[a],points[otherDiag[0]],points[b],points[otherDiag[1]]};
      int currentDif=countIntersections(edgeL[j][i],edgeLIncrX,points)-
        countIntersections(otherDiag,edgeLIncrX,points);

      if (canBeFlipped(quad)
        && currentDif>0) {
        sortedEdges.push_back({currentDif,i});
        }
    }
  }

  std::sort(sortedEdges.rbegin(),sortedEdges.rend());
  for (i2 curEdge : sortedEdges) {
    int i=curEdge[1];
    int a=edgeL[j][i][0],b=edgeL[j][i][1];
    if (a>b)
      std::swap(a,b);
    int edgeUid=a*n+b;

    i2 otherDiag=tJ.otherDiagonal(i);

    if (flipped.count(edgeUid)==0){
      int c=otherDiag[0],d=otherDiag[1];
      flipped.insert(edgeUid);
      flipped.insert(std::min(a,c)*n+std::max(a,c));
      flipped.insert(std::min(b,c)*n+std::max(b,c));
      flipped.insert(std::min(b,d)*n+std::max(b,d));
      flipped.insert(std::min(a,d)*n+std::max(a,d));

      edgeL[j][i]=otherDiag;

      if (!toFlip) {
        toFlip=true;
      }

      if (j==0)
        tJ.flip(i);

      if (j==0)
        flips.push_back({a,b});
      else
        flips.push_back({c,d});

      flipsIdx.push_back(i);
      flipsThisRound++;
    }
  }

  return countIntersections(edgeL[0],edgeL[1],points);
}

int flipToReduceIntesectionsBiDir(facesDS& t, facesDS& t1, size_t n, const vP& points, std::vector<vi2>& flipsPerRound) {

  int F=0;
  //stop when no edge has been flipped in the previous round
  bool toFlip=true;

  int totalFlips=0;

  std::array<vi2,2> edgeL={t.getEdges(),t1.getEdges()};

  // std::array<vi2,2> edgeLSorted={edgeL[0],edgeL[1]};

  // std::array<vi2,2> edgeLIncrX={edgeL[0],edgeL[1]};

  // for (int j=0;j<2;j++) {
  //   std::sort(edgeL[j].begin(),edgeL[j].end());
  //   std::sort(edgeLIncrX[j].begin(),edgeLIncrX[j].end(),[points](i2 edgeA, i2 edgeB) {
  //     return std::min(points[edgeA[0]].x(),points[edgeA[1]].x())<std::min(points[edgeB[0]].x(),points[edgeB[1]].x());
  //   });
  // }

  std::vector<vi2> flipsPerRound1;

  // std::cout<<"Intersections in the same triangulation (should be 0): "<<countIntersections(edgeL1,edgeL1,points)<<'\n';
  // std::cout<<"Intersections between t0 and t1: "<<countIntersections(edgeL,edgeL1,points)<<'\n';

  while (toFlip && totalFlips<(n*n)) {

    std::array<facesDS,2> tTmp={facesDS(t),facesDS(t1)};
    std::array<vi2,2> flipsTmp;
    std::array<vi,2> flipsIdx;
    std::array<std::array<vi2,2>,2> edgeLTmp={std::array<vi2,2>{tTmp[0].getEdges(),edgeL[1]}, std::array<vi2,2>{edgeL[0],tTmp[1].getEdges()}};
    std::array<int,2> cntIts;
    for (int j=0;j<2;j++)
      cntIts[j]=parallelFlip(j,tTmp[j],edgeLTmp[j],n,points,flipsTmp[j],flipsIdx[j]);
    // int cntItsFlip1=parallelFlip(1,t1Tmp,edgeL,n,points,flipsTmp[1]);
    if (cntIts[0]==0 || cntIts[1]==0)
      toFlip=false;
    // else {
    if (cntIts[0]<cntIts[1]) {
      // t=facesDS(tTmp[0]);
      // edgeL=edgeLTmp[0];
      // edgeL[0]=t.getEdges();
      for (int i : flipsIdx[0]) {
        t.flip(i);
        edgeL[0][i]=t.otherDiagonal(i);
      }

      flipsPerRound.push_back(flipsTmp[0]);
      totalFlips+=flipsTmp[0].size();

    }
    else{
      // t1=facesDS(tTmp[1]);
      // edgeL=edgeLTmp[1];
      // edgeL[1]=t1.getEdges();
      for (int i : flipsIdx[1]) {
        t1.flip(i);
        edgeL[1][i]=t1.otherDiagonal(i);
      }
      flipsPerRound1.push_back(flipsTmp[1]);
      totalFlips+=flipsTmp[1].size();
    }
    // }
    // std::array<vi2,2> edgeLSorted={edgeL[0],edgeL[1]};
    //
    // std::array<vi2,2> edgeLIncrX={edgeL[0],edgeL[1]};
    //
    // for (int j=0;j<2;j++) {
    //   std::sort(edgeLSorted[j].begin(),edgeLSorted[j].end());
    //   std::sort(edgeLIncrX[j].begin(),edgeLIncrX[j].end(),[points](i2 edgeA, i2 edgeB) {
    //     return std::min(points[edgeA[0]].x(),points[edgeA[1]].x())<std::min(points[edgeB[0]].x(),points[edgeB[1]].x());
    //   });
    // }
    //
    // // std::array<vi2,2> roundFlipsTmp;
    // int flipsThisRound=0;
    // //todo choose best option rather than alternating
    // int j=F%2;
    // // for (int j=0;j<2;j++) {
    //   toFlip=false;
    //   std::set<int> flipped;
    //
    //   Graph G0(3*n);
    //   int cntFlippable=0;
    //   //todo always flip isolated edges
    //   for (int i=0;i<edgeL[j].size();i++) {
    //     int a=edgeL[j][i][0],b=edgeL[j][i][1];
    //
    //     i2 otherDiag;
    //     if (j==0)
    //      otherDiag=t.otherDiagonal(i);
    //     else
    //       otherDiag=t1.otherDiagonal(i);
    //
    //     if (otherDiag[0]!=-1
    //       && std::binary_search(edgeLSorted[(j+1)%2].begin(),edgeLSorted[(j+1)%2].end(),otherDiag)){
    //
    //       std::array<P,4> quad={points[a],points[otherDiag[0]],points[b],points[otherDiag[1]]};
    //       if (canBeFlipped(quad)){
    //         i2 uv;
    //         if (j==0)
    //          uv=t.incidentFaces(i);
    //         else
    //           uv=t1.incidentFaces(i);
    //         boost::add_edge(uv[0],uv[1],G0);
    //         G0[boost::edge(uv[0],uv[1],G0).first].idx=i;
    //         cntFlippable++;
    //       }
    //
    //       }
    //   }
    //
    //   if (cntFlippable>0) {
    //     size_t colors = boost::edge_coloring(G0, get(&EdgeProperty::color, G0));
    //
    //     vi colorFreq(colors,0);
    //
    //     boost::graph_traits<Graph>::edge_iterator ei, ei_end;
    //     for (boost::tie(ei, ei_end) = boost::edges(G0); ei != ei_end; ++ei)
    //       colorFreq[G0[boost::edge(boost::source(*ei, G0),boost::target(*ei,G0),G0).first].color]++;
    //
    //     int maxFreq=colorFreq[0],bestColor=0;
    //     for (int i=1;i<colors;i++)
    //       if (colorFreq[i]>maxFreq) {
    //         maxFreq=colorFreq[i];
    //         bestColor=i;
    //       }
    //     for (boost::tie(ei, ei_end) = boost::edges(G0); ei != ei_end; ++ei) {
    //       int u=boost::source(*ei, G0), v=boost::target(*ei,G0);
    //       if (G0[boost::edge(u,v,G0).first].color==bestColor /*|| (boost::out_degree(u,G0)<=1 && boost::out_degree(v,G0)<=1)*/) {
    //         int i=G0[boost::edge(u,v,G0).first].idx;
    //         int a=edgeL[j][i][0],b=edgeL[j][i][1];
    //         i2 otherDiag;
    //         if (j==0)
    //           otherDiag=t.otherDiagonal(i);
    //         else
    //           otherDiag=t1.otherDiagonal(i);
    //         int c=otherDiag[0],d=otherDiag[1];
    //         if (a>b)
    //           std::swap(a,b);
    //
    //         flipped.insert(a*n+b);
    //         flipped.insert(std::min(a,c)*n+std::max(a,c));
    //         flipped.insert(std::min(b,c)*n+std::max(b,c));
    //         flipped.insert(std::min(b,d)*n+std::max(b,d));
    //         flipped.insert(std::min(a,d)*n+std::max(a,d));
    //
    //         edgeL[j][i]=otherDiag;
    //         assert(c<d);
    //
    //         if (!toFlip) {
    //           toFlip=true;
    //           if (j==0)
    //             flipsPerRound.push_back(vi2());
    //           else
    //             flipsPerRound1.push_back(vi2());
    //         }
    //
    //         if (j==0)
    //           t.flip(i);
    //         else
    //           t1.flip(i);
    //
    //         if (j==0)
    //           flipsPerRound[flipsPerRound.size()-1].push_back({a,b});
    //         else
    //           flipsPerRound1[flipsPerRound1.size()-1].push_back({c,d});
    //
    //         // if (j==1 && c==77 && (d==89 || d==115))
    //         //   std::cout<<"a b in wrong flip: "<<a<<' '<<b<<'\n';
    //         flipsThisRound++;
    //       }
    //     }
    //   }
    //
    //   vi2 sortedEdges;
    //   for (int i=0;i<edgeL[j].size();i++) {
    //     int a=edgeL[j][i][0],b=edgeL[j][i][1];
    //     if (a>b)
    //       std::swap(a,b);
    //     int edgeUid=a*n+b;
    //
    //     i2 otherDiag;
    //     if (j==0)
    //       otherDiag=t.otherDiagonal(i);
    //     else
    //       otherDiag=t1.otherDiagonal(i);
    //
    //     if (flipped.count(edgeUid)==0 && otherDiag[0]!=-1) {
    //       std::array<P,4> quad={points[a],points[otherDiag[0]],points[b],points[otherDiag[1]]};
    //       int currentDif=countIntersections(edgeL[j][i],edgeLIncrX[(j+1)%2],points)-
    //         countIntersections(otherDiag,edgeLIncrX[(j+1)%2],points);
    //
    //       if (canBeFlipped(quad)
    //         && currentDif>0) {
    //         sortedEdges.push_back({currentDif,i});
    //         }
    //     }
    //   }
    //
    //   std::sort(sortedEdges.rbegin(),sortedEdges.rend());
    //   for (i2 curEdge : sortedEdges) {
    //     int i=curEdge[1];
    //     int a=edgeL[j][i][0],b=edgeL[j][i][1];
    //     if (a>b)
    //       std::swap(a,b);
    //     int edgeUid=a*n+b;
    //
    //     i2 otherDiag;
    //     if (j==0)
    //       otherDiag=t.otherDiagonal(i);
    //     else
    //       otherDiag=t1.otherDiagonal(i);
    //
    //     if (flipped.count(edgeUid)==0){
    //       int c=otherDiag[0],d=otherDiag[1];
    //       flipped.insert(edgeUid);
    //       flipped.insert(std::min(a,c)*n+std::max(a,c));
    //       flipped.insert(std::min(b,c)*n+std::max(b,c));
    //       flipped.insert(std::min(b,d)*n+std::max(b,d));
    //       flipped.insert(std::min(a,d)*n+std::max(a,d));
    //
    //       edgeL[j][i]=otherDiag;
    //
    //       if (!toFlip) {
    //         toFlip=true;
    //         if (j==0)
    //           flipsPerRound.push_back(vi2());
    //         else
    //           flipsPerRound1.push_back(vi2());
    //       }
    //
    //       if (j==0)
    //         t.flip(i);
    //       else
    //         t1.flip(i);
    //
    //       if (j==0)
    //         flipsPerRound[flipsPerRound.size()-1].push_back({a,b});
    //       else
    //         flipsPerRound1[flipsPerRound1.size()-1].push_back({c,d});
    //
    //       // if (j==1 && c==77 && (d==89 || d==115))
    //       //   std::cout<<"a b in wrong flip: "<<a<<' '<<b<<'\n';
    //       flipsThisRound++;
    //     }
    //   }
      //
      // if (cntFlippable>0)
      //   std::cout<<"Flipped "<<flipsThisRound-prev<<" out of "<<cntFlippable<<'\n';
    // }
    // if (toFlip)
      F++;
  }

  for (int i=flipsPerRound1.size()-1; i>=0; i--)
    flipsPerRound.push_back(flipsPerRound1[i]);
  // if (flipsPerRound.size()!=F)
    // std::cout<<"Flips per round: "<<flipsPerRound.size()<<" F: "<<F<<'\n';
  return F;
}

int flipToReduceIntesections(facesDS& t, facesDS& t1, size_t n, const vP& points, std::vector<vi2>& flipsPerRound, int dir) {

  int F=0;
  //stop when no edge has been flipped in the previous round
  bool toFlip=true;

  int totalFlips=0;

  std::array<vi2,2> edgeL={t.getEdges(),t1.getEdges()};

  // std::array<vi2,2> edgeLSorted={edgeL[0],edgeL[1]};

  // std::array<vi2,2> edgeLIncrX={edgeL[0],edgeL[1]};

  // for (int j=0;j<2;j++) {
  //   std::sort(edgeL[j].begin(),edgeL[j].end());
  //   std::sort(edgeLIncrX[j].begin(),edgeLIncrX[j].end(),[points](i2 edgeA, i2 edgeB) {
  //     return std::min(points[edgeA[0]].x(),points[edgeA[1]].x())<std::min(points[edgeB[0]].x(),points[edgeB[1]].x());
  //   });
  // }

  std::vector<vi2> flipsPerRound1;

  // std::cout<<"Intersections in the same triangulation (should be 0): "<<countIntersections(edgeL1,edgeL1,points)<<'\n';
  // std::cout<<"Intersections between t0 and t1: "<<countIntersections(edgeL,edgeL1,points)<<'\n';

  while (toFlip && totalFlips<(n*n)) {

    std::array<facesDS,2> tTmp={facesDS(t),facesDS(t1)};
    std::array<vi2,2> flipsTmp;
    std::array<vi,2> flipsIdx;
    std::array<std::array<vi2,2>,2> edgeLTmp={std::array<vi2,2>{tTmp[0].getEdges(),edgeL[1]}, std::array<vi2,2>{edgeL[0],tTmp[1].getEdges()}};
    std::array<int,2> cntIts;
    int j=dir;
    cntIts[j]=parallelFlip(j,tTmp[j],edgeLTmp[j],n,points,flipsTmp[j],flipsIdx[j]);
    if (cntIts[dir]==0)
      toFlip=false;
    if (dir==0) {
      for (int i : flipsIdx[0]) {
        t.flip(i);
        edgeL[0][i]=t.otherDiagonal(i);
      }

      flipsPerRound.push_back(flipsTmp[0]);
      totalFlips+=flipsTmp[0].size();

    }
    else{
      for (int i : flipsIdx[1]) {
        t1.flip(i);
        edgeL[1][i]=t1.otherDiagonal(i);
      }
      flipsPerRound1.push_back(flipsTmp[1]);
      totalFlips+=flipsTmp[1].size();
    }
      F++;
  }

  for (int i=flipsPerRound1.size()-1; i>=0; i--)
    flipsPerRound.push_back(flipsPerRound1[i]);

  return F;
}

int flipToReduceIntesections(facesDS& t, facesDS& t1, size_t n, const vP& points, std::vector<vi2>& flipsPerRound, mii& intersectionsWithTarget) {
  // flipsPerRound.clear();
  t.init();
  t1.init();
  int F=0;
  //stop when no edge has been flipped in the previous round
  bool toFlip=true;

  int totalFlips=0;

  vi2 edgeL=t.getEdges();
  vi2 edgeL1=t1.getEdges();
  vi2 edgeL1IncrX(edgeL1.size());

  for (size_t i=0;i<edgeL1.size();i++){
    if (edgeL[i][0]>edgeL[i][1])
      std::swap(edgeL[i][0],edgeL[i][1]);
    // assert(edgeL[i][0]<edgeL[i][1]);
    assert(t.otherDiagonal(i)[0]<=t.otherDiagonal(i)[1]);
    if (edgeL1[i][0]>edgeL1[i][1])
      std::swap(edgeL1[i][0],edgeL1[i][1]);
    edgeL1IncrX[i]=edgeL1[i];
  }
  std::sort(edgeL1.begin(),edgeL1.end());
  std::sort(edgeL1IncrX.begin(),edgeL1IncrX.end(),[points](i2 edgeA, i2 edgeB) {
    return std::min(points[edgeA[0]].x(),points[edgeA[1]].x())<std::min(points[edgeB[0]].x(),points[edgeB[1]].x());
  });

  // std::cout<<"Intersections in the same triangulation (should be 0): "<<countIntersections(edgeL1,edgeL1,points)<<'\n';
  // std::cout<<"Intersections between t0 and t1: "<<countIntersections(edgeL,edgeL1,points)<<'\n';

  while (toFlip && totalFlips<(n*n)) {
    toFlip=false;
    std::set<int> flipped;

    int flipsThisRound=0;


    Graph G0(3*n);
    // Graph G0(n);
    int cntFlippable=0;
    //done always flip isolated edges
    for (int i=0;i<edgeL.size();i++) {
      int a=edgeL[i][0],b=edgeL[i][1];

      i2 otherDiag=t.otherDiagonal(i);
      if (otherDiag[0]!=-1
        && std::binary_search(edgeL1.begin(),edgeL1.end(),otherDiag)){

        std::array<P,4> quad={points[a],points[otherDiag[0]],points[b],points[otherDiag[1]]};
        if (canBeFlipped(quad)){
          i2 uv=t.incidentFaces(i);
          boost::add_edge(uv[0],uv[1],G0);
          G0[boost::edge(uv[0],uv[1],G0).first].idx=i;
          // boost::add_edge(a,b,G0);
          // G0[boost::edge(a,b,G0).first].idx=i;
          cntFlippable++;
        }

        }
    }

    if (cntFlippable>0) {
      size_t colors = boost::edge_coloring(G0, get(&EdgeProperty::color, G0));

      vi colorFreq(colors,0);

      boost::graph_traits<Graph>::edge_iterator ei, ei_end;
      for (boost::tie(ei, ei_end) = boost::edges(G0); ei != ei_end; ++ei)
        colorFreq[G0[boost::edge(boost::source(*ei, G0),boost::target(*ei,G0),G0).first].color]++;

      int maxFreq=colorFreq[0],bestColor=0;
      for (int i=1;i<colors;i++)
        if (colorFreq[i]>maxFreq) {
          maxFreq=colorFreq[i];
          bestColor=i;
        }
      for (boost::tie(ei, ei_end) = boost::edges(G0); ei != ei_end; ++ei) {
        int u=boost::source(*ei, G0), v=boost::target(*ei,G0);
        if (G0[boost::edge(u,v,G0).first].color==bestColor || (boost::out_degree(u,G0)<=1 && boost::out_degree(v,G0)<=1)) {
          int i=G0[boost::edge(u,v,G0).first].idx;
          int a=edgeL[i][0],b=edgeL[i][1];
          i2 otherDiag=t.otherDiagonal(i);
          int c=otherDiag[0],d=otherDiag[1];
          if (a>b)
            std::swap(a,b);

          flipped.insert(a*n+b);
          flipped.insert(std::min(a,c)*n+std::max(a,c));
          flipped.insert(std::min(b,c)*n+std::max(b,c));
          flipped.insert(std::min(b,d)*n+std::max(b,d));
          flipped.insert(std::min(a,d)*n+std::max(a,d));

          edgeL[i]=otherDiag;
          assert(c<d);

          if (!toFlip) {
            toFlip=true;
            flipsPerRound.push_back(vi2());
          }
          t.flip(i);

          flipsPerRound[F].push_back({a,b});
          flipsThisRound++;
        }
      }
    }

    vi2 sortedEdges;
    for (int i=0;i<edgeL.size();i++) {
      int a=edgeL[i][0],b=edgeL[i][1];
      int edgeUid=a*n+b;

      i2 otherDiag=t.otherDiagonal(i);

      if (flipped.count(edgeUid)==0 && otherDiag[0]!=-1) {
        std::array<P,4> quad={points[a],points[otherDiag[0]],points[b],points[otherDiag[1]]};
        int currentDif=countIntersections(edgeL[i],edgeL1IncrX,points,intersectionsWithTarget)-
          countIntersections(otherDiag,edgeL1IncrX,points,intersectionsWithTarget);

        if (canBeFlipped(quad)
          && currentDif>0) {
          sortedEdges.push_back({currentDif,i});
          }
      }
    }

    std::sort(sortedEdges.rbegin(),sortedEdges.rend());
    for (i2 curEdge : sortedEdges) {
      int i=curEdge[1];
      int a=edgeL[i][0],b=edgeL[i][1];
      int edgeUid=a*n+b;

      i2 otherDiag=t.otherDiagonal(i);

      if (flipped.count(edgeUid)==0){
        int c=otherDiag[0],d=otherDiag[1];
        flipped.insert(edgeUid);
        flipped.insert(std::min(a,c)*n+std::max(a,c));
        flipped.insert(std::min(b,c)*n+std::max(b,c));
        flipped.insert(std::min(b,d)*n+std::max(b,d));
        flipped.insert(std::min(a,d)*n+std::max(a,d));

        edgeL[i]=otherDiag;

        if (!toFlip) {
          toFlip=true;
          flipsPerRound.push_back(vi2());
        }
        t.flip(i);

        flipsPerRound[F].push_back({a,b});
        flipsThisRound++;
      }
    }
    //
    // if (cntFlippable>0)
    //   std::cout<<"Flipped "<<flipsThisRound-prev<<" out of "<<cntFlippable<<'\n';

    if (toFlip)
      F++;
    totalFlips+=flipsThisRound;
  }

  return F;
}


int flipToReduceIntesections(facesDS& t, facesDS& t1, size_t n, const vP& points, std::vector<vi2>& flipsPerRound, mii& intersectionsWithTarget, int cntTtoT1) {
  // flipsPerRound.clear();
  t.init();
  t1.init();
  int F=0;
  //stop when no edge has been flipped in the previous round
  bool toFlip=true;

  int totalFlips=0;

  vi2 edgeL=t.getEdges();
  vi2 edgeL1=t1.getEdges();
  vi2 edgeL1IncrX(edgeL1.size());

  for (size_t i=0;i<edgeL1.size();i++){
    if (edgeL[i][0]>edgeL[i][1])
      std::swap(edgeL[i][0],edgeL[i][1]);
    assert(t.otherDiagonal(i)[0]<=t.otherDiagonal(i)[1]);
    if (edgeL1[i][0]>edgeL1[i][1])
      std::swap(edgeL1[i][0],edgeL1[i][1]);
    edgeL1IncrX[i]=edgeL1[i];
  }
  std::sort(edgeL1.begin(),edgeL1.end());
  std::sort(edgeL1IncrX.begin(),edgeL1IncrX.end(),[points](i2 edgeA, i2 edgeB) {
    return std::min(points[edgeA[0]].x(),points[edgeA[1]].x())<std::min(points[edgeB[0]].x(),points[edgeB[1]].x());
  });

  // std::cout<<"Intersections in the same triangulation (should be 0): "<<countIntersections(edgeL1,edgeL1,points)<<'\n';
  // std::cout<<"Intersections between t0 and t1: "<<countIntersections(edgeL,edgeL1,points)<<'\n';

  while (toFlip && F<cntTtoT1) {
    toFlip=false;
    std::set<int> flipped;

    int flipsThisRound=0;


    Graph G0(3*n);
    int cntFlippable=0;

    for (int i=0;i<edgeL.size();i++) {
      int a=edgeL[i][0],b=edgeL[i][1];

      i2 otherDiag=t.otherDiagonal(i);
      if (otherDiag[0]!=-1
        && std::binary_search(edgeL1.begin(),edgeL1.end(),otherDiag)){

        std::array<P,4> quad={points[a],points[otherDiag[0]],points[b],points[otherDiag[1]]};
        if (canBeFlipped(quad)){
          i2 uv=t.incidentFaces(i);
          boost::add_edge(uv[0],uv[1],G0);
          G0[boost::edge(uv[0],uv[1],G0).first].idx=i;
          cntFlippable++;
        }

        }
    }

    if (cntFlippable>0) {
      size_t colors = boost::edge_coloring(G0, get(&EdgeProperty::color, G0));

      vi colorFreq(colors,0);

      boost::graph_traits<Graph>::edge_iterator ei, ei_end;
      for (boost::tie(ei, ei_end) = boost::edges(G0); ei != ei_end; ++ei)
        colorFreq[G0[boost::edge(boost::source(*ei, G0),boost::target(*ei,G0),G0).first].color]++;

      int maxFreq=colorFreq[0],bestColor=0;
      for (int i=1;i<colors;i++)
        if (colorFreq[i]>maxFreq) {
          maxFreq=colorFreq[i];
          bestColor=i;
        }
      for (boost::tie(ei, ei_end) = boost::edges(G0); ei != ei_end; ++ei) {
        int u=boost::source(*ei, G0), v=boost::target(*ei,G0);
        if (G0[boost::edge(u,v,G0).first].color==bestColor || (boost::out_degree(u,G0)<=1 && boost::out_degree(v,G0)<=1)) {
          int i=G0[boost::edge(u,v,G0).first].idx;
          int a=edgeL[i][0],b=edgeL[i][1];
          i2 otherDiag=t.otherDiagonal(i);
          int c=otherDiag[0],d=otherDiag[1];
          if (a>b)
            std::swap(a,b);

          flipped.insert(a*n+b);
          flipped.insert(std::min(a,c)*n+std::max(a,c));
          flipped.insert(std::min(b,c)*n+std::max(b,c));
          flipped.insert(std::min(b,d)*n+std::max(b,d));
          flipped.insert(std::min(a,d)*n+std::max(a,d));

          edgeL[i]=otherDiag;
          assert(c<d);

          if (!toFlip) {
            toFlip=true;
            flipsPerRound.push_back(vi2());
          }
          t.flip(i);

          flipsPerRound[F].push_back({a,b});
          flipsThisRound++;
        }
      }
    }

    vi2 sortedEdges;
    for (int i=0;i<edgeL.size();i++) {
      int a=edgeL[i][0],b=edgeL[i][1];
      int edgeUid=a*n+b;

      i2 otherDiag=t.otherDiagonal(i);

      if (flipped.count(edgeUid)==0 && otherDiag[0]!=-1) {
        std::array<P,4> quad={points[a],points[otherDiag[0]],points[b],points[otherDiag[1]]};
        int currentDif=countIntersections(edgeL[i],edgeL1IncrX,points,intersectionsWithTarget)-
          countIntersections(otherDiag,edgeL1IncrX,points,intersectionsWithTarget);

        if (canBeFlipped(quad)
          && currentDif>0) {
          sortedEdges.push_back({currentDif,i});
          }
      }
    }

    std::sort(sortedEdges.rbegin(),sortedEdges.rend());
    for (i2 curEdge : sortedEdges) {
      int i=curEdge[1];
      int a=edgeL[i][0],b=edgeL[i][1];
      int edgeUid=a*n+b;

      i2 otherDiag=t.otherDiagonal(i);

      if (flipped.count(edgeUid)==0){
        int c=otherDiag[0],d=otherDiag[1];
        flipped.insert(edgeUid);
        flipped.insert(std::min(a,c)*n+std::max(a,c));
        flipped.insert(std::min(b,c)*n+std::max(b,c));
        flipped.insert(std::min(b,d)*n+std::max(b,d));
        flipped.insert(std::min(a,d)*n+std::max(a,d));

        edgeL[i]=otherDiag;

        if (!toFlip) {
          toFlip=true;
          flipsPerRound.push_back(vi2());
        }
        t.flip(i);

        flipsPerRound[F].push_back({a,b});
        flipsThisRound++;
      }
    }
    //
    // if (cntFlippable>0)
    //   std::cout<<"Flipped "<<flipsThisRound-prev<<" out of "<<cntFlippable<<'\n';

    if (toFlip)
      F++;
    totalFlips+=flipsThisRound;
  }

  mii intersectionsWithTarget1;
  std::vector<vi2> flipsPerRound1;
  F+=flipToReduceIntesectionsRevDir(t1,t,n,points,flipsPerRound1,intersectionsWithTarget1);
  for (const auto& flips1 : flipsPerRound1)
    flipsPerRound.push_back(flips1);
  return F;
}
int flipToReduceIntesectionsRevDir(facesDS& t, facesDS& t1, size_t n, const vP& points, std::vector<vi2>& flipsPerRound, mii& intersectionsWithTarget) {
  // flipsPerRound.clear();
  int F=0;
  //stop when no edge has been flipped in the previous round
  bool toFlip=true;

  int totalFlips=0;

  vi2 edgeL=t.getEdges();
  vi2 edgeL1=t1.getEdges();
  vi2 edgeL1IncrX(edgeL1.size());

  for (size_t i=0;i<edgeL1.size();i++){
    if (edgeL[i][0]>edgeL[i][1])
      std::swap(edgeL[i][0],edgeL[i][1]);
    // assert(edgeL[i][0]<edgeL[i][1]);
    // assert(t.otherDiagonal(i)[0]<=t.otherDiagonal(i)[1]);
    if (edgeL1[i][0]>edgeL1[i][1])
      std::swap(edgeL1[i][0],edgeL1[i][1]);
    edgeL1IncrX[i]=edgeL1[i];
  }
  std::sort(edgeL1.begin(),edgeL1.end());
  std::sort(edgeL1IncrX.begin(),edgeL1IncrX.end(),[points](i2 edgeA, i2 edgeB) {
    return std::min(points[edgeA[0]].x(),points[edgeA[1]].x())<std::min(points[edgeB[0]].x(),points[edgeB[1]].x());
  });

  while (toFlip && totalFlips<(n*n)) {
    toFlip=false;
    std::set<int> flipped;

    int flipsThisRound=0;


    Graph G0(3*n);
    // Graph G0(n);
    int cntFlippable=0;
    //todo always flip isolated edges
    for (int i=0;i<edgeL.size();i++) {
      int a=edgeL[i][0],b=edgeL[i][1];

      i2 otherDiag=t.otherDiagonal(i);
      if (otherDiag[0]!=-1
        && std::binary_search(edgeL1.begin(),edgeL1.end(),otherDiag)){

        std::array<P,4> quad={points[a],points[otherDiag[0]],points[b],points[otherDiag[1]]};
        if (canBeFlipped(quad)){
          i2 uv=t.incidentFaces(i);
          boost::add_edge(uv[0],uv[1],G0);
          G0[boost::edge(uv[0],uv[1],G0).first].idx=i;
          // boost::add_edge(a,b,G0);
          // G0[boost::edge(a,b,G0).first].idx=i;
          cntFlippable++;
        }

        }
    }

    if (cntFlippable>0) {
      size_t colors = boost::edge_coloring(G0, get(&EdgeProperty::color, G0));

      vi colorFreq(colors,0);

      boost::graph_traits<Graph>::edge_iterator ei, ei_end;
      for (boost::tie(ei, ei_end) = boost::edges(G0); ei != ei_end; ++ei)
        colorFreq[G0[boost::edge(boost::source(*ei, G0),boost::target(*ei,G0),G0).first].color]++;

      int maxFreq=colorFreq[0],bestColor=0;
      for (int i=1;i<colors;i++)
        if (colorFreq[i]>maxFreq) {
          maxFreq=colorFreq[i];
          bestColor=i;
        }
      for (boost::tie(ei, ei_end) = boost::edges(G0); ei != ei_end; ++ei) {
        int u=boost::source(*ei, G0), v=boost::target(*ei,G0);
        if (G0[boost::edge(u,v,G0).first].color==bestColor || (boost::out_degree(u,G0)<=1 && boost::out_degree(v,G0)<=1)) {
          int i=G0[boost::edge(u,v,G0).first].idx;
          int a=edgeL[i][0],b=edgeL[i][1];
          i2 otherDiag=t.otherDiagonal(i);
          int c=otherDiag[0],d=otherDiag[1];
          if (a>b)
            std::swap(a,b);

          flipped.insert(a*n+b);
          flipped.insert(std::min(a,c)*n+std::max(a,c));
          flipped.insert(std::min(b,c)*n+std::max(b,c));
          flipped.insert(std::min(b,d)*n+std::max(b,d));
          flipped.insert(std::min(a,d)*n+std::max(a,d));

          edgeL[i]=otherDiag;
          assert(c<d);

          if (!toFlip) {
            toFlip=true;
            flipsPerRound.push_back(vi2());
          }
          t.flip(i);

          flipsPerRound[F].push_back({c,d});
          flipsThisRound++;
        }
      }
    }


    vi2 sortedEdges;
    for (int i=0;i<edgeL.size();i++) {
      int a=edgeL[i][0],b=edgeL[i][1];
      int edgeUid=a*n+b;

      i2 otherDiag=t.otherDiagonal(i);

      if (flipped.count(edgeUid)==0 && otherDiag[0]!=-1) {
        std::array<P,4> quad={points[a],points[otherDiag[0]],points[b],points[otherDiag[1]]};
        int currentDif=countIntersections(edgeL[i],edgeL1IncrX,points,intersectionsWithTarget)-
          countIntersections(otherDiag,edgeL1IncrX,points,intersectionsWithTarget);

        if (canBeFlipped(quad)
          && currentDif>0) {
          sortedEdges.push_back({currentDif,i});
          }
      }
    }

    std::sort(sortedEdges.rbegin(),sortedEdges.rend());
    for (i2 curEdge : sortedEdges) {
      int i=curEdge[1];
      int a=edgeL[i][0],b=edgeL[i][1];
      int edgeUid=a*n+b;

      i2 otherDiag=t.otherDiagonal(i);

      if (flipped.count(edgeUid)==0){
        int c=otherDiag[0],d=otherDiag[1];
        flipped.insert(edgeUid);
        flipped.insert(std::min(a,c)*n+std::max(a,c));
        flipped.insert(std::min(b,c)*n+std::max(b,c));
        flipped.insert(std::min(b,d)*n+std::max(b,d));
        flipped.insert(std::min(a,d)*n+std::max(a,d));

        edgeL[i]=otherDiag;

        if (!toFlip) {
          toFlip=true;
          flipsPerRound.push_back(vi2());
        }
        t.flip(i);

        flipsPerRound[F].push_back({c,d});
        flipsThisRound++;
      }
    }

    if (toFlip)
      F++;
    totalFlips+=flipsThisRound;
  }

  std::reverse(flipsPerRound.begin(),flipsPerRound.end());
  return F;
}

std::vector<facesDS> flipToReduceIntesectionsWithHistory(facesDS& t, facesDS& t1, size_t n, const vP& points, std::vector<vi2>& flipsPerRound, mii& intersectionsWithTarget) {
  t.init();
  t1.init();
  std::vector<facesDS> intermediateTriangulations;
  int F=0;
  //stop when no edge has been flipped in the previous round
  bool toFlip=true;

  int totalFlips=0;

  vi2 edgeL=t.getEdges();
  vi2 edgeL1=t1.getEdges();
  vi2 edgeL1IncrX(edgeL1.size());

  for (size_t i=0;i<edgeL1.size();i++){
    // assert(edgeL[i][0]<edgeL[i][1]);
    if (edgeL[i][0]>edgeL[i][1])
      std::swap(edgeL[i][0],edgeL[i][1]);
    assert(t.otherDiagonal(i)[0]<=t.otherDiagonal(i)[1]);
    if (edgeL1[i][0]>edgeL1[i][1])
      std::swap(edgeL1[i][0],edgeL1[i][1]);
    edgeL1IncrX[i]=edgeL1[i];
  }
  std::sort(edgeL1.begin(),edgeL1.end());
  std::sort(edgeL1IncrX.begin(),edgeL1IncrX.end(),[points](i2 edgeA, i2 edgeB) {
    return std::min(points[edgeA[0]].x(),points[edgeA[1]].x())<std::min(points[edgeB[0]].x(),points[edgeB[1]].x());
  });

  // std::cout<<"Intersections in the same triangulation (should be 0): "<<countIntersections(edgeL1,edgeL1,points)<<'\n';
  // std::cout<<"Intersections between t0 and t1: "<<countIntersections(edgeL,edgeL1,points)<<'\n';

  while (toFlip && totalFlips<(n*n)) {
    toFlip=false;
    std::set<int> flipped;

    int flipsThisRound=0;


    Graph G0(3*n);
    int cntFlippable=0;
    for (int i=0;i<edgeL.size();i++) {
      int a=edgeL[i][0],b=edgeL[i][1];

      i2 otherDiag=t.otherDiagonal(i);
      if (otherDiag[0]!=-1
        && std::binary_search(edgeL1.begin(),edgeL1.end(),otherDiag)){

        std::array<P,4> quad={points[a],points[otherDiag[0]],points[b],points[otherDiag[1]]};
        if (canBeFlipped(quad)){
          i2 uv=t.incidentFaces(i);
          boost::add_edge(uv[0],uv[1],G0);
          G0[boost::edge(uv[0],uv[1],G0).first].idx=i;
          // boost::add_edge(a,b,G0);
          // G0[boost::edge(a,b,G0).first].idx=i;
          /* also, this version of the heuristic was kind of abandoned because usually improving the centers gives better results. So there's something I can try to improve in it.
           */

          cntFlippable++;
        }

        }
    }

    if (cntFlippable>0) {
      size_t colors = boost::edge_coloring(G0, get(&EdgeProperty::color, G0));

      vi colorFreq(colors,0);

      boost::graph_traits<Graph>::edge_iterator ei, ei_end;
      for (boost::tie(ei, ei_end) = boost::edges(G0); ei != ei_end; ++ei)
        colorFreq[G0[boost::edge(boost::source(*ei, G0),boost::target(*ei,G0),G0).first].color]++;

      int maxFreq=colorFreq[0],bestColor=0;
      for (int i=1;i<colors;i++)
        if (colorFreq[i]>maxFreq) {
          maxFreq=colorFreq[i];
          bestColor=i;
        }
      for (boost::tie(ei, ei_end) = boost::edges(G0); ei != ei_end; ++ei) {
        int u=boost::source(*ei, G0), v=boost::target(*ei,G0);
        if (G0[boost::edge(u,v,G0).first].color==bestColor || (boost::out_degree(u,G0)<=1 && boost::out_degree(v,G0)<=1)) {
          int i=G0[boost::edge(u,v,G0).first].idx;
          int a=edgeL[i][0],b=edgeL[i][1];
        // int a=boost::source(*ei, G0), b=boost::target(*ei,G0);
        // if (G0[boost::edge(a,b,G0).first].color==bestColor || (boost::out_degree(a,G0)==1 && boost::out_degree(b,G0)==1)) {
        //   int i=G0[boost::edge(a,b,G0).first].idx;
          i2 otherDiag=t.otherDiagonal(i);
          int c=otherDiag[0],d=otherDiag[1];
          if (a>b)
            std::swap(a,b);

          flipped.insert(a*n+b);
          flipped.insert(std::min(a,c)*n+std::max(a,c));
          flipped.insert(std::min(b,c)*n+std::max(b,c));
          flipped.insert(std::min(b,d)*n+std::max(b,d));
          flipped.insert(std::min(a,d)*n+std::max(a,d));

          edgeL[i]=otherDiag;
          assert(c<d);

          if (!toFlip) {
            toFlip=true;
            flipsPerRound.push_back(vi2());
          }
          t.flip(i);

          flipsPerRound[F].push_back({a,b});
          flipsThisRound++;
        }
      }
    }

    vi2 sortedEdges;
    for (int i=0;i<edgeL.size();i++) {
      int a=edgeL[i][0],b=edgeL[i][1];
      int edgeUid=a*n+b;

      i2 otherDiag=t.otherDiagonal(i);

      if (flipped.count(edgeUid)==0 && otherDiag[0]!=-1) {
        std::array<P,4> quad={points[a],points[otherDiag[0]],points[b],points[otherDiag[1]]};
        int currentDif=countIntersections(edgeL[i],edgeL1IncrX,points,intersectionsWithTarget)-
          countIntersections(otherDiag,edgeL1IncrX,points,intersectionsWithTarget);

        if (canBeFlipped(quad)
          && currentDif>0) {
          sortedEdges.push_back({currentDif,i});
          }
      }
    }

    std::sort(sortedEdges.rbegin(),sortedEdges.rend());
    for (i2 curEdge : sortedEdges) {
      int i=curEdge[1];
      int a=edgeL[i][0],b=edgeL[i][1];
      int edgeUid=a*n+b;

      i2 otherDiag=t.otherDiagonal(i);

      if (flipped.count(edgeUid)==0){
        int c=otherDiag[0],d=otherDiag[1];
        flipped.insert(edgeUid);
        flipped.insert(std::min(a,c)*n+std::max(a,c));
        flipped.insert(std::min(b,c)*n+std::max(b,c));
        flipped.insert(std::min(b,d)*n+std::max(b,d));
        flipped.insert(std::min(a,d)*n+std::max(a,d));

        edgeL[i]=otherDiag;

        if (!toFlip) {
          toFlip=true;
          flipsPerRound.push_back(vi2());
        }
        t.flip(i);

        flipsPerRound[F].push_back({a,b});
        flipsThisRound++;
      }
    }

    if (toFlip) {
      F++;
      intermediateTriangulations.push_back(facesDS(t));
    }
    totalFlips+=flipsThisRound;
  }

  return intermediateTriangulations;
}

std::vector<simpleDS> flipToReduceIntesectionsWithHistory(flippingDS& t, flippingDS& t1, size_t n, const vP& points, std::vector<vi2>& flipsPerRound, mii& intersectionsWithTarget) {
  t.init();
  t1.init();
  std::vector<simpleDS> intermediateTriangulations;
  int F=0;
  //stop when no edge has been flipped in the previous round
  bool toFlip=true;

  int totalFlips=0;

  vi2 edgeL=t.getEdges();
  vi2 edgeL1=t1.getEdges();
  vi2 edgeL1IncrX(edgeL1.size());

  for (size_t i=0;i<edgeL1.size();i++){
    assert(edgeL[i][0]<edgeL[i][1]);
    assert(t.otherDiagonal(i)[0]<=t.otherDiagonal(i)[1]);
    if (edgeL1[i][0]>edgeL1[i][1])
      std::swap(edgeL1[i][0],edgeL1[i][1]);
    edgeL1IncrX[i]=edgeL1[i];
  }
  std::sort(edgeL1.begin(),edgeL1.end());
  std::sort(edgeL1IncrX.begin(),edgeL1IncrX.end(),[points](i2 edgeA, i2 edgeB) {
    return std::min(points[edgeA[0]].x(),points[edgeA[1]].x())<std::min(points[edgeB[0]].x(),points[edgeB[1]].x());
  });

  // std::cout<<"Intersections in the same triangulation (should be 0): "<<countIntersections(edgeL1,edgeL1,points)<<'\n';
  // std::cout<<"Intersections between t0 and t1: "<<countIntersections(edgeL,edgeL1,points)<<'\n';

  while (toFlip && totalFlips<(n*n)) {
    toFlip=false;
    std::set<int> flipped;

    int flipsThisRound=0;


    Graph G0(n);
    int cntFlippable=0;
    for (int i=0;i<edgeL.size();i++) {
      int a=edgeL[i][0],b=edgeL[i][1];

      i2 otherDiag=t.otherDiagonal(i);
      if (otherDiag[0]!=-1
        && std::binary_search(edgeL1.begin(),edgeL1.end(),otherDiag)){

        std::array<P,4> quad={points[a],points[otherDiag[0]],points[b],points[otherDiag[1]]};
        if (canBeFlipped(quad)){
          boost::add_edge(a,b,G0);
          G0[boost::edge(a,b,G0).first].idx=i;
          cntFlippable++;
        }

        }
    }

    if (cntFlippable>0) {
      size_t colors = boost::edge_coloring(G0, get(&EdgeProperty::color, G0));

      vi colorFreq(colors,0);

      boost::graph_traits<Graph>::edge_iterator ei, ei_end;
      for (boost::tie(ei, ei_end) = boost::edges(G0); ei != ei_end; ++ei)
        colorFreq[G0[boost::edge(boost::source(*ei, G0),boost::target(*ei,G0),G0).first].color]++;

      int maxFreq=colorFreq[0],bestColor=0;
      for (int i=1;i<colors;i++)
        if (colorFreq[i]>maxFreq) {
          maxFreq=colorFreq[i];
          bestColor=i;
        }
      for (boost::tie(ei, ei_end) = boost::edges(G0); ei != ei_end; ++ei) {
        int a=boost::source(*ei, G0), b=boost::target(*ei,G0);
        if (G0[boost::edge(a,b,G0).first].color==bestColor || (boost::out_degree(a,G0)==1 && boost::out_degree(b,G0)==1)) {
          int i=G0[boost::edge(a,b,G0).first].idx;
          i2 otherDiag=t.otherDiagonal(i);
          int c=otherDiag[0],d=otherDiag[1];
          if (a>b)
            std::swap(a,b);

          flipped.insert(a*n+b);
          flipped.insert(std::min(a,c)*n+std::max(a,c));
          flipped.insert(std::min(b,c)*n+std::max(b,c));
          flipped.insert(std::min(b,d)*n+std::max(b,d));
          flipped.insert(std::min(a,d)*n+std::max(a,d));

          edgeL[i]=otherDiag;
          assert(c<d);

          if (!toFlip) {
            toFlip=true;
            flipsPerRound.push_back(vi2());
          }
          t.flip(i);

          flipsPerRound[F].push_back({a,b});
          flipsThisRound++;
        }
      }
    }

    for (int i=0;i<edgeL.size();i++) {
      int a=edgeL[i][0],b=edgeL[i][1];
      // assert(a<b);
      int edgeUid=a*n+b;

      i2 otherDiag=t.otherDiagonal(i);

      if (flipped.count(edgeUid)==0 && otherDiag[0]!=-1
        // && !std::binary_search(edgeL1.begin(),edgeL1.end(),otherDiag)
        // && countIntersections(otherDiag,edgeL1IncrX,points,intersectionsWithTarget)<countIntersections(edgeL[i],edgeL1IncrX,points,intersectionsWithTarget)
        ){
        std::array<P,4> quad={points[a],points[otherDiag[0]],points[b],points[otherDiag[1]]};

        if (canBeFlipped(quad)
          && countIntersections(otherDiag,edgeL1IncrX,points,intersectionsWithTarget)<countIntersections(edgeL[i],edgeL1IncrX,points,intersectionsWithTarget)) {
          int c=otherDiag[0],d=otherDiag[1];
          flipped.insert(edgeUid);
          flipped.insert(std::min(a,c)*n+std::max(a,c));
          flipped.insert(std::min(b,c)*n+std::max(b,c));
          flipped.insert(std::min(b,d)*n+std::max(b,d));
          flipped.insert(std::min(a,d)*n+std::max(a,d));

          edgeL[i]=otherDiag;

          if (!toFlip) {
            toFlip=true;
            flipsPerRound.push_back(vi2());
          }
          t.flip(i);

          flipsPerRound[F].push_back({a,b});
          flipsThisRound++;
        }

      }
    }

    if (toFlip) {
      F++;
      intermediateTriangulations.push_back(simpleDS(t));
    }
    totalFlips+=flipsThisRound;
  }

  return intermediateTriangulations;
}

int improveCenter(facesDS& center, std::vector<facesDS>& tVec, size_t n, const vP& points, std::vector<vi2>& flipsPerRound) {
  mii intersectionsWithTarget;
  std::vector<mii> intersectionsWithTargetOfEachTri(tVec.size());
  int F=0;
  //stop when no edge has been flipped in the previous round
  bool toFlip=true;

  int totalFlips=0;

  vi2 edgeL=center.getEdges();
  std::vector<vi2> edgeVecOther(tVec.size());
  for (size_t i=0;i<edgeVecOther.size();i++) {
    vi2 edgeL1=tVec[i].getEdges();
    for (size_t j=0;j<edgeL1.size();j++)
      if (edgeL1[i][0]>edgeL1[i][1])
        edgeL1[i]={edgeL1[i][1],edgeL1[i][0]};
    std::sort(edgeL1.begin(),edgeL1.end(),[points](i2 edgeA, i2 edgeB) {
    return std::min(points[edgeA[0]].x(),points[edgeA[1]].x())<std::min(points[edgeB[0]].x(),points[edgeB[1]].x());
  });
    edgeVecOther[i]=edgeL1;
  }

  for (size_t i=0;i<edgeL.size();i++){
    // assert(edgeL[i][0]<edgeL[i][1]);
    if (edgeL[i][0]<edgeL[i][1])
      edgeL[i]={edgeL[i][1],edgeL[i][0]};
    assert(center.otherDiagonal(i)[0]<=center.otherDiagonal(i)[1]);
  }
  // std::sort(edgeL.begin(),edgeL.end());


  // std::cout<<"Intersections in the same triangulation (should be 0): "<<countIntersections(edgeL1,edgeL1,points)<<'\n';
  // std::cout<<"Intersections between t0 and t1: "<<countIntersections(edgeL,edgeL1,points)<<'\n';
  // toFlip=true;

  while (toFlip && totalFlips<(n*n)) {
    toFlip=false;
    std::set<int> flipped;

    int flipsThisRound=0;


    for (int i=0;i<edgeL.size();i++) {
      int a=edgeL[i][0],b=edgeL[i][1];
      int edgeUid=a*n+b;

      i2 otherDiag=center.otherDiagonal(i);

      if (flipped.count(edgeUid)==0 && otherDiag[0]!=-1){
        std::array<P,4> quad={points[a],points[otherDiag[0]],points[b],points[otherDiag[1]]};

        if (canBeFlipped(quad)
          // && compareItsReductionCardinality(edgeL[i],otherDiag,edgeVecOther,points,intersectionsWithTargetOfEachTri)
          && countIntersections(otherDiag,edgeVecOther,points,intersectionsWithTarget)<countIntersections(edgeL[i],edgeVecOther,points,intersectionsWithTarget)
          ) {
          int c=otherDiag[0],d=otherDiag[1];
          flipped.insert(edgeUid);
          flipped.insert(std::min(a,c)*n+std::max(a,c));
          flipped.insert(std::min(b,c)*n+std::max(b,c));
          flipped.insert(std::min(b,d)*n+std::max(b,d));
          flipped.insert(std::min(a,d)*n+std::max(a,d));

          edgeL[i]=otherDiag;

          if (!toFlip) {
            toFlip=true;
            flipsPerRound.push_back(vi2());
          }
          center.flip(i);

          flipsPerRound[F].push_back({a,b});
          flipsThisRound++;
        }

      }
    }
    if (toFlip)
      F++;
    totalFlips+=flipsThisRound;
  }

  toFlip=false;

  while (toFlip && totalFlips<(n*n)) {
    toFlip=false;
    std::set<int> flipped;

    int flipsThisRound=0;

    for (int i=0;i<edgeL.size();i++) {
      int a=edgeL[i][0],b=edgeL[i][1];
      int edgeUid=a*n+b;

      i2 otherDiag=center.otherDiagonal(i);

      if (flipped.count(edgeUid)==0 && otherDiag[0]!=-1){
        std::array<P,4> quad={points[a],points[otherDiag[0]],points[b],points[otherDiag[1]]};

        if (canBeFlipped(quad)
          // && compareItsReductionCardinality(edgeL[i],otherDiag,edgeVecOther,points,intersectionsWithTargetOfEachTri)
          && countIntersections(otherDiag,edgeVecOther,points,intersectionsWithTarget)<countIntersections(edgeL[i],edgeVecOther,points,intersectionsWithTarget)
          ) {
          int c=otherDiag[0],d=otherDiag[1];
          flipped.insert(edgeUid);
          flipped.insert(std::min(a,c)*n+std::max(a,c));
          flipped.insert(std::min(b,c)*n+std::max(b,c));
          flipped.insert(std::min(b,d)*n+std::max(b,d));
          flipped.insert(std::min(a,d)*n+std::max(a,d));

          edgeL[i]=otherDiag;

          if (!toFlip) {
            toFlip=true;
            flipsPerRound.push_back(vi2());
          }
          center.flip(i);

          flipsPerRound[F].push_back({a,b});
          flipsThisRound++;
          }

      }
    }

    for (int i=0;i<edgeL.size();i++) {
      int a=edgeL[i][0],b=edgeL[i][1];
      int edgeUid=a*n+b;

      i2 otherDiag=center.otherDiagonal(i);

      if (flipped.count(edgeUid)==0 && otherDiag[0]!=-1){
        std::array<P,4> quad={points[a],points[otherDiag[0]],points[b],points[otherDiag[1]]};

        if (canBeFlipped(quad)
          && compareItsReductionCardinality(edgeL[i],otherDiag,edgeVecOther,points,intersectionsWithTargetOfEachTri)
          // && countIntersections(otherDiag,edgeVecOther,points,intersectionsWithTarget)<countIntersections(edgeL[i],edgeVecOther,points,intersectionsWithTarget)
          ) {
          int c=otherDiag[0],d=otherDiag[1];
          flipped.insert(edgeUid);
          flipped.insert(std::min(a,c)*n+std::max(a,c));
          flipped.insert(std::min(b,c)*n+std::max(b,c));
          flipped.insert(std::min(b,d)*n+std::max(b,d));
          flipped.insert(std::min(a,d)*n+std::max(a,d));

          edgeL[i]=otherDiag;

          if (!toFlip) {
            toFlip=true;
            flipsPerRound.push_back(vi2());
          }
          center.flip(i);

          flipsPerRound[F].push_back({a,b});
          flipsThisRound++;
          }

      }
    }
    if (toFlip)
      F++;
    totalFlips+=flipsThisRound;
  }

  return F;
}

int omniDirFlips(int centerIdx, std::vector<facesDS>& tVec, size_t n, const vP& points, std::vector<std::vector<vi2>>& flipsPerRound) {
  int F=0;
  //stop when no edge has been flipped in the previous round
  bool toFlip=true;

  int totalFlips=0;

  std::vector<vi2> edgeListVector(tVec.size());
  for (size_t i=0;i<edgeListVector.size();i++) {
    vi2 edgeL1=tVec[i].getEdges();
  //   for (size_t j=0;j<edgeL1.size();j++)
  //     if (edgeL1[i][0]>edgeL1[i][1])
  //       edgeL1[i]={edgeL1[i][1],edgeL1[i][0]};
  //   std::sort(edgeL1.begin(),edgeL1.end(),[points](i2 edgeA, i2 edgeB) {
  //   return std::min(points[edgeA[0]].x(),points[edgeA[1]].x())<std::min(points[edgeB[0]].x(),points[edgeB[1]].x());
  // });
    edgeListVector[i]=edgeL1;
  }

  // for (size_t i=0;i<edgeL.size();i++){
  //   assert(edgeL[i][0]<edgeL[i][1]);
  //   assert(center.otherDiagonal(i)[0]<=center.otherDiagonal(i)[1]);
  // }
  // std::sort(edgeL.begin(),edgeL.end());


  // std::cout<<"Intersections in the same triangulation (should be 0): "<<countIntersections(edgeL1,edgeL1,points)<<'\n';
  // std::cout<<"Intersections between t0 and t1: "<<countIntersections(edgeL,edgeL1,points)<<'\n';

  while (toFlip && totalFlips<(n*n*n*tVec.size())) {

    int flipsThisRound=0;

    while (toFlip) {
      std::set<int> flipped;

      toFlip=false;
      for (int i=0;i<edgeListVector[centerIdx].size();i++) {
        int a=edgeListVector[centerIdx][i][0],b=edgeListVector[centerIdx][i][1];
        int edgeUid=a*n+b;

        i2 otherDiag=tVec[centerIdx].otherDiagonal(i);

        if (flipped.count(edgeUid)==0 && otherDiag[0]!=-1){
          std::array<P,4> quad={points[a],points[otherDiag[0]],points[b],points[otherDiag[1]]};

          if (canBeFlipped(quad)
            && countIntersections(otherDiag,edgeListVector,points)<countIntersections(edgeListVector[centerIdx][i],edgeListVector,points)) {
            int c=otherDiag[0],d=otherDiag[1];
            flipped.insert(edgeUid);
            flipped.insert(std::min(a,c)*n+std::max(a,c));
            flipped.insert(std::min(b,c)*n+std::max(b,c));
            flipped.insert(std::min(b,d)*n+std::max(b,d));
            flipped.insert(std::min(a,d)*n+std::max(a,d));

            edgeListVector[centerIdx][i]=otherDiag;

            if (!toFlip) {
              toFlip=true;
              flipsPerRound[centerIdx].push_back(vi2());
            }
            tVec[centerIdx].flip(i);

            flipsPerRound[centerIdx][flipsPerRound[centerIdx].size()-1].push_back({a,b});
            flipsThisRound++;
            }

        }
      }
      if (toFlip)
        F++;
      break;
    }

    toFlip=false;

    for (size_t triIdx=0;triIdx<tVec.size();triIdx++) {
      if (triIdx!=centerIdx) {
        bool thisToFlip=false;
        std::set<int> flipped;

        for (int i=0;i<edgeListVector[triIdx].size();i++) {
          int a=edgeListVector[triIdx][i][0],b=edgeListVector[triIdx][i][1];
          int edgeUid=a*n+b;

          i2 otherDiag=tVec[triIdx].otherDiagonal(i);

          if (flipped.count(edgeUid)==0 && otherDiag[0]!=-1){
            std::array<P,4> quad={points[a],points[otherDiag[0]],points[b],points[otherDiag[1]]};

            if (canBeFlipped(quad)
              && countIntersections(otherDiag,edgeListVector[centerIdx],points)==0 &&
              countIntersections(edgeListVector[triIdx][i],edgeListVector[centerIdx],points)==1) {
              int c=otherDiag[0],d=otherDiag[1];
              flipped.insert(edgeUid);
              flipped.insert(std::min(a,c)*n+std::max(a,c));
              flipped.insert(std::min(b,c)*n+std::max(b,c));
              flipped.insert(std::min(b,d)*n+std::max(b,d));
              flipped.insert(std::min(a,d)*n+std::max(a,d));

              edgeListVector[triIdx][i]=otherDiag;

              if (!thisToFlip) {
                toFlip=true;
                thisToFlip=true;
                F++;
                flipsPerRound[triIdx].push_back(vi2());
              }
              tVec[triIdx].flip(i);

              flipsPerRound[triIdx][flipsPerRound[triIdx].size()-1].push_back({a,b});
              flipsThisRound++;
              }

          }
        }

        vi2 sortedEdges;

        for (int i=0;i<edgeListVector[triIdx].size();i++) {
          int a=edgeListVector[triIdx][i][0],b=edgeListVector[triIdx][i][1];
          int edgeUid=a*n+b;

          i2 otherDiag=tVec[triIdx].otherDiagonal(i);

          if (flipped.count(edgeUid)==0 && otherDiag[0]!=-1) {
            std::array<P,4> quad={points[a],points[otherDiag[0]],points[b],points[otherDiag[1]]};
            if (canBeFlipped(quad)) {
              int newIts=countIntersections(otherDiag,edgeListVector[centerIdx],points);
              int currentIts=countIntersections(edgeListVector[triIdx][i],edgeListVector[centerIdx],points);
              int itsDiff=currentIts-newIts;

              if (itsDiff>0) {
                sortedEdges.push_back({itsDiff,i});
              }
            }
          }
        }

        std::sort(sortedEdges.rbegin(),sortedEdges.rend());

        for (i2 curEdge : sortedEdges) {
          int i=curEdge[1];
          int a=edgeListVector[triIdx][i][0],b=edgeListVector[triIdx][i][1];
          int edgeUid=a*n+b;

          i2 otherDiag=tVec[triIdx].otherDiagonal(i);

          if (flipped.count(edgeUid)==0){
            int c=otherDiag[0],d=otherDiag[1];
            flipped.insert(edgeUid);
            flipped.insert(std::min(a,c)*n+std::max(a,c));
            flipped.insert(std::min(b,c)*n+std::max(b,c));
            flipped.insert(std::min(b,d)*n+std::max(b,d));
            flipped.insert(std::min(a,d)*n+std::max(a,d));

            edgeListVector[triIdx][i]=otherDiag;

            if (!thisToFlip) {
              toFlip=true;
              thisToFlip=true;
              F++;
              flipsPerRound[triIdx].push_back(vi2());
            }
            tVec[triIdx].flip(i);

            flipsPerRound[triIdx][flipsPerRound[triIdx].size()-1].push_back({a,b});
            flipsThisRound++;
            }

        }


        // for (int i=0;i<edgeListVector[triIdx].size();i++) {
        //   int a=edgeListVector[triIdx][i][0],b=edgeListVector[triIdx][i][1];
        //   int edgeUid=a*n+b;
        //
        //   i2 otherDiag=tVec[triIdx].otherDiagonal(i);
        //
        //   if (flipped.count(edgeUid)==0 && otherDiag[0]!=-1){
        //     std::array<P,4> quad={points[a],points[otherDiag[0]],points[b],points[otherDiag[1]]};
        //
        //     if (canBeFlipped(quad)
        //       && countIntersections(otherDiag,edgeListVector[centerIdx],points)<countIntersections(edgeListVector[triIdx][i],edgeListVector[centerIdx],points)) {
        //       int c=otherDiag[0],d=otherDiag[1];
        //       flipped.insert(edgeUid);
        //       flipped.insert(std::min(a,c)*n+std::max(a,c));
        //       flipped.insert(std::min(b,c)*n+std::max(b,c));
        //       flipped.insert(std::min(b,d)*n+std::max(b,d));
        //       flipped.insert(std::min(a,d)*n+std::max(a,d));
        //
        //       edgeListVector[triIdx][i]=otherDiag;
        //
        //       if (!thisToFlip) {
        //         toFlip=true;
        //         thisToFlip=true;
        //         F++;
        //         flipsPerRound[triIdx].push_back(vi2());
        //       }
        //       tVec[triIdx].flip(i);
        //
        //       flipsPerRound[triIdx][flipsPerRound[triIdx].size()-1].push_back({a,b});
        //       flipsThisRound++;
        //       }
        //
        //   }
        // }
      }
    }

   //  bool thisToFlip=false;
   // while (thisToFlip) {
   //    std::set<int> flipped;
   //
   //    thisToFlip=false;
   //    for (int i=0;i<edgeListVector[centerIdx].size();i++) {
   //      int a=edgeListVector[centerIdx][i][0],b=edgeListVector[centerIdx][i][1];
   //      int edgeUid=a*n+b;
   //
   //      i2 otherDiag=tVec[centerIdx].otherDiagonal(i);
   //
   //      if (flipped.count(edgeUid)==0 && otherDiag[0]!=-1){
   //        std::array<P,4> quad={points[a],points[otherDiag[0]],points[b],points[otherDiag[1]]};
   //
   //        if (canBeFlipped(quad)
   //          && countIntersections(otherDiag,edgeListVector,points)<countIntersections(edgeListVector[centerIdx][i],edgeListVector,points)) {
   //          int c=otherDiag[0],d=otherDiag[1];
   //          flipped.insert(edgeUid);
   //          flipped.insert(std::min(a,c)*n+std::max(a,c));
   //          flipped.insert(std::min(b,c)*n+std::max(b,c));
   //          flipped.insert(std::min(b,d)*n+std::max(b,d));
   //          flipped.insert(std::min(a,d)*n+std::max(a,d));
   //
   //          edgeListVector[centerIdx][i]=otherDiag;
   //
   //          if (!thisToFlip) {
   //            toFlip=true;
   //            thisToFlip=true;
   //            flipsPerRound[centerIdx].push_back(vi2());
   //          }
   //          tVec[centerIdx].flip(i);
   //
   //          flipsPerRound[centerIdx][flipsPerRound[centerIdx].size()-1].push_back({a,b});
   //          flipsThisRound++;
   //          }
   //
   //      }
   //    }
   //    if (thisToFlip)
   //      F++;
   //    // break;
   //  }

    totalFlips+=flipsThisRound;
    // std::cout<<F<<std::endl;

  }
  // std::cout<<F<<std::endl;

  return F;
}
// int flipToReduceIntesections(flippingDS& t, flippingDS& t1, size_t n, const vP& points, std::vector<vi2>& flipsPerRound) {
//   t.init();
//   t1.init();
//   int F=0;
//   //stop when no edge has been flipped in the previous round
//   bool toFlip=true;
//
//   int totalFlips=0;
//
//   vi2 edgeL=t.getEdges();
//   vi2 edgeL1=t1.getEdges();
//   vi2 edgeL1IncrX(edgeL1.size());
//
//   for (size_t i=0;i<edgeL1.size();i++){
//     assert(edgeL[i][0]<edgeL[i][1]);
//     assert(t.otherDiagonal(i)[0]<=t.otherDiagonal(i)[1]);
//     if (edgeL1[i][0]>edgeL1[i][1])
//       std::swap(edgeL1[i][0],edgeL1[i][1]);
//     edgeL1IncrX[i]=edgeL1[i];
//   }
//   std::sort(edgeL1.begin(),edgeL1.end());
//   std::sort(edgeL1IncrX.begin(),edgeL1IncrX.end(),[points](i2 edgeA, i2 edgeB) {
//     return std::min(points[edgeA[0]].x(),points[edgeA[1]].x())<std::min(points[edgeB[0]].x(),points[edgeB[1]].x());
//   });
//
//   // std::cout<<"Intersections in the same triangulation (should be 0): "<<countIntersections(edgeL1,edgeL1,points)<<'\n';
//   // std::cout<<"Intersections between t0 and t1: "<<countIntersections(edgeL,edgeL1,points)<<'\n';
//
//   while (toFlip && totalFlips<(n*n)) {
//     toFlip=false;
//     std::set<int> flipped;
//
//     int flipsThisRound=0;
//
//     for (int i=0;i<edgeL.size();i++) {
//       int a=edgeL[i][0],b=edgeL[i][1];
//       int edgeUid=a*n+b;
//
//       i2 otherDiag=t.otherDiagonal(i);
//       if (flipped.count(edgeUid)==0 && otherDiag[0]!=-1 && !std::binary_search(edgeL1.begin(),edgeL1.end(),edgeL[i])) {
//         if (std::binary_search(edgeL1.begin(),edgeL1.end(),otherDiag)){
//           std::array<P,4> quad={points[a],points[otherDiag[0]],points[b],points[otherDiag[1]]};
//
//           if (canBeFlipped(quad)) {
//             int c=otherDiag[0],d=otherDiag[1];
//             flipped.insert(edgeUid);
//             flipped.insert(std::min(a,c)*n+std::max(a,c));
//             flipped.insert(std::min(b,c)*n+std::max(b,c));
//             flipped.insert(std::min(b,d)*n+std::max(b,d));
//             flipped.insert(std::min(a,d)*n+std::max(a,d));
//
//             edgeL[i]=otherDiag;
//             assert(c<d);
//
//             if (!toFlip) {
//               toFlip=true;
//               flipsPerRound.push_back(vi2());
//             }
//             t.flip(i);
//
//             flipsPerRound[F].push_back({a,b});
//             flipsThisRound++;
//           }
//         }
//       }
//     }
//
//     for (int i=0;i<edgeL.size();i++) {
//       int a=edgeL[i][0],b=edgeL[i][1];
//       assert(a<b);
//       int edgeUid=a*n+b;
//
//       i2 otherDiag=t.otherDiagonal(i);
//
//       if (flipped.count(edgeUid)==0 && otherDiag[0]!=-1
//         && countIntersections(otherDiag,edgeL1IncrX,points)<countIntersections(edgeL[i],edgeL1IncrX,points)){
//         std::array<P,4> quad={points[a],points[otherDiag[0]],points[b],points[otherDiag[1]]};
//
//         if (canBeFlipped(quad)) {
//           int c=otherDiag[0],d=otherDiag[1];
//           flipped.insert(edgeUid);
//           flipped.insert(std::min(a,c)*n+std::max(a,c));
//           flipped.insert(std::min(b,c)*n+std::max(b,c));
//           flipped.insert(std::min(b,d)*n+std::max(b,d));
//           flipped.insert(std::min(a,d)*n+std::max(a,d));
//
//           edgeL[i]=otherDiag;
//
//           if (!toFlip) {
//             toFlip=true;
//             flipsPerRound.push_back(vi2());
//           }
//           t.flip(i);
//
//           flipsPerRound[F].push_back({a,b});
//           flipsThisRound++;
//         }
//
//       }
//     }
//
//     if (toFlip)
//       F++;
//     totalFlips+=flipsThisRound;
//   }
//
//   return F;
// }