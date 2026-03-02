#include "../include/quadrilateralGraph.hpp"

#include <CGAL/Polygon_2.h>
#include <CGAL/convex_hull_2.h>
typedef CGAL::Polygon_2<IK> Polygon_2;
typedef CGAL::Triangle_2<IK> Triangle_2;

typedef IK::Segment_2 S;

vvi crossingEdges(vP points) {
  int n=points.size();
  vvi crossingEdgesVec(n*n,vi());
  int cntCrossings=0;
  for (int i=0;i<n;i++)
    for (int j=i+1;j<n;j++)
      for (int k=j+1;k<n;k++)
        for (int l=k+1;l<n;l++) {
          S e0(points[i],points[j]);
          S e1(points[k],points[l]);
          if (CGAL::do_intersect(e0,e1)) {
            crossingEdgesVec[i*n+j].push_back(k*n+l);
            cntCrossings++;
          }
          S e2(points[i],points[k]);
          S e3(points[j],points[l]);
          if (CGAL::do_intersect(e2,e3)) {
            crossingEdgesVec[i*n+k].push_back(j*n+l);
            cntCrossings++;
          }
          S e4(points[i],points[l]);
          S e5(points[j],points[k]);
          if (CGAL::do_intersect(e4,e5)) {
            crossingEdgesVec[i*n+l].push_back(j*n+k);
            cntCrossings++;
          }
        }
  std::cout<<"There are "<<cntCrossings<<" crossing edges"<<std::endl;
  return crossingEdgesVec;
}
void buildEmptyTriangle(const vP& points, int n, vvvb& emptyTriangle) {
  for (int i=0;i<n;i++)
    for (int j=i+1;j<n;j++)
      for (int k=j+1;k<n;k++) {
        Triangle_2 t(points[i],points[j],points[k]);
        //todo you can speed up this part with some sorting (remember to store the initial indices)
        for (int e=0;e<n && emptyTriangle[i][j][k];e++) {
          if (t.bounded_side(points[e])==CGAL::ON_BOUNDED_SIDE) {
            emptyTriangle[i][j][k]=false;
          }
        }
        emptyTriangle[i][k][j]=emptyTriangle[i][j][k];
        emptyTriangle[j][i][k]=emptyTriangle[i][j][k];
        emptyTriangle[j][k][i]=emptyTriangle[i][j][k];
        emptyTriangle[k][i][j]=emptyTriangle[i][j][k];
        emptyTriangle[k][j][i]=emptyTriangle[i][j][k];
      }
}
vvi buildQGPlusNonConvex(vP points, const std::string& fileName, vvvb& emptyTriangle) {
  std::ifstream input(fileName);
  if (input.is_open()) {
    int n=points.size();
    vvi qg(n*n);
    for (int i=0;i<n*n;i++) {
      int degree;
      input>>degree;
      qg[i]=vi(degree);
      for (int j=0;j<degree;j++) {
        input>>qg[i][j];
      }
    }

    bool tmp;
    input>>tmp;
    if (input.eof()) {
      buildEmptyTriangle(points,n,emptyTriangle);
      std::ofstream output(fileName,std::ios::app);
      for (int i=0;i<n;i++)
        for (int j=0;j<n;j++) {
          for (int k=0;k<n;k++) {
            output<<emptyTriangle[i][j][k]<<' ';
          }
          output<<'\n';
        }
      output.close();
    }
    else {
      for (int i=0;i<n;i++)
        for (int j=0;j<n;j++)
          for (int k=0;k<n;k++) {
            emptyTriangle[i][j][k]=tmp;
            input>>tmp;
          }
    }
    input.close();
    return qg;
  }

  int n=points.size();
  vvi qg(n*n,vi());
  //to speed up emptyness check, precompute triangles with points inside them
  //(you need to check the right triangles for non convex quadrilaterals)
  // vvvb emptyTriangle(n,vvb(n,vb(n,true)));
  vvb emptySegment(n,vb(n,true));
  for (int i=0;i<n;i++)
    for (int j=i+1;j<n;j++) {
      S s0(points[i],points[j]);
      for (int k=0;k<n && emptySegment[i][j];k++) {
        if (k!=i && k!=j && s0.has_on(points[k])) {
          emptySegment[i][j]=false;
          emptySegment[j][i]=false;
        }
      }
    }
  buildEmptyTriangle(points,n,emptyTriangle);
  int cntEmptyQuad=0;
  for (int i=0;i<n;i++)
    for (int j=i+1;j<n;j++)
      for (int k=j+1;k<n;k++)
        //think about speeding this up by using emptyTriangle (be careful with non convex quadr)
        for (int l=k+1;l<n;l++) {

          //sorts j,k,l so that i,j,k,l are in convex order (tries all 6 permutations until the right one is found
          int indices[3]={j,k,l};
          do {
            Polygon_2 quadr;
            quadr.push_back(points[i]);
            for (int idx : indices)
              quadr.push_back(points[idx]);

            if (quadr.is_simple() && quadr.is_counterclockwise_oriented()) {
              //compute the diagonals (to linearize indices, the vertex indices must be in increasing order)
              int d00,d01,d10,d11;
              if (i<indices[1]) {
                d00=i;
                d01=indices[1];
              }
              else {
                d00=indices[1];
                d01=i;
              }
              if (indices[0]<indices[2]) {
                d10=indices[0];
                d11=indices[2];
              }
              else {
                d10=indices[2];
                d11=indices[0];
              }
              bool d0MidpointBoundedSide=quadr.bounded_side(CGAL::midpoint(points[d00],points[d01]))==CGAL::ON_BOUNDED_SIDE;
              bool d1MidpointBoundedSide=quadr.bounded_side(CGAL::midpoint(points[d10],points[d11]))==CGAL::ON_BOUNDED_SIDE;

              std::array<int,3> t0={d00,d01,d10},t1={d00,d01,d11};
              if (!d0MidpointBoundedSide) {
                t0={d10,d11,d00};
                t1={d10,d11,d01};
              }
              bool empty=emptyTriangle[t0[0]][t0[1]][t0[2]] && emptyTriangle[t1[0]][t1[1]][t1[2]] &&
                emptySegment[t0[0]][t0[1]];
              // bool empty=true;
              // for (int e=0;e<n && empty;e++) {
              //   if (quadr.bounded_side(points[e])==CGAL::ON_BOUNDED_SIDE) {
              //     empty=false;
              //   }
              // }

              // if (i==1 && j==3 && k==11 && l==18)
              //   std::cout<<empty<<'\n';
              if (empty) {
                // std::cout<<"quadr: "<<i<<' ';
                // for (int idx : indices)
                //   std::cout<<idx<<' ';
                // std::cout<<'\n';
                cntEmptyQuad++;
                // if (i==2 && j==10 && k==11 && l==16)

                // assert(d00<n);
                // assert(d01<n);
                // assert(d10<n);
                // assert(d11<n);
                int d0=d00*n+d01;
                int d1=d10*n+d11;
                // assert(d0<n*n);
                // assert(d1<n*n);
                // if (i==0 && j==25 && k==42 && l==48) {
                //   std::cout<<d00<<' '<<d01<<' '<<d10<<' '<<d11<<'\n';
                //   std::cout<<quadr.is_convex()<<'\n';
                //   std::cout<<d0MidpointBoundedSide<<' '<<
                //     d1MidpointBoundedSide<<'\n';
                // }
                if (quadr.is_convex() &&
                  d0MidpointBoundedSide &&
                  d1MidpointBoundedSide
                  ) {
                  qg[d0].push_back(d1);
                  qg[d1].push_back(d0);
                }
                else {
                  // std::cout<<d00<<' '<<d01<<' '<<d10<<' '<<d11<<'\n';
                  if (d0MidpointBoundedSide)
                    qg[d0].push_back(-d1-1);
                  if (d1MidpointBoundedSide)
                    qg[d1].push_back(-d0-1);

                }
              }
              // else {
              //   std::cout<<"Point inside ";
              //   std::cout<<"quadr: "<<i<<' ';
              //   for (int idx : indices)
              //     std::cout<<idx<<' ';
              //   std::cout<<'\n';
              // }
            }


          }while (std::next_permutation(indices,indices+3));



          // if (i==1 && j==3 && k==11 && l==18) {
          //   std::cout<<"1 3 11 18 quadrilateral: "<<i<<' ';
          //   for (int idx : indices)
          //     std::cout<<idx<<' ';
          //   std::cout<<'\n';
          // }

        }
  std::cout<<"There are "<<cntEmptyQuad<<" empty (not necessarily convex) quadrilaterals in the point set\n";
  std::ofstream output(fileName);
  if (output.is_open()) {
    for (int i=0;i<n*n;i++) {
      output<<qg[i].size()<<' ';
      for (int otherDiag : qg[i])
        output<<otherDiag<<' ';
      output<<'\n';
    }
    for (int i=0;i<n;i++)
      for (int j=0;j<n;j++) {
        for (int k=0;k<n;k++) {
          output<<emptyTriangle[i][j][k]<<' ';
        }
        output<<'\n';
      }
    output.close();
  }
  return qg;
}

vvi buildQG(vP points) {
  int n=points.size();
  vvi qg(n*n,vi());
  int cntEmptyQuad=0;
  for (int i=0;i<n;i++)
    for (int j=i+1;j<n;j++)
      for (int k=j+1;k<n;k++)
        for (int l=k+1;l<n;l++) {

          //sorts j,k,l so that i,j,k,l are in convex order (tries all 6 permutations until the right one is found
          int indices[3]={j,k,l};
          bool found=false;
          do {
            //todo is the third check necessary? from what I tried it shouldn't be
            if (CGAL::left_turn(points[i],points[indices[0]],points[indices[1]]) &&
              CGAL::left_turn(points[indices[0]],points[indices[1]],points[indices[2]])
              && CGAL::left_turn(points[indices[1]],points[indices[2]],points[i])
              )
              found=true;
            if (i==1 && j==3 && k==11 && l==18) {
              std::cout<<"1 3 11 18 current permutation: "<<i<<' ';
              for (int idx : indices)
                std::cout<<idx<<' ';
              std::cout<<'\n';
              std::cout<<CGAL::orientation(points[i],points[indices[0]],points[indices[1]])<<'\n';
              std::cout<<CGAL::orientation(points[indices[0]],points[indices[1]],points[indices[2]])<<'\n';
              std::cout<<CGAL::orientation(points[indices[1]],points[indices[2]],points[i])<<'\n';
              if (found)
                std::cout<<"found\n";
            }

          }while (!found && std::next_permutation(indices,indices+3));

          Polygon_2 quadr;
          quadr.push_back(points[i]);
          for (int idx : indices)
            quadr.push_back(points[idx]);

          if (i==1 && j==3 && k==11 && l==18) {
            std::cout<<"1 3 11 18 quadrilateral: "<<i<<' ';
            for (int idx : indices)
              std::cout<<idx<<' ';
            std::cout<<'\n';
          }
          // if (i==2 && j==10 && k==11 && l==16) {
          //   //   indices[0]=10;
          //   // indices[1]=11;
          //   // indices[2]=16;
          //     std::cout<<"2 10 11 16 quadrilateral: "<<i<<' ';
          //     for (int idx : indices)
          //       std::cout<<idx<<' ';
          //     std::cout<<'\n';
          //   // quadr.clear();
          //   // quadr.push_back(points[i]);
          //   // for (int idx : indices)
          //   //   quadr.push_back(points[idx]);
          //   if (quadr.is_simple() && quadr.is_convex() && (quadr.is_counterclockwise_oriented() || quadr.is_clockwise_oriented()))
          //     std::cout<<"It's simple and convex and cw or ccw\n";
          //   else
          //     std::cout<<quadr.is_simple()<<" "<<quadr.is_convex()<<'\n';
          // }
          //todo there might still be a bug here, it should always be ccw but it isn't
          if (quadr.is_simple() /*&& quadr.is_convex()*/ && (quadr.is_counterclockwise_oriented() || quadr.is_clockwise_oriented())) {
            //check that what I did above is right
            // std::cout<<"quadr: "<<i<<' ';
            // for (int idx : indices)
            //   std::cout<<idx<<' ';
            // std::cout<<'\n';
            // std::cout<<quadr.orientation()<<'\n';
            // assert(quadr.is_counterclockwise_oriented());
            // if (!quadr.is_clockwise_oriented()) {
            //   std::cout<<"Not cw: ";
            //   std::cout<<quadr.orientation()<<' ';
            //   std::cout<<i<<' ';
            //   for (int idx : indices)
            //     std::cout<<idx<<' ';
            //   std::cout<<'\n';
            // }
            bool empty=true;
            //todo you can speed up this part with some sorting (remember to store the initial indices)
            for (int e=0;e<n && empty;e++) {
              if (quadr.bounded_side(points[e])==CGAL::ON_BOUNDED_SIDE) {
                empty=false;
                // if (i==2 && j==10 && k==11 && l==16)
                //   std::cout<<e<<' '<<points[e].x()<<' '<<points[e].y()<< '\n';
              }
            }

            // if (i==2 && j==10 && k==11 && l==16)
            //   std::cout<<empty<<'\n';
            if (i==1 && j==3 && k==11 && l==18)
              std::cout<<empty<<'\n';
            if (empty) {
              // std::cout<<"quadr: "<<i<<' ';
              // for (int idx : indices)
              //   std::cout<<idx<<' ';
              // std::cout<<'\n';
              cntEmptyQuad++;

              //compute the diagonals (to linearize indices, the vertex indices must be in increasing order)
              int d00,d01,d10,d11;
              if (i<indices[1]) {
                d00=i;
                d01=indices[1];
              }
              else {
                d00=indices[1];
                d01=i;
              }
              if (indices[0]<indices[2]) {
                d10=indices[0];
                d11=indices[2];
              }
              else {
                d10=indices[2];
                d11=indices[0];
              }
              // if (i==2 && j==10 && k==11 && l==16)
              //   std::cout<<d00<<' '<<d01<<' '<<d10<<' '<<d11<<'\n';
              int d0=d00*n+d01;
              int d1=d10*n+d11;
              if (quadr.is_convex()) {
                qg[d0].push_back(d1);
                qg[d1].push_back(d0);
              }
              else {
                std::cout<<d00<<' '<<d01<<' '<<d10<<' '<<d11<<'\n';
                if (quadr.bounded_side(CGAL::midpoint(points[d00],points[d01]))==CGAL::ON_BOUNDED_SIDE)
                  qg[d0].push_back(-d1-1);
                if (quadr.bounded_side(CGAL::midpoint(points[d10],points[d11]))==CGAL::ON_BOUNDED_SIDE)
                  qg[d1].push_back(-d0-1);

              }
            }
            // else {
            //   std::cout<<"Point inside ";
            //   std::cout<<"quadr: "<<i<<' ';
            //   for (int idx : indices)
            //     std::cout<<idx<<' ';
            //   std::cout<<'\n';
            // }
          }
          // P quadrPoints[4]={points[i],points[j],points[k],points[l]};
          // P result[4];
          // P *ptr=CGAL::convex_hull_2(quadrPoints,quadrPoints+4,result);
          // //check that the four points are in convex position
          // if (ptr-result==4) {
          //   Polygon_2 quadr(result,result+4);
          //   bool empty=true;
          //   //todo you can speed up this part with some sorting
          //   for (int e=0;e<n && empty;e++) {
          //     if (quadr.bounded_side(points[e])==CGAL::ON_BOUNDED_SIDE)
          //       empty=false;
          //   }

            // if (empty) {
            //   int d00,d01,d10,d11;
            //
            // }

          // }
        }
  std::cout<<"There are "<<cntEmptyQuad<<" empty (not necessarily convex) quadrilaterals in the point set\n";
  return qg;
}

vi bfs(vvi& qg, int src) {
  int len=qg.size();
  vi dist(len,-1);
  dist[src]=0;
  std::queue<int> q;
  q.push(src);
  while (!q.empty()) {
    int u=q.front(); q.pop();
    // std::cout<<src<<' '<<u<<' '<<dist[u]<<'\n';
    for (int v : qg[u]) {
      if (dist[v]==-1) {
        dist[v]=dist[u]+1;
        q.push(v);
      }
    }
  }
  return dist;
}
int longestSP( flippingDS& t0,  flippingDS& t1, vvi& qg, int nPoints) {
  int maxMinDistance=0;
  vi2 edges0=t0.getEdges(), edges1=t1.getEdges();
  for (int k=0;k<edges0.size();k++) {
    i2 e0=edges0[k];
    if (t0.otherDiagonal(k)[0]!=-1) {
      // std::cout<<"Current e0: "<<e0[0]<<' '<<e0[1]<<'\n';
      int i=e0[0]*nPoints+e0[1];
      int minDist=qg.size();
      vi dist=bfs(qg,i);
      for (i2 e1 : edges1) {
        int j=e1[0]*nPoints+e1[1];

        if (dist[j]>=0) {
          // std::cout<<e1[0]<<' '<<e1[1]<<' '<<dist[j]<<'\n';
          minDist=std::min(minDist,dist[j]);
        }
      }


      maxMinDistance=std::max(maxMinDistance,minDist);
    }
  }
  return maxMinDistance;
}