#include<generateCNF.hpp>
#include<geometryUtils.hpp>
#include<sstream>
typedef std::vector<vi2> vvi2;

typedef std::array<int,4> i4;
typedef std::vector<i4> vi4;
typedef std::map<i4,int> mi4;

typedef std::vector<bool> vb;
typedef std::vector<vb> vvb;
typedef std::vector<vvb> vvvb;

int toIdx(int isPresentVar, int edgeIdx, int step, int triIdx, vi4& idxToVar, mi4& varToIdx) {
  // int isPresentInt=isPresentVar ? 1 : 0;
  i4 tmp={isPresentVar,edgeIdx,step,triIdx};
  if (varToIdx.count(tmp)==0) {
    varToIdx[tmp]=idxToVar.size();
    idxToVar.push_back(tmp);
  }
  return varToIdx[tmp]+1;
}

void printIdx(int idx, int n, const vi4& idxToVar) {
  std::cout<<"Idx: "<<idx<<' ';
  if (idx<0) {
    // std::cout<<"Negated ";
    idx*=-1;
  }
  idx--;
  bool isPresentVar=idxToVar[idx][0]==1;
  int u=idxToVar[idx][1]/n;
  int v=idxToVar[idx][1]%n;
  int step=idxToVar[idx][2];
  int triIdx=idxToVar[idx][3];
  std::cout<<"isPresent: "<<isPresentVar<<" u v: "<<u<<' '<<v<<" step: "<<step<<" "<<triIdx<<'\n';
}

void fromIdx(int idx, bool& isPresentVar, int& u, int& v, int& step, int& triIdx, int n, const vi4& idxToVar) {
  idx--;
  isPresentVar=idxToVar[idx][0]==1;
  u=idxToVar[idx][1]/n;
  v=idxToVar[idx][1]%n;
  step=idxToVar[idx][2];
  triIdx=idxToVar[idx][3];
  // std::cout<<"isPresent: "<<isPresentVar<<" u v: "<<u<<' '<<v<<" step: "<<step<<" "<<triIdx<<'\n';
}

int pointIndicesToEdgeIdx(int n, int a, int b) {
  if (a<b)
    return a*n+b;
  return b*n+a;
}

int toIdxPrint(int n, int distCenter, int m, bool isPresentVar, int edgeIdx, int step, int triIdx) {
  //we want to predicate that after the next step we have non intersecting edges, so step is in [0,distCenter]
  distCenter++;
  std::cout<<isPresentVar<<' '<<edgeIdx<<' '<<step<<' '<<triIdx<<'\n';

  int idx=triIdx*n*n*distCenter+step*n*n+edgeIdx;
  if (isPresentVar)
    idx+=m*n*n*distCenter;
  //+1 is there to ensure that all clauses are >=1
  return idx+1;
}
int toIdx(int n, int distCenter, int m, bool isPresentVar, int edgeIdx, int step, int triIdx) {
  //we want to predicate that after the next step we have non intersecting edges, so step is in [0,distCenter]
  distCenter++;

  int idx=triIdx*n*n*distCenter+step*n*n+edgeIdx;
  if (isPresentVar)
    idx+=m*n*n*distCenter;
  //+1 is there to ensure that all clauses are >=1
  return idx+1;
}

void fromIdxOutput(int idx, int n, int distCenter, int m) {
  bool isPresentVar;
  int u, v, step, triIdx;
  //we want to predicate that after the next step we have non intersecting edges, so step is in [0,distCenter]
  distCenter++;

  idx-=1;
  if (idx>=m*n*n*distCenter) {
    isPresentVar=true;
    idx-=m*n*n*distCenter;
  }
  else
    isPresentVar=false;

  triIdx=idx/(n*n*distCenter);
  idx-=triIdx*n*n*distCenter;
  step=idx/(n*n);
  idx-=step*n*n;
  u=idx/n;
  idx-=u*n;
  v=idx;
  std::cout<<"Present u v step triIdx: ";
  std::cout<<isPresentVar<<' '<<u<<' '<<v<<' '<<step<<' '<<triIdx<<'\n';
}
void fromIdx(int idx, int n, int distCenter, int m, bool& isPresentVar, int& u, int& v, int& step, int& triIdx) {
  //we want to predicate that after the next step we have non intersecting edges, so step is in [0,distCenter]
  distCenter++;

  idx-=1;
  if (idx>=m*n*n*distCenter) {
    isPresentVar=true;
    idx-=m*n*n*distCenter;
  }
  else
    isPresentVar=false;

  triIdx=idx/(n*n*distCenter);
  idx-=triIdx*n*n*distCenter;
  step=idx/(n*n);
  idx-=step*n*n;
  u=idx/n;
  idx-=u*n;
  v=idx;
}

// The formulation below is incomplete and wrong
// void genCNFMoreVar(const vP& points, const vvi& qg, const vvvb& emptyTriangle, vvi& crossingEdgesVec, std::vector<facesDS>& tVec, int n, vi& distCenter, const std::string& outFileName, const std::string& idxFileName) {
//   int m=tVec.size();
//   int nSquared=n*n;
//   vvi2 edgeL(tVec.size());
//   for (int i=0;i<tVec.size();i++)
//     edgeL[i]=tVec[i].getEdges();
//   vi4 idxToVar;
//   mi4 varToIdx;
//   vvi clauses;
//
//   int maxDistCenter=distCenter[0];
//
//   for (int d : distCenter)
//     maxDistCenter=std::max(d,maxDistCenter);
//   //this stores for each triangulation, for each edge (pair of points) the earliest step in which it could be present
//   vvi mightBePresentFromStep(m,vi(n*n,2*(maxDistCenter+1)));
//
//   int firstNonNeg=-1;
//   for (int t=0;t<tVec.size();t++)
//     if (distCenter[t]>=0 && firstNonNeg==-1)
//       firstNonNeg=t;
//
//   //step 0 clauses (here we number steps from 0, unlike in the pdf formulation, just so that it's then easier
//   //to convert variables to flipList)
//   for (int t=0;t<tVec.size();t++) {
//     if (distCenter[t]==0) {
//       for (int i=0;i<edgeL[t].size();i++) {
//         int a=edgeL[t][i][0],b=edgeL[t][i][1];
//         int edgeIdx=a*n+b;
//         clauses.push_back({toIdx(1,edgeIdx,distCenter[firstNonNeg],firstNonNeg,idxToVar,varToIdx)});
//         mightBePresentFromStep[t][edgeIdx]=0;
//       }
//     }
//     else if (distCenter[t]>0) {
//       std::vector<bool> couldBeFlipped(edgeL[t].size(),false);
//       std::vector<bool> isFlippable(n*n,false);
//
//       for (int i=0;i<edgeL[t].size();i++) {
//         int a=edgeL[t][i][0],b=edgeL[t][i][1];
//         int edgeIdx=a*n+b;
//         mightBePresentFromStep[t][edgeIdx]=0;
//         i2 otherDiag=tVec[t].otherDiagonal(i);
//
//         if (otherDiag[0]!=-1) {
//           std::array<P,4> quad={points[a],points[otherDiag[0]],points[b],points[otherDiag[1]]};
//           if (canBeFlipped(quad)){
//
//             isFlippable[edgeIdx]=true;
//             int otherDiagIdx=otherDiag[0]*n+otherDiag[1];
//             mightBePresentFromStep[t][otherDiagIdx]=std::min(mightBePresentFromStep[t][otherDiagIdx],1);
//             int flipEdge=toIdx(0,edgeIdx,0,t,idxToVar,varToIdx);
//             int presentNextStep;
//             if (distCenter[t]==1)
//               presentNextStep=toIdx(1,edgeIdx,distCenter[firstNonNeg],firstNonNeg,idxToVar,varToIdx);
//             else
//               presentNextStep=toIdx(1,edgeIdx,1,t,idxToVar,varToIdx);
//             int presentOtherDiagNextStep;
//             if (distCenter[t]==1)
//               presentOtherDiagNextStep=toIdx(1,otherDiagIdx,distCenter[firstNonNeg],firstNonNeg,idxToVar,varToIdx);
//             else
//               presentOtherDiagNextStep=toIdx(1,otherDiagIdx,1,t,idxToVar,varToIdx);
//
//             clauses.push_back({flipEdge,presentNextStep});
//             clauses.push_back({-flipEdge,presentOtherDiagNextStep});
//             clauses.push_back({-presentNextStep,-presentOtherDiagNextStep});
//
//             couldBeFlipped[i]=true;
//           }
//         }
//       }
//
//       for (int i=0;i<nSquared;i++) {
//
//         if (!isFlippable[i]) {
//
//           if (mightBePresentFromStep[t][i]==0) {
//             //todo try removing these clauses and removing the corresponding term in clauses at the next step where they are not negated
//             int presentEdgeNext;
//             if (distCenter[t]==1) //todo might be redundant...
//               presentEdgeNext=toIdx(1,i,distCenter[firstNonNeg],firstNonNeg,idxToVar,varToIdx);
//             else
//               presentEdgeNext=toIdx(1,i,1,t,idxToVar,varToIdx);
//             clauses.push_back({presentEdgeNext});
//           }
//
//         }
//       }
//
//       //independence clauses
//       for (int i=0;i<edgeL[t].size();i++)
//         for (int j=i+1;j<edgeL[t].size();j++)
//           if (couldBeFlipped[i] && couldBeFlipped[j] && tVec[t].shareTriangle(i,j)) {
//             int a0=edgeL[t][i][0],b0=edgeL[t][i][1];
//             int a1=edgeL[t][j][0],b1=edgeL[t][j][1];
//             int edgeIdx0=a0*n+b0;
//             int edgeIdx1=a1*n+b1;
//             int flipEdge0=toIdx(0,edgeIdx0,0,t,idxToVar,varToIdx);
//             int flipEdge1=toIdx(0,edgeIdx1,0,t,idxToVar,varToIdx);
//             clauses.push_back({-flipEdge0,-flipEdge1});
//           }
//     }
//   }
//
//
//   //fill mightBeAtStep for steps [2,distCenter] (this is essentially some variant of bfs in the QG)
//   for (int t=0;t<tVec.size();t++) {
//     for (int s=1;s<distCenter[t];s++) {
//       for (int e0=0;e0<nSquared;e0++) {
//         if (mightBePresentFromStep[t][e0]<=s) {
//           int a=e0/n;
//           int c=e0-a*n;
//           for(int e1 : qg[e0]) {
//             //if mightBePresentFromStep[t][e1]>s, it actually has to be the initial value (2*distCenter)
//             if (e1>=0 && mightBePresentFromStep[t][e1]>s) {
//               int b=e1/n;
//               int d=e1-b*n;
//               vi quadrVertices={a,b,c,d};
//               bool quadrAtStepS=true;
//               for (int j=0;j<quadrVertices.size() && quadrAtStepS;j++) {
//                 int edgeIdx=pointIndicesToEdgeIdx(n,quadrVertices[j],quadrVertices[(j+1)%quadrVertices.size()]);
//                 quadrAtStepS=(mightBePresentFromStep[t][edgeIdx]<=s);
//               }
//               if (quadrAtStepS) {
//                 mightBePresentFromStep[t][e1]=s+1;
//               }
//             }
//           }
//         }
//       }
//     }
//   }
//
//   //this array stored the following information:
//   //if an edge might be in the center, ie mightBePresentFromStep is <= distCenter for all triangulations, then it's 0 (which means that this edge might be in the center)
//   //the other values are filled with multi source bfs from the edges at value 0
//   vi minDistanceFromCenter(n*n,2*(maxDistCenter+1));
//
//   for (int e=0;e<nSquared;e++) {
//     bool mightBeInCenter=true;
//     for (int t=0;t<tVec.size() && mightBeInCenter;t++)
//       if (distCenter[t]>=0)
//         mightBeInCenter=mightBePresentFromStep[t][e]<=distCenter[t];
//     if (mightBeInCenter) {
//       minDistanceFromCenter[e]=0;
//     }
//   }
//
//   for (int s=0;s<maxDistCenter;s++) {
//     for (int e0=0;e0<nSquared;e0++) {
//       if (minDistanceFromCenter[e0]<=s) {
//         int a=e0/n;
//         int c=e0-a*n;
//         for (int e1 : qg[e0]) {
//           if (e1>=0 && minDistanceFromCenter[e1]>s) {
//             int b=e1/n;
//             int d=e1-b*n;
//             vi quadrVertices={a,b,c,d};
//             bool quadrAtStepS=true;
//             for (int j=0;j<quadrVertices.size() && quadrAtStepS;j++) {
//               int edgeIdx=pointIndicesToEdgeIdx(n,quadrVertices[j],quadrVertices[(j+1)%quadrVertices.size()]);
//               quadrAtStepS=(minDistanceFromCenter[edgeIdx]<=s);
//             }
//             if (quadrAtStepS) {
//               minDistanceFromCenter[e1]=s+1;
//             }
//           }
//         }
//       }
//     }
//   }
//
//   //steps 1 to distCenter-1
//   std::cout<<"Computed helper distance bounds arrays"<<std::endl;
//   vvvb isFlippable(m);
//
//
//   for (int i=0;i<nSquared;i++)
//     //if the edge is too far from any potential center edge, it can't be present
//     if (minDistanceFromCenter[i]>0) {
//       int presentEdge=toIdx(1,i,distCenter[firstNonNeg],firstNonNeg,idxToVar,varToIdx);
//       clauses.push_back({-presentEdge});
//     }
//
//   for (int t=0;t<tVec.size();t++) {
//     isFlippable[t]=vvb(distCenter[t]+1,vb(n*n,false));
//     for (int s=1;s<distCenter[t];s++) {
//
//       for (int e0=0;e0<nSquared;e0++) {
//         if (mightBePresentFromStep[t][e0]<=s) {
//           int a=e0/n;
//           int c=e0-a*n;
//
//           for(int e1 : qg[e0]) {
//             if (e1>=0 && minDistanceFromCenter[e1]<distCenter[t]-s) {
//               //the quadrilateral is {a,b,c,d} and the diagonal are ac and bd (a,b,c,d might not be in
//               //convex order, but it doesn't matter); since the graph is undirected, the symmetric clause
//               //for bd being the current diagonal will be inserted by the symmetric edge
//               int b=e1/n;
//               int d=e1-b*n;
//
//               vi quadrVertices={a,b,c,d};
//               bool quadrAtStepS=true;
//               for (int j=0;j<quadrVertices.size() && quadrAtStepS;j++) {
//                 int edgeIdx=pointIndicesToEdgeIdx(n,quadrVertices[j],quadrVertices[(j+1)%quadrVertices.size()]);
//                 quadrAtStepS=(mightBePresentFromStep[t][edgeIdx]<=s);
//               }
//
//               if (quadrAtStepS) {
//                 vi presentVar;
//                 for (int j=0;j<quadrVertices.size();j++) {
//                   int edgeIdx=pointIndicesToEdgeIdx(n,quadrVertices[j],quadrVertices[(j+1)%quadrVertices.size()]);
//                   presentVar.push_back(toIdx(1,edgeIdx,s,t,idxToVar,varToIdx));
//                 }
//
//                 vi clauseTmp;
//                 int presentAc=toIdx(1,e0,s,t,idxToVar,varToIdx);
//                 int flipAc=toIdx(0,e0,s,t,idxToVar,varToIdx);
//
//                 isFlippable[t][s][e0]=true;
//                 int flipToBd;
//                 if (s+1==distCenter[t])
//                   flipToBd=toIdx(2,e1,distCenter[firstNonNeg]-1,firstNonNeg,idxToVar,varToIdx);
//                 else
//                   flipToBd=toIdx(2,e1,s,t,idxToVar,varToIdx);
//                 clauseTmp={-presentAc,-flipAc,flipToBd};
//
//                 for (int var : presentVar)
//                   clauseTmp.push_back(-var);
//                 clauses.push_back(clauseTmp);
//
//               }
//
//             }
//           }
//           if (isFlippable[t][s][e0]) {
//             for(int e1 : qg[e0]){
//               if (e1<0 || minDistanceFromCenter[e1]>=distCenter[t]-s) {
//                 //if the quadrilateral is not convex, we add the following clause:
//                 //the quadrilateral is missing, or ac is not flipped at this step
//                 int actualIdx1=-(e1+1);
//                 if (e1>=0)
//                   actualIdx1=e1;
//                 int b=actualIdx1/n;
//                 int d=actualIdx1-b*n;
//                 vi quadrVertices={a,b,c,d};
//                 //if e0 cannot be present at this step, it will not be flipped because of the clause p_{e_0} or not f_{e_0}
//                 bool quadrAtStepS=mightBePresentFromStep[t][e0]<=s;
//                 for (int j=0;j<quadrVertices.size() && quadrAtStepS;j++) {
//                   int edgeIdx=pointIndicesToEdgeIdx(n,quadrVertices[j],quadrVertices[(j+1)%quadrVertices.size()]);
//                   quadrAtStepS=(mightBePresentFromStep[t][edgeIdx]<=s);
//                 }
//                 if (quadrAtStepS) {
//                   vi presentVar;
//                   for (int j=0;j<quadrVertices.size();j++) {
//                     int edgeIdx=pointIndicesToEdgeIdx(n,quadrVertices[j],quadrVertices[(j+1)%quadrVertices.size()]);
//                     presentVar.push_back(toIdx(1,edgeIdx,s,t,idxToVar,varToIdx));
//                   }
//
//                   int flipAc=toIdx(0,e0,s,t,idxToVar,varToIdx);
//                   vi clauseTmp={-flipAc};
//                   for (int var : presentVar)
//                     clauseTmp.push_back(-var);
//                   clauses.push_back(clauseTmp);
//                 }
//               }
//             }
//           }
//
//
//           //independence clauses
//           if (isFlippable[t][s][e0]) {
//             i2 e0Points={a,c};
//             for (int i=0;i<2;i++) {
//               int u=e0Points[i],v=e0Points[(i+1)%2];
//               for (int y=0;y<n;y++) {
//                 //todo you could speed this up by precomputing the empty triangles from a given pair of points, but this shouldn't be the bottleneck anyway
//                 if (emptyTriangle[u][v][y]) {
//                   int uyIdx=pointIndicesToEdgeIdx(n,u,y);
//                   int vyIdx=pointIndicesToEdgeIdx(n,v,y);
//                   if (e0>uyIdx && isFlippable[t][s][uyIdx]
//                     && mightBePresentFromStep[t][vyIdx]<=s) {
//                     int flipUv=toIdx(0,e0,s,t,idxToVar,varToIdx),
//                     flipUy=toIdx(0,uyIdx,s,t,idxToVar,varToIdx),
//                     presentVy=toIdx(1,vyIdx,s,t,idxToVar,varToIdx);
//                     clauses.push_back({-presentVy,-flipUv,-flipUy});
//                     }
//                 }
//               }
//             }
//           }
//         }
//       }
//
//
//       for (int i=0;i<nSquared;i++){
//         //if the edge is too far from any potential center edge, it can't be present
//         //(I do this for step s+1 because I need to impose this for steps [2,distCenter[t]))
//         if (s+1<distCenter[t] && minDistanceFromCenter[i]>=distCenter[t]-s) {
//           int presentEdgeNext=toIdx(1,i,s+1,t,idxToVar,varToIdx);
//           clauses.push_back({-presentEdgeNext});
//         }
//         //If the edge might be present at step s:
//         //if it's flippable, add 2 clauses: (present at step s and not flipped) => (present at step s+1); not present => not flipped
//         //if it's not, add 1 clause:  present at step s => present at step s+1
//         if (mightBePresentFromStep[t][i]<=s) {
//           int presentEdge=toIdx(1,i,s,t,idxToVar,varToIdx);
//           int presentEdgeNext;
//           if (s+1==distCenter[t])
//             presentEdgeNext=toIdx(1,i,distCenter[firstNonNeg],firstNonNeg,idxToVar,varToIdx);
//           else
//             presentEdgeNext=toIdx(1,i,s+1,t,idxToVar,varToIdx);
//           int flipToEdge=toIdx(2,i,s,t,idxToVar,varToIdx);
//
//           clauses.push_back({-flipToEdge,presentEdgeNext});
//           clauses.push_back({flipToEdge,presentEdge,-presentEdgeNext});
//
//           if (isFlippable[t][s][i]) {
//             int flipEdge=toIdx(0,i,s,t,idxToVar,varToIdx);
//             clauses.push_back({-presentEdge,flipEdge,presentEdgeNext});
//             clauses.push_back({presentEdge,-flipEdge});
//             clauses.push_back({flipToEdge,-presentEdge,-flipEdge,-presentEdgeNext});
//           }
//           else
//             clauses.push_back({-presentEdge,presentEdgeNext});
//
//         }
//
//       }
//
//     }
//   }
//
//   std::cout<<"In total, there are "<<clauses.size()<<" clauses"<<std::endl;
//
//   std::ofstream output(outFileName);
//   output<<"p cnf "<<idxToVar.size()<<' '<<clauses.size()<<'\n';
//   for (auto & clause : clauses) {
//     for (int literal : clause)
//       output<<literal<<' ';
//     output<<"0\n";
//   }
//   output.close();
//
//   //output the idxToVar vector, using append mode so that you can first store the uid in the main
//   std::ofstream outIdx{idxFileName,std::ios::app};
//   outIdx<<n<<' '<<m<<' '<<idxToVar.size()<<'\n';
//   for (int d : distCenter) {
//     outIdx<<d<<' ';
//   }
//   outIdx<<'\n';
//   for (i4 var : idxToVar) {
//     for (int i=0;i<4;i++)
//       outIdx<<var[i]<<' ';
//     outIdx<<'\n';
//   }
//   outIdx.close();
// }

//Computes the array that stores for each triangulation, for each edge (pair of points) the earliest step in which it could be present
vvi computeMightBePresent(const vP& points, const vvi& qg, std::vector<facesDS>& tVec, int n, vi& distCenter, vvi2& edgeL) {
  int m=edgeL.size();
  int nSquared=n*n;
  int maxDistCenter=distCenter[0];
  for (int d : distCenter)
    maxDistCenter=std::max(d,maxDistCenter);

  vvi mightBePresentFromStep(m,vi(n*n,2*(maxDistCenter+1)));
  for (int t=0;t<m;t++)
    for (int i=0;i<edgeL[t].size();i++) {
      int a=edgeL[t][i][0],b=edgeL[t][i][1];
      int edgeIdx=a*n+b;
      mightBePresentFromStep[t][edgeIdx]=0;
      if (distCenter[t]>0) {
        i2 otherDiag=tVec[t].otherDiagonal(i);
        if (otherDiag[0]!=-1) {
          std::array<P,4> quad={points[a],points[otherDiag[0]],points[b],points[otherDiag[1]]};
          if (canBeFlipped(quad)) {
            int otherDiagIdx=otherDiag[0]*n+otherDiag[1];
            mightBePresentFromStep[t][otherDiagIdx]=std::min(mightBePresentFromStep[t][otherDiagIdx],1);
          }
        }
      }
    }

  //fill mightBeAtStep for steps [2,distCenter] (this is essentially some variant of bfs in the QG)
  for (int t=0;t<m;t++) {

    for (int s=1;s<distCenter[t];s++) {
      for (int e0=0;e0<nSquared;e0++) {
        if (mightBePresentFromStep[t][e0]<=s) {
          int a=e0/n;
          int c=e0-a*n;
          for(int e1 : qg[e0]) {
            //if mightBePresentFromStep[t][e1]>s, it actually has to be the initial value (2*distCenter)
            if (e1>=0 && mightBePresentFromStep[t][e1]>s) {
              int b=e1/n;
              int d=e1-b*n;
              vi quadrVertices={a,b,c,d};
              bool quadrAtStepS=true;
              for (int j=0;j<quadrVertices.size() && quadrAtStepS;j++) {
                int edgeIdx=pointIndicesToEdgeIdx(n,quadrVertices[j],quadrVertices[(j+1)%quadrVertices.size()]);
                quadrAtStepS=(mightBePresentFromStep[t][edgeIdx]<=s);
              }
              if (quadrAtStepS) {
                mightBePresentFromStep[t][e1]=s+1;
              }
            }
          }
        }
      }
    }
  }
  return mightBePresentFromStep;
}

//computes the array that stores the following information:
//if an edge might be in the center, ie mightBePresentFromStep is <= distCenter for all triangulations, then it's 0 (which means that this edge might be in the center)
//the other values are filled with graph exploration from the edges at value 0 (considering the quadrilateral graph)
vi computeMinDistanceFromCenter(int n, const vi& distCenter, const vvi& mightBePresentFromStep, const vvi& qg) {
  int m=distCenter.size();
  int maxDistCenter=distCenter[0];
  for (int d : distCenter)
    maxDistCenter=std::max(d,maxDistCenter);
  int nSquared=n*n;
  vi minDistanceFromCenter(nSquared,2*(maxDistCenter+1));

  for (int e=0;e<nSquared;e++) {
    bool mightBeInCenter=true;
    for (int t=0;t<m && mightBeInCenter;t++)
      if (distCenter[t]>=0)
        mightBeInCenter=mightBePresentFromStep[t][e]<=distCenter[t];
    if (mightBeInCenter) {
      // std::cout<<e<<" might be in center"<<'\n';
      minDistanceFromCenter[e]=0;
    }
  }

  for (int s=0;s<maxDistCenter;s++) {
    for (int e0=0;e0<nSquared;e0++) {
      if (minDistanceFromCenter[e0]<=s) {
        int a=e0/n;
        int c=e0-a*n;
        for (int e1 : qg[e0]) {
          if (e1>=0 && minDistanceFromCenter[e1]>s) {
            int b=e1/n;
            int d=e1-b*n;
            vi quadrVertices={a,b,c,d};
            bool quadrAtStepS=true;
            for (int j=0;j<quadrVertices.size() && quadrAtStepS;j++) {
              int edgeIdx=pointIndicesToEdgeIdx(n,quadrVertices[j],quadrVertices[(j+1)%quadrVertices.size()]);
              quadrAtStepS=(minDistanceFromCenter[edgeIdx]<=s);
            }
            if (quadrAtStepS) {
              minDistanceFromCenter[e1]=s+1;
            }
          }
        }
      }
    }
  }
  return minDistanceFromCenter;
}
void genCNF(const vP& points, const vvi& qg, const vvvb& emptyTriangle, vvi& crossingEdgesVec, std::vector<facesDS>& tVec, int n, vi& distCenter, const std::string& outFileName, const std::string& idxFileName) {
  // genCNFMoreVar(points,qg,emptyTriangle,crossingEdgesVec,tVec,n,distCenter,outFileName,idxFileName);
  // return;
  int m=tVec.size();
  int nSquared=n*n;
  vvi2 edgeL(tVec.size());
  for (int i=0;i<tVec.size();i++)
    edgeL[i]=tVec[i].getEdges();
  vi4 idxToVar;
  mi4 varToIdx;
  vvi clauses;

  int maxDistCenter=distCenter[0];
  // for (int d : distCenter)
  //   std::cout<<d<<std::endl;
  for (int d : distCenter)
    maxDistCenter=std::max(d,maxDistCenter);
  //this stores for each triangulation, for each edge (pair of points) the earliest step in which it could be present
  // vvi mightBePresentFromStep(m,vi(n*n,2*(maxDistCenter+1)));
  vvi mightBePresentFromStep=computeMightBePresent(points,qg,tVec,n,distCenter,edgeL);
  int firstNonNeg=-1;
  for (int t=0;t<tVec.size();t++)
    if (distCenter[t]>=0 && firstNonNeg==-1)
      firstNonNeg=t;
  // std::cout<<firstNonNeg<<std::endl;
  //step 0 clauses (here we number steps from 0, unlike in the pdf formulation, just so that it's then easier
  //to convert variables to flipList)
  for (int t=0;t<tVec.size();t++) {
    if (distCenter[t]==0) {
      for (int i=0;i<edgeL[t].size();i++) {
        int a=edgeL[t][i][0],b=edgeL[t][i][1];
        int edgeIdx=a*n+b;
        clauses.push_back({toIdx(1,edgeIdx,distCenter[firstNonNeg],firstNonNeg,idxToVar,varToIdx)});
        // mightBePresentFromStep[t][edgeIdx]=0;
      }
      // for (int i=0;i<nSquared;i++)
      //   if (mightBePresentFromStep[t][i]>0)
      //     clauses.push_back({-toIdx(true,i,0,t,idxToVar,varToIdx)});
    }
    else if (distCenter[t]>0) {
      std::vector<bool> couldBeFlipped(edgeL[t].size(),false);
      std::vector<bool> isFlippable(n*n,false);

      for (int i=0;i<edgeL[t].size();i++) {
        int a=edgeL[t][i][0],b=edgeL[t][i][1];
        int edgeIdx=a*n+b;
        // mightBePresentFromStep[t][edgeIdx]=0;
        i2 otherDiag=tVec[t].otherDiagonal(i);

        if (otherDiag[0]!=-1) {
          std::array<P,4> quad={points[a],points[otherDiag[0]],points[b],points[otherDiag[1]]};
          if (canBeFlipped(quad)){

            isFlippable[edgeIdx]=true;
            int otherDiagIdx=otherDiag[0]*n+otherDiag[1];
            // mightBePresentFromStep[t][otherDiagIdx]=std::min(mightBePresentFromStep[t][otherDiagIdx],1);
            int flipEdge=toIdx(0,edgeIdx,0,t,idxToVar,varToIdx);
            int presentNextStep;
            if (distCenter[t]==1)
              presentNextStep=toIdx(1,edgeIdx,distCenter[firstNonNeg],firstNonNeg,idxToVar,varToIdx);
            else
              presentNextStep=toIdx(1,edgeIdx,1,t,idxToVar,varToIdx);
            int presentOtherDiagNextStep;
            if (distCenter[t]==1)
              presentOtherDiagNextStep=toIdx(1,otherDiagIdx,distCenter[firstNonNeg],firstNonNeg,idxToVar,varToIdx);
            else
              presentOtherDiagNextStep=toIdx(1,otherDiagIdx,1,t,idxToVar,varToIdx);
            // std::cout<<"Clause :"<<clauses.size()<<'\n';
            // printIdx(flipEdge,n,idxToVar);
            // printIdx(presentNextStep,n,idxToVar);
            clauses.push_back({flipEdge,presentNextStep});
            clauses.push_back({-flipEdge,presentOtherDiagNextStep});
            clauses.push_back({-presentNextStep,-presentOtherDiagNextStep});

            couldBeFlipped[i]=true;
          }
        }
      }

      for (int i=0;i<nSquared;i++) {
        //todo triple check that these are useless. If you add them, put them later, and check that mightBePresentFromStep[t][i]<=distCenter[t]
        // if (mightBePresentFromStep[t][i]>1) {
        //   int presentEdgeNext=toIdx(1,i,1,t,idxToVar,varToIdx);
        //   clauses.push_back({-presentEdgeNext});
        // }
        if (!isFlippable[i]) {
          // int flipEdge=toIdx(0,i,0,t,idxToVar,varToIdx);
          // clauses.push_back({-flipEdge});
          if (mightBePresentFromStep[t][i]==0) {
            //todo try removing these clauses and removing the corresponding term in clauses at the next step where they are not negated
            int presentEdgeNext;
            if (distCenter[t]==1) //todo might be redundant...
              presentEdgeNext=toIdx(1,i,distCenter[firstNonNeg],firstNonNeg,idxToVar,varToIdx);
            else
              presentEdgeNext=toIdx(1,i,1,t,idxToVar,varToIdx);
            clauses.push_back({presentEdgeNext});
          }

        }
      }

      //independence clauses
      for (int i=0;i<edgeL[t].size();i++)
        for (int j=i+1;j<edgeL[t].size();j++)
          if (couldBeFlipped[i] && couldBeFlipped[j] && tVec[t].shareTriangle(i,j)) {
            int a0=edgeL[t][i][0],b0=edgeL[t][i][1];
            int a1=edgeL[t][j][0],b1=edgeL[t][j][1];
            int edgeIdx0=a0*n+b0;
            int edgeIdx1=a1*n+b1;
            int flipEdge0=toIdx(0,edgeIdx0,0,t,idxToVar,varToIdx);
            int flipEdge1=toIdx(0,edgeIdx1,0,t,idxToVar,varToIdx);
            clauses.push_back({-flipEdge0,-flipEdge1});
          }
    }
  }

  //this array stores the following information:
  //if an edge might be in the center, ie mightBePresentFromStep is <= distCenter for all triangulations, then it's 0 (which means that this edge might be in the center)
  //the other values are filled with multi source bfs from the edges at value 0
  // vi minDistanceFromCenter(n*n,2*(maxDistCenter+1));

  vi minDistanceFromCenter=computeMinDistanceFromCenter(n,distCenter,mightBePresentFromStep,qg);
  // for (int e=0;e<nSquared;e++)
  //   std::cout<<minDistanceFromCenter[e]<<' ';
  // std::cout<<'\n';
  //steps 1 to distCenter-1
  std::cout<<"Computed helper distance bounds arrays"<<std::endl;
  vvvb isFlippable(m);


  for (int i=0;i<nSquared;i++)
    //if the edge is too far from any potential center edge, it can't be present
    if (minDistanceFromCenter[i]>0) {
      int presentEdge=toIdx(1,i,distCenter[firstNonNeg],firstNonNeg,idxToVar,varToIdx);
      clauses.push_back({-presentEdge});
    }

  for (int t=0;t<tVec.size();t++) {
    isFlippable[t]=vvb(distCenter[t]+1,vb(n*n,false));
    for (int s=1;s<distCenter[t];s++) {
      // std::cout<<clauses.size()<<std::endl;
      // for (int e0=0;e0<nSquared;e0++)
      //   if (mightBePresentFromStep[t][e0]>s) {
      //     clauses.push_back({-toIdx(1,e0,s,t,idxToVar,varToIdx)});
      //   }


      for (int e0=0;e0<nSquared;e0++) {
        if (mightBePresentFromStep[t][e0]<=s) {
          int a=e0/n;
          int c=e0-a*n;

          for(int e1 : qg[e0]) {
            if (e1>=0 && minDistanceFromCenter[e1]<distCenter[t]-s) {
              //the quadrilateral is {a,b,c,d} and the diagonal are ac and bd (a,b,c,d might not be in
              //convex order, but it doesn't matter); since the graph is undirected, the symmetric clause
              //for bd being the current diagonal will be inserted by the symmetric edge
              int b=e1/n;
              int d=e1-b*n;

              vi quadrVertices={a,b,c,d};
              bool quadrAtStepS=true;
              for (int j=0;j<quadrVertices.size() && quadrAtStepS;j++) {
                int edgeIdx=pointIndicesToEdgeIdx(n,quadrVertices[j],quadrVertices[(j+1)%quadrVertices.size()]);
                quadrAtStepS=(mightBePresentFromStep[t][edgeIdx]<=s);
              }

              if (quadrAtStepS) {
                vi presentVar;
                // std::cout<<"Triangulation "<<t<<", step: "<<s<<", edge: "<<a<<" "<<c<<' '<<b<<' '<<d<<'\n';
                for (int j=0;j<quadrVertices.size();j++) {
                  int edgeIdx=pointIndicesToEdgeIdx(n,quadrVertices[j],quadrVertices[(j+1)%quadrVertices.size()]);
                  presentVar.push_back(toIdx(1,edgeIdx,s,t,idxToVar,varToIdx));
                }

                vi clauseTmp;
                int presentAc=toIdx(1,e0,s,t,idxToVar,varToIdx);
                int flipAc=toIdx(0,e0,s,t,idxToVar,varToIdx);

                isFlippable[t][s][e0]=true;
                int presentBdNext;
                if (s+1==distCenter[t])
                  presentBdNext=toIdx(1,e1,distCenter[firstNonNeg],firstNonNeg,idxToVar,varToIdx);
                else
                  presentBdNext=toIdx(1,e1,s+1,t,idxToVar,varToIdx);
                clauseTmp={-presentAc,-flipAc,presentBdNext};

                for (int var : presentVar)
                  clauseTmp.push_back(-var);
                clauses.push_back(clauseTmp);

              }

            }
          }
          if (isFlippable[t][s][e0]) {
            for(int e1 : qg[e0]){
              if (e1<0 || minDistanceFromCenter[e1]>=distCenter[t]-s) {
                //if the quadrilateral is not convex, we add the following clause:
                //the quadrilateral is missing, or ac is not flipped at this step
                int actualIdx1=-(e1+1);
                if (e1>=0)
                  actualIdx1=e1;
                int b=actualIdx1/n;
                int d=actualIdx1-b*n;
                vi quadrVertices={a,b,c,d};
                //if e0 cannot be present at this step, it will not be flipped because of the clause p_{e_0} or not f_{e_0}
                bool quadrAtStepS=mightBePresentFromStep[t][e0]<=s;
                for (int j=0;j<quadrVertices.size() && quadrAtStepS;j++) {
                  int edgeIdx=pointIndicesToEdgeIdx(n,quadrVertices[j],quadrVertices[(j+1)%quadrVertices.size()]);
                  quadrAtStepS=(mightBePresentFromStep[t][edgeIdx]<=s);
                }
                if (quadrAtStepS) {
                  vi presentVar;
                  // std::cout<<"Triangulation "<<t<<", step: "<<s<<", edge: "<<a<<" "<<c<<' '<<b<<' '<<d<<'\n';
                  for (int j=0;j<quadrVertices.size();j++) {
                    int edgeIdx=pointIndicesToEdgeIdx(n,quadrVertices[j],quadrVertices[(j+1)%quadrVertices.size()]);
                    presentVar.push_back(toIdx(1,edgeIdx,s,t,idxToVar,varToIdx));
                  }

                  int flipAc=toIdx(0,e0,s,t,idxToVar,varToIdx);
                  vi clauseTmp={-flipAc};
                  for (int var : presentVar)
                    clauseTmp.push_back(-var);
                  clauses.push_back(clauseTmp);
                }
              }
            }
          }


          //independence clauses
          //TODO check for duplicates (e0>uyIdx should remove them, but double check that there is no duplicate left)
          if (isFlippable[t][s][e0]) {
            i2 e0Points={a,c};
            for (int i=0;i<2;i++) {
              int u=e0Points[i],v=e0Points[(i+1)%2];
              for (int y=0;y<n;y++) {
                //todo you could speed this up by precomputing the empty triangles from a given pair of points, but this shouldn't be the bottleneck anyway
                if (emptyTriangle[u][v][y]) {
                  int uyIdx=pointIndicesToEdgeIdx(n,u,y);
                  int vyIdx=pointIndicesToEdgeIdx(n,v,y);
                  //!qg[e0].empty() && !qg[uyIdx].empty() && mightBePresentFromStep[t][e0]<=s && mightBePresentFromStep[t][uyIdx]<=s
                  if (e0>uyIdx && isFlippable[t][s][uyIdx]
                    && mightBePresentFromStep[t][vyIdx]<=s) {
                    int flipUv=toIdx(0,e0,s,t,idxToVar,varToIdx),
                    flipUy=toIdx(0,uyIdx,s,t,idxToVar,varToIdx),
                    presentVy=toIdx(1,vyIdx,s,t,idxToVar,varToIdx);
                    clauses.push_back({-presentVy,-flipUv,-flipUy});
                    }
                }
              }
            }
          }
        }
      }


      for (int i=0;i<nSquared;i++){
        //if the edge is too far from any potential center edge, it can't be present
        //(I do this for step s+1 because I need to impose this for steps [2,distCenter[t]))
        if (s+1<distCenter[t] && minDistanceFromCenter[i]>=distCenter[t]-s) {
          int presentEdgeNext=toIdx(1,i,s+1,t,idxToVar,varToIdx);
          clauses.push_back({-presentEdgeNext});
        }
        //If the edge might be present at step s:
        //if it's flippable, add 2 clauses: (present at step s and not flipped) => (present at step s+1); not present => not flipped
        //if it's not, add 1 clause:  present at step s => present at step s+1
        if (mightBePresentFromStep[t][i]<=s) {
          int presentEdge=toIdx(1,i,s,t,idxToVar,varToIdx);
          int presentEdgeNext;//=toIdx(1,i,s+1,t,idxToVar,varToIdx);
          if (s+1==distCenter[t])
            presentEdgeNext=toIdx(1,i,distCenter[firstNonNeg],firstNonNeg,idxToVar,varToIdx);
          else
            presentEdgeNext=toIdx(1,i,s+1,t,idxToVar,varToIdx);

          if (isFlippable[t][s][i]) {
            int flipEdge=toIdx(0,i,s,t,idxToVar,varToIdx);
            clauses.push_back({-presentEdge,flipEdge,presentEdgeNext});
            clauses.push_back({presentEdge,-flipEdge});
          }
          else
            clauses.push_back({-presentEdge,presentEdgeNext});

        }

      }

    }
  }


  //no crossing edges are allowed within each triangulation at any step
  //for the last step, it is not necessary to impose this
  std::cout<<"Computing no crossing clauses, current number of clauses is "<<clauses.size()<<std::endl;
  for (int t=0;t<tVec.size();t++) {
    int maxS=distCenter[t]-1;
    if (t==firstNonNeg)
      maxS++;

    for (int e0=0;e0<nSquared;e0++) {
      if (mightBePresentFromStep[t][e0]<=maxS) {
        for (int e1 : crossingEdgesVec[e0]) {
          int minS=std::max(mightBePresentFromStep[t][e0],mightBePresentFromStep[t][e1]);
          int maxSConsideringDistToCenter=std::min(maxS,distCenter[t]-std::max(minDistanceFromCenter[e0],minDistanceFromCenter[e1]));
          for (int s=minS;s<=maxSConsideringDistToCenter;s++){
            int presentE0=toIdx(1,e0,s,t,idxToVar,varToIdx);
            int presentE1=toIdx(1,e1,s,t,idxToVar,varToIdx);
            clauses.push_back({-presentE0,-presentE1});
          }
        }
      }
    }

  }
  // std::cout<<"Computing last step clauses, current number of clauses is "<<clauses.size()<<std::endl;

  // for (int t1=firstNonNeg+1;t1<tVec.size();t1++) {
  //   int t0=firstNonNeg;
  //   if (distCenter[t0]>=0 && distCenter[t1]>=0)
  //     for (int e=0;e<nSquared;e++) {
  //       if (minDistanceFromCenter[e]==0) {
  //         int presentET0=toIdx(1,e,distCenter[t0],t0,idxToVar,varToIdx),
  //         presentET1=toIdx(1,e,distCenter[t1],t1,idxToVar,varToIdx);
  //         clauses.push_back({-presentET0,presentET1});
  //         clauses.push_back({presentET0,-presentET1});
  //       }
  //       else {
  //         int presentET0=toIdx(1,e,distCenter[t0],t0,idxToVar,varToIdx),
  //         presentET1=toIdx(1,e,distCenter[t1],t1,idxToVar,varToIdx);
  //         clauses.push_back({-presentET0});
  //         clauses.push_back({-presentET1});
  //       }
  //     }
  // }
  // for (int t0=0;t0<tVec.size();t0++) {
  //   for (int t1=t0+1;t1<tVec.size();t1++) {
  //     if (distCenter[t0]>=0 && distCenter[t1]>=0)
  //       for (int e0=0;e0<nSquared;e0++) {
  //         for (int e1 : qg[e0]) {
  //           if (e1>=0){// && mightBePresentFromStep[t0][e0]<=distCenter[t0] && mightBePresentFromStep[t1][e1]<=distCenter[t1]) {
  //             int presentAcT0=toIdx(1,e0,distCenter[t0],t0,idxToVar,varToIdx),
  //             presentBdT1=toIdx(1,e1,distCenter[t1],t1,idxToVar,varToIdx);
  //             clauses.push_back({-presentAcT0,-presentBdT1});
  //           }
  //         }
  //       }
  //   }
  // }
  std::cout<<"In total, there are "<<clauses.size()<<" clauses"<<std::endl;

  std::ofstream output(outFileName);
  output<<"p cnf "<<idxToVar.size()<<' '<<clauses.size()<<'\n';
  int cnt=0;
  for (auto & clause : clauses) {
    // if (distCenter[0]==2 && distCenter[1]==-1 && distCenter[2]==1) {
    //   std::cout<<"Clause :"<<cnt<<'\n';
    //   for (int literal : clause)
    //     printIdx(literal,n,idxToVar);
    //   cnt++;
    // }
    for (int literal : clause)
      output<<literal<<' ';
    output<<"0\n";
  }
  output.close();

  //output the idxToVar vector, using append mode so that you can first store the uid in the main
  std::ofstream outIdx{idxFileName,std::ios::app};
  outIdx<<n<<' '<<m<<' '<<idxToVar.size()<<'\n';
  for (int d : distCenter) {
    outIdx<<d<<' ';
  }
  outIdx<<'\n';
  for (i4 var : idxToVar) {
    for (int i=0;i<4;i++)
      outIdx<<var[i]<<' ';
    outIdx<<'\n';
  }
  outIdx.close();
}

void genCNF(const vP& points, const vvi& qg, vvi& crossingEdgesVec, std::vector<facesDS>& tVec, int n, int distCenter, const std::string& outFileName) {
  int m=tVec.size();
  int numberVar=2*m*n*n*(distCenter+1)+1;
  int nSquared=n*n;

  // std::cout<<toIdx(n,distCenter,m,true,16*n+33,1,0)<<'\n';
  vvi clauses;
  vvi2 edgeL(tVec.size());
  // std::vector<bool> testAssignment(numberVar,false);
  for (int i=0;i<tVec.size();i++) {
    edgeL[i]=tVec[i].getEdges();
  }

  // std::cout<<qg[397][0]<<'\n';
  // assert(qg[397][0]==2798);

  //todo consider removing isPresent and mightBeThereAtNextStep (this is a generalization of those)
  //this stores for each triangulation, for each edge (pair of points) the earliest step in which it could be present
  vvi mightBePresentFromStep(m,vi(n*n,2*distCenter));

  //step 0 clauses (here we number steps from 0, unlike in the pdf formulation, just so that it's then easier
  //to convert variables to flipList)
  for (int t=0;t<tVec.size();t++) {
    std::vector<bool> couldBeFlipped(edgeL[t].size(),false);
    std::vector<bool> isFlippable(n*n,false);
    std::vector<bool> isPresent(n*n,false);
    std::vector<bool> mightBeThereAtNextStep(n*n,false);

    for (int i=0;i<edgeL[t].size();i++) {
      int a=edgeL[t][i][0],b=edgeL[t][i][1];
      int edgeIdx=a*n+b;
      isPresent[edgeIdx]=true;
      mightBePresentFromStep[t][edgeIdx]=0;
      // for (int s=0;s<=distCenter;s++)
      //   mightBeAtStep[t][s][edgeIdx]=true;
      mightBeThereAtNextStep[edgeIdx]=true;
      i2 otherDiag=tVec[t].otherDiagonal(i);
      // testAssignment[toIdx(n,distCenter,m,true,edgeIdx,1,t)]=true;
      // testAssignment[toIdx(n,distCenter,m,true,edgeIdx,2,t)]=true;

      // std::cout<<"Triangulation "<<t<<" ";
      // std::cout<<a<<' '<<b<<" other diag: "<<otherDiag[0]<<' '<<otherDiag[1]<<'\n';
      if (otherDiag[0]!=-1) {
        std::array<P,4> quad={points[a],points[otherDiag[0]],points[b],points[otherDiag[1]]};
        // std::cout<<a<<' '<<b<<" other diag: "<<otherDiag[0]<<' '<<otherDiag[1]<<'\n';
        // std::cout<<"AAA\n";
        if (canBeFlipped(quad)){

          isFlippable[edgeIdx]=true;
          int otherDiagIdx=otherDiag[0]*n+otherDiag[1];
          mightBeThereAtNextStep[otherDiagIdx]=true;
          mightBePresentFromStep[t][otherDiagIdx]=std::min(mightBePresentFromStep[t][otherDiagIdx],1);
          // for (int s=1;s<=distCenter;s++)
          //   mightBeAtStep[t][s][otherDiagIdx]=true;
          int flipEdge=toIdx(n,distCenter,m,false,edgeIdx,0,t);
          int presentNextStep=toIdx(n,distCenter,m,true,edgeIdx,1,t);
          int presentOtherDiagNextStep=toIdx(n,distCenter,m,true,otherDiagIdx,1,t);
          // std::cout<<a<<' '<<b<<" other diag: "<<otherDiag[0]<<' '<<otherDiag[1]<<'\n';
          clauses.push_back({flipEdge,presentNextStep});
          clauses.push_back({-flipEdge,presentOtherDiagNextStep});
          clauses.push_back({-presentNextStep,-presentOtherDiagNextStep});

          couldBeFlipped[i]=true;
        }
      }

    }

    for (int i=0;i<nSquared;i++) {
      if (!mightBeThereAtNextStep[i]) {
        int presentEdgeNext=toIdx(n,distCenter,m,true,i,1,t);
        clauses.push_back({-presentEdgeNext});
      }
      if (!isFlippable[i]) {
        int flipEdge=toIdx(n,distCenter,m,false,i,0,t);
        // std::cout<<"Triangulation "<<t<<" edge "<<i/n<<' '<<i%n<<'\n';
        // std::cout<<"variable: "<<flipEdge<<"\n";
        clauses.push_back({-flipEdge});
        if (isPresent[i]) {
          int presentEdgeNext=toIdx(n,distCenter,m,true,i,1,t);
          clauses.push_back({presentEdgeNext});
        }

      }
    }

    //independence clauses
    for (int i=0;i<edgeL[t].size();i++)
      for (int j=i+1;j<edgeL[t].size();j++)
        if (couldBeFlipped[i] && couldBeFlipped[j] && tVec[t].shareTriangle(i,j)) {
          int a0=edgeL[t][i][0],b0=edgeL[t][i][1];
          int a1=edgeL[t][j][0],b1=edgeL[t][j][1];
          int edgeIdx0=a0*n+b0;
          int edgeIdx1=a1*n+b1;
          int flipEdge0=toIdx(n,distCenter,m,false,edgeIdx0,0,t);
          int flipEdge1=toIdx(n,distCenter,m,false,edgeIdx1,0,t);
          clauses.push_back({-flipEdge0,-flipEdge1});
        }

  }

  //fill mightBeAtStep for steps [2,distCenter] (this is essentially some variant of bfs in the QG)
  for (int s=1;s<distCenter;s++) {
    for (int t=0;t<tVec.size();t++) {
      // std::cout<<clauses.size()<<std::endl;
      std::vector<bool> isFlippable(n*n,false);
      for (int e0=0;e0<nSquared;e0++) {
        if (mightBePresentFromStep[t][e0]<=s) {
          int a=e0/n;
          int c=e0-a*n;
          for(int e1 : qg[e0]) {
            //if mightBePresentFromStep[t][e1]>s, it actually has to be the initial value (2*distCenter)
            if (e1>=0 && mightBePresentFromStep[t][e1]>s) {
              int b=e1/n;
              int d=e1-b*n;
              vi quadrVertices={a,b,c,d};
              bool quadrAtStepS=true;
              for (int j=0;j<quadrVertices.size() && quadrAtStepS;j++) {
                int edgeIdx=pointIndicesToEdgeIdx(n,quadrVertices[j],quadrVertices[(j+1)%quadrVertices.size()]);
                quadrAtStepS=(mightBePresentFromStep[t][edgeIdx]<=s);
              }
              if (quadrAtStepS) {
                mightBePresentFromStep[t][e1]=s+1;
                // for (int s1=s+1;s1<=distCenter;s1++)
                //   mightBeAtStep[t][s1][e1]=true;
              }
            }
          }
        }
      }
    }
  }

  //steps 1 to distCenter-1

  // std::cout<<mightBePresentFromStep[0][6*n+45]<<'\n';
  // for (int e1 : crossingEdgesVec[6*n+45])
  //   std::cout<<e1/n<<' '<<e1%n<<'\n';
  for (int s=1;s<distCenter;s++) {
    for (int t=0;t<tVec.size();t++) {
      // std::cout<<clauses.size()<<std::endl;
      for (int e0=0;e0<nSquared;e0++)
        if (mightBePresentFromStep[t][e0]>s) {
          clauses.push_back({-toIdx(n,distCenter,m,true,e0,s,t)});
        }

      std::vector<bool> isFlippable(n*n,false);
      for (int e0=0;e0<nSquared;e0++) {
        if (mightBePresentFromStep[t][e0]<=s) {
          int a=e0/n;
          int c=e0-a*n;
          int presentAc=toIdx(n,distCenter,m,true,e0,s,t),
              presentAcNext=toIdx(n,distCenter,m,true,e0,s+1,t),
              flipAc=toIdx(n,distCenter,m,false,e0,s,t);

          //ac is present at step s+1 only if it was present at step s, or one of the edges that intersect it was flipped at step s
          //TODO this is not enough, you need the right quadrilateral when you flip some other edge. Current solution I'm thinking of:
          //for all pairs of crossing edges (not only those that form empty quadrilaterals), they can't both be present at step s
          std::array<vi,5> causalityClauses;
          for (int i=0;i<causalityClauses.size();i++)
            causalityClauses[i]={presentAc,-presentAcNext};
          for(int e1 : qg[e0]){
            if (e1>=0) {
              //the quadrilateral is {a,b,c,d} and the diagonal are ac and bd (a,b,c,d might not be in
              //convex order, but it doesn't matter); since the graph is undirected, the symmetric clause
              //for bd being the current diagonal will be inserted by the symmetric edge
              // if (e1>=n*n) {
              //   std::cout<<"Something off in QG: "<<a<<' '<<c<<' '<<e1<<" triangulation: "<<t<<'\n';
              // }
              int b=e1/n;
              int d=e1-b*n;
              isFlippable[e0]=true;

              int presentBdNext=toIdx(n,distCenter,m,true,e1,s+1,t);

              vi quadrVertices={a,b,c,d};
              bool quadrAtStepS=true;
              for (int j=0;j<quadrVertices.size() && quadrAtStepS;j++) {
                int edgeIdx=pointIndicesToEdgeIdx(n,quadrVertices[j],quadrVertices[(j+1)%quadrVertices.size()]);
                quadrAtStepS=(mightBePresentFromStep[t][edgeIdx]<=s);
              }

              if (quadrAtStepS) {
                int flipBd=toIdx(n,distCenter,m,false,e1,s,t);
                // if (flipBd>numberVar) {
                //   std::cout<<"Too large var "<<flipBd<<" something went wrong\n";
                //   toIdxPrint(n,distCenter,m,false,e1,s,t);
                // }
                causalityClauses[0].push_back(flipBd);
                vi presentVar;
                // std::cout<<"Triangulation "<<t<<", step: "<<s<<", edge: "<<a<<" "<<c<<' '<<b<<' '<<d<<'\n';
                for (int j=0;j<quadrVertices.size();j++) {
                  int edgeIdx=pointIndicesToEdgeIdx(n,quadrVertices[j],quadrVertices[(j+1)%quadrVertices.size()]);
                  presentVar.push_back(toIdx(n,distCenter,m,true,edgeIdx,s,t));
                  causalityClauses[j+1].push_back(toIdx(n,distCenter,m,true,edgeIdx,s,t));
                }

                vi clauseTmp;

                clauseTmp={-presentAc,-flipAc,presentBdNext};
                for (int var : presentVar)
                  clauseTmp.push_back(-var);
                clauses.push_back(clauseTmp);
              }

              //this is the same regardless of what the current diagonal is, so don't insert it twice
              //included in the more general non crossing clauses
              // if (e0<e1) {
              //   clauses.push_back({-presentAcNext,-presentBdNext});
              // }
            }
            else {
              //if the quadrilateral is not convex, we add the following clause:
              //the quadrilateral is missing, or ac is not flipped at this step
              int actualIdx1=-(e1+1);
              int b=actualIdx1/n;
              int d=actualIdx1-b*n;
              vi quadrVertices={a,b,c,d};
              bool quadrAtStepS=true;
              for (int j=0;j<quadrVertices.size() && quadrAtStepS;j++) {
                int edgeIdx=pointIndicesToEdgeIdx(n,quadrVertices[j],quadrVertices[(j+1)%quadrVertices.size()]);
                quadrAtStepS=(mightBePresentFromStep[t][edgeIdx]<=s);
              }
              if (quadrAtStepS) {
                vi presentVar;
                // std::cout<<"Triangulation "<<t<<", step: "<<s<<", edge: "<<a<<" "<<c<<' '<<b<<' '<<d<<'\n';
                for (int j=0;j<quadrVertices.size();j++) {
                  int edgeIdx=pointIndicesToEdgeIdx(n,quadrVertices[j],quadrVertices[(j+1)%quadrVertices.size()]);
                  presentVar.push_back(toIdx(n,distCenter,m,true,edgeIdx,s,t));
                }

                vi clauseTmp={-flipAc};
                for (int var : presentVar)
                  clauseTmp.push_back(-var);
                clauses.push_back(clauseTmp);
              }
            }
          }

          for (const vi& causalityClause : causalityClauses)
            clauses.push_back(causalityClause);

          //independence clauses
          //TODO check for duplicates
          i2 e0Points={a,c};
          for (int i=0;i<2;i++) {
            int u=e0Points[i],v=e0Points[(i+1)%2];
            for (int y=0;y<n;y++) {
              int uyIdx=pointIndicesToEdgeIdx(n,u,y);
              int vyIdx=pointIndicesToEdgeIdx(n,v,y);
              //TODO you should check if qg[e0] and qg[uyIdx] have any entry >=0 to reduce the number of clauses (while maintaining correctness)
              if (!qg[e0].empty() && !qg[uyIdx].empty() &&
                mightBePresentFromStep[t][e0]<=s && mightBePresentFromStep[t][uyIdx]<=s && mightBePresentFromStep[t][vyIdx]<=s) {
                int flipUv=toIdx(n,distCenter,m,false,e0,s,t),
                flipUy=toIdx(n,distCenter,m,false,uyIdx,s,t),
                presentVy=toIdx(n,distCenter,m,true,vyIdx,s,t);
                clauses.push_back({-presentVy,-flipUv,-flipUy});
              }
            }
          }
        }
      }

      //TODO use mightBePresentFromStep here too to reduce the number of clauses
      for (int i=0;i<nSquared;i++){
        int flipEdge=toIdx(n,distCenter,m,false,i,s,t);
        int presentEdge=toIdx(n,distCenter,m,true,i,s,t);
        int presentEdgeNext=toIdx(n,distCenter,m,true,i,s+1,t);
        // std::cout<<"var: "<<flipEdge<<", triangulation "<<t<<" step "<<s<<" edge "<<i/n<<' '<<i%n<<'\n';

        clauses.push_back({-presentEdge,flipEdge,presentEdgeNext});

        clauses.push_back({presentEdge,-flipEdge});
      }

      //edges that are part of no convex quadrilateral can never be flipped
      for (int i=0;i<nSquared;i++) {
        if (!isFlippable[i]) {
          int flipEdge=toIdx(n,distCenter,m,false,i,s,t);
          clauses.push_back({-flipEdge});
        }
      }
    }
  }

  for (int s=1;s<=distCenter;s++) {
    for (int t=0;t<tVec.size();t++) {
      for (int e0=0;e0<nSquared;e0++) {
        if (mightBePresentFromStep[t][e0]<=s) {
          for (int e1 : crossingEdgesVec[e0]) {
            if (mightBePresentFromStep[t][e1]<=s) {
              // if (e1==6*n+45)
              //   std::cout<<"Edge that crosses 6,45: "<<e0/n<<' '<<e0%n<<'\n';
              int presentE0=toIdx(n,distCenter,m,true,e0,s,t);
              int presentE1=toIdx(n,distCenter,m,true,e1,s,t);
              clauses.push_back({-presentE0,-presentE1});
            }
          }
        }
      }
    }
  }

  // std::cout<<"Before crossing clauses we have "<<clauses.size()<<" clauses"<<std::endl;
  //clauses to say that at the end no crossing edge must be present
  for (int t0=0;t0<tVec.size();t0++) {
    for (int t1=t0+1;t1<tVec.size();t1++) {
      for (int e0=0;e0<nSquared;e0++) {
        for (int e1 : qg[e0]) {
          if (e1>=0) {
            int presentAcT0=toIdx(n,distCenter,m,true,e0,distCenter,t0),
            presentBdT1=toIdx(n,distCenter,m,true,e1,distCenter,t1);
            clauses.push_back({-presentAcT0,-presentBdT1});
          }
        }
      }
    }
  }
  std::cout<<"After crossing clauses we have "<<clauses.size()<<" clauses"<<std::endl;
  // for (auto & clause : clauses) {
  //   bool isSat=false;
  //   for (int literal : clause) {
  //     int var=literal>0 ? literal : -literal;
  //     bool isPresentVar;
  //     int u,v,step,triIdx;
  //
  //     fromIdx(var,n,distCenter,m,isPresentVar,u,v,step,triIdx);
  //     if (!isPresentVar)
  //       isSat=true;
  //     else {
  //       if ((literal>0 && testAssignment[var]) || (literal<0 && !testAssignment[var]))
  //         isSat=true;
  //     }
  //   }
  //   if (!isSat) {
  //     std::cout<<"Clause: \n";
  //     for (int literal : clause) {
  //       int var=literal>0 ? literal : -literal;
  //       bool isPresentVar;
  //       int u,v,step,triIdx;
  //
  //       fromIdx(var,n,distCenter,m,isPresentVar,u,v,step,triIdx);
  //
  //
  //       std::cout<<"IsPresent, triangulation: "<<isPresentVar<<' '<<triIdx<<" step "<<step<<" edge: "<<u<<' '<<v<<std::endl;
  //       std::cout<<"variable "<<literal<<"\n";
  //     }
  //     break;
  //   }
  // }
  std::ofstream output(outFileName);
  output<<"p cnf "<<numberVar<<' '<<clauses.size()<<'\n';
  for (auto & clause : clauses) {
    for (int literal : clause)
      output<<literal<<' ';
    output<<"0 ";
  }
  output.close();
}

void cnfToFlipList(const std::string& inFileName, const std::string& idxFileName, std::vector<std::vector<vi2>>& flipsPerRound) {
  std::ifstream input(inFileName);
  std::ifstream inIdx(idxFileName);

  std::string tmp;
  input>>tmp; input>>tmp;
  //this is useless here, but it's present in the file, so I read it
  std::string instanceUid;
  inIdx>>instanceUid;
  // std::cout<<"uid: "<<instanceUid<<'\n';
  int n,m,varCnt;
  inIdx>>n>>m>>varCnt;
  vi distCenter(m);
  for (int i=0;i<m;i++)
    inIdx>>distCenter[i];
  vi4 idxToVar(varCnt);
  for (int i=0;i<varCnt;i++)
    for (int j=0;j<4;j++)
      inIdx>>idxToVar[i][j];
  inIdx.close();
  if (tmp=="SATISFIABLE") {
    input>>tmp;
    flipsPerRound.resize(m);
    for (int i=0;i<m;i++)
      if (distCenter[i]>0)
        flipsPerRound[i]=std::vector<vi2>(distCenter[i]);
      else
        flipsPerRound[i]=std::vector<vi2>(0);
    //mallob format: a single string of comma separated values, enclosed in squared brackets
    if (tmp!="v") {
      std::stringstream ss(tmp);
      char bracket;
      ss>>bracket;
      for (int var; ss >> var;) {
        if (var>0) {
          bool isPresentVar;
          int u,v,step,triIdx;

          fromIdx(var,isPresentVar,u,v,step,triIdx,n,idxToVar);

          if (!isPresentVar && step<distCenter[triIdx]) {
            flipsPerRound[triIdx][step].push_back({u,v});
          }
        }
        if (ss.peek() == ',' || ss.peek() == ']')
          ss.ignore();
      }
    }
    else
      while (tmp!="0") {
        if (tmp!="v") {
          int var=atoi(tmp.c_str());
          if (var>0) {
            bool isPresentVar;
            int u,v,step,triIdx;

            fromIdx(var,isPresentVar,u,v,step,triIdx,n,idxToVar);

            if (!isPresentVar && step<distCenter[triIdx]) {
              flipsPerRound[triIdx][step].push_back({u,v});
            }
          }
        }
        input>>tmp;
      }
  }
  // for (int i=0;i<m;i++)
  //   for (int j=0;j<distCenter[i];j++)
  //     for (i2 flip : flipsPerRound[i][j])
  //       std::cout<<i<<' '<<j<<' '<<flip[0]<<' '<<flip[1]<<'\n';
  inIdx.close();
  input.close();
}

void cnfToFlipList(const std::string& inFileName, size_t n, size_t m, int distCenter, std::vector<std::vector<vi2>>& flipsPerRound) {
  std::ifstream input(inFileName);
  std::string tmp;
  input>>tmp; input>>tmp;
  if (tmp=="SATISFIABLE") {
    input>>tmp;
    for (int i=0;i<m;i++)
      flipsPerRound[i]=std::vector<vi2>(distCenter);
    while (tmp!="0") {
      if (tmp!="v") {
        int var=atoi(tmp.c_str());
        if (var>0) {
          bool isPresentVar;
          int u,v,step,triIdx;

          fromIdx(var,n,distCenter,m,isPresentVar,u,v,step,triIdx);
          // std::cout<<u<<' '<<v<<'\n';
          // if (isPresentVar && step==1) {
          //   // std::cout<<var<<std::endl;
          //   // std::cout<<isPresentVar<<' '<<triIdx<<' '<<step<<' '<<u<<' '<<v<<std::endl;
          //   std::cout<<"This is a present variable. Triangulation, step, edge: "<<triIdx<<' '<<step<<' '<<u<<' '<<v<<std::endl;
          //
          // }
          if (!isPresentVar && step<distCenter) {
            // std::cout<<"Triangulation, step, edge: "<<triIdx<<' '<<step<<' '<<u<<' '<<v<<std::endl;
            // std::cout<<"variable "<<var<<"\n";
            flipsPerRound[triIdx][step].push_back({u,v});
          }
        }
      }
      input>>tmp;
    }
  }
  input.close();
}
