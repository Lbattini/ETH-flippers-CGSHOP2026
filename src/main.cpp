#include "geometryUtils.hpp"
#include "../include/korderDelaunay.hpp"
#include "../include/inOut.hpp"
#include "../include/differentDiagonals.hpp"
#include "../include/flippingDS.hpp"
#include "../include/quadrilateralGraph.hpp"
#include "../include/generateCNF.hpp"

//include the header below if needed
// #include "../include/genRandomTriangulations.hpp"
#define INF 100000000

void updateBestCenter(int m, const vi& totDistanceVec, const std::vector<std::vector<std::vector<vi2>>>& curFlips, int& totDistance, std::vector<std::vector<vi2>>& flipsPerRound, bool verbose, bool improveCenterVar, int baseCenter) {
  for (int i=0;i<totDistanceVec.size();i++)
    if ((totDistanceVec[i]<totDistance)) {
      totDistance=totDistanceVec[i];
      flipsPerRound=curFlips[i];
      if (verbose) {
        int center=baseCenter+i;
        if(improveCenterVar)
          std::cout<<"With improved center "<<center<<" the pfd is "<<totDistanceVec[i]<<std::endl;
        else
          std::cout<<"With center "<<center<<" the pfd is "<<totDistanceVec[i]<<std::endl;
      }
    }
}

int reduceItsFromFixedCenter(int n, int m, vP& points, bool verbose, std::vector<facesDS>& triVectors, facesDS& triCenter, std::vector<std::vector<vi2>>& curFlips, int center, int previousBestTotDistance, const std::vector<vi2>& flipsCenter) {
  // std::vector<vi2> flipsCenter=curFlips[center];
  bool allEqual=true;
  int totDistance=0;

  mii intersectionsWithCenter;

  for (int i=0;i<m;i++) {
    //if (i!=center) {
    if (totDistance>previousBestTotDistance)
      return INF;

    bool curEqual;
    facesDS triVectorsICopy(triVectors[i]), centerCopy(triCenter);
    mii intersectionsWithTriVectorsI;

    // int curDistance=flipToReduceIntesectionsBiDir(triVectors[i],centerCopy,n,points,curFlips[i]);
    int curDistance=flipToReduceIntesections(triVectors[i],triCenter,n,points,curFlips[i],intersectionsWithCenter);
    std::vector<vi2> curFlipsTmp;
    // int curDistance1=flipToReduceIntesections(triVectorsICopy,centerCopy,n,points,curFlipsTmp,1);
    int curDistance1=flipToReduceIntesectionsRevDir(centerCopy,triVectorsICopy,n,points,curFlipsTmp,intersectionsWithTriVectorsI);
    if (triVectorsICopy.equal(centerCopy) && curDistance1<curDistance) {
      curDistance=curDistance1;
      curFlips[i]=curFlipsTmp;
    }
    if(i==center && flipsCenter.size()<curDistance){
      curDistance=flipsCenter.size();
      curFlips[i]=flipsCenter;
    }
    totDistance+=curDistance;
    // curEqual=triVectors[i].equal(centerCopy);
    curEqual=triVectors[i].equal(triCenter);
    if (curEqual){
      if (verbose)
        std::cout<<center<<' '<<i<<" are equal and at distance "<<curDistance<<"\n";
    }
    else {
      if (verbose)
        std::cout<<center<<' '<<i<<" are different\n";
      allEqual=false;
    }

    //}
  }
  if (allEqual)
    return totDistance;
  return INF;
}

int reduceItsFromFixedCenter(int n, int m, vP& points, bool verbose, std::vector<facesDS> triVectors, facesDS& triCenter, std::vector<std::vector<vi2>>& curFlips, int center, int previousBestTotDistance, bool tryAllPrefForPairwise=false) {
  bool allEqual=true;
  int totDistance=0;

  mii intersectionsWithCenter;

  for (int i=0;i<m;i++) {
    if (totDistance>previousBestTotDistance)
      return INF;

    bool curEqual;
    facesDS initialTriVectorI(triVectors[i]), initialCenter(triCenter);
    facesDS triVectorsICopy(triVectors[i]), centerCopy(triCenter);
    mii intersectionsWithTriVectorsI;

    int curDistance=flipToReduceIntesections(triVectors[i],triCenter,n,points,curFlips[i],intersectionsWithCenter);
    std::vector<vi2> curFlipsTmp;
    int curDistance1=flipToReduceIntesectionsRevDir(centerCopy,triVectorsICopy,n,points,curFlipsTmp,intersectionsWithTriVectorsI);
    if (triVectorsICopy.equal(centerCopy) && curDistance1<curDistance) {
      curDistance=curDistance1;
      curFlips[i]=curFlipsTmp;
    }

    if (tryAllPrefForPairwise) {
      for (int j=1;j<curDistance;j++) {
        std::vector<vi2> curFlipsTmp1;
        facesDS triVectorsICopy1(initialTriVectorI), centerCopy1(initialCenter);
        curDistance1=flipToReduceIntesections(triVectorsICopy1,centerCopy1,n,points,curFlipsTmp1,intersectionsWithCenter,j);
        if (triVectorsICopy1.equal(centerCopy1) && curDistance1<curDistance) {
          std::cout<<curDistance<<' '<<curDistance1<<'\n';
          curDistance=curDistance1;
          curFlips[i]=curFlipsTmp1;
        }
      }
    }

    totDistance+=curDistance;
    curEqual=triVectors[i].equal(triCenter);
    if (curEqual){
      if (verbose)
        std::cout<<center<<' '<<i<<" are equal and at distance "<<curDistance<<"\n";
    }
    else {
      if (verbose)
        std::cout<<center<<' '<<i<<" are different\n";
      allEqual=false;
    }
  }
  if (allEqual)
    return totDistance;
  return INF;
}


bool parallelReduceItsLessRam(vConstrainedTriangulation& triangulations, pointToIdx& idxMap, int n, vP& points, bool verbose, std::vector<std::vector<vi2>>& flipsPerRound, int& totDistance, bool improveCenterVar) {
  int m=triangulations.size();
  int step=16;
  int firstI=0;
  if (m==200) {
    firstI=64;
    std::vector<facesDS> triVectors;
    std::vector<std::vector<vi2>> curFlips(m);
    for (int k=0;k<m;k++)
      triVectors.push_back(facesDS(triangulations[k],idxMap));

    int center=10;
    std::vector<vi2> flipsCenter;
    facesDS triCenter(triVectors[center]);
    if (improveCenterVar) {
      totDistance=improveCenter(triCenter,triVectors,n,points,flipsCenter);
      if (verbose)
        std::cout<<"To improve center "<<center<<", "<<totDistance<<" flips were performed\n";
    }
    totDistance=reduceItsFromFixedCenter(n,m,points,verbose,triVectors,triCenter,curFlips,center,INF,flipsCenter);

    flipsPerRound=curFlips;
  }

  for (int i=firstI;i<m;i+=step) {
    int supCenter=std::min(i+step,m);
    int nCenters=supCenter-i;
    std::vector<int> totDistanceVec(nCenters,0);
    std::vector<std::vector<facesDS>> triVectors(nCenters);
    std::vector<std::vector<std::vector<vi2>>> curFlips(nCenters,std::vector<std::vector<vi2>>(m));
    for (int j=0;j<nCenters;j++)
      for (int k=0;k<m;k++)
        triVectors[j].push_back(facesDS(triangulations[k],idxMap));


    #pragma omp parallel for
    for (int center=i;center<supCenter; center++) {
      facesDS triCenter(triVectors[center-i][center]);
      std::vector<vi2> flipsCenter;
      if (improveCenterVar) {
        //facesDS centerBeforeImpr(triVectors[center-i][center]);
        totDistanceVec[center-i]=improveCenter(triCenter,triVectors[center-i],n,points,flipsCenter);
        //facesDS centerAfterImpr(triVectors[center-i][center]), centerCopy(triCenter);
        if (verbose)
          std::cout<<"To improve center "<<center<<", "<<totDistanceVec[center-i]<<" flips were performed\n";
      }
      totDistanceVec[center-i]=reduceItsFromFixedCenter(n,m,points,verbose,triVectors[center-i],triCenter,curFlips[center-i],center, totDistance, flipsCenter);
    }

    updateBestCenter(m,totDistanceVec,curFlips,totDistance,flipsPerRound,verbose,improveCenterVar,i);

  }

  return totDistance<INF;
}

bool parallelReduceItsLessRam(vConstrainedTriangulation& triangulations, pointToIdx& idxMap, int n, vP& points, bool verbose, std::vector<std::vector<vi2>>& flipsPerRound, int& totDistance) {
  totDistance=INF;
  return parallelReduceItsLessRam(triangulations,idxMap,n,points,verbose,flipsPerRound,totDistance,true) ||
    parallelReduceItsLessRam(triangulations,idxMap,n,points,verbose,flipsPerRound,totDistance,false);
}

bool parallelReduceItsLessRam(vConstrainedTriangulation& triangulations, pointToIdx& idxMap, int n, vP& points, bool verbose, std::vector<std::vector<vi2>>& flipsPerRound, int& totDistance, std::vector<facesDS>& candidateCenters, const std::vector<std::pair<int,int>>& centerIndices, int nCandidateCenters) {
  int m=triangulations.size();
  int step=32;
  int firstI=0;
  int lastI=firstI+nCandidateCenters;

  if (n==6000 && m==20) {
    // int center=179;
    // facesDS triCenter(candidateCenters[centerIndices[center].second]);
    //
    //
    // std::vector<facesDS> triVectors;
    // for (int k=0;k<m;k++)
    //   triVectors.push_back(facesDS(triangulations[k],idxMap));
    //
    // std::vector<std::vector<std::vector<vi2>>> curFlips(1,std::vector<std::vector<vi2>>(m));
    // totDistance=reduceItsFromFixedCenter(n,m,points,verbose,triVectors,triCenter,curFlips[0],center, totDistance,true);
    //
    // std::vector<int> totDistanceVec(1,totDistance);
    // updateBestCenter(m,totDistanceVec,curFlips,totDistance,flipsPerRound,verbose,false,0);
    // std::cout<<"With center "<<center<<" the pfd is: "<<totDistance<<std::endl;
    step=16;
    firstI=174;
    lastI=firstI+step;
  }
  for (int i=firstI;i<lastI;i+=step) {
    std::cout<<i<<"/"<<lastI-1<<std::endl;
    int supCenter=std::min(i+step,lastI);
    int nCenters=supCenter-i;
    std::vector<int> totDistanceVec(nCenters,0);
    std::vector<std::vector<facesDS>> triVectors(nCenters);
    std::vector<std::vector<std::vector<vi2>>> curFlips(nCenters,std::vector<std::vector<vi2>>(m));
    for (int j=0;j<nCenters;j++)
      for (int k=0;k<m;k++)
        triVectors[j].push_back(facesDS(triangulations[k],idxMap));


    #pragma omp parallel for
    for (int center=i;center<supCenter; center++) {
      facesDS triCenter(candidateCenters[centerIndices[center].second]);


      totDistanceVec[center-i]=reduceItsFromFixedCenter(n,m,points,verbose,triVectors[center-i],triCenter,curFlips[center-i],center, totDistance);

      if (totDistanceVec[center-i]<=totDistance+m) {
        curFlips[center-i]=std::vector<std::vector<vi2>>(m);
        totDistanceVec[center-i]=reduceItsFromFixedCenter(n,m,points,verbose,triVectors[center-i],triCenter,curFlips[center-i],center, totDistance,true);
      }
    }

    updateBestCenter(m,totDistanceVec,curFlips,totDistance,flipsPerRound,verbose,false,i);

  }

  return totDistance<INF;
}
bool parallelReduceIts(vConstrainedTriangulation& triangulations, pointToIdx& idxMap, int n, vP& points, bool verbose, std::vector<std::vector<vi2>>& flipsPerRound, int& totDistance) {
  int m=triangulations.size();
  std::vector<bool> allEqual(m,false);
  totDistance=INF;
  std::vector<int> totDistanceVec(m,INF);
  std::vector<std::vector<facesDS>> triVectors(m);
  std::vector<std::vector<std::vector<vi2>>> curFlips(m,std::vector<std::vector<vi2>>(m));

  for (int j=0;j<m;j++)
    for (int i=0;i<m;i++)
      triVectors[j].push_back(facesDS(triangulations[i],idxMap));

  #pragma omp parallel for
  for (int center=0;center<m;center++) {
    allEqual[center]=true;
    totDistanceVec[center]=0;

    mii intersectionsWithCenter;

    for (int i=0;i<m;i++) {
      if (i!=center) {
        bool curEqual;

        int curDistance=flipToReduceIntesections(triVectors[center][i],triVectors[center][center],n,points,curFlips[center][i],intersectionsWithCenter);
        // int curDistance=flipToReduceIntesections(triVector[i],triVector[center],n,points,curFlips[i]);

        // std::cout<<curDistance<<'\n';
        totDistanceVec[center]+=curDistance;
        curEqual=triVectors[center][i].equal(triVectors[center][center]);
        if (curEqual){
          if (verbose)
            std::cout<<center<<' '<<i<<" are equal and at distance "<<curDistance<<"\n";
        }
        else {
          if (verbose)
            std::cout<<center<<' '<<i<<" are different\n";
          allEqual[center]=false;
        }

      }
    }

    // if (n>2000)
    //   break;
  }

  for (int center=0;center<m;center++)
    if (allEqual[center] && (totDistanceVec[center]<totDistance)) {
      totDistance=totDistanceVec[center];
      flipsPerRound=curFlips[center];
      if (verbose)
        // if(center>=m)
        //   std::cout<<"With improved center "<<center-m<<" the pfd is "<<totDistanceVec[center]<<'\n';
        // else
          std::cout<<"With center "<<center<<" the pfd is "<<totDistanceVec[center]<<'\n';
    }
  curFlips=std::vector<std::vector<std::vector<vi2>>>(m,std::vector<std::vector<vi2>>(m));
  
  #pragma omp parallel for
  for (int center=0;center<m;center++) {
    allEqual[center]=true;
    totDistanceVec[center]=improveCenter(triVectors[center][center],triVectors[center],n,points,curFlips[center][center]);

    if (verbose)
      std::cout<<"To improve center "<<center<<", "<<totDistanceVec[center]<<" flips were performed\n";

    mii intersectionsWithCenter;

    for (int i=0;i<m;i++) {
      if (i!=center) {
        bool curEqual;

        int curDistance=flipToReduceIntesections(triVectors[center][i],triVectors[center][center],n,points,curFlips[center][i],intersectionsWithCenter);

        totDistanceVec[center]+=curDistance;
        curEqual=triVectors[center][i].equal(triVectors[center][center]);
        if (curEqual){
          if (verbose)
            std::cout<<center<<' '<<i<<" are equal and at distance "<<curDistance<<"\n";
        }
        else {
          if (verbose)
            std::cout<<center<<' '<<i<<" are different\n";
          allEqual[center]=false;
        }

      }
    }

    // if (n>2000)
    //   break;
  }

  for (int center=0;center<m;center++)
    if (allEqual[center] && (totDistanceVec[center]<totDistance)) {
      totDistance=totDistanceVec[center];
      flipsPerRound=curFlips[center];
      if (verbose)
        // if(center>=m)
          std::cout<<"With improved center "<<center-m<<" the pfd is "<<totDistanceVec[center]<<'\n';
        // else
        //   std::cout<<"With center "<<center<<" the pfd is "<<totDistanceVec[center]<<'\n';
    }
  return totDistance<INF;
}

int reduceItsFromCenterWithHistory(size_t n, bool verbose, std::vector<std::vector<vi2>> &curFlips, bool &allEqual, const vP& points, std::vector<facesDS> triVector, std::vector<facesDS>& candidateCenters, std::vector<std::pair<int, int>>& candidateCentersWithTotDist, int &maxDistance, int &centerMinIts, int& cntMinIts, const std::vector<facesDS>& initialTriVector, facesDS triCenter, int center) {
  allEqual=true;
  int m=initialTriVector.size();
  // std::vector<std::vector<vi2>> curFlips(m);
  int totDistanceFromCenter=0;
  mii intersectionsWithCenter;

  facesDS improvedCenter(triCenter);

  if (center>=0 && center<m)
    improveCenter(improvedCenter,triVector,n,points,curFlips[center]);
  else {
    std::vector<vi2> curFlipsTmp;
    improveCenter(improvedCenter,triVector,n,points,curFlipsTmp);
  }


  for (int i=0;i<m;i++) {
    bool curEqual;

    std::vector<facesDS> tHistory=flipToReduceIntesectionsWithHistory(triVector[i],improvedCenter,n,points,curFlips[i],intersectionsWithCenter);

    int curDistance=tHistory.size();

    if (curDistance>maxDistance) {
      maxDistance=curDistance;
    }

    for (int ii=0;ii<curDistance-1;ii++) {
      int cntTotIntersections=0;
      for (int t=0;t<initialTriVector.size(); t++)
        cntTotIntersections+=tHistory[ii].cntIntersections(initialTriVector[t]);
      if (cntTotIntersections<cntMinIts) {
        std::cout<<"Previous cntMinIts "<<cntMinIts<<", current: "<<cntTotIntersections<<'\n';
        cntMinIts=cntTotIntersections;
        centerMinIts=candidateCenters.size();
      }
      candidateCentersWithTotDist.push_back({cntTotIntersections,candidateCenters.size()});
      candidateCenters.push_back(tHistory[ii]);
    }

    totDistanceFromCenter+=curDistance;
    curEqual=triVector[i].equal(improvedCenter);
    if (curEqual){
      if (verbose)
        std::cout<<center<<' '<<i<<" are equal and at distance "<<curDistance<<"\n";
    }
    else {
      if (verbose)
        std::cout<<center<<' '<<i<<" are different\n";
      allEqual=false;
    }
  }

  // if (allEqual && (totDistanceFromCenter<totDistance)) {
  //   totDistance=totDistanceFromCenter;
  //   flipsPerRound=curFlips;
  //   if (verbose)
  //     std::cout<<"With center "<<center<<" the pfd is "<<totDistanceFromCenter<<'\n';
  // }

  for (int i=0;i<m;i++)
    triVector[i]=facesDS(initialTriVector[i]);
  return totDistanceFromCenter;
}

int main(int argc, char * argv[]) {
  // int n,minCoord,maxCoord;
  // std::cin >> n >> minCoord >> maxCoord;
  // std::vector<IPoint> points=genRandomPoints(n,minCoord,maxCoord);
  // Triangulation t=scanTriangulation(points);
  // // printEdges(t);
  // CGAL::draw(t);
  // randomFlips(t,1,1);
  // // printEdges(t);
  // CGAL::draw(t);

  size_t n;
  pointToIdx idxMap;
  std::string inputFile;
  if (argc>=2)
    inputFile=argv[1];
  else
    inputFile="../../example_instances_cgshop2026/examples/example_ps_20_nt2_pfd5_random.json";

  std::cout<<"Computing a solution for file "<<inputFile<<'\n';
  std::string outputFile;
  if (argc>=3)
    outputFile=argv[2];
  else
    //todo only use the actual file name and change folder
      outputFile="./outputFiles/output0.json";

  std::string strategy="default";
  if (argc>=4)
    strategy=argv[3];

  bool verbose=true;

  if (argc>=5 && std::string(argv[4])=="nonVerbose") {
    verbose=false;
  }

  std::vector<std::vector<vi2>> flipsPerRound;
  std::string instance_uid;
  int totDistance=0;
  bool allEqual=true;
  if (strategy=="assignmentToFlipList") {
    int distCenter=1;
    if (argc>=6)
      distCenter=atoi(argv[5]);
    n=atoi(argv[6]);
    int m=2;
    if (argc>=8)
      m=atoi(argv[7]);
    instance_uid=argv[8];
    // flipsPerRound=std::vector<std::vector<vi2>>(m);
    flipsPerRound.resize(m);
    std::cout<<"n m distance to center "<<n<<' '<<m<<' '<<distCenter<<std::endl;
    cnfToFlipList(inputFile,n,m,distCenter,flipsPerRound);
  }
  else {
    if (strategy=="assignmentToFlipListNU") {
      std::cout<<"Converting assignment to flip list\n";
      std::string idxFileName=argv[5];
      std::ifstream inIdx(idxFileName);
      inIdx>>instance_uid;
      // std::cout<<idxFileName<<'\n';
      std::cout<<"The instance_uid is: "<<instance_uid<<'\n';
      inIdx.close();
      cnfToFlipList(inputFile,idxFileName,flipsPerRound);
    }
    else{
      vP points;
      vConstrainedTriangulation triangulations=parseInput(inputFile,n,idxMap,points,instance_uid);

      flipsPerRound=std::vector<std::vector<vi2>>(triangulations.size());
      // std::vector<CGALconstrainedTriangulation> triCGAL;
      std::vector<facesDS> triVector;
      for (int i=0;i<triangulations.size();i++) {
        // triCGAL.push_back(CGALconstrainedTriangulation(triangulations[i],idxMap));
        triVector.push_back(facesDS(triangulations[i],idxMap));
      }

      // std::cout<<"There are "<<n<<" points in the input\n";
      // vvi qg=buildQG(points);
      // for (int i=0;i<triangulations.size();i++) {
      //   for (int j=i+1;j<triangulations.size();j++) {
      //     int lowerBound=longestSP(triVector[i],triVector[j],qg,points.size());
      //     std::cout<<"Lower bound between "<<i<<" and "<<j<<" is: "<<lowerBound<<std::endl;
      //   }
      // }


      // switch (strategy) {
      //   case "flipDiff":;
      //   default: ;
      // }
      if (strategy=="generateCNF") {
        std::cout<<"There are "<<n<<" points in the input\n";
        std::string emptyQuadrFileName="./emptyQuadr/"+instance_uid+".txt";
        vvvb emptyTriangle(n,vvb(n,vb(n,true)));
        vvi qg=buildQGPlusNonConvex(points,emptyQuadrFileName,emptyTriangle);
        vvi crossingEdgesVec=crossingEdges(points);
        int distCenter=1;
        if (argc>=6)
          distCenter=atoi(argv[5]);
        std::cout<<"The distance to center is: "<<distCenter<<'\n';
        genCNF(points,qg,crossingEdgesVec,triVector,n,distCenter,outputFile);
      }
      else {
        if (strategy=="generateCNFNU") {
          std::cout<<"There are "<<n<<" points in the input"<<std::endl;
          std::string idxFileName=argv[5];
          std::ofstream outIdx{idxFileName};
          outIdx<<instance_uid<<std::endl;
          outIdx.close();
          int m=triangulations.size();
          vi distCenter(m,1);
          if (argc==7)
            for (int i=0;i<m;i++)
              distCenter[i]=atoi(argv[6]);
          else
            for (int i=0;i<m;i++)
              distCenter[i]=atoi(argv[6+i]);
          std::cout<<"Distances to center: ";
          for (int i=0;i<m;i++)
            std::cout<<distCenter[i]<<' ';
          std::cout<<std::endl;
          std::string emptyQuadrFileName="./emptyQuadr/"+instance_uid+".txt";
          vvvb emptyTriangle(n,vvb(n,vb(n,true)));
          vvi qg=buildQGPlusNonConvex(points,emptyQuadrFileName,emptyTriangle);
          std::cout<<"Computed QG"<<std::endl;
          vvi crossingEdgesVec=crossingEdges(points);
          std::cout<<"Computed crossing edges"<<std::endl;
          genCNF(points,qg,emptyTriangle,crossingEdgesVec,triVector,n,distCenter,outputFile,idxFileName);
        }
        else{
          if (strategy=="flipDiff") {
            allEqual=false;
            for (int center=0;center<triangulations.size() && !allEqual;center++) {
              // flipsPerRound.clear();
              // std::vector<std::vector<vi2>> curFlips;
              allEqual=true;
              for (int i=0;i<triangulations.size();i++) {
                if (i!=center) {
                  bool curEqual;
                  int curDistance=flipDifferentDiagonals(triVector[i],triVector[center],n,points,flipsPerRound[i]);
                  if (verbose)
                    std::cout<<curDistance<<'\n';
                  totDistance+=curDistance;
                  // CGAL::draw(triangulations[i]);
                  curEqual=triVector[i].equal(triVector[center]);
                  if (curEqual){
                    if (verbose)
                      std::cout<<center<<' '<<i<<" are equal\n";
                  }
                  else {
                    if (verbose)
                      std::cout<<center<<' '<<i<<" are different\n";
                    allEqual=false;
                  }
                }
              }
            }
          }
          else {
            if (strategy=="reduceIts") {
              allEqual=false;
              totDistance=INF;

              // std::cout<<"Here"<<std::endl;
              allEqual=parallelReduceItsLessRam(triangulations,idxMap,n,points,verbose,flipsPerRound,totDistance);

            }
            else {
              if (strategy=="reduceItsCenters") {
                totDistance=INF;
                std::vector<facesDS> candidateCenters;
                std::vector<std::pair<int,int>> candidateCentersPfdIdx;
                int maxDistance=0;
                int centerMinIts=0,cntMinIts=INF;
                std::vector<facesDS> initialTriVector;
                for (int i=0;i<triangulations.size();i++) {
                  initialTriVector.push_back(facesDS(triangulations[i],idxMap));
                }

                std::vector<std::vector<std::vector<vi2>>> curFlips(triangulations.size(),std::vector<std::vector<vi2>>(triangulations.size()));
                std::vector<int> totDistanceVec(triangulations.size());
                std::vector<int> allEqualVec(triangulations.size());
                for (int center=0;center<triangulations.size();center++) {
                  bool allEqualTmp;
                  totDistanceVec[center]=reduceItsFromCenterWithHistory(n, verbose, curFlips[center], allEqualTmp, points,
                                                 triVector, candidateCenters, candidateCentersPfdIdx, maxDistance, centerMinIts,
                                                 cntMinIts, initialTriVector, triVector[center], center);
                  allEqualVec[center]=allEqualTmp;
                }
                for (int center=0;center<triangulations.size();center++) {
                  if (allEqualVec[center] && (totDistanceVec[center]<totDistance)) {
                    totDistance=totDistanceVec[center];
                    flipsPerRound=curFlips[center];
                    if (verbose)
                      std::cout<<"With center "<<center<<" the pfd is "<<totDistance<<'\n';
                  }
                }

                std::vector<std::vector<vi2>> curFlipsMinIts(triangulations.size());
                int totDistMinIts=reduceItsFromFixedCenter(n,triangulations.size(),points,verbose,triVector,candidateCenters[centerMinIts],curFlipsMinIts,-1,totDistance);
                if (allEqual && (totDistMinIts<totDistance)) {
                  totDistance=totDistMinIts;
                  flipsPerRound=curFlipsMinIts;
                  if (verbose)
                    std::cout<<"With center with minimum intersections, "<<" the pfd is "<<totDistMinIts<<'\n';
                }

                std::sort(candidateCentersPfdIdx.begin(),candidateCentersPfdIdx.end());

                parallelReduceItsLessRam(triangulations,idxMap,n,points,verbose,flipsPerRound,totDistance,candidateCenters, candidateCentersPfdIdx, 128);

              }
              else {
                if (strategy=="omniDir") {
                  allEqual=false;
                  totDistance=INF;

                  for (int center=0;center<triangulations.size();center++) {
                    allEqual=true;
                    std::vector<std::vector<vi2>> curFlips(triangulations.size());
                    int totDistanceFromCenter=0;


                    bool curEqual=true;

                    int curDistance=omniDirFlips(center,triVector,n,points,curFlips);

                    // std::cout<<curDistance<<'\n';
                    totDistanceFromCenter+=curDistance;
                    for (int i=0;i<triangulations.size();i++) {
                      curEqual=curEqual && triVector[i].equal(triVector[center]);
                      if (!curEqual) {
                        if (verbose)
                          std::cout<<center<<' '<<i<<" are different\n";
                        allEqual=false;
                      }
                    }

                    if (allEqual && (totDistanceFromCenter<totDistance)) {
                      totDistance=totDistanceFromCenter;
                      flipsPerRound=curFlips;
                      if (verbose)
                        std::cout<<"With center "<<center<<" the pfd is "<<totDistanceFromCenter<<'\n';
                    }
                    for (int i=0;i<triangulations.size();i++)
                      triVector[i]=facesDS(triangulations[i],idxMap);
                    if (n>2000)
                      break;
                  }
                }
                else {
                  int k;
                  std::cin>>k;
                  for (int i=0;i<triangulations.size();i++) {
                    int curDistance=korderDelaunay(triangulations[i],k,n,idxMap);
                    if (verbose)
                      std::cout<<curDistance<<'\n';
                    totDistance+=curDistance;
                    if (i>0 && !sameEdges(triangulations[i],triangulations[i-1]))
                      allEqual=false;
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if (strategy!="generateCNF" && strategy!="generateCNFNU")
    writeOutput(flipsPerRound,outputFile,instance_uid,strategy);

  if (strategy!="generateCNF" && strategy!="generateCNFNU" &&
    strategy!="assignmentToFlipList" && strategy!="assignmentToFlipListNU") {
    if (allEqual)
      std::cout<<"All the triangulations coincide\n";
    else
      std::cout<<"Not all the triangulations coincide\n";
    std::cout<<"The pfd is: "<<totDistance<<'\n';
  }

  // edgeVectorToTriangulation();
}
