#include "../include/genRandomTriangulations.hpp"
#include "../include/geometryUtils.hpp"

//For initial testing, something similar will probably be inside the utils we'll be given
std::vector<IPoint> genRandomPoints(int n, double minCoord, double maxCoord){
  std::random_device rd;
  std::mt19937 rng(rd());
  std::uniform_real_distribution<double> distribution(minCoord,maxCoord);

  std::vector<IPoint> points(n);
  for(int i=0;i<n;i++) {
    double x,y;
    x = distribution(rng);
    y = distribution(rng);
    std::cout<<x<<" "<<y<<std::endl;
    points[i]={P(x,y),i};
  }

  return points;
}

//genarate random points with integer coordinates (for debugging/testing)
std::vector<IPoint> genRandomPoints(int n, int minCoord, int maxCoord){
  std::random_device rd;
  std::mt19937 rng(rd());
  std::uniform_int_distribution<> distribution(minCoord,maxCoord);

  std::vector<IPoint> points(n);
  for(int i=0;i<n;i++) {
    double x,y;
    x = distribution(rng);
    y = distribution(rng);
    std::cout<<x<<" "<<y<<std::endl;
    points[i]={P(x,y),i};
  }

  return points;
}

//This is also for initial testing. Scans all the edges and
//if the edge can be flipped, sample uar in [0,1], if >=(1-flipProbability) flips the edge
//This is repeated for nIterations
void randomFlips(Triangulation& t, int nIterations, double flipProbability) {
  double doNotFlipProbability = 1-flipProbability;
  std::random_device rd;
  std::mt19937 rng(rd());
  std::uniform_real_distribution<double> distribution(0,1);
  for (int i=0;i<nIterations;i++) {
    for (auto e = t.finite_edges_begin(); e != t.finite_edges_end(); ++e) {
      if (distribution(rng)>=doNotFlipProbability && canBeFlipped(t,e)) {
        t.flip(e->first,e->second);
      }
    }
  }
}
