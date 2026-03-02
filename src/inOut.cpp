#include "../include/inOut.hpp"
#include <CGAL/draw_triangulation_2.h>


void edgeVectorToTriangulation(const vi2& edgeVector, const vP& points, size_t nEdges, ConstrainedTriangulation& ct) {
  // CGAL::Triangulation_2<IK,CGAL::Default> t;
  // CGAL::Default& triangulationDS=t.tds();

  for (int i=0;i<nEdges;i++)
    ct.insert_constraint(points[edgeVector[i][0]],points[edgeVector[i][1]]);

  for (auto e = ct.finite_edges_begin(); e != ct.finite_edges_end(); ++e) {
    ct.remove_constrained_edge(e->first,e->second);
  }

  assert(ct.is_valid());
  // CGAL::draw(ct);
  // triangulationDS=ct.tds();
  // printEdges(ct);
}

vConstrainedTriangulation parseInput(const std::string& filename, size_t& n, pointToIdx& idxMap, vP& points, std::string& instance_uid) {
  std::ifstream input(filename);
  nlohmann::json data=nlohmann::json::parse(input);
  instance_uid=data["instance_uid"];
  // int firstCoordinate=data["points_x"][0];
  n=data["points_x"].size();
  // vP points(n);
  points.resize(n);
  //I use a constrained triangulation since I am given an edge list. But unless I'm mistaken,
  //it does not allow vertices with information (index)
  //So I store the mapping from point to index in a separate unordered map
  for (int i=0;i<n;i++) {
    points[i]={data["points_x"][i].template get<int>(),data["points_y"][i].template get<int>()};
    idxMap[points[i]]=i;
  }
  size_t m=data["triangulations"].size();

  //the number of edges is constant for plane straight line triangulations of n fixed vertices
  size_t nEdges=data["triangulations"][0].size();
  vConstrainedTriangulation triangulations(m);
  for (int i=0;i<m;i++) {
    assert(data["triangulations"][i].size()==nEdges);
    vi2 edgeVector(nEdges);
    for (int j=0;j<nEdges;j++) {
      edgeVector[j]={data["triangulations"][i][j][0].template get<int>(),data["triangulations"][i][j][1].template get<int>()};
    }
    edgeVectorToTriangulation(edgeVector,points,nEdges,triangulations[i]);
  }

  // std::cout<<firstCoordinate<<std::endl;
  input.close();
  return triangulations;
}

void writeOutput(const std::vector<std::vector<vi2>>& flipsPerRound, const std::string& filename, std::string instance_uid, std::string strategy) {
  std::ofstream output(filename);
  nlohmann::json data;
  data["content_type"]="CGSHOP2026_Solution";
  data["instance_uid"]=instance_uid;
  // output<<std::setw(4)<<data["content_type"];
  // output<<std::setw(4)<<data["instance_uid"];

  for (int i=0;i<flipsPerRound.size();i++) {
    // std::cout<<flipsPerRound[i].size()<<"\n";
    if (flipsPerRound[i].empty()) {
      data["flips"][i]=nlohmann::json::array();
    }
    else {
      int skipEmpty=0;
      for (int j=0;j<flipsPerRound[i].size();j++) {
        if (!flipsPerRound[i][j].empty())
          for (int k=0;k<flipsPerRound[i][j].size();k++) {
            data["flips"][i][j-skipEmpty][k]={flipsPerRound[i][j][k][0],flipsPerRound[i][j][k][1]};
          }
        else
          skipEmpty++;
      }
    }

  }
  // output<<data["flips"];
  data["meta"]["author"]="ETHFlippers";
  data["meta"]["algorithm"]=strategy;
  output<<data;
  // output<<std::setw(4)<<data["meta"];
  // output<<std::setw(2)<<data;
  // output<<data.dump(0);
  // std::cout<<data;
  output.close();
}