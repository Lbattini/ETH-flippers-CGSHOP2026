import os
import argparse
import cgshop2026_pyutils
from cgshop2026_pyutils.schemas import CGSHOP2026Instance, CGSHOP2026Solution
from cgshop2026_pyutils.geometry import (    FlippableTriangulation,
    draw_flips,
    Point,
    expand_edges_by_convex_hull_edges,
    is_triangulation,
)
from cgshop2026_pyutils.verify import check_for_errors

from cgshop2026_pyutils.io import read_instance, read_solution

from matplotlib import pyplot as plt
from cgshop2026_pyutils.visualize import create_instance_plot

parser = argparse.ArgumentParser()
parser.add_argument("fileIn")
parser.add_argument("fileOut")
parser.add_argument("--visualize",action="store_true")
args = parser.parse_args()

#fileIn=input("Input file: ")
#fileOut=input("Output file: ")
fileIn=args.fileIn
fileOut=args.fileOut
instance=read_instance(fileIn)
solution=read_solution(fileOut)
print("Objective value: ",solution.objective_value)

errors = check_for_errors(instance, solution)
print("Errors:", errors or "None ✔")

triIdx=0
points=[]
for i in range(len(instance.points_x)):
  points.append(Point(instance.points_x[i],instance.points_y[i]))

for flipSequence in solution.flips:
    flippableTriang=FlippableTriangulation.from_points_edges(
      points,instance.triangulations[triIdx] )
    flipIdx=0
    print(len(flipSequence))
    #print(triIdx)
    triIdx+=1

solFolder="outputFilesFromScript/"
for fileSol in os.listdir(solFolder):
  if False and fileSol.startswith("woc-185-tsplib-b8dd6b77"):
    print(fileSol)
    solution=read_solution(solFolder+fileSol)
    print("Objective value: ",solution.objective_value)

    errors = check_for_errors(instance, solution)
    print("Errors:", errors or "None ✔")
if args.visualize:
  triIdx=0
  for flipSequence in solution.flips:
    flippableTriang=FlippableTriangulation.from_points_edges(
      points,instance.triangulations[triIdx] )
    flipIdx=0
    print(len(flipSequence))
    print(triIdx)
    for parallelFlip in flipSequence:
      draw_flips(
      flippableTriang, show_indices=True, title="Before committing the flip."
      )
      plt.show()
      print(flipIdx)
      flipIdx+=1
      for flip in parallelFlip:
        print(flip)
        flippableTriang.add_flip(flip)
      flippableTriang.commit()
      draw_flips(
      flippableTriang, show_indices=True, title="After committing the flip."
      )
      plt.show()
  
    triIdx+=1

