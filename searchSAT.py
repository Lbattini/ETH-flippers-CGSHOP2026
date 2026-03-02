import os
import subprocess
import json
from json2table import convert
from multiprocessing import Pool
from functools import partial
from itertools import repeat
#CGSHOP pyutils includes
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

def lexSmaller(distA, distB,j):
  for i in range(j):
    if distA[i]<distB[i]:
      return True
    elif distA[i]>distB[i]:
      return False
  return False

def executeAndVerify2(satSolverPath,executablePath,cnfDir,outDir,idxDir,assignmentDir,distance,inFile):
  cnfFile=cnfDir+os.path.splitext(inFile)[0]+".cnf"
  idxFile=idxDir+os.path.splitext(inFile)[0]+".out"
  assignmentFile=assignmentDir+os.path.splitext(inFile)[0]+".out"
  outFile=outDir+os.path.splitext(inFile)[0]+".solution.json"
  args=[executablePath,inDir+inFile,cnfFile,"generateCNF","nonVerbose",idxFile]
  # for distance in distances:
  args.append(str(distance))
  subprocess.run(args,capture_output=True)
  with open(assignmentFile,"w") as asgnFile:
    subprocess.run([satSolverPath,cnfFile,"-q"],stdout=asgnFile)
  line0=""
  with open(assignmentFile,"r") as asgnFile:
    line0=asgnFile.readline()
    line0=line0.strip()
  # print(line0)
  if line0=="s SATISFIABLE":
    # print(line0)
    # print(outFile)
    subprocess.run([executablePath,assignmentFile,outFile,"assignmentToFlipListNU","nonVerbose",idxFile],capture_output=True)
    solution=read_solution(outFile)
    triangulationsPosDist=[]
    flipsPosDist=[]
    allPos=True
    #Can I do this with some variant of [triangulation for triangulation in instance.triangulations]?
    for i in range(len(instance.triangulations)):
      if distances[i]>=0:
        triangulationsPosDist.append(instance.triangulations[i])
        flipsPosDist.append(solution.flips[i])
      else:
        allPos=False

    instancePosDist = CGSHOP2026Instance(
      instance_uid=instance.instance_uid,
      points_x=instance.points_x,
      points_y=instance.points_y,
      triangulations=triangulationsPosDist
    )
    solutionPosDist=CGSHOP2026Solution(
      instance_uid=solution.instance_uid,
      flips=flipsPosDist
    )


    errors = check_for_errors(instancePosDist, solutionPosDist)
    print("Errors:", errors or "None ✔")
    if not errors:
      print("Objective value: ",solutionPosDist.objective_value)

      if allPos:
        with open(executionStatsFile,"a") as exStatsFile:
          exStatsFile.write("Input file: "+inFile+" pfd: "+str(solution.objective_value)+"\n")
      return True
  return False

def executeAndVerify(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inFile, inDir, distances):
  instanceUidAndDistances=os.path.splitext(inFile)[0]

  for distance in distances:
    instanceUidAndDistances+="_"
    instanceUidAndDistances+=str(distance)

  cnfFile=cnfDir+instanceUidAndDistances+".cnf"
  idxFile=idxDir+instanceUidAndDistances+".out"
  assignmentFile=assignmentDir+instanceUidAndDistances+".out"
  outFile=outDir+instanceUidAndDistances+".solution.json"

  args=[executablePath,inDir+inFile,cnfFile,"generateCNFNU","nonVerbose",idxFile]
  for distance in distances:
    args.append(str(distance))
  subprocess.run(args,capture_output=True)
  if satSolverPath=="./build/mallob":
    os.chdir("../mallob/")
    # pathPrefix="../cgshop2026/"
    subprocess.run([satSolverPath,"-mono="+cnfFile,"-q","-t=12","-satsolver=k","-s2f="+assignmentFile],capture_output=True)
    os.chdir("../cgshop2026/")
  else:
    with open(assignmentFile,"w") as asgnFile:
      subprocess.run([satSolverPath,cnfFile,"-q"],stdout=asgnFile)
  line0=""
  with open(assignmentFile,"r") as asgnFile:
    line0=asgnFile.readline()
    line0=line0.strip()
  # print(line0)
  if line0=="s SATISFIABLE":
    # print(line0)
    # print(outFile)
    subprocess.run([executablePath,assignmentFile,outFile,"assignmentToFlipListNU","nonVerbose",idxFile],capture_output=True)
    solution=read_solution(outFile)
    triangulationsPosDist=[]
    flipsPosDist=[]
    allPos=True
    instance=read_instance(inDir+inFile)
    #Can I do this with some variant of [triangulation for triangulation in instance.triangulations]?
    for i in range(len(instance.triangulations)):
      if distances[i]>=0:
        #print(i)
        triangulationsPosDist.append(instance.triangulations[i])
        flipsPosDist.append(solution.flips[i])
      else:
        allPos=False

    instancePosDist = CGSHOP2026Instance(
      instance_uid=instance.instance_uid,
      points_x=instance.points_x,
      points_y=instance.points_y,
      triangulations=triangulationsPosDist
    )
    solutionPosDist=CGSHOP2026Solution(
      instance_uid=solution.instance_uid,
      flips=flipsPosDist
    )


    errors = check_for_errors(instancePosDist, solutionPosDist)
    # print(distances)
    if allPos or errors:
      print("Errors:", errors or "None ✔")
    # points=[]
    # for i in range(len(instance.points_x)):
    #   points.append(Point(instance.points_x[i],instance.points_y[i]))
    # flippableTriang=FlippableTriangulation.from_points_edges(
    #   points,instancePosDist.triangulations[0] )
    #
    # draw_flips(
    #   flippableTriang, show_indices=True, title="Initially."
    # )
    # if errors:
    #   triIdx=0
    #   for flipSequence in solutionPosDist.flips:
    #     flippableTriang=FlippableTriangulation.from_points_edges(
    #       points,instancePosDist.triangulations[triIdx] )
    #     flipIdx=0
    #     print("Number of parallel flips: "+str(len(flipSequence)))
    #     print("Triangulation: "+str(triIdx))
    #     if len(flipSequence)==0:
    #       draw_flips(
    #         flippableTriang, show_indices=True, title="Initially."
    #       )
    #     for parallelFlip in flipSequence:
    #       draw_flips(
    #         flippableTriang, show_indices=True, title="Before committing the flip."
    #       )
    #       plt.show()
    #       print(flipIdx)
    #       flipIdx+=1
    #       for flip in parallelFlip:
    #         #print(flip)
    #         flippableTriang.add_flip(flip)
    #       flippableTriang.commit()
    #       draw_flips(
    #         flippableTriang, show_indices=True, title="After committing the flip."
    #       )
    #       plt.show()
    #
    #     triIdx+=1
    if not errors:
      if allPos:
        print("Objective value: ",solutionPosDist.objective_value)

      if allPos:
        with open(executionStatsFile,"a") as exStatsFile:
          exStatsFile.write("Input file: "+inFile+" pfd: "+str(solution.objective_value)+"\n")
        subprocess.run(["cp",outFile,outDir+instance.instance_uid+".solution.json"])
      else:
        subprocess.run(["rm",idxFile])
        subprocess.run(["rm",cnfFile])
      return True
  else:
    subprocess.run(["rm",idxFile])
    subprocess.run(["rm",cnfFile])

  return False

# def completeSearch(satSolverPath,executablePath,cnfDir,outDir,idxDir,assignmentDir,distances,inFile,i,m,totDist):
#   if i==m-1:
#     distances[i]=totDist-sum(distances)+distances[i]
#     isSat=executeAndVerify(satSolverPath,executablePath,cnfDir,outDir,idxDir,assignmentDir,distances,inFile)
#     if isSat:
#       return True
#     else:
#       return False
#   curSol=False
#   if(sum(distances)<totDist):
#     j=1
#     while(sum(distances)<=totDist):
#       distances[i]=j
#       distances2=[1]*m
#       for k in range(m):
#         distances2[k]=distances[k]
#       curSol=completeSearch(satSolverPath,executablePath,cnfDir,outDir,idxDir,assignmentDir,distances2,inFile,i+1,m,totDist)
#       if curSol:
#         return True
#       j+=1
#       # distances[i]-=1
#   # if not curSol:
#   #   curSol=completeSearch(satSolverPath,executablePath,cnfDir,outDir,idxDir,assignmentDir,distances,inFile,i+1,m,totDist)
#   # distances[i]=1
#   return curSol

def pairwiseDistance(satSolverPath,executablePath,cnfDir,outDir,idxDir,assignmentDir,inDir, inFile, m,t0,t1):
  print(str(t0)+"  "+str(t1))
  distances=[-1]*m
  distances[t0]=0
  distances[t1]=0
  distT0T1=0
  foundSol=False
  while not foundSol:
    print(distances)
    foundSol= executeAndVerify(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inFile, inDir,
                               distances)
    if not foundSol:
      if distT0T1%2==0:
        distances[t0]+=1
      else:
        distances[t1]+=1
      distT0T1+=1
  print("Distance between "+str(t0)+" and "+str(t1)+" is "+str(distT0T1))
  return distT0T1

def pairwiseDistances(satSolverPath,executablePath,cnfDir,outDir,idxDir,assignmentDir,inDir, inFile,m,pairwiseDistFileName):
  if os.path.isfile(pairwiseDistFileName):
    with open(pairwiseDistFileName,"r") as inFile:
      pairwiseDist=[[int(dist) for dist in line.split()] for line in inFile]
    return pairwiseDist

  pairwiseDist=[[]]
  # print("Here")
  for t0 in range(m):
    pairwiseDist.append([])
    for t1 in range(t0):
      pairwiseDist[t0].append(pairwiseDist[t1][t0])
    pairwiseDist[t0].append(0)
    # def wrapperFunction(t1):
    #   return pairwiseDistance(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inDir, inFile, m, t0, t1)
    # res=map(wrapperFunction, range(t0+1,m))
    # print(list(res))
    res=(list(map(partial(pairwiseDistance,satSolverPath,executablePath,cnfDir,outDir,idxDir,assignmentDir,inDir,inFile,m,t0),range(t0+1,m))))
    for dist in res:
      pairwiseDist[t0].append(dist)

    # with Pool(6) as p:
    #   # res=(list(p.map(wrapperFunction, range(t0+1,m))))
    #   # res=(list(p.map(pairwiseDistance, [range(t0+1,m) ])))
    #   # res=(list(map(partial(pairwiseDistance,satSolverPath,executablePath,cnfDir,outDir,idxDir,assignmentDir,inDir, inFile,m,t0),range(t0+1,m))))
    #
    #   res=(list(p.map(partial(pairwiseDistance,satSolverPath,executablePath,cnfDir,outDir,idxDir,assignmentDir,inDir, inFile,m,t0),range(t0+1,m))))
    #   for dist in res:
    #     pairwiseDist[t0].append(dist)


    # print(pairwiseDist)


  with open(pairwiseDistFileName,"w") as outFile:
    for t0 in range(m):
      for t1 in range(m):
        if t1==m-1:
          outFile.write(str(pairwiseDist[t0][t1])+'\n')
        else:
          outFile.write(str(pairwiseDist[t0][t1])+' ')
  return pairwiseDist

# TODO the commented code below is wrong, but the idea is worth thinking about a bit more:
# can you stop searching a bit earlier if the distance combination does not respect pairwise distances?
# def completeSearchPrune(satSolverPath,executablePath,cnfDir,outDir,idxDir,assignmentDir,inFile,totDist,pairwiseDist,distances,m,i):
#   if i==m-1:
#     distances[i]=totDist
#     for j in range(i):
#       distances[i]-=distances[j]
#     if distances[i]<0:
#       return 2
#     for j in range(i):
#       if distances[i]+distances[j]<pairwiseDist[i][j]:
#         return 2
#     print(distances)
#     if executeAndVerify(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, distances, inFile):
#       return 1
#     return 0
#   lowerBoundDistI=0
#   upperBoundDistI=totDist
#   for j in range(i):
#     lowerBoundDistI=max(lowerBoundDistI,pairwiseDist[i][j]-distances[j])
#     upperBoundDistI-=distances[j]
#   if lowerBoundDistI>upperBoundDistI:
#     return 2
#   for distI in range(lowerBoundDistI,upperBoundDistI+1):
#     distances2=[0]*m
#     for j in range(i):
#       distances2[j]=distances[j]
#     distances2[i]=distI
#     tmp=completeSearch(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inFile, totDist, pairwiseDist, distances2, m, i+1)
#     if tmp>=1:
#       return tmp
#
#   return 0

def completeSearch(satSolverPath,executablePath,cnfDir,outDir,idxDir,assignmentDir,inDir,inFile,totDist,pairwiseDist,distances,m,i,minFrom):
  if i>0 and totDist==64 and inFile=="woc-120-tsplib-1ac0c01d.json" and lexSmaller(distances,[10, 7, 13, 12, 7, 7, 8],i):
    #print(distances)
    return False
  if i==m-1:
    distances[i]=totDist
    for j in range(i):
      if distances[j]>=0:
        distances[i]-=distances[j]
    if distances[i]<0:
      return False
    for j in range(i):
      if distances[j]>=0 and distances[i]+distances[j]<pairwiseDist[i][j]:
        return False
    print(distances)
    return executeAndVerify(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inFile, inDir, distances)
  lowerBoundDistI=0
  upperBoundDistI=totDist
  for j in range(i):
    if distances[j]>=0:
      lowerBoundDistI=max(lowerBoundDistI,pairwiseDist[i][j]-distances[j])
      upperBoundDistI-=distances[j]
  if i<=m//2 and upperBoundDistI<minFrom[m//2]:
    return False
  if i>m//2 and upperBoundDistI<minFrom[i]:
    return False
  if lowerBoundDistI<=upperBoundDistI and (i==10 or i==15) and distances[0]>=0:#(i==10 or i==15):
    # print("min from 10: "+str(minFrom10)+" current ub: "+str(upperBoundDistI))
    distances2=[-1]*m
    # distances3=[-1]*m
    for j in range(i):
      distances2[j]=distances[j]
    # print(distances2)
    if not executeAndVerify(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inFile, inDir, distances2): #\
            # or not completeSearch(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inDir, inFile, upperBoundDistI, pairwiseDist, distances3, m, i):
      return False
    print("So far it's sat")
    print(distances2)

  for distI in range(lowerBoundDistI,upperBoundDistI+1):
    distances2=[0]*m
    for j in range(i):
      distances2[j]=distances[j]
    distances2[i]=distI
    if completeSearch(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inDir, inFile, totDist, pairwiseDist, distances2, m, i+1,minFrom):
      return True

  return False

# def completeSearch20(satSolverPath,executablePath,cnfDir,outDir,idxDir,assignmentDir,inDir, inFile,totDist,pairwiseDist):

def binSearchSuffix(satSolverPath,executablePath,cnfDir,outDir,idxDir,assignmentDir,inDir, inFile,pairwiseDist,low,high,m,suffixStart):
  minFrom=[0]*m
  while low<high:
    mid=(low+high)//2
    print("Searching at pfd: "+str(mid))
    foundSol=completeSearch(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inDir, inFile,mid,pairwiseDist,[-1]*m,m,suffixStart,minFrom)
    if foundSol:
      high=mid
    else:
      low=mid+1
  print(low)
  return low
def completeSearch3(satSolverPath,executablePath,cnfDir,outDir,idxDir,assignmentDir,inDir, inFile,totDist,pairwiseDist):
  for d0 in range(0,totDist+1):
    for d1 in range(0,totDist+1-d0):
      d2=totDist-d0-d1
      if d0+d1>=pairwiseDist[0][1] and d0+d2>=pairwiseDist[0][2] and d1+d2>=pairwiseDist[1][2]:
        distances=[d0,d1,d2]
        # print("Here")
        print(distances)
        isSat= executeAndVerify(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inFile, inDir,
                                distances)
        if isSat:
          return True
  return False

def completeSearch4(satSolverPath,executablePath,cnfDir,outDir,idxDir,assignmentDir,inDir, inFile,totDist,pairwiseDist):
  #larger distances increase the number of clauses, so search starting from the distances with the lowest maximum distance
  candidateDist=[]
  for d0 in range(0,totDist+1):
    for d1 in range(max(0,pairwiseDist[0][1]-d0),totDist+1-d0):
      for d2 in range(max(0,pairwiseDist[0][2]-d0,pairwiseDist[1][2]-d1),totDist+1-d0-d1):
        d3=totDist-d0-d1-d2
        if (d0+d3>=pairwiseDist[0][3] and d1+d3>=pairwiseDist[1][3]
              and d2+d3>=pairwiseDist[2][3]):
            candidateDist.append([max(d0,d1,d2,d3),d0,d1,d2,d3])

  candidateDist.sort(key=lambda distances : distances[0])
  for distWithKey in candidateDist:
    distances=distWithKey[1:]
    print(distances)
    isSat= executeAndVerify(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inFile, inDir, distances)
    if isSat:
      return True
  return False
  # for d0 in range(0,totDist+1):
  #   for d1 in range(max(0,pairwiseDist[0][1]-d0),totDist+1-d0):
  #     for d2 in range(max(0,pairwiseDist[0][2]-d0,pairwiseDist[1][2]-d1),totDist+1-d0-d1):
  #       d3=totDist-d0-d1-d2
  #       if (d0+d3>=pairwiseDist[0][3] and d1+d3>=pairwiseDist[1][3]
  #               and d2+d3>=pairwiseDist[2][3]):
  #         distances=[d0,d1,d2,d3]
  #         print(distances)
  #         isSat= executeAndVerify(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inFile, inDir,
  #                                 distances)
  #         if isSat:
  #           return True
  # return False
# def completeSearch5(satSolverPath,executablePath,cnfDir,outDir,idxDir,assignmentDir,inFile,totDist,pairwiseDist):
#   for d0 in range(0,totDist+1):
#     for d1 in range(max(0,pairwiseDist[0][1]-d0),totDist+1-d0):
#       for d2 in range(max(0,pairwiseDist[0][2]-d0,pairwiseDist[1][2]-d1),totDist+1-d0-d1):
#         for d3 in range(max(0,pairwiseDist[0][3]-d0,pairwiseDist[1][3]-d1,pairwiseDist[2][3]-d2),totDist+1-d0-d1-d2):
#           d4=totDist-d0-d1-d2-d3
#           if (d0+d4>=pairwiseDist[0][4] and d1+d4>=pairwiseDist[1][4]
#                   and d2+d4>=pairwiseDist[2][4] and d3+d4>=pairwiseDist[3][4]):
#             distances=[d0,d1,d2,d3,d4]
#             print(distances)
#             isSat=executeAndVerify(satSolverPath,executablePath,cnfDir, outDir, idxDir, assignmentDir, distances, inFile)
#             if isSat:
#               return True
#   return False
def completeSearch5(satSolverPath,executablePath,cnfDir,outDir,idxDir,assignmentDir,inDir, inFile,totDist,pairwiseDist):
  #larger distances increase the number of clauses, so search starting from the distances with the lowest maximum distance
  candidateDist=[]
  for d0 in range(0,totDist+1):
    for d1 in range(max(0,pairwiseDist[0][1]-d0),totDist+1-d0):
      for d2 in range(max(0,pairwiseDist[0][2]-d0,pairwiseDist[1][2]-d1),totDist+1-d0-d1):
        for d3 in range(max(0,pairwiseDist[0][3]-d0,pairwiseDist[1][3]-d1,pairwiseDist[2][3]-d2),totDist+1-d0-d1-d2):
            d4=totDist-d0-d1-d2-d3
            if (d0+d4>=pairwiseDist[0][4] and d1+d4>=pairwiseDist[1][4]
                    and d2+d4>=pairwiseDist[2][4] and d3+d4>=pairwiseDist[3][4]):
              candidateDist.append([max(d0,d1,d2,d3,d4),d0,d1,d2,d3,d4])

  candidateDist.sort(key=lambda distances : distances[0])
  for distWithKey in candidateDist:
    distances=distWithKey[1:]
    print(distances)
    isSat= executeAndVerify(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inFile, inDir, distances)
    if isSat:
      return True
  return False

def completeSearch6(satSolverPath,executablePath,cnfDir,outDir,idxDir,assignmentDir,inDir, inFile,totDist,pairwiseDist):
  #larger distances increase the number of clauses, so search starting from the distances with the lowest maximum distance
  candidateDist=[]
  for d0 in range(0,totDist+1):
    for d1 in range(max(0,pairwiseDist[0][1]-d0),totDist+1-d0):
      for d2 in range(max(0,pairwiseDist[0][2]-d0,pairwiseDist[1][2]-d1),totDist+1-d0-d1):
        for d3 in range(max(0,pairwiseDist[0][3]-d0,pairwiseDist[1][3]-d1,pairwiseDist[2][3]-d2),totDist+1-d0-d1-d2):
          for d4 in range(max(0,pairwiseDist[0][4]-d0,pairwiseDist[1][4]-d1,pairwiseDist[2][4]-d2,pairwiseDist[3][4]-d3),totDist+1-d0-d1-d2):
            d5=totDist-d0-d1-d2-d3-d4
            if (d0+d5>=pairwiseDist[0][5] and d1+d5>=pairwiseDist[1][5]
                    and d2+d5>=pairwiseDist[2][5] and d3+d5>=pairwiseDist[3][5] and d4+d5>=pairwiseDist[4][5]):
              candidateDist.append([max(d0,d1,d2,d3,d4,d5),d0,d1,d2,d3,d4,d5])

  candidateDist.sort(key=lambda distances : distances[0])
  j=0
  # for distWithKey in candidateDist:
  # for i in range(0,len(candidateDist)-len(candidateDist)%8,8):
  #   j=i+8
  #   distancesL=[candidateDist[i][1:],candidateDist[i+1][1:],candidateDist[i+2][1:],candidateDist[i+3][1:],
  #               candidateDist[i+4][1:],candidateDist[i+5][1:],candidateDist[i+6][1:],candidateDist[i+7][1:]]
  #   print(str(i)+"/"+str(len(candidateDist)-1))
  #   print(distancesL)
  #   # if not (inFile=="woc-155-tsplib-16bd1a72.json" and totDist==41 and i<1872):
  #   with Pool(8) as p:
  #
  #     partialFunc=partial(executeAndVerify,satSolverPath,
  #                         executablePath,cnfDir, outDir, idxDir, assignmentDir,
  #                         inFile, inDir)
  #     # partialFunc(distances=distances)
  #     res=(list(p.map(partialFunc,distancesL)))
  #     for isSat in res:
  #       if isSat:
  #         return True
      # isSat=executeAndVerify(satSolverPath,executablePath,cnfDir, outDir, idxDir, assignmentDir, inFile, inDir, distances)
      # if isSat:
      #   return True
    # i+=1
  for distWithKey in candidateDist:
    distances=distWithKey[1:]
    # print(str(j)+"/"+str(len(candidateDist)-1))
    print(distances)
    isSat=executeAndVerify(satSolverPath,executablePath,cnfDir, outDir, idxDir, assignmentDir, inFile, inDir, distances)
    if isSat:
      return True

  # for d0 in range(0,totDist+1):
  #   for d1 in range(max(0,pairwiseDist[0][1]-d0),totDist+1-d0):
  #     for d2 in range(max(0,pairwiseDist[0][2]-d0,pairwiseDist[1][2]-d1),totDist+1-d0-d1):
  #       for d3 in range(max(0,pairwiseDist[0][3]-d0,pairwiseDist[1][3]-d1,pairwiseDist[2][3]-d2),totDist+1-d0-d1-d2):
  #         for d4 in range(max(0,pairwiseDist[0][4]-d0,pairwiseDist[1][4]-d1,pairwiseDist[2][4]-d2,pairwiseDist[3][4]-d3),totDist+1-d0-d1-d2):
  #           d5=totDist-d0-d1-d2-d3-d4
  #           if (d0+d5>=pairwiseDist[0][5] and d1+d5>=pairwiseDist[1][5]
  #                   and d2+d5>=pairwiseDist[2][5] and d3+d5>=pairwiseDist[3][5] and d4+d5>=pairwiseDist[4][5]):
  #             distances=[d0,d1,d2,d3,d4,d5]
  #             print(distances)
  #             isSat=executeAndVerify(satSolverPath,executablePath,cnfDir, outDir, idxDir, assignmentDir, distances, inFile)
  #             if isSat:
  #               return True
  return False

def readCurrentBest(lowerBoundDict,upperBoundDict,initialPathDict,optValuesFileName):
  # dataFromJson=[]
  with open(optValuesFileName,mode="r", encoding="utf-8") as optValuesFile:
  # with open("./scriptOutputs/currentBest/instancesStatsTestShort.json",mode="r", encoding="utf-8") as optValuesFile:
    dataFromJson=json.load(optValuesFile)

  # print(dataFromJson['instances'])
  for instanceData in dataFromJson['instances']:

    instance_uid=instanceData["instance_uid"]
    if (not instance_uid in lowerBoundDict) or lowerBoundDict[instance_uid]<instanceData["lower_bound"]:
      lowerBoundDict[instance_uid]=instanceData["lower_bound"]

    if (not instance_uid in upperBoundDict) or upperBoundDict[instance_uid]>instanceData["upper_bound"]:
      upperBoundDict[instance_uid]=instanceData["upper_bound"]
      initialPathDict[instance_uid]=instanceData["initial_path"]
    # print(instanceData["instance_uid"])

def copyCurrentBest(optValuesFileName,currentBestDir):
  with open(optValuesFileName,mode="r", encoding="utf-8") as optValuesFile:
    dataFromJson=json.load(optValuesFile)

  for instanceData in dataFromJson['instances']:
    instance_uid=instanceData["instance_uid"]
    subprocess.run(["cp",instanceData["initial_path"],currentBestDir+instance_uid+".solution.json"])
def writeCurrentBest(lowerBoundDict,upperBoundDict,initialPathDict,optValuesFileName):
  readCurrentBest(lowerBoundDict, upperBoundDict, initialPathDict, optValuesFileName)
  data=[]
  for instance_uid in upperBoundDict:
    data.append({"instance_uid":instance_uid,
                 "lower_bound":lowerBoundDict[instance_uid],
                 "upper_bound":upperBoundDict[instance_uid],
                 "initial_path":initialPathDict[instance_uid]})

  with open(optValuesFileName,mode="w", encoding="utf-8") as optValuesFile:
    optValuesFile.write('{"instances":')
  with open(optValuesFileName,mode="a", encoding="utf-8") as optValuesFile:
    json.dump(data,optValuesFile,indent=2)
  with open(optValuesFileName,mode="a", encoding="utf-8") as optValuesFile:
    optValuesFile.write('}')

  with open(optValuesFileName,mode="r", encoding="utf-8") as optValuesFile:
    dataFromJson=json.load(optValuesFile)

  html=convert(dataFromJson)
  with open('./scriptOutputs/currentBest/instancesStats.html',mode="w", encoding="utf-8") as tableFile:
    tableFile.write(html)

def fillOptValues(solutionsDir):
  currentBestDir="./scriptOutputs/currentBest/"
  optValuesFileName="./scriptOutputs/currentBest/instancesStats.json"
  instanceDir="../benchmark_instances_rev1/benchmark_instances/"
  # lowerBoundDict=dict()
  #store upper bound on the pfd of each instance (indexed by uid)
  lowerBoundDict=dict()
  upperBoundDict=dict()
  #this is to store the folder where the current best solution came from (just to store some extra information, it's not really necessary)
  initialPathDict=dict()
  readCurrentBest(lowerBoundDict, upperBoundDict, initialPathDict, optValuesFileName)
  for inFile in os.listdir(solutionsDir):
    if inFile.endswith("json") and inFile!="instancesStats.json":
      # try:
        solution=read_solution(solutionsDir+inFile)
        instance_uid=solution.instance_uid
        print(inFile+" "+instance_uid)

        if inFile.startswith("example"):
          instance=read_instance("../example_instances_cgshop2026/examples/"+instance_uid+".json")
        else:
          instance=read_instance(instanceDir+instance_uid+".json")
        errors = check_for_errors(instance, solution)
        print("Errors:", errors or "None ✔")
        if not errors:
          if (not (instance_uid in upperBoundDict) or upperBoundDict[instance_uid]>solution.objective_value):
            if not (instance_uid in upperBoundDict):
              lowerBoundDict[instance_uid]=0
            upperBoundDict[instance_uid]=solution.objective_value
            initialPathDict[instance_uid]=solutionsDir+inFile
            subprocess.run(["cp",solutionsDir+inFile,currentBestDir+instance_uid+".solution.json"])
            # print("Objective value: ",solution.objective_value)
      # except:
      #   print(inFile+" has some syntax error")

  writeCurrentBest(lowerBoundDict,upperBoundDict,initialPathDict,optValuesFileName)
  # print("Reached this point")
  # data=[]
  # for instance_uid in upperBoundDict:
  #   data.append({"instance_uid":instance_uid,
  #         "lower_bound":0,
  #         "upper_bound":upperBoundDict[instance_uid],
  #         "initial_path":initialPathDict[instance_uid]})
  #
  # with open(optValuesFileName,mode="w", encoding="utf-8") as optValuesFile:
  #   optValuesFile.write('{"instances":')
  # with open(optValuesFileName,mode="a", encoding="utf-8") as optValuesFile:
  #   json.dump(data,optValuesFile)
  # with open(optValuesFileName,mode="a", encoding="utf-8") as optValuesFile:
  #   optValuesFile.write('}')
    # if inFile.startswith("values") and inFile != "values":
      # with open(executionStatsFile,"r") as exStatsFile:
      #   line = exStatsFile.readline()
      #   while line:
      #     instanceUid=line.strip()
      #     instanceUid=instanceUid.removeprefix("Computing a solution for file ../benchmark_instances_rev1/benchmark_instances/")
      #     instanceUid=instanceUid.removesuffix(".json")
      #     line=exStatsFile.readline()
      #     validFlipSequence=line.strip()
      #     line=exStatsFile.readline()

      # subprocess.run(["cp",])


if __name__ == '__main__':
  optValuesFileName="./scriptOutputs/currentBest/instancesStats.json"
  solutionsDir="./scriptOutputs/satSolutions_114/"
  currentBestDir="./scriptOutputs/currentBest/"
  # fillOptValues(solutionsDir)
  # copyCurrentBest(optValuesFileName,"./scriptOutputs/currentBest/")
  lowerBoundDict=dict()
  upperBoundDict=dict()
  initialPathDict=dict()
  readCurrentBest(lowerBoundDict,upperBoundDict,initialPathDict,optValuesFileName)

  inDir="../benchmark_instances_rev1/benchmark_instances/"
  pairwiseDistDir="./scriptOutputs/pairwiseDistDir1/"
  dirCnt=0
  cnfDir="/fastStorage/scriptOutputs/cnfDir_"+str(dirCnt)+"/"
  while os.path.isdir(cnfDir):
    dirCnt+=1
    cnfDir="/fastStorage/scriptOutputs/cnfDir_"+str(dirCnt)+"/"

  outDir="/fastStorage/scriptOutputs/satSolutions_"+str(dirCnt)+"/"
  idxDir="/fastStorage/scriptOutputs/idxDir_"+str(dirCnt)+"/"
  assignmentDir="/fastStorage/scriptOutputs/assignments_"+str(dirCnt)+"/"
  os.mkdir(cnfDir)
  os.mkdir(idxDir)
  os.mkdir(outDir)
  os.mkdir(assignmentDir)
  executablePath="cmake-build-release/centralTriangulationClion"
  #satSolverPath="../kissat/build/kissat"
  satSolverPath="./build/mallob"
  cnt=0
  executionStatsFile=outDir+"values.out"

  #this is to avoid executing too large files
  maxN=10
  maxM=10
  cntSolved=0
  # pairwiseDistFileName=pairwiseDistDir+"random_instance_110_15_3_test3"+".out"
  # pairwiseDistances(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, "random_instance_110_15_3.json", 3,pairwiseDistFileName)
  #
  # executeAndVerify(satSolverPath,executablePath,cnfDir, outDir, idxDir, assignmentDir, [4,3], "random_instance_482_80_2.json")
  # pairwiseDistFileName=pairwiseDistDir+"random_instance_840_15_3_test6"+".out"
  # pairwiseDist=pairwiseDistances(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, "random_instance_840_15_3.json", 3,pairwiseDistFileName)
  #
  # completeSearch3(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, "random_instance_840_15_3.json",8,pairwiseDist)
  for inFile in os.listdir(inDir):
    if (not inFile.startswith("rirs")) and (inFile.endswith(".json")):
    #if inFile.endswith(".json"):
      instance=read_instance(inDir+inFile)
      n=len(instance.points_x)
      m=len(instance.triangulations)
      instance_uid=instance.instance_uid
      #and m<=maxM:
      # if n<=60 and instance_uid in lowerBoundDict and lowerBoundDict[instance_uid]==upperBoundDict[instance_uid]:
      #   print("Solving: "+instance_uid)
      #   if m==2:
      #     pfd=lowerBoundDict[instance_uid]
      #     print("Searching at pfd: "+str(pfd))
      #     distances=[pfd//2,pfd-pfd//2]
      #     foundSol=executeAndVerify(satSolverPath,executablePath,cnfDir, outDir, idxDir, assignmentDir, distances, inFile)
      #
      #     if foundSol:
      #       print("Fixed solution file")
      #       upperBoundDict[instance_uid]=pfd
      #       initialPathDict[instance_uid]=outDir+instance_uid+".solution.json"
      #       subprocess.run(["cp",initialPathDict[instance_uid],currentBestDir+instance_uid+".solution.json"])
      #   else:
      #     pairwiseDistFileName=pairwiseDistDir+instance_uid+".out"
      #     pairwiseDist=pairwiseDistances(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inFile, m,pairwiseDistFileName)
      #     print(pairwiseDist)
      #     pfd=lowerBoundDict[instance_uid]
      #
      #     print("Searching at pfd: "+str(pfd))
      #
      #     foundSol=completeSearch(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inFile,pfd,pairwiseDist,[0]*m,m,0)
      #
      #     if foundSol:
      #       print("Fixed solution file")
      #       upperBoundDict[instance_uid]=pfd
      #       initialPathDict[instance_uid]=outDir+instance_uid+".solution.json"
      #       subprocess.run(["cp",initialPathDict[instance_uid],currentBestDir+instance_uid+".solution.json"])
      #   writeCurrentBest(lowerBoundDict, upperBoundDict, initialPathDict, optValuesFileName)
      if n<=maxN and m<=maxM: #and (instance_uid in lowerBoundDict and lowerBoundDict[instance_uid]==upperBoundDict[instance_uid]):#(n<=320 and m==2) or (n<=15 and m<=10) or (n<=160 and m<=3) or (n<=60 and m<=10):# or upperBoundDict[instance_uid]<20:
        #lowerBoundDict[instance_uid]=0
        if (not instance_uid in lowerBoundDict) or lowerBoundDict[instance_uid]<upperBoundDict[instance_uid]:
          print("Solving: "+instance_uid)
          if m==2:
            pfd=0
            distances=[0]*m
            foundSol=False
            low=lowerBoundDict[instance_uid]
            high=upperBoundDict[instance_uid]
            while low<high:
              mid=(low+high)//2
              print("Searching at pfd: "+str(mid))
              distances[0]=mid//2
              distances[1]=mid-distances[0]
              print(distances)
              foundSol= executeAndVerify(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inFile,
                                         inDir, distances)
              if foundSol:
                high=mid
              else:
                low=mid+1

            if not foundSol:
              print("Searching at pfd: "+str(low))
              distances[0]=low//2
              distances[1]=low-distances[0]
              foundSol= executeAndVerify(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inFile,
                                         inDir, distances)

            lowerBoundDict[instance_uid]=low
            if low<upperBoundDict[instance_uid]:
              print("Improved solution from pfd "+str(upperBoundDict[instance_uid])+" to pfd "+str(low))
              upperBoundDict[instance_uid]=low
              initialPathDict[instance_uid]=outDir+instance_uid+".solution.json"
              subprocess.run(["cp",initialPathDict[instance_uid],currentBestDir+instance_uid+".solution.json"])
          else:
            pairwiseDistFileName=pairwiseDistDir+instance_uid+".out"
            pairwiseDist=pairwiseDistances(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inDir, inFile, m,pairwiseDistFileName)
            print(pairwiseDist)
            minFrom=[0, 0, 0, 27, 16, 9, 0]#[0]*m
            # if m==20:
            for suffixStart in range(m//2,0):
              low1=lowerBoundDict[instance_uid]
              high1=upperBoundDict[instance_uid]
              if suffixStart>m//2:
                high1=minFrom[suffixStart-1]
              minFrom[suffixStart]=binSearchSuffix(satSolverPath,executablePath,cnfDir,outDir,idxDir,assignmentDir,inDir,inFile,pairwiseDist,low1,high1,m,suffixStart)

            print("MIN FROM: "+str(minFrom))
            low=0 #max(minFrom[m//2],lowerBoundDict[instance_uid])
            high=upperBoundDict[instance_uid]
            while low<high:
              mid=(low+high)//2
              print("Searching at pfd: "+str(mid)+"["+str(low)+","+str(high)+"]")
              foundSol=False
              if m==3:
                foundSol=completeSearch3(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inDir, inFile,mid,pairwiseDist)
              elif m==4:
                foundSol=completeSearch4(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inDir, inFile,mid,pairwiseDist)
              elif m==5:
                foundSol=completeSearch5(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inDir, inFile,mid,pairwiseDist)
              elif m==6:
                foundSol=completeSearch6(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inDir, inFile,mid,pairwiseDist)
              # elif m==20:
              #   foundSol=completeSearch(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inDir, inFile,mid,pairwiseDist,[0]*m,m,0,minFrom10=minFrom10)
              else:
                foundSol=completeSearch(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inDir, inFile,mid,pairwiseDist,[0]*m,m,0,minFrom)
              if foundSol:
                high=mid
              else:
                low=mid+1
            lowerBoundDict[instance_uid]=low
            # foundSol=completeSearch(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inDir, inFile,low,pairwiseDist,[0]*m,m,0)
            if low<upperBoundDict[instance_uid]:
              print("Improved solution from pfd "+str(upperBoundDict[instance_uid])+" to pfd "+str(low))
              upperBoundDict[instance_uid]=low
              instanceUidAndDistances=os.path.splitext(inFile)[0]
              # for distance in distances:
              #   instanceUidAndDistances+="_"
              #   instanceUidAndDistances+=str(distance)
              initialPathDict[instance_uid]=outDir+instance_uid+".solution.json"
              subprocess.run(["cp",initialPathDict[instance_uid],currentBestDir+instance_uid+".solution.json"])

          writeCurrentBest(lowerBoundDict, upperBoundDict, initialPathDict, optValuesFileName)

          # distances=[-1]*m
            # distances[0]=0
            # distances[1]=0
            # foundSol=False
            # while not foundSol:
            #   foundSol=executeAndVerify(satSolverPath,executablePath,cnfDir, outDir, idxDir, assignmentDir, distances, inFile)
            #   distances[0]+=1
            # tmpDist=distances[0]
            # distances[0]=tmpDist/2
            # distances[1]=tmpDist-distances[0]
            # for i in range(2,m):
            #   distances[i]=0
            #   foundSol=False
            #   while not foundSol:
            #     foundSol=executeAndVerify(satSolverPath,executablePath,cnfDir, outDir, idxDir, assignmentDir, distances, inFile)
            #     distances[i]+=1

            # uniformDistance=0
            # foundSol=False
            # i=0
            # distances=[0]*m
            # while not foundSol:
            #   # distances=[uniformDistance]*m
            #   foundSol=executeAndVerify(satSolverPath,executablePath,cnfDir, outDir, idxDir, assignmentDir, distances, inFile)
            #   # uniformDistance+=1
            #   distances[i%m]+=1
            #   i+=1

      if lowerBoundDict[instance_uid]==upperBoundDict[instance_uid]:
        cntSolved+=1

        # args=[executablePath,cnfFile,
        #print(inFile)
        # if True: #not os.path.isfile(outFile):
          #print(outFile)
          # cnt+=1

          #if(cnt>10):
          #  break
          #with open(executionStatsFile,"a") as exStatsFile:
          #  subprocess.run([executablePath,inFolder+inFile,outFile,strategy,"nonVerbose"],stdout=exStatsFile)
  # writeCurrentBest(lowerBoundDict, upperBoundDict, initialPathDict, optValuesFileName)
  print("Solved "+str(cntSolved)+" instances")
