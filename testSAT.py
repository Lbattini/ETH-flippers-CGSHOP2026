import os
import subprocess
import json
from json2table import convert
from multiprocessing import Pool
from functools import partial
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

def executeAndVerify(satSolverPath,executablePath,cnfDir,outDir,idxDir,assignmentDir,distances,inFile):
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
        print(i)
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

def pairwiseDistance(satSolverPath,executablePath,cnfDir,outDir,idxDir,assignmentDir,inFile,m,t0,t1):
  print(str(t0)+"  "+str(t1))
  distances=[-1]*m
  distances[t0]=0
  distances[t1]=0
  distT0T1=0
  foundSol=False
  while not foundSol:
    print(distances)
    foundSol=executeAndVerify(satSolverPath,executablePath,cnfDir,outDir,idxDir,assignmentDir,distances,inFile)
    if not foundSol:
      if distT0T1%2==0:
        distances[t0]+=1
      else:
        distances[t1]+=1
      distT0T1+=1
  print("Distance between "+str(t0)+" and "+str(t1)+" is "+str(distT0T1))
  return distT0T1

def pairwiseDistances(satSolverPath,executablePath,cnfDir,outDir,idxDir,assignmentDir,inFile,m,pairwiseDistFileName):
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
    #   return pairwiseDistance(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inFile, m, t0, t1)
    # res=map(wrapperFunction, range(t0+1,m))
    # print(list(res))
    res=(list(map(partial(pairwiseDistance,satSolverPath,executablePath,cnfDir,outDir,idxDir,assignmentDir,inFile,m,t0),range(t0+1,m))))
    for dist in res:
      pairwiseDist[t0].append(dist)

    # with Pool(6) as p:
    #   # res=(list(p.map(wrapperFunction, range(t0+1,m))))
    #   # res=(list(p.map(pairwiseDistance, [range(t0+1,m) ])))
    #   # res=(list(map(partial(pairwiseDistance,satSolverPath,executablePath,cnfDir,outDir,idxDir,assignmentDir,inFile,m,t0),range(t0+1,m))))
    #
    #   res=(list(p.map(partial(pairwiseDistance,satSolverPath,executablePath,cnfDir,outDir,idxDir,assignmentDir,inFile,m,t0),range(t0+1,m))))
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



def completeSearch(satSolverPath,executablePath,cnfDir,outDir,idxDir,assignmentDir,inFile,totDist,pairwiseDist,distances,m,i):
  if i==m-1:
    distances[i]=totDist
    for j in range(i):
      distances[i]-=distances[j]
    if distances[i]<0:
      return False
    for j in range(i):
      if distances[i]+distances[j]<pairwiseDist[i][j]:
        return False
    print(distances)
    return executeAndVerify(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, distances, inFile)
  lowerBoundDistI=0
  upperBoundDistI=totDist
  for j in range(i):
    lowerBoundDistI=max(lowerBoundDistI,pairwiseDist[i][j]-distances[j])
    upperBoundDistI-=distances[j]
  for distI in range(lowerBoundDistI,upperBoundDistI+1):
    distances2=[0]*m
    for j in range(i):
      distances2[j]=distances[j]
    distances2[i]=distI
    if completeSearch(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inFile, totDist, pairwiseDist, distances2, m, i+1):
      return True

  return False


def completeSearch3(satSolverPath,executablePath,cnfDir,outDir,idxDir,assignmentDir,inFile,totDist,pairwiseDist):
  for d0 in range(0,totDist+1):
    for d1 in range(0,totDist+1-d0):
      d2=totDist-d0-d1
      if d0+d1>=pairwiseDist[0][1] and d0+d2>=pairwiseDist[0][2] and d1+d2>=pairwiseDist[1][2]:
        distances=[d0,d1,d2]
        # print("Here")
        print(distances)
        isSat=executeAndVerify(satSolverPath,executablePath,cnfDir, outDir, idxDir, assignmentDir, distances, inFile)
        if isSat:
          return True
  return False

def completeSearch4(satSolverPath,executablePath,cnfDir,outDir,idxDir,assignmentDir,inFile,totDist,pairwiseDist):
  for d0 in range(0,totDist+1):
    for d1 in range(max(0,pairwiseDist[0][1]-d0),totDist+1-d0):
      for d2 in range(max(0,pairwiseDist[0][2]-d0,pairwiseDist[1][2]-d1),totDist+1-d0-d1):
        d3=totDist-d0-d1-d2
        if (d0+d3>=pairwiseDist[0][3] and d1+d3>=pairwiseDist[1][3]
                and d2+d3>=pairwiseDist[2][3]):
          distances=[d0,d1,d2,d3]
          print(distances)
          isSat=executeAndVerify(satSolverPath,executablePath,cnfDir, outDir, idxDir, assignmentDir, distances, inFile)
          if isSat:
            return True
  return False

def completeSearch5(satSolverPath,executablePath,cnfDir,outDir,idxDir,assignmentDir,inFile,totDist,pairwiseDist):
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
    isSat=executeAndVerify(satSolverPath,executablePath,cnfDir, outDir, idxDir, assignmentDir, distances, inFile)
    if isSat:
      return True
  return False

def completeSearch6(satSolverPath,executablePath,cnfDir,outDir,idxDir,assignmentDir,inFile,totDist,pairwiseDist):
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
  for distWithKey in candidateDist:
    distances=distWithKey[1:]
    print(distances)
    isSat=executeAndVerify(satSolverPath,executablePath,cnfDir, outDir, idxDir, assignmentDir, distances, inFile)
    if isSat:
      return True

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


if __name__ == '__main__':
  optValuesFileName="./scriptOutputs/currentBest/instancesStats.json"
  currentBestDir="./scriptOutputs/currentBest/"

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
  # satSolverPath="../kissat/build/kissat"
  satSolverPath="./build/mallob"
  cnt=0
  executionStatsFile=outDir+"values.out"

  #this is to avoid executing too large files
  maxN=80
  maxM=5
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

      if n<=maxN and m<=maxM and (instance_uid in lowerBoundDict and lowerBoundDict[instance_uid]==upperBoundDict[instance_uid]):
        if (not instance_uid in lowerBoundDict) or lowerBoundDict[instance_uid]<=upperBoundDict[instance_uid]:
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
              foundSol=executeAndVerify(satSolverPath,executablePath,cnfDir, outDir, idxDir, assignmentDir, distances, inFile)
              if foundSol:
                high=mid
              else:
                low=mid+1

            if not foundSol:
              print("Searching at pfd: "+str(low))
              distances[0]=low//2
              distances[1]=low-distances[0]
              foundSol=executeAndVerify(satSolverPath,executablePath,cnfDir, outDir, idxDir, assignmentDir, distances, inFile)

            lowerBoundDict[instance_uid]=low
            if low<upperBoundDict[instance_uid]:
              print("Improved solution from pfd "+str(upperBoundDict[instance_uid])+" to pfd "+str(low))
              upperBoundDict[instance_uid]=low
              initialPathDict[instance_uid]=outDir+instance_uid+".solution.json"
              subprocess.run(["cp",initialPathDict[instance_uid],currentBestDir+instance_uid+".solution.json"])
          else:
            pairwiseDistFileName=pairwiseDistDir+instance_uid+".out"
            pairwiseDist=pairwiseDistances(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inFile, m,pairwiseDistFileName)
            print(pairwiseDist)
            low=lowerBoundDict[instance_uid]
            high=upperBoundDict[instance_uid]
            while low<high:
              mid=(low+high)//2
              print("Searching at pfd: "+str(mid))
              foundSol=False
              if m==3:
                foundSol=completeSearch3(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inFile,mid,pairwiseDist)
              elif m==4:
                foundSol=completeSearch4(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inFile,mid,pairwiseDist)
              elif m==5:
                foundSol=completeSearch5(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inFile,mid,pairwiseDist)
              elif m==6:
                foundSol=completeSearch6(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inFile,mid,pairwiseDist)
              else:
                foundSol=completeSearch(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inFile,mid,pairwiseDist,[0]*m,m,0)
              if foundSol:
                high=mid
              else:
                low=mid+1
            lowerBoundDict[instance_uid]=low
            foundSol=completeSearch(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inFile,low,pairwiseDist,[0]*m,m,0)
            if low<upperBoundDict[instance_uid]:
              print("Improved solution from pfd "+str(upperBoundDict[instance_uid])+" to pfd "+str(low))
              upperBoundDict[instance_uid]=low

              initialPathDict[instance_uid]=outDir+instance_uid+".solution.json"
              # subprocess.run(["cp",initialPathDict[instance_uid],currentBestDir+instance_uid+".solution.json"])

          # writeCurrentBest(lowerBoundDict, upperBoundDict, initialPathDict, optValuesFileName)


      if lowerBoundDict[instance_uid]==upperBoundDict[instance_uid]:
        cntSolved+=1

  print("Solved "+str(cntSolved)+" instances")
