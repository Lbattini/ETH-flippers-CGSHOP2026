#This is the script used to run heuristics on multiple input files (at the moment they require a single execution of the C executable per instance)
import os
import subprocess
import json
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

def custom_sort_key(s):
    return [(chr(0) if c == '_' else c) for c in s]
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

  # html=convert(dataFromJson)
  # with open('./scriptOutputs/currentBest/instancesStats.html',mode="w", encoding="utf-8") as tableFile:
  #   tableFile.write(html)

optValuesFileName="./scriptOutputs/currentBest/instancesStats.json"
currentBestDir="./scriptOutputs/currentBest/"

baseDir="/fastStorage"
dirCnt=0
outDir=baseDir+"/scriptOutputs/heuristicSolutions_"+str(dirCnt)+"/"

while os.path.isdir(outDir):
  dirCnt+=1
  outDir=baseDir+"/scriptOutputs/heuristicSolutions_"+str(dirCnt)+"/"

os.mkdir(outDir)

lowerBoundDict=dict()
upperBoundDict=dict()
initialPathDict=dict()
readCurrentBest(lowerBoundDict,upperBoundDict,initialPathDict,optValuesFileName)
inDir="../benchmark_instances_rev1/benchmark_instances/"

executablePath="./centralTriangulation"
i=0
strategy="reduceItsCenters"
cnt=0
executionStatsFile=outDir+"values.out"
# executionStatsFile=outDir+"values_"+str(i)+".out"
# while os.path.isfile(executionStatsFile):
#   i+=1
#   executionStatsFile=out+"values_"+str(i)+".out"
#executionStatsFiles=outFolder+"values"+str(i)+a
filesToSolveWithKey=[]
maxN=200
minN=10
for inFile in sorted(os.listdir(inDir),key=custom_sort_key):
  if(inFile.endswith(".json")):
  #if (inFile.startswith("rirs")) and (inFile.endswith(".json")):
  #if inFile.endswith(".json"): inFile=="woc-185-tsplib-b8dd6b77.json"
    with open(inDir+inFile, 'r') as file:
      data = json.load(file)
    n=len(data["points_x"])
    m=len(data["triangulations"])
    instance_uid=data["instance_uid"] 
    #print(instance_uid)
    #print(upperBoundDict[instance_uid])
    if ((not instance_uid in lowerBoundDict) or lowerBoundDict[instance_uid]<upperBoundDict[instance_uid]):
      #print(upperBoundDict[instance_uid])
      #print(instance_uid+" n: "+str(n)+" m: "+str(m)+"lower bound: "+str(lowerBoundDict[instance_uid])+" upper bound: "+str(upperBoundDict[instance_uid])+" "+str(upperBoundDict[instance_uid]/m))
      if n<=maxN: #(n>=minN and m==20) or (n>=3000 and m==50) or (n>2500 and m==75): #n==12500 and m==200: #n<=maxN and m==20
        filesToSolveWithKey.append([n*m,inFile])
filesToSolveWithKey.sort(key=lambda fileWithKey : fileWithKey[0])
cnt=0
for fileWithKey in filesToSolveWithKey[cnt:]:
  inFile=fileWithKey[1]
  # print(cnt)
  cnt+=1
  with open(inDir+inFile, 'r') as file:
    data = json.load(file)
  n=len(data["points_x"])
  m=len(data["triangulations"])
  instance_uid=data["instance_uid"] #os.path.splitext(inFile)[0]
  # outFile=outDir+instance_uid+"_out"+str(i)+".solution.json"
  outFile=outDir+instance_uid+".solution.json"
  print(instance_uid+" n: "+str(n)+" m: "+str(m)+"lower bound: "+str(lowerBoundDict[instance_uid])+" upper bound: "+str(upperBoundDict[instance_uid])+" "+str(upperBoundDict[instance_uid]/m))

  #print(inFile) rirs-12500-100-0c057cab"
  if (not instance_uid in lowerBoundDict) or lowerBoundDict[instance_uid]<upperBoundDict[instance_uid]: #(not instance_uid in lowerBoundDict) and#or instance_uid in ["rirs-1000-20-ee506ff1","rirs-2000-20-fc840c19","rirs-1000-75-50e9d715"]: #or lowerBoundDict[instance_uid]<upperBoundDict[instance_uid]: #not os.path.isfile(outFile):
    #print(outFile)
    print("Solving: "+instance_uid)
    #cnt+=1
    #if(cnt>10):
    #  break
    with open(executionStatsFile,"a") as exStatsFile:
      subprocess.run([executablePath,inDir+inFile,outFile,strategy,"verbose"],stdout=exStatsFile)
    instance=read_instance(inDir+inFile)
    solution=read_solution(outFile)

    if n<=500:
      errors = check_for_errors(instance, solution)
      print("Errors:", errors or "None ✔")
    if True: #not errors: skipped error checking because it's very slow for rirs
      print("Objective value: ",solution.objective_value)
      if (not instance_uid in upperBoundDict) or solution.objective_value<upperBoundDict[instance_uid]:
        if instance_uid in upperBoundDict:
          print("Improved solution from pfd "+str(upperBoundDict[instance_uid])+" to pfd "+str(solution.objective_value))
        lowerBoundDict[instance_uid]=0
        upperBoundDict[instance_uid]=solution.objective_value
        initialPathDict[instance_uid]=outFile
        subprocess.run(["cp",initialPathDict[instance_uid],currentBestDir+instance_uid+".solution.json"])
    writeCurrentBest(lowerBoundDict, upperBoundDict, initialPathDict, optValuesFileName)


