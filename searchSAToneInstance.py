import os
import subprocess
import json
from multiprocessing import Pool
from multiprocessing import Process, Queue, current_process, freeze_support
from functools import partial
from itertools import repeat
import argparse

maxThreads=17

def lexSmaller(distA, distB,j):
    for i in range(j):
        if distA[i]<distB[i]:
            return True
        elif distA[i]>distB[i]:
            return False
    return False

#
# Function run by worker processes (from docs.python.org)
#

def worker(input, output):
  for func, args in iter(input.get, 'STOP'):
    result = calculate(func, args)
    output.put(result)

#
# Function used to calculate result
#

def calculate(func, args):
  result = func(*args)
  return '%s says that %s%s = %s' % \
    (current_process().name, func.__name__, args, result)

def executeAndVerify(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inFile, inDir, distances):
  instance_uid=os.path.splitext(inFile)[0]
  instanceUidAndDistances=instance_uid

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
  # print(inFile)
  # print(args)
  subprocess.run(args,capture_output=True)
  if satSolverPath=="./build/mallob":
    os.chdir("../mallob/")
    # pathPrefix="../cgshop2026/"
    subprocess.run([satSolverPath,"-mono="+cnfFile,"-q","-t=12","-satsolver=k","-s2f="+assignmentFile],capture_output=True)
    os.chdir("../cgshop26/")
  else:
    with open(assignmentFile,"w") as asgnFile:
      subprocess.run([satSolverPath,cnfFile,"-q"],stdout=asgnFile)
  line0=""
  with open(assignmentFile,"r") as asgnFile:
    line0=asgnFile.readline()
    line0=line0.strip()
  # print(line0)
  allPos=True
  for d in distances:
    if d<0:
        allPos=False
  if not allPos:
    subprocess.run(["rm",idxFile])
    subprocess.run(["rm",cnfFile])
    #subprocess.run(["rm",assignmentFile])
  if line0=="s SATISFIABLE":
    subprocess.run([executablePath,assignmentFile,outFile,"assignmentToFlipListNU","nonVerbose",idxFile],capture_output=True)
    if allPos:
        #subprocess.run([executablePath,assignmentFile,outFile,"assignmentToFlipListNU","nonVerbose",idxFile],capture_output=True)
        subprocess.run(["cp",outFile,outDir+instance_uid+".solution.json"])
    return True
  if allPos:
    subprocess.run(["rm",idxFile])
    subprocess.run(["rm",cnfFile])
    subprocess.run(["rm",assignmentFile])
  return False
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
  pairwiseDist=[[]]
  if os.path.isfile(pairwiseDistFileName):
    with open(pairwiseDistFileName,"r") as readFile:
      pairwiseDist=[[int(dist) for dist in line.split()] for line in readFile]
    if len(pairwiseDist)==m and len(pairwiseDist[m-1])==m:
      return pairwiseDist
    # else: restart from where you ended
    #   subprocess.run(["rm",pairwiseDistFileName])

  # print("Here")
  startT0=len(pairwiseDist)
  if startT0>0:
    t0=startT0-1
    startT1=len(pairwiseDist[t0])
    if startT1<t0:
      with open(pairwiseDistFileName,"a") as outFile:
        for t1 in range(startT1,t0):
          pairwiseDist[t0].append(pairwiseDist[t1][t0])
          outFile.write(str(pairwiseDist[t0][t1])+' ')
    if startT1<=t0:
      pairwiseDist[t0].append(0)
      with open(pairwiseDistFileName,"a") as outFile:
        if t0==m-1:
          outFile.write("0\n")
        else:
          outFile.write("0 ")

    for t1 in range(max(startT1,t0+1),m):
      pairwiseDist[t0].append(pairwiseDistance(satSolverPath,executablePath,cnfDir,outDir,idxDir,assignmentDir,inDir,inFile,m,t0,t1))
      with open(pairwiseDistFileName,"a") as outFile:
        if t1==m-1:
          outFile.write(str(pairwiseDist[t0][t1])+'\n')
        else:
          outFile.write(str(pairwiseDist[t0][t1])+' ')
  for t0 in range(startT0,m):
    pairwiseDist.append([])
    with open(pairwiseDistFileName,"a") as outFile:
      for t1 in range(t0):
        pairwiseDist[t0].append(pairwiseDist[t1][t0])
        outFile.write(str(pairwiseDist[t0][t1])+' ')
    pairwiseDist[t0].append(0)
    with open(pairwiseDistFileName,"a") as outFile:
      if t0==m-1:
        outFile.write("0\n")
      else:
        outFile.write("0 ")
    nThreads=min(maxThreads,m-t0)
    with Pool(nThreads) as p:
      partialPairwiseDist=partial(pairwiseDistance,satSolverPath,executablePath,cnfDir,outDir,idxDir,assignmentDir,inDir,inFile,m,t0)
      res=list(p.map(partialPairwiseDist,range(t0+1,m)))
      for dist in res:
        pairwiseDist[t0].append(dist)
    for t1 in range(t0+1,m):
      #pairwiseDist[t0].append(pairwiseDistance(satSolverPath,executablePath,cnfDir,outDir,idxDir,assignmentDir,inDir,inFile,m,t0,t1))
      with open(pairwiseDistFileName,"a") as outFile:
        if t1==m-1:
          outFile.write(str(pairwiseDist[t0][t1])+'\n')
        else:
          outFile.write(str(pairwiseDist[t0][t1])+' ')

  return pairwiseDist

  #for t0 in range(m):
  #  pairwiseDist.append([])
  #  for t1 in range(t0):
  #    pairwiseDist[t0].append(pairwiseDist[t1][t0])
  #  pairwiseDist[t0].append(0)
  #  res=(list(map(partial(pairwiseDistance,satSolverPath,executablePath,cnfDir,outDir,idxDir,assignmentDir,inDir,inFile,m,t0),range(t0+1,m))))
  #  for dist in res:
  #    pairwiseDist[t0].append(dist)




  # with open(pairwiseDistFileName,"w") as outFile:
  #   for t0 in range(m):
  #     for t1 in range(m):
  #       if t1==m-1:
  #         outFile.write(str(pairwiseDist[t0][t1])+'\n')
  #       else:
  #         outFile.write(str(pairwiseDist[t0][t1])+' ')
  # return pairwiseDist

def completeSearch(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inDir, inFile, totDist,
                   pairwiseDist, m, minFrom, i, distances):
  # if i>0 and totDist==66 and inFile=="random_instance_982_320_10.json" and lexSmaller(distances,[6,6,6,9,6,6,8,6,7]):
  #   return False
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
  firstNonNeg=i
  for j in range(i):
    if distances[j]>=0:
      firstNonNeg=min(j,firstNonNeg)
      lowerBoundDistI=max(lowerBoundDistI,pairwiseDist[i][j]-distances[j])
      upperBoundDistI-=distances[j]

  if upperBoundDistI<minFrom[i]:
    return False

  for k in range(i+1,m):
    minForK=0
    for j in range(i):
      if distances[j]>=0:
        minForK=max(minForK,pairwiseDist[j][k]-distances[j])
    upperBoundDistI-=minForK

  if lowerBoundDistI<=upperBoundDistI and ( i==(m+firstNonNeg)//2+5) and distances[0]>=0:#(i==10 or i==15): i==(m+firstNonNeg)//2 or
    distances2=[-1]*m
    for j in range(i):
      distances2[j]=distances[j]
    if not executeAndVerify(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inFile, inDir, distances2):
      return False
    print("So far it's sat")
    print(distances2)

  for distI in range(lowerBoundDistI,upperBoundDistI+1):
    distances2=[0]*m
    for j in range(i):
      distances2[j]=distances[j]
    distances2[i]=distI
    if completeSearch(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inDir, inFile, totDist,
                      pairwiseDist, m, minFrom, i + 1, distances2):
      return True

  return False

def generateDistancesList(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inDir, inFile, totDist,
                   pairwiseDist, m, minFrom, i, maxI, distances):
  if i==maxI:
    distances[i]=totDist
    for j in range(i):
      if distances[j]>=0:
        distances[i]-=distances[j]
    if distances[i]<0:
      return []
    for j in range(i):
      if distances[j]>=0 and distances[i]+distances[j]<pairwiseDist[i][j]:
        return []
    # print(distances)
    return [distances]
  lowerBoundDistI=0
  upperBoundDistI=totDist
  firstNonNeg=i
  for j in range(i):
    if distances[j]>=0:
      firstNonNeg=min(j,firstNonNeg)
      lowerBoundDistI=max(lowerBoundDistI,pairwiseDist[i][j]-distances[j])
      upperBoundDistI-=distances[j]

  if upperBoundDistI<minFrom[i]:
    return []

  for k in range(i+1,m):
    minForK=0
    for j in range(i):
      if distances[j]>=0:
        minForK=max(minForK,pairwiseDist[j][k]-distances[j])
    upperBoundDistI-=minForK

  #if lowerBoundDistI<=upperBoundDistI and (i==(m+firstNonNeg)//2 or i==(m+firstNonNeg)//2+5):# and distances[0]>=0:#(i==10 or i==15):
  #   distances2=[-1]*m
  #   for j in range(i):
  #     distances2[j]=distances[j]
  #   if not executeAndVerify(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inFile, inDir, distances2):
  #     return []
  #   print("So far it's sat")
  #   print(distances2)

  distanceL=[]
  for distI in range(lowerBoundDistI,upperBoundDistI+1):
    distances2=[0]*m
    for j in range(i):
      distances2[j]=distances[j]
    distances2[i]=distI
    tmpList=generateDistancesList(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inDir, inFile, totDist,
                           pairwiseDist, m, minFrom, i + 1, m-1, distances2)
    if len(tmpList)>0:
      distanceL+=tmpList


  return distanceL

def searchDistanceList(satSolverPath,executablePath,cnfDir, outDir, idxDir, assignmentDir, inFile, inDir,distanceList,parallelEx=False):
  distListWithKey=[]
  #print(distanceList)
  m=len(distanceList[0])
  firstNonNeg=0
  while distanceList[0][firstNonNeg]==-1:
    firstNonNeg+=1
  halfIdx=(m+firstNonNeg)//2
  distanceList.sort()
  i=0
  while i<len(distanceList):
    distances=distanceList[i]
    distances2=[-1]*m
    for j in range(halfIdx):
      distances2[j]=distances[j]
    print(distances)
    if not executeAndVerify(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inFile, inDir, distances2):
      print("First half not sat: "+str(distances2))
      k=i
      while i<len(distanceList) and distanceList[k][:halfIdx]==distanceList[i][:halfIdx]:
        i+=1
    else:
      k=i
      print("First half sat: "+str(distances2))
      while i<len(distanceList) and distanceList[k][:halfIdx-1]==distanceList[i][:halfIdx-1]: #I can continue when halfIdx-1 changes because it's in lexicographic order, so it can only increase
        distances=distanceList[i]
        avgDist=sum(distances)/len(distances)
        sumAbsDiff=0
        for d in distances:
          sumAbsDiff+=abs(d-avgDist)
        distListWithKey.append([sumAbsDiff]+distances)
        i+=1
  #for distance in distanceList:
  #  avgDist=sum(distance)/len(distance)
  #  sumAbsDiff=0
  #  for d in distance:
  #    sumAbsDiff+=abs(d-avgDist)
  #  distListWithKey.append([sumAbsDiff]+distance)
  distListWithKey.sort(key=lambda distances : distances[0])
  j=0
  nThreads=min(1+len(distListWithKey),maxThreads)
  sortedDistList=[]
  foundLastEx=False
  for distWithKey in distListWithKey:
    distances=distWithKey[1:]
    if inFile=="random_instance_915_320_20.json" and distances[0:12]==[-1]*12 and sum(distances[13:20])==44:

      if distances[13:20]==[4,5,4,4,4,6,17]: # 4,3,4,3,4,5,4,5,4,4,13 [4,5,4,3,6,7,5,4,6,4,5]: # 4,4,5,6,4,6,5,4,4,5,6]:
        foundLastEx=True
    else:
      if inFile=="random_instance_552_320_20.json" and distances[0:10]==[-1]*10 and sum(distances[10:20])==58:
        if distances[10:20]==[4,6,5,6,6,7,5,6,6,7]:
          foundLastEx=True
      else:
        foundLastEx=True
    if foundLastEx:
      sortedDistList.append(distances)
      #print("To search: "+str(j)+"/"+str(len(distListWithKey)-1))
    j+=1
  print("To search: "+str(len(sortedDistList)))
  if parallelEx:
    partialEx=partial(executeAndVerify,satSolverPath,executablePath,cnfDir, outDir, idxDir, assignmentDir, inFile, inDir)
    nThreads=maxThreads
    with Pool(nThreads) as p:
      resIt=(p.imap_unordered(partialEx,sortedDistList))
      for r in resIt:
        if r:
          return True
    return False
  
  j=0
  for distances in sortedDistList:
    #distances=distWithKey[1:]
    print(str(j)+"/"+str(len(sortedDistList)-1))
    j+=1
    print(distances)
    isSat=executeAndVerify(satSolverPath,executablePath,cnfDir, outDir, idxDir, assignmentDir, inFile, inDir, distances)
    if isSat:
      return True

  return False

def binSearchSuffix(satSolverPath,executablePath,cnfDir,outDir,idxDir,assignmentDir,inDir, inFile,pairwiseDist,low,high,m,minFrom,suffixStart):
  print("Suffix start: "+str(suffixStart))
  while low<high:
    mid=(low+high)//2
    print("Searching at pfd: "+str(mid))

    distanceList=generateDistancesList(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inDir, inFile, mid,
                             pairwiseDist, m, minFrom, suffixStart, m-1, [-1] * m)
    print("Generated distance list")
    foundSol=searchDistanceList(satSolverPath,executablePath,cnfDir, outDir, idxDir, assignmentDir, inFile, inDir, distanceList,True)
    if foundSol:
      high=mid
    else:
      low=mid+1
  print(low)
  return low

def linSearchSuffix(satSolverPath,executablePath,cnfDir,outDir,idxDir,assignmentDir,inDir, inFile,pairwiseDist,low,m,minFrom,suffixStart):
  print("Suffix start: "+str(suffixStart))
  curDist=low
  while True:
    print("Searching at pfd: "+str(curDist))

    distanceList=generateDistancesList(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inDir, inFile, curDist,
                             pairwiseDist, m, minFrom, suffixStart, m-1, [-1] * m)
    print("Generated distance list")
    foundSol=searchDistanceList(satSolverPath,executablePath,cnfDir, outDir, idxDir, assignmentDir, inFile, inDir, distanceList,True)
    if foundSol:
      return curDist
      print("Suffix distance is: "+curDist)

def binSearchSuffixFindHigh(satSolverPath,executablePath,cnfDir,outDir,idxDir,assignmentDir,inDir, inFile,pairwiseDist,low,m,minFrom,suffixStart):
  high=max(low,1)
  foundSol= completeSearch(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inDir, inFile, high,
                           pairwiseDist, m, minFrom, suffixStart, [-1] * m)

  while not foundSol:
    high*=2
    foundSol= completeSearch(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inDir, inFile, high,
                             pairwiseDist, m, minFrom, suffixStart, [-1] * m)

  while low<high:
    mid=(low+high)//2
    print("Searching at pfd: "+str(mid)+"["+str(low)+","+str(high)+"]")
    foundSol= completeSearch(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inDir, inFile, mid,
                             pairwiseDist, m, minFrom, suffixStart, [-1] * m)
    if foundSol:
      high=mid
    else:
      low=mid+1
  print(low)
  return low
def completeSearch3(satSolverPath,executablePath,cnfDir,outDir,idxDir,assignmentDir,inDir, inFile,totDist,pairwiseDist):
  #larger distances increase the number of clauses, so search starting from the distances with the lowest maximum distance
  candidateDist=[]
  for d0 in range(0,totDist+1):
    for d1 in range(max(0,pairwiseDist[0][1]-d0),totDist+1-d0):
      d2=totDist-d0-d1
      if d0+d1>=pairwiseDist[0][1] and d0+d2>=pairwiseDist[0][2] and d1+d2>=pairwiseDist[1][2]:
        candidateDist.append([max(d0,d1,d2),d0,d1,d2])

  candidateDist.sort(key=lambda distances : distances[0])
  for distWithKey in candidateDist:
    distances=distWithKey[1:]
    print(distances)
    isSat= executeAndVerify(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inFile, inDir, distances)
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
          for d4 in range(max(0,pairwiseDist[0][4]-d0,pairwiseDist[1][4]-d1,pairwiseDist[2][4]-d2,pairwiseDist[3][4]-d3),totDist+1-d0-d1-d2-d3):
            d5=totDist-d0-d1-d2-d3-d4
            if (d0+d5>=pairwiseDist[0][5] and d1+d5>=pairwiseDist[1][5]
                    and d2+d5>=pairwiseDist[2][5] and d3+d5>=pairwiseDist[3][5] and d4+d5>=pairwiseDist[4][5]):
              candidateDist.append([max(d0,d1,d2,d3,d4,d5),d0,d1,d2,d3,d4,d5])

  candidateDist.sort(key=lambda distances : distances[0])
  # j=0
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
  j=0
  for distWithKey in candidateDist:
    distances=distWithKey[1:]
    print(str(j)+"/"+str(len(candidateDist)-1))
    j+=1
    print(distances)
    isSat=executeAndVerify(satSolverPath,executablePath,cnfDir, outDir, idxDir, assignmentDir, inFile, inDir, distances)
    if isSat:
      return True

  return False

def completeSearch7(satSolverPath,executablePath,cnfDir,outDir,idxDir,assignmentDir,inDir, inFile,totDist,pairwiseDist,minFrom):
  #larger distances increase the number of clauses, so search starting from the distances with the lowest maximum distance
  candidateDist=[]
  for d0 in range(0,totDist+1):
    for d1 in range(max(0,pairwiseDist[0][1]-d0),totDist+1-d0):
      for d2 in range(max(0,pairwiseDist[0][2]-d0,pairwiseDist[1][2]-d1),totDist+1-d0-d1):
        for d3 in range(max(0,pairwiseDist[0][3]-d0,pairwiseDist[1][3]-d1,pairwiseDist[2][3]-d2),totDist+1-d0-d1-d2):
          if totDist-d0-d1-d2-d3>=minFrom[4]:
            for d4 in range(max(0,pairwiseDist[0][4]-d0,pairwiseDist[1][4]-d1,pairwiseDist[2][4]-d2,pairwiseDist[3][4]-d3),totDist+1-d0-d1-d2-d3):
              if totDist-d0-d1-d2-d3-d4>=minFrom[5]:
                for d5 in range(max(0,pairwiseDist[0][5]-d0,pairwiseDist[1][5]-d1,pairwiseDist[2][5]-d2,pairwiseDist[3][5]-d3,pairwiseDist[4][5]-d4),totDist+1-d0-d1-d2-d3-d4):
                  d6=totDist-d0-d1-d2-d3-d4-d5
                  if (d0+d6>=pairwiseDist[0][6] and d1+d6>=pairwiseDist[1][6]
                          and d2+d6>=pairwiseDist[2][6] and d3+d6>=pairwiseDist[3][6] and d4+d6>=pairwiseDist[4][6] and d5+d6>=pairwiseDist[5][6]):
                    candidateDist.append([max(d0,d1,d2,d3,d4,d5),d0,d1,d2,d3,d4,d5,d6])

  candidateDist.sort(key=lambda distances : distances[0])
  j=0
  for distWithKey in candidateDist:
    distances=distWithKey[1:]
    print(str(j)+"/"+str(len(candidateDist)-1))
    j+=1
    print(distances)
    isSat=executeAndVerify(satSolverPath,executablePath,cnfDir, outDir, idxDir, assignmentDir, inFile, inDir, distances)
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

  # html=convert(dataFromJson)
  # with open('./scriptOutputs/currentBest/instancesStats.html',mode="w", encoding="utf-8") as tableFile:
  #   tableFile.write(html)

def parallelSearch1(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir,
                           inDir, inFile, m, pairwiseDist, minFrom, mid):
  prefixLen=9
  distanceLPref=generateDistancesList(satSolverPath,executablePath,cnfDir,outDir,idxDir,assignmentDir,inDir,inFile,mid,pairwiseDist,m,minFrom,0,prefixLen-1,[-1]*m)
  distListWithKey=[]
  for distances in distanceLPref:
    avgDist=sum(distances)/len(distances)
    sumAbsDiff=0
    for d in distances:
      sumAbsDiff+=abs(d-avgDist)
    distListWithKey.append([sumAbsDiff]+distances)
  distListWithKey.sort(key=lambda distances : distances[0])
  lenBuffer=maxThreads*10
  distListBuffer=[]
  partialCompleteSearch=partial(completeSearch,satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inDir, inFile,
                                mid, pairwiseDist, m, minFrom, prefixLen)
  for distWithKey in distListWithKey:
    distances=distWithKey[1:]
    i=prefixLen
    lowerBoundDistI=0
    upperBoundDistI=mid
    for j in range(i):
      if distances[j]>=0:
        lowerBoundDistI=max(lowerBoundDistI,pairwiseDist[i][j]-distances[j])
        upperBoundDistI-=distances[j]

    if upperBoundDistI>=minFrom[i]:

      for k in range(i+1,m):
        minForK=0
        for j in range(i):
          if distances[j]>=0:
            minForK=max(minForK,pairwiseDist[j][k]-distances[j])
        upperBoundDistI-=minForK

      if lowerBoundDistI<=upperBoundDistI:
        distances2=[-1]*m
        for j in range(i):
          distances2[j]=distances[j]
        if executeAndVerify(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inFile, inDir, distances2):
          print("added to buffer: "+str(distances))
          distListBuffer.append(distances)

    if len(distListBuffer)==lenBuffer:
      print("Searching buffer")
      with Pool(nThreads) as p:
        # print(len(distancesL[len(distancesL)-len(distancesL)%nThreads:]))
        resIt=(p.imap_unordered(partialCompleteSearch,distListBuffer))#[len(distancesL)-len(distancesL)%nThreads:]))
        for r in resIt:
          if r:
            return True
      distListBuffer.clear()
  return False

def parallelSearch(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir,
                   inDir, inFile, m, pairwiseDist, minFrom, mid):
  nThreads=maxThreads

  # j=0
  # while j<=mid:
  distancesL=[]
  for d0 in range(0,mid+1):
    maxD1=mid+1-d0
    for k in range(2,m):
      minForK=max(0,pairwiseDist[0][k]-d0)
      maxD1-=minForK
    for d1 in range(max(0,pairwiseDist[0][1]-d0),maxD1):
      maxD2=mid+1-d0-d1
      for k in range(3,m):
        minForK=max(0,pairwiseDist[0][k]-d0,pairwiseDist[1][k]-d1)
        maxD2-=minForK
      minD2=max(0,pairwiseDist[0][2]-d0,pairwiseDist[1][2]-d1)
      for d2 in range(minD2,maxD2):
        maxD3=mid+1-d0-d1-d2
        for k in range(4,m):
          minForK=max(0,pairwiseDist[0][k]-d0,pairwiseDist[1][k]-d1,pairwiseDist[2][k]-d2)
          maxD3-=minForK
        for d3 in range(max(0,pairwiseDist[0][3]-d0,pairwiseDist[1][3]-d1,pairwiseDist[2][3]-d2),maxD3):
          maxD4=mid+1-d0-d1-d2-d3
          for k in range(5,m):
            minForK=max(0,pairwiseDist[0][k]-d0,pairwiseDist[1][k]-d1,pairwiseDist[2][k]-d2,pairwiseDist[3][k]-d3)
            maxD4-=minForK
          for d4 in range(max(0,pairwiseDist[0][4]-d0,pairwiseDist[1][4]-d1,pairwiseDist[2][4]-d2,pairwiseDist[3][4]-d3),maxD4):
            if mid-d0-d1-d2-d3-d4>=minFrom[5]:
              maxD5=mid+1-d0-d1-d2-d3-d4
              for k in range(6,m):
                minForK=max(0,pairwiseDist[0][k]-d0,pairwiseDist[1][k]-d1,pairwiseDist[2][k]-d2,pairwiseDist[3][k]-d3,pairwiseDist[4][k]-d4)
                maxD5-=minForK
              for d5 in range(max(0,pairwiseDist[0][5]-d0,pairwiseDist[1][5]-d1,pairwiseDist[2][5]-d2,pairwiseDist[3][5]-d3,pairwiseDist[4][5]-d4),maxD5):
                distances2=[-1]*m
                distances2[0]=d0
                distances2[1]=d1
                distances2[2]=d2
                distances2[3]=d3
                distances2[4]=d4
                distances2[5]=d5
                distancesL.append(distances2)
  #print(distancesL)

  partialCompleteSearch=partial(completeSearch,satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inDir, inFile,
                                mid, pairwiseDist, m, minFrom, 6)


  # for i in range(0):#,len(distancesL),nThreads):
  #   with Pool(nThreads) as p:
  #     for distances in distancesL[i:i+nThreads]:
  #       print(str(distances[0])+" "+str(distances[1]))
  #     # print((distancesL[i:i+nThreads]))
  #     resIt=(p.imap_unordered(partialCompleteSearch,distancesL[i:i+nThreads]))
  #     for r in resIt:
  #       if r:
  #         return True
  with Pool(nThreads) as p:
    # print(len(distancesL[len(distancesL)-len(distancesL)%nThreads:]))
    resIt=(p.imap_unordered(partialCompleteSearch,distancesL))#[len(distancesL)-len(distancesL)%nThreads:]))
    for r in resIt:
      if r:
        return True
  return False
def solveInstance(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir,
                  inDir, pairwiseDistDir, lowerBoundDict,upperBoundDict,inFile):
  with open(inDir+inFile, 'r') as file:
    data = json.load(file)
  n=len(data["points_x"])
  m=len(data["triangulations"])
  instance_uid=data["instance_uid"]

  print("Solving: "+instance_uid)
  if m==2:
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
    minFrom=[0]*m
    minFrom[m-2]=pairwiseDist[m-2][m-1]
    # startMinFrom=m-2
    # if m==20:
    startMinFrom=m//2
    if m==20:
      startMinFrom=m//2
    minFromFileName="./scriptOutputs/minFromDir/"+instance_uid+"_minFrom.out"
    if os.path.isfile(minFromFileName):
      with open(minFromFileName,"r") as readFile:
        minFrom=[[int(dist) for dist in line.split()] for line in readFile][0]
    else:
      if instance_uid=="random_instance_213_80_20":
        minFrom=[45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 40, 36, 32, 29, 24, 20, 16, 12, 8, 0]
      else:
        #solution=read_solution("./scriptOutputs/currentBest/"+instance_uid+".solution.json")
        curBestPath="./scriptOutputs/currentBest/"+instance_uid+".solution.json"
        with open(curBestPath, 'r') as file:
          dataCurBest = json.load(file)
        distancesCurBest=[]
        for flipSequence in dataCurBest["flips"]:
          print(len(flipSequence))
          distancesCurBest.append(len(flipSequence))
        
        #high1=upperBoundDict[instance_uid]
        #for i in range(m):
        #  curSum=0
        #  for j in range(startMinFrom,m):
        #    curSum+=pairwiseDist[i][j]
        #  high1=min(high1,curSum)
        #for suffixStart in range(startMinFrom,m-2):
        #  low1=0
          #high1=upperBoundDict[instance_uid]
        #  if suffixStart>startMinFrom:
        #    high1=minFrom[suffixStart-1]
        #    for i in range(m):
        #      curSum=0
        #      for j in range(suffixStart,m):
        #        curSum+=pairwiseDist[i][j]
        #      high1=min(high1,curSum)
        #  minFrom[suffixStart]=binSearchSuffix(satSolverPath,executablePath,cnfDir,outDir,idxDir,assignmentDir,inDir,inFile,pairwiseDist,low1,high1,m,minFrom, suffixStart)
        endMinFrom=m-3
        if instance_uid=="random_instance_915_320_20":
          endMinFrom=13
          minFrom[17]=15
          minFrom[16]=22
          minFrom[15]=30
          minFrom[14]=36
        for suffixStart in range(endMinFrom,startMinFrom-1,-1):
        #   # print("Suffix start: "+str(suffixStart))
           low1=minFrom[suffixStart+1]
           if instance_uid=="random_instance_915_320_20" and suffixStart==14:
            minFrom[suffixStart]=binSearchSuffix(satSolverPath,executablePath,cnfDir,outDir,idxDir,assignmentDir,inDir,inFile,pairwiseDist
,low1,37,m,minFrom, suffixStart)
           else: 
            minFrom[suffixStart]=binSearchSuffix(satSolverPath,executablePath,cnfDir,outDir,idxDir,assignmentDir,inDir,inFile,pairwiseDist
,low1,sum(distancesCurBest[suffixStart:]),m,minFrom, suffixStart)
        for suffixStart in range(startMinFrom):
          minFrom[suffixStart]=minFrom[startMinFrom]
    # if instance_uid=="random_instance_603_40_20":
    #   minFrom=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 35, 31, 27, 22, 18, 15, 11, 7, 0, 4]
        with open(minFromFileName,"a") as outFile:
          for minFromValue in minFrom:
            outFile.write(str(minFromValue)+' ')
    print(instance_uid+" MIN FROM: "+str(minFrom))
    low=max(minFrom[startMinFrom],lowerBoundDict[instance_uid])
    high=upperBoundDict[instance_uid]
    while low<high:
      mid=(low+high)//2
      print(instance_uid+" searching at pfd: "+str(mid)+"["+str(low)+","+str(high)+"]")
      foundSol=False
      if m<20:
        distanceList=generateDistancesList(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inDir, inFile, mid,
                                           pairwiseDist, m, minFrom, 0,m-1,[-1] * m)
        print("Generated distance list")
        foundSol=searchDistanceList(satSolverPath,executablePath,cnfDir, outDir, idxDir, assignmentDir, inFile, inDir, distanceList)
      else:
        foundSol=parallelSearch1(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inDir, inFile, m, pairwiseDist, minFrom, mid)
        # foundSol= completeSearch(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inDir, inFile,
        #                          mid, pairwiseDist, m, minFrom, 0, [0] * m)
      # foundSol=False
      # if m==3:
      #   foundSol=completeSearch3(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inDir, inFile,mid,pairwiseDist)
      # elif m==4:
      #   foundSol=completeSearch4(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inDir, inFile,mid,pairwiseDist)
      # elif m==5:
      #   foundSol=completeSearch5(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inDir, inFile,mid,pairwiseDist)
      # elif m==6:
      #   foundSol=completeSearch6(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inDir, inFile,mid,pairwiseDist)
      # elif m==7:
      #   foundSol=completeSearch7(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inDir, inFile,mid,pairwiseDist,minFrom)
      # else:
      #   # foundSol=parallelSearch(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inDir, inFile, m, pairwiseDist, minFrom, mid)
      #   foundSol= completeSearch(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inDir, inFile,
      #                            mid, pairwiseDist, m, minFrom, 0, [0] * m)
      if foundSol:
        high=mid
      else:
        low=mid+1
      #lowerBoundDict[instance_uid]=low
      #upperBoundDict[instance_uid]=high
      #writeCurrentBest(lowerBoundDict, upperBoundDict, initialPathDict, optValuesFileName)
    lowerBoundDict[instance_uid]=low
    # foundSol=completeSearch(satSolverPath, executablePath, cnfDir, outDir, idxDir, assignmentDir, inDir, inFile,low,pairwiseDist,[0]*m,m,0)
    if low<upperBoundDict[instance_uid]:
      print("Improved solution from pfd "+str(upperBoundDict[instance_uid])+" to pfd "+str(low))
      upperBoundDict[instance_uid]=low
      initialPathDict[instance_uid]=outDir+instance_uid+".solution.json"
      subprocess.run(["cp",initialPathDict[instance_uid],currentBestDir+instance_uid+".solution.json"])

  writeCurrentBest(lowerBoundDict, upperBoundDict, initialPathDict, optValuesFileName)

if __name__ == '__main__':
  # print("Test")
  freeze_support()
  optValuesFileName="./scriptOutputs/currentBest/instancesStats.json"
  solutionsDir="./scriptOutputs/satSolutions_114/"
  currentBestDir="./scriptOutputs/currentBest/"

  lowerBoundDict=dict()
  upperBoundDict=dict()
  initialPathDict=dict()
  readCurrentBest(lowerBoundDict,upperBoundDict,initialPathDict,optValuesFileName)

  baseDir="/fastStorage"
  #baseDir="/cluster/scratch/lbattini"
  inDir="../benchmark_instances_rev1/benchmark_instances/"
  #inDir="../benchmark_instances/"
  pairwiseDistDir="./scriptOutputs/pairwiseDistDir1/"
  dirCnt=0
  cnfDir=baseDir+"/scriptOutputs/cnfDir_"+str(dirCnt)+"/"
  while os.path.isdir(cnfDir):
    dirCnt+=1
    cnfDir=baseDir+"/scriptOutputs/cnfDir_"+str(dirCnt)+"/"

  outDir=baseDir+"/scriptOutputs/satSolutions_"+str(dirCnt)+"/"
  idxDir=baseDir+"/scriptOutputs/idxDir_"+str(dirCnt)+"/"
  assignmentDir=baseDir+"/scriptOutputs/assignments_"+str(dirCnt)+"/"
  os.mkdir(cnfDir)
  os.mkdir(idxDir)
  os.mkdir(outDir)
  os.mkdir(assignmentDir)
  executablePath="cmake-build-release/centralTriangulationClion"
  #executablePath="./centralTriangulationClion"
  satSolverPath="../kissat/build/kissat"
  #satSolverPath="./build/mallob"
  cnt=0
  executionStatsFile=outDir+"values.out"
  
  parser = argparse.ArgumentParser()
  parser.add_argument("fileIn")
  args = parser.parse_args()
  fileIn=args.fileIn

  solveInstance(satSolverPath,executablePath,cnfDir,outDir,idxDir,assignmentDir,inDir,pairwiseDistDir,lowerBoundDict,upperBoundDict,fileIn)
  print("Solved 1 instance")
