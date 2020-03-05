
"""
Preamble
# Program: localAlignment50000.py
# Author: Jonathan Myers
# Date: Thu Jan 16 13:33:28 2020
# Purpose: Package of methods for performing local alignment on Holder format 
DNA and RNA sequence data.
# Arguments: None.
# Loads: None.
# Calls: None.
# Returns: None.
"""

import random as r
import multiprocessing as mp
import numpy as np
import time

def dataFileOpen(stg):
    file = open(stg,"r")
    filedata = file.readlines()
    data = [element.split()[0] for element in filedata]
    file.close()
    return data

def randomizeDeltaVector(window,data):
    deltaVector = [0]*len(data)
    for i in range(len(deltaVector)):
        deltaVector[i] = r.randint(0,(len(data[i])-window))
    return deltaVector

def getProbabilities(window,deltaVector,seqIndex,data):
    tempData = [data[i] for i in range(len(data)) if i != seqIndex]
    tempDeltaVector = [deltaVector[i] for i in range(len(deltaVector)) if i != seqIndex]
    dataLength = len(tempData)
    probabilities = [[]]*window # [[A,T,G,C],[A,T,G,C],...]
    for i in range(window):
        count = [0,0,0,0] # [A,T,G,C]
        for j in range(dataLength):
            residueIndex = tempDeltaVector[j]+i
            residue = tempData[j][residueIndex].capitalize()
            if residue == "A":
                count[0] += 1
            elif residue == "T":
                count[1] += 1
            elif residue == "G":
                count[2] += 1
            else: # residue == "C"
                count[3] += 1   
        probabilities[i] = [element/dataLength for element in count]
    return probabilities

def getPropensities(window,deltaVector,seqIndex,data):
    possibleIndices = len(data[seqIndex])-window
    propensitiesVector = [0]*possibleIndices
    probabilities = getProbabilities(window,deltaVector,seqIndex,data)
    for i in range(possibleIndices):
        propensity = 1
        for j in range(window):
            residue = data[seqIndex][i+j].capitalize()
            if residue == "A":
                probability = probabilities[j][0]
            elif residue == "T":
                probability = probabilities[j][1]
            elif residue == "G":
                probability = probabilities[j][2]
            else: # residue == "C":
                probability = probabilities[j][3]
            propensity = propensity*probability
        propensitiesVector[i] = propensity
    return propensitiesVector

def getIdealResidueIndex(propensities):
    length = len(propensities)
    total = sum(propensities)
    if total == 0:
        return None
    PDFlist = [element/total for element in propensities]
    CDFlist = [PDFlist[0]]*length
    for i in range(1,length):
        CDFlist[i] = CDFlist[i-1] + PDFlist[i]
    ranNum = r.uniform(0,1)
    for i in range(length):
        if ranNum < CDFlist[i]:
            idealResidueIndex = i
            break
    return idealResidueIndex

def updateDeltaTrain(seqIndex,idealResidueIndex,deltaVector,window,data):
    if idealResidueIndex == None:
        newSeqIndex = r.randint(0,len(data)-1)
        return newSeqIndex,deltaVector
    deltaVector[seqIndex] = idealResidueIndex
    newSeqIndex = r.randint(0,len(data)-1)
    return newSeqIndex,deltaVector

def updateDeltaRecord(seqIndex,idealResidueIndex,deltaVector,trackingVector,
                trackingVectorNumber,window,data):
    if idealResidueIndex == None:
        newSeqIndex = r.randint(0,len(data)-1)
        return newSeqIndex,deltaVector,trackingVector,trackingVectorNumber
    deltaVector[seqIndex] = idealResidueIndex
    if deltaVector in trackingVector:
        location = trackingVector.index(deltaVector)
        trackingVectorNumber[location] += 1
    else:
        trackingVector.append([element for element in deltaVector])
        trackingVectorNumber.append(1)
    newSeqIndex = r.randint(0,len(data)-1)
    return newSeqIndex,deltaVector,trackingVector,trackingVectorNumber

def localAlign(window,data):
    # stg = "DPalign1.txt"
    # data = dataFileOpen(stg)
    trackingVector = []
    trackingVectorNumber = []
    seqIndex = r.randint(0,len(data)-1)
    deltaVector = randomizeDeltaVector(window,data)
    for i in range(50000):
        propensities = getPropensities(window,deltaVector,seqIndex,data)
        idealResidueIndex = getIdealResidueIndex(propensities)
        seqIndex,deltaVector = updateDeltaTrain(seqIndex,idealResidueIndex,
                                                deltaVector,window,data)
    for i in range(50000):
        propensities = getPropensities(window,deltaVector,seqIndex,data)
        idealResidueIndex = getIdealResidueIndex(propensities)
        seqIndex,deltaVector,trackingVector,trackingVectorNumber = updateDeltaRecord(
            seqIndex,idealResidueIndex,deltaVector,trackingVector,trackingVectorNumber,
            window,data)
    if trackingVectorNumber == []:
        return [-1]*len(data)
    maxIndex = trackingVectorNumber.index(max(trackingVectorNumber))
    bestLocalSeq = trackingVector[maxIndex]
    return bestLocalSeq

def test(stg):
    t1 = time.perf_counter()
    data = dataFileOpen(stg)
    window1 = 20
    results1 = localAlign(window1,data)
    window2 = 40
    results2 = localAlign(window2,data)
    t2 = time.perf_counter()
    print(t2-t1)
    np.savetxt("AlignmentNumbers1.csv",results1,fmt="%1i",delimiter=",")
    np.savetxt("AlignmentNumbers2.csv",results2,fmt="%1i",delimiter=",")

def main(poolsize,stg,windowMin,windowMax):
    t1 = time.perf_counter()
    pool = mp.Pool(poolsize)
    data = dataFileOpen(stg)
    windowSize = np.arange(windowMin,windowMax+1,5)
    results = np.array(pool.starmap(localAlign,[(window,data) for window in windowSize]))
    pool.close()
    t2 = time.perf_counter()
    timeRun = [int(t2-t1),windowMin,windowMax]+[-1]*(len(results[0])-3)
    resultsTime = np.append(results,[timeRun],axis=0)
    np.savetxt("AlignmentNumbers.csv",resultsTime,fmt="%1i",delimiter=",")

main(16,"DPalign1.txt",10,30)
#test()

# M02 End Program