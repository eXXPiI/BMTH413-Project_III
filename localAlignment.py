
"""
Preamble
# Program: localAlignment.py
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

def getMotif(window,data,bestLocalSeqDelta):
    dataLength = len(data)
    totals = [[]]*window # [[A,T,G,C],[A,T,G,C],...]
    residues = [""]*window
    for i in range(window):
        count = [0,0,0,0] # [A,T,G,C]
        for j in range(dataLength):
            residueIndex = bestLocalSeqDelta[j]+i
            residue = data[j][residueIndex].capitalize()
            if residue == "A":
                count[0] += 1
            elif residue == "T":
                count[1] += 1
            elif residue == "G":
                count[2] += 1
            else: # residue == "C"
                count[3] += 1   
        totals[i] = [element for element in count]
    maxIndices = [counts.index(max(counts)) for counts in totals]
    for i in range(len(maxIndices)):
        num = maxIndices[i]
        if num == 0:
            residues[i] = "A"
        elif num == 1:
            residues[i] = "T"
        elif num == 2:
            residues[i] = "G"
        else:
            residues[i] = "C"
    return residues

def localAlign(window,data):
    # stg = "DPalign1.txt"
    # data = dataFileOpen(stg)
    trackingVector = []
    trackingVectorNumber = []
    seqIndex = r.randint(0,len(data)-1)
    deltaVector = randomizeDeltaVector(window,data)
    for i in range(200000):
        propensities = getPropensities(window,deltaVector,seqIndex,data)
        idealResidueIndex = getIdealResidueIndex(propensities)
        seqIndex,deltaVector = updateDeltaTrain(seqIndex,idealResidueIndex,
                                                deltaVector,window,data)
    for i in range(100000):
        propensities = getPropensities(window,deltaVector,seqIndex,data)
        idealResidueIndex = getIdealResidueIndex(propensities)
        seqIndex,deltaVector,trackingVector,trackingVectorNumber = updateDeltaRecord(
            seqIndex,idealResidueIndex,deltaVector,trackingVector,trackingVectorNumber,
            window,data)
    if trackingVectorNumber == []:
        return [""]*len(data)
    maxIndex = trackingVectorNumber.index(max(trackingVectorNumber))
    bestLocalSeqDelta = trackingVector[maxIndex]
    bestLocalSeq = getMotif(window,data,bestLocalSeqDelta)
    return bestLocalSeq

def test(stg):
    data = dataFileOpen(stg)
    window1 = 3
    results1 = localAlign(window1,data)
    window2 = 4
    results2 = localAlign(window2,data)
    file = open("AlignmentMotifs.csv","w")
    results = [results1,results2]
    for line in results:
        for i in line:
            file.write(str(i))
        file.write("\n")
    file.close()

def main(poolsize,stg,windowMin,windowMax):
    t1 = time.perf_counter()
    pool = mp.Pool(poolsize)
    data = dataFileOpen(stg)
    windowSize = np.arange(windowMin,windowMax+1,1)
    results = pool.starmap(localAlign,[(window,data) for window in windowSize])
    pool.close()
    t2 = time.perf_counter()
    timeRun = [int(t2-t1),windowMin,windowMax]
    results.append(timeRun)
    file = open("AlignmentMotifs.csv","w")
    for line in results:
        for i in line:
            file.write(str(i))
        file.write("\n")
    file.close()

main(42,"DPalign1.txt",5,55)
#test("DPalign3.txt")

# M02 End Program