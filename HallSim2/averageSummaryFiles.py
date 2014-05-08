# -*- coding: utf-8 -*-
import os
import csv
import numpy
import matplotlib
import matplotlib.pylab as pylab
import matplotlib.pyplot
import pdb
from numpy import array, linspace
from scipy.interpolate import spline
import re



def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

os.chdir("C:\Users\Jenna\Documents\Visual Studio 2012\Projects\HallSim2\Output\Working3")

#Edit these each run
numRuns = 1 #This is the number of iterations (i in visual studio)
numReps = 37 #This is the number of steps the program took (so if we go to 14800, it is that divided by 4, which is 37) (9600=24)
import sys
import os.path

cancerCells = []
regularCells = []
data = []
totalCells = []
sgCells = []
igiCells = []
aaCells = []
itCells = []
aCells = []
drCells = []
aiCells = []
guCells = []
aiCells = []
letterFinal = 'x'

folder = "C:\Users\Jenna\Documents\Visual Studio 2012\Projects\HallSim2\Output\Working3"
for filename in os.listdir (folder):
    #add the time (from file name?) and the number of cancer cells to cancerCell list
    #add the time and the number of regular cells
    #Plot the cancer cells and regular cells on a graph against time
    
    '''    
    “NUMPY F {'descr': '|S94', 'fortran_order': False, 'shape': (), }             
    Total 1338
    Cancer 539
    Noncancer 799
    all 112
    sg 168
    igi 332
    aa 137
    it 156
    gu 145
    '''
    #Load in data
    data[:] = []
    with open(filename) as f:
        next(f)
        for line in f:
            data.append(line)
    
    #currentStep = re.sub("summary_Cells_[a-z]_11_18_1_", "", filename)
    currentStep = re.sub("summary_Cells_[a-z]_01_2_\d+_", "", filename)
    #currentStep = filename.replace("summary_Cells_[a-z]_11_13_2_", "")
    currentStep = re.sub("_it_\d+.txt", "", currentStep)
    currentStep = (int)(currentStep)/400
    
    currentLetter = filename.replace("summary_Cells_", "")
    #currentLetter = currentLetter.replace("_11_13_2_","")
    letter = re.sub(r"_01_2_\d+_\d+_it_\d+.txt", "", currentLetter)
    letterFinal = letter
    #[Total x, Cancer x, noncancer x, all x, sg x, igi x, aa x, it x, gu x]
    for piece in data:

        temp = piece.split()
        if temp[0] == "Total":
            value = [currentStep, temp[1], letter]
            totalCells.append(value)
            #total.append(temp[1])
        if temp[0] == "Cancer":
            value =[currentStep, temp[1], letter]
            cancerCells.append(value)
            #cancerCells.append(temp[1])
        if temp[0] == "sg":
            value = [currentStep, temp[1], letter]
            sgCells.append(value)
        if temp[0] == "igi":
            value = [currentStep, temp[1], letter]
            igiCells.append(value)
        if temp[0] == "aa":
            value = [currentStep, temp[1], letter]
            aaCells.append(value)
        if temp[0] == "it":
            value = [currentStep, temp[1], letter]
            itCells.append(value)
        if temp[0] == "a":
            value = [currentStep, temp[1], letter]
            aCells.append(value)
        #if temp[0] == "dr":
         #   value = [currentStep, temp[1], letter]
          #  drCells.append(value)
        #if temp[0] == "ai":
         #   value = [currentStep, temp[1], letter]
          #  aiCells.append(value)
        if temp[0] == "gu":
            value = [currentStep, temp[1], letter]
            guCells.append(value)
        if temp[0] == "ai":
            value = [currentStep, temp[1], letter]
            aiCells.append(value)            
        if temp[0] == "HealthyAlive":
            value = [currentStep, temp[1], letter]
            regularCells.append(value)


##Get a list of the final cancer values to use for stats test
finalCancerNumbers = []
for i in cancerCells:
    if i[0] == numReps:
        finalCancerNumbers.append(i[1])

print finalCancerNumbers

outString = ", ".join(str(e) for e in finalCancerNumbers) 
print finalCancerNumbers
outfileName = "..\\Stats2\\" + letterFinal + ".txt"
#Print both the final cancer values for all of the runs, and the average cancer values to use in the big graph
output = numpy.save(open(outfileName, 'w'), outString)


#Now total cells has what we want to plot.. so does Cancer cells
totalCells = numpy.array(totalCells)
cancerCells = numpy.array(cancerCells)
sgCells = numpy.array(sgCells)
igiCells = numpy.array(igiCells)
aaCells = numpy.array(aaCells)
itCells = numpy.array(itCells)
aCells = numpy.array(aCells)
guCells = numpy.array(guCells)
aiCells = numpy.array(aiCells)
regularCells = numpy.array(regularCells)

#pdb.set_trace()


totalCellCount = [0]*numReps
cancerCellCount =  [0]*numReps
sgCellCount =  [0]*numReps
igiCellCount =  [0]*numReps
aaCellCount =  [0]*numReps
itCellCount =  [0]*numReps
aCellCount =   [0]*numReps
guCellCount =  [0]*numReps
aiCellCount =   [0]*numReps
regularCellCount =  [0]*numReps

cancerCellNumbers = []

#Now get the totals for each run to get averages
for i in range(1, numReps):
    for row in totalCells:
        #print row[0]
        if row[0] == str(i): 
            #HERE
            totalCellCount[i] = totalCellCount[i] + int(row[1])
    for row in cancerCells:
        #print row[0]
        if row[0] == str(i): 
            #HERE
            cancerCellCount[i] = cancerCellCount[i] + int(row[1])
            if (i==15):
                cancerCellNumbers.append(row[1])
    for row in sgCells:
        #print row[0]
        if row[0] == str(i): 
            #HERE
            sgCellCount[i] = sgCellCount[i] + int(row[1])
    for row in igiCells:
        #print row[0]
        if row[0] == str(i): 
            #HERE
            igiCellCount[i] = igiCellCount[i] + int(row[1])
    for row in aaCells:
       # print row[0]
        if row[0] == str(i): 
            #HERE
            aaCellCount[i] = aaCellCount[i] + int(row[1])
    for row in itCells:
        #print row[0]
        if row[0] == str(i): 
            #HERE
            itCellCount[i] = itCellCount[i] + int(row[1])
    for row in aCells:
        #print row[0]
        if row[0] == str(i): 
            #HERE
            aCellCount[i] = aCellCount[i] + int(row[1])            
    for row in guCells:
        #print row[0]
        if row[0] == str(i): 
            #HERE
            guCellCount[i] = guCellCount[i] + int(row[1])
    for row in aiCells:
        #print row[0]
        if row[0] == str(i): 
            #HERE
            aiCellCount[i] = aiCellCount[i] + int(row[1])    
    for row in regularCells:
        #print row[0]
        if row[0] == str(i): 
            #HERE
            regularCellCount[i] = regularCellCount[i] + int(row[1])

totalCellCount = numpy.array(totalCellCount)
cancerCellCount = numpy.array(cancerCellCount)
sgCellCount = numpy.array(sgCellCount)
igiCellCount = numpy.array(igiCellCount)
aaCellCount = numpy.array(aaCellCount)
itCellCount = numpy.array(itCellCount)
aCellCount = numpy.array(aCellCount)
guCellCount = numpy.array(guCellCount)
aiCellCount = numpy.array(aiCellCount)
regularCellCount = numpy.array(regularCellCount)

#Calculate averages
avTotal = (totalCellCount/numRuns)
avCancer = (cancerCellCount/numRuns)
avSG = (sgCellCount/numRuns)
avIGI = (igiCellCount/numRuns)
avAA = (aaCellCount/numRuns)
avIT = (itCellCount/numRuns)
avA = (aCellCount/numRuns)
avGU = (guCellCount/numRuns)
avAI = (aiCellCount/numRuns)
avReg = (regularCellCount/numRuns)

#print avCancer

#Convert back into numpy arrays and delete first column
totalCells = numpy.array(avTotal)
totalCells = totalCells[1:]
cancerCells = numpy.array(avCancer)
cancerCells = cancerCells[1:]
sgCells = numpy.array(avSG)
sgCells = sgCells[1:]
igiCells = numpy.array(avIGI)
igiCells = igiCells[1:]
aaCells = numpy.array(avAA)
aaCells = aaCells[1:]
itCells = numpy.array(avIT)
itCells = itCells[1:]
aCells = numpy.array(avA)
aCells = aCells[1:]
guCells = numpy.array(avGU)
guCells = guCells[1:]
aiCells = numpy.array(avAI)
aiCells = aiCells[1:]
regularCells = numpy.array(avReg)
regularCells = regularCells[1:]

#pdb.set_trace()
#Jenna - I think we have the correct numbers at this point... now gotta fix the pictures

#pylab.ion()
#pylab.draw()
#pylab.show()
#matplotlib.pyplot.scatter(totalCells[:,0], totalCells[:,1], color="blue")
#Save the regular and cancer cell plot
runs = range(1, numReps)
runs = numpy.array(runs)
pylab.clf()
pylab.show()
pylab.draw()
c = matplotlib.pyplot.scatter(runs, cancerCells, color="red")
n = matplotlib.pyplot.scatter(runs, regularCells, color="blue")
pylab.legend((c, n), ('Cancer', 'Regular'), scatterpoints=1,loc='upper left', fontsize=22)
pylab.xlabel("Simulation step (step size = 400)", fontsize=22)
pylab.ylabel("Cell count", fontsize=22)
pylab.ylim(-100,3000)
pylab.xlim(0, 40)
pictureFileName = "..\\Stats2\\" + letter + "_cancerAndNon.png"
pylab.savefig(pictureFileName, dpi=800, figsize=(4.0, 3.0))

#plot all types
pylab.clf()
can1 = pylab.scatter(runs, cancerCells, color="red")
sg1 = pylab.scatter(runs, sgCells, color="darkgoldenrod", marker="+")
igi1 = matplotlib.pyplot.scatter(runs, igiCells, color="green", marker="^")
aa1 = pylab.scatter(runs, aaCells, color="magenta", marker="v")
it1 = pylab.scatter(runs, itCells, color="purple", marker="<")
a1 = pylab.scatter(runs, aCells, color="darkblue", marker=">")
#dr1 = pylab.scatter(drCells[:,0], drCells[:,1], color="magenta")
gu1 = pylab.scatter(runs, guCells, color="black", marker="8")
ai1 = pylab.scatter(runs, aiCells, color="saddlebrown", marker="s")
pylab.legend((can1, sg1, igi1, aa1, it1, a1, gu1, ai1), ( 'All cancer cells', 'Self growth', 'Ignores growth inhibition','Avoids apoptosis', 'Ignores telomeres', 
    'Angiogenesis', 'Genetical unstable', 'Avoids immunity'), loc='upper left', fontsize=18)
pylab.xlabel("Simulation step (step size = 400)", fontsize=22)
pylab.ylabel("Number of cells with \n hallmark activated", fontsize=22)
pylab.ylim(-100,3000)
pylab.xlim(0, 40)
pictureFileName2 = "..\\Stats2\\" + letter + "_allMut.png"
pylab.savefig(pictureFileName2, dpi=800, figsize=(3.0, 1.4))


outfileName = "..\\Stats2\\" + letterFinal + "_averagePerStep.txt"
outString = numpy.array_str(cancerCells)
print outString
#Print both the final cancer values for all of the runs, and the average cancer values to use in the big graph
output = numpy.save(open(outfileName, 'w'), outString)