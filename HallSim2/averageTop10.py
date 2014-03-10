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

os.chdir("C:\Users\Jenna\Documents\Visual Studio 2012\Projects\HallSim2\Output\Working")

#Edit these each run
numRuns = 5
numReps = 15

import sys
import os.path

#pdb.set_trace()

topPhenos = {}

cancerCells = []
regularCells = []
data = []
totalCells = []
sgCells = []
igiCells = []
aaCells = []
itCells = []
drCells = []
aiCells = []
guCells = []
letterFinal = 'x'

folder = "C:\Users\Jenna\Documents\Visual Studio 2012\Projects\HallSim2\Output\Working"
for filename in os.listdir (folder):


    #Load in data
    data[:] = []
    with open(filename) as f:
        next(f)
        for line in f:
            data.append(line)
    
   
    #currentStep = re.sub("summary_Cells_[a-z]_11_18_1_", "", filename)
    currentStep = re.sub("Top10_Cells_[a-z]_01_2_\d+_", "", filename)
    pdb.set_trace()
    
    #currentStep = filename.replace("summary_Cells_[a-z]_11_13_2_", "")
    currentStep = re.sub("_it_\d+.txt", "", currentStep)
    currentStep = (int)(currentStep)/400
    
    currentLetter = filename.replace("Top10_Cells_", "")
    #currentLetter = currentLetter.replace("_11_13_2_","")
    letter = re.sub(r"_01_2_\d+_\d+_it_\d+.txt", "", currentLetter)
    letterFinal = letter
    #[Total x, Cancer x, noncancer x, all x, sg x, igi x, aa x, it x, gu x]
    for piece in data:
        pdb.set_trace()
        #pdb.set_trace()
             
        #Add each of the top 10 phenos to a dictionary 
        #If it has been seen in another file, we add the total number as well as the number of files it has been seen in
        #It looks like this:
        #key: phenoInt value: tuple (total number of times, total number of files)
        #so 128, (2134, 3)
        #means the phenotype 128 has been in the top 10 of 3 runs and has appeared a total of 2134 times
        temp = piece.split()
        if (temp[0] in topPhenos):
            currentNumber = topPhenos[temp[0]][1]
            topPhenos[temp[0]]+=(temp[1], currentNumber+1)
        else:
            topPhenos[temp[0]] = (temp[1], 0)
