#Script to count the number of occurances of each type of mutation for last file of every run
import os
import csv
import numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as pylab
import matplotlib.pyplot as plt
import pdb
import datetime
import scipy.spatial.distance
import numpy as np


def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1



def get_bit(byteval,idx):
    return ((byteval&(1<<idx))!=0)
    

def decode_phenos(phenoNumbers):
    output = ""
    for x in phenoNumbers:
        output = output + str(x) + " " + str(phenoNumbers[x]) + "\n"
    return output
        

    
#get current date time for the time-stamped folder
now = datetime.datetime.now()

#Sharcnet files
#newLocation = "/work/jcamer7/HallSim/Output" + str(now.month) + "_" + str(now.day) + "_" + str(now.year) + "_" + str(now.hour)+ "_" + str(now.minute)
#os.makedirs(newLocation)
#newLocation = newLocation + "/"
#os.chdir("/work/jcamer7/HallSim/Output")
#folder = "/work/jcamer7/HallSim/Output"

#Local home files
fileTime = "PhenoCount_" + str(now.month) + "_" + str(now.day) + "_" + str(now.year) + "_" + str(now.hour)+ "_" + str(now.minute)
newLocation = "C:\\Users\\Jenna\\Documents\\Visual Studio 2012\\Projects\\HallSim2\\Output\\" + fileTime
#newLocation = "C:\\Users\\Jenna\\Documents\\Visual Studio 2012\\Projects\\HallSim2\\Output\\" + str(now.month) + "_" + str(now.day) + "_" + str(now.year) + "_" + str(now.hour)+ "_" + str(now.minute)
os.chdir("C:\Users\Jenna\Documents\Visual Studio 2012\Projects\HallSim2\Output\Working") #Changes to the director with the input files
os.makedirs(newLocation) #Make the directory for the output
newLocation = newLocation + "\\"
folder = "C:\Users\Jenna\Documents\Visual Studio 2012\Projects\HallSim2\Output\Working"

#Local work files
#newLocation = "/Users/jenna/Research/HallSim2/" + str(now.month) + "_" + str(now.day) + "_" + str(now.year) + "_" + str(now.hour)+ "_" + str(now.minute)
#os.makedirs(newLocation)
#newLocation = newLocation + "/"
#os.chdir("/Users/jenna/Research/HallSim2/Output")
#folder = "/Users/jenna/Research/HallSim2/Output")

import sys
import os.path

#pdb.set_trace()

for filename in os.listdir (folder):
    total = file_len(filename)-1;

    #Create file name
    #eg infilename = 'Cells_g_11_13_2_end.txt'
    infilename = filename
    #outfilename = 'summary_11_' + run + '.txt'
    outfilename = newLocation + 'phenoCounts_' + filename;
    #picturefilename = 'tumour_' + run + '.jpg'

    
    #Load in data
    i, j, state, mut, sg, igi, aa, it, a, gu, ai = numpy.loadtxt(infilename, dtype=int, unpack=True)
   # pdb.set_trace()
    # pylab.show()
    # pylab.ion()
    
    #Variables
    allMutations = 0
    sg_count = 0
    igi_count = 0
    aa_count = 0
    it_count = 0
    a_count = 0
    #dr_count = 0
    #ai_count = 0
    gu_count = 0
    ai_count = 0
    cancer = 0;
    cancerCellsi = []
    cancerCellsj = []
    regularCellsi = []
    regularCellsj = []
    aliveCellsi = []
    aliveCellsj = []
    deadCellsi = []
    deadCellsj = []
    angioCellsi = []
    angioCellsj = []
    igiCellsi = []
    igiCellsj = []
    aaCellsi = []
    aaCellsj = []
    a_igi_Cellsi = []
    a_igi_Cellsj = []
    a_aa_Cellsi = []
    a_aa_Cellsj = []
    alive = 0
    dead = 0
    mutInts = []
    phenotypes = {}
    phenoIJ = {}
    
    
    for x in range(0, total):
        # get the binary value of the mutations
        mutValue =  "0b" + str(state[x]) + str(sg[x]) + str(igi[x]) + str(aa[x]) + str(it[x]) + str(a[x]) + str(gu[x]) + str(ai[x])
        #Convert into integer
        mutInt = int(mutValue, 2)
        #Keep track of how many of each phenotype are present
        if (mutInt in phenotypes):
            phenotypes[mutInt]+=1
        else:
            phenotypes[mutInt] = 1

            
        #Attempting to get a list of (i,j) points for every hallmark phenotype, that I can then do spatial analysis on
        if (mutInt not in phenoIJ):
            phenoIJ[mutInt] = [(i[x], j[x])]
        if (mutInt in phenoIJ):
            phenoIJ[mutInt].append((i[x], j[x]))
            
        #mutInt is a unique integer for every phenotype...
        #Keep track of each cells phenotype
        mutInts.append(str(mutInt))

      
    import re        
    currentStep = re.sub("Cells_[a-z]_01_2_1_", "", filename)
    
    #Cells_a_01_2_1_4000_it_2

    currentStep = re.sub("_it_\d+.txt", "", currentStep)
    currentStep = (int)(currentStep)/2000
    
    phenoCounts = decode_phenos(phenotypes)
    output = numpy.save(open(outfilename, 'w'), phenoCounts);