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


#Order of phenotypes in binary lowest to highest bit
one = "avoids immunity"
two = "genome unstable"
three = "angiogenesis"
four = "ignores telomere"
five = "avoids apoptosis"
six = "ignores growth inhibition"
seven = "self grows"
eight = "dead" 

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1
    
def make_r(x):
    x = int(x)
    if (x >= 128): #Cell is dead so colour black
        return 0
    if (x == 0): 
        return 100
    if (x*4 < 255):
        return x*4
    if (x*3 < 255):
        return x*3
    return x*2    
    
def make_g(x):
    x = int(x) 
    if (x >= 128): #Cell is dead so colour black
        return 0
    if (x == 0): #cell is healthy and alive 
        return 149
    else:
        g = x*3
        if (g > 255):
            g = g/4
            return g
        return g
        
def make_b(x):
    x = int(x)
    if (x >= 128): #Cell is dead so colour black
        return 0
    if (x == 0): #cell is healthy and alive 
        return 237
    return x/2
            
    
#def makeColourMap(mutInts):
#    colours={}
#    for x in range(0, len(mutInts)):
#         colours[str(mutInts[x])] = rgb_to_hex(x)
#    return colours

def get_bit(byteval,idx):
    return ((byteval&(1<<idx))!=0)
    

def convert_to_hallmarks(mutBin):
    '''Function to convert a binary number to a list of hallmarks
    The hallsmarks are written out in this order:
        
    selfGrows(), ignoresGrowthInhibition, avoidsApoptosis, ignoresTelomere, angiogenesis, genomeUnstable, avoidsImmunity
    I will check each bit to see if set'''
    pheno = "Mutations: "
    if (type(mutBin) != int):
        mutBin = int(mutBin, 2)
    if (get_bit(mutBin, 0) == 1):
        pheno = pheno + one + ", "
    if (get_bit(mutBin, 1) == 1):
        pheno = pheno + two + ", "
    if (get_bit(mutBin, 2) == 1):
        pheno = pheno + three + ", "
    if (get_bit(mutBin, 3) == 1):
        pheno = pheno + four + ", "
    if (get_bit(mutBin, 4) == 1):
        pheno = pheno + five + ", "
    if (get_bit(mutBin, 5) == 1):
        pheno = pheno + six + ", "
    if (get_bit(mutBin, 6) == 1):
        pheno = pheno + seven + ", "
    if (get_bit(mutBin, 7) == 1):
        pheno = pheno + eight + ", "
    if (pheno == "Mutations: dead"):
        return "Dead"
    if (pheno == "Mutations: "):
        return "Normal"
    else:
        return pheno
        
def decode_phenos(phenoNumbers):
    output = "Phenotypes by number: \n"
    for x in phenoNumbers:
        output = output + "Pheno number: " + x + " is " + convert_to_hallmarks(x) + " and occurs " + str(phenoNumbers[x]) + " times \n"
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
fileTime = "AveragePhenos" + str(now.month) + "_" + str(now.day) + "_" + str(now.year) + "_" + str(now.hour)+ "_" + str(now.minute)
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

mutTotals = {}

pdb.set_trace()
for filename in os.listdir (folder):
    total = file_len(filename)-1;

    #Create file name
    #eg infilename = 'Cells_g_11_13_2_end.txt'
    infilename = filename
    #outfilename = 'summary_11_' + run + '.txt'
    outfilename = newLocation + 'summary_' + filename;
    #picturefilename = 'tumour_' + run + '.jpg'

    
    #Load in data as just a list of mutation numbers and the number of times they occur
    mut, occurance = numpy.loadtxt(infilename, dtype=int, unpack=True)
   # pdb.set_trace()
    # pylab.show()
    # pylab.ion()
    
    #Go through all of the mutations for this particular file
    for x in range(0, total+1):
        #Get the current mutation int
        currentMut = mut[x]
        #if we have already seen it, just add the number of occurances to the list
        if currentMut in mutTotals:
            mutTotals[currentMut].append(occurance[x])
        # if we haven't seen it, start a new list for that number
        else:
            mutTotals[currentMut] = [(occurance[x])]
            
#once we have been through every file, mutTotals will have lists for containing the number of occurances for each mutation for each run of the simulator
# We want to get the averages
pdb.set_trace()
allkeys = mutTotals.keys()
averages = ""
for mut in allkeys:
    length = len(mutTotals[mut])
    sumOfOccurances = sum(mutTotals[mut])
    average = sumOfOccurances/length
    averages = averages + "mut " + str(mut) + + " which is: " + convert_to_hallmarks(mut) + " occurs " + str(length) + " times out of the runs and is present in " + str(average) + " number of cells each run on average \n" 

#write out the averages string
outfilename = newLocation + 'averagePhenos' + filename;
output = numpy.save(open(outfilename, 'w'), averages);