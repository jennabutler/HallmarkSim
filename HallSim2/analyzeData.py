import os
import csv
import numpy
import matplotlib
import matplotlib.pylab as pylab
import matplotlib.pyplot
import pdb
import datetime

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1
    
def make_colour(x):
    x = int(x)   
    r=0
    g=0
    b=0
    if (x >= 128): #Cells 128> are dead, so colour them black
        return "#000000"
    if (x == 0): #Cells with all 0 are alive and noncancerous, so they are a nice bright blue
        return "#6495ed"
    if (x == 1):
        return "#FF00FF"
    if (x == 2):
        return "#800000"
    if (x == 3):
        return "#BA55D3"
    r = x
    if ((x*2) < 255):
        g = x*2
    b = x/2
    if (b < 16):
        b = 17
    if (r < 16):
        r = r * 10
    if (g < 16):
        g = g * 8 + (b) + (r/5)     
    r = hex(r)
    g = hex(g)
    b = hex(b)
    newColour = "#" + r[2:] + g[2:] + b[2:]
    return newColour
    
def makeColourMap(mutInts):
    colours={}
    for x in range(0, len(mutInts)):
        colours[str(mutInts[x])] = make_colour(mutInts[x])
    return colours
         

#get current date time for the time-stamped folder
now = datetime.datetime.now()

#Sharcnet files
#newLocation = "/work/jcamer7/HallSim/Output" + str(now.month) + "_" + str(now.day) + "_" + str(now.year) + "_" + str(now.hour)+ "_" + str(now.minute)
#os.makedirs(newLocation)
#newLocation = newLocation + "/"
#os.chdir("/work/jcamer7/HallSim/Output")
#folder = "/work/jcamer7/HallSim/Output"

#Local home files
newLocation = "C:\\Users\\Jenna\\Documents\\Visual Studio 2012\\Projects\\HallSim2\\Output\\" + str(now.month) + "_" + str(now.day) + "_" + str(now.year) + "_" + str(now.hour)+ "_" + str(now.minute)
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
    #infilename = 'Cells_g_11_13_2_end.txt'
    infilename = filename
    #outfilename = 'summary_11_' + run + '.txt'
    outfilename = newLocation + 'summary_' + filename;
    #picturefilename = 'tumour_' + run + '.jpg'

    
    #Load in data
    #istream << "i, j, state, mut, sg, igi, aa, it, gu \n";
    i, j, state, mut, sg, igi, aa, it, a, gu, ai = numpy.loadtxt(infilename, dtype=int, unpack=True)
   # pdb.set_trace()
    # pylab.show()
    # pylab.ion()
    
    #Draw the whole tumour
    #matplotlib.pyplot.scatter(i,j)
    
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

    
    
    
    for x in range(0, total):
        mutValue =  "0b" + str(state[x]) + str(sg[x]) + str(igi[x]) + str(aa[x]) + str(it[x]) + str(a[x]) + str(gu[x]) + str(ai[x])
        mutInt = int(mutValue, 2)
        #mutInt is a unique integer for every phenotype...
        #mutInt 0 is an alive normal cell
        #mutInt of 1 -> 127 means it is an alive cancer cell
        #mutInt of 128->255 is dead
        mutInts.append(str(mutInt))
#        
#        
        if a[x] == 1 and igi[x] == 1 and state[x] == 0:
            a_igi_Cellsi.append(i[x])
            a_igi_Cellsj.append(j[x])
        if a[x] == 1 and aa[x] == 1 and state[x] == 0:
            a_aa_Cellsi.append(i[x])
            a_aa_Cellsj.append(j[x])
        if sg[x] == 1 and state[x] == 0:
            sg_count += 1
        if sg[x] == 1 and igi[x] == 1 and aa[x] ==1 and it[x] == 1 and gu[x] == 1 and state[x] == 0:
            allMutations += 1
        if igi[x] == 1 and state[x] == 0:
            igi_count += 1
            igiCellsi.append(i[x])
            igiCellsj.append(j[x])
        if aa[x] == 1 and state[x] == 0:
            aa_count +=1
            aaCellsi.append(i[x])
            aaCellsj.append(j[x])
        if it[x] == 1 and state[x] == 0:
            it_count += 1
        if a[x] == 1 and state[x] == 0:
            a_count += 1
            angioCellsi.append(i[x]);
            angioCellsj.append(j[x]);
        #if dr[x] == 1:
         #   dr_count+=1
        #if ai[x] == 1 and state[x] == 0:
         #   ai_count += 1
        if gu[x] == 1 and state[x] == 0:
            gu_count += 1
        if ai[x] == 1 and state[x] == 0:
            ai_count += 1
        if mut[x] == 1 and state[x] == 0:
            cancer += 1
            cancerCellsi.append(i[x])
            cancerCellsj.append(j[x])
        if mut[x] == 0 and state[x] == 0:
            regularCellsi.append(i[x])
            regularCellsj.append(j[x])
        if state[x] == 0:
            aliveCellsi.append(i[x])
            aliveCellsj.append(j[x])
            alive += 1
        if state[x] == 1:
            deadCellsi.append(i[x])
            deadCellsj.append(j[x])
            dead += 1

    
    import re        
    currentStep = re.sub("Cells_[a-z]_01_2_1_", "", filename)
    
    #Cells_a_01_2_1_4000_it_2

    currentStep = re.sub("_it_\d+.txt", "", currentStep)
    currentStep = (int)(currentStep)/2000
    
    matplotlib.pyplot.clf()
    #matplotlib.pyplot.scatter(cancerCellsi, cancerCellsj, color="red", alpha=0.6, s=9)
    #matplotlib.pyplot.scatter(regularCellsi, regularCellsj, color="blue", alpha=0.6, s=9)
    #matplotlib.pyplot.scatter(deadCellsi, deadCellsj, color="black", alpha=0.6, s=9)
    #matplotlib.pyplot.scatter(angioCellsi, angioCellsj, color="green", alpha=0.6, s=9)
    #matplotlib.pyplot.scatter(a_aa_Cellsi, a_aa_Cellsj, color="yellow", alpha=0.6, s=9)
    #matplotlib.pyplot.scatter(a_igi_Cellsi, a_igi_Cellsj, color="cyan", alpha=0.6, s=9)
 
    pylab.xlim(0, 300) #2400,2600
    pylab.ylim(0, 300)
    
    ##Scatter where every phenotype is coloured
    ##Mark do you like this? How can I make it more informative?
    matplotlib.pyplot.clf()
    use_colours = makeColourMap(mutInts) #I wanted to pass in a set to save computation time but they aren't indexible?
    #use_colours = {'5': "MediumBlue", "9":"LightSalmon", "6": "Fuchsia", "32": "red", "64": "green", "8": "yellow", "0": "MistyRose", "3": "DarkGreen", "4": "MediumOrchid", "52": "cyan", "1": "red", "2": "Gray", "16": "Brown", "17": "Coral", "2": "Cornsilk", "20":"Crimson", "18":"DarkGoldenRod", "12":"LightGray"}
    #MARK here is where it gets bizare
    #If you are tracing it, you can print set(mutInts) and see all of the unique phenotypes
    #You can then print use_colours and see that every phenotype has an associated colour mapped to it
    #Then if you continue you will get an error for key 3 or key 5, neither of which are phenotypes
    pdb.set_trace()
    matplotlib.pyplot.scatter(i,j,c=[use_colours[x[0]] for x in mutInts],s=9, lw=0)
    matplotlib.pyplot.scatter(deadCellsi, deadCellsj, color="black", s=9, lw=0)

    picturefilename = newLocation + 'picture_' + filename[:-4] + '.jpg'
    matplotlib.pyplot.savefig(picturefilename, dpi=150)
        
            

   # 
    regCells = total - cancer
   #         
    total = "Total " + (str)(total) + "\n"
    cancer = "Cancer " + (str)(cancer) + "\n"
    reg = "Noncancer " + (str)(regCells) + "\n"
    aliveString = "Alive " + (str)(alive) + "\n"
    deadString = "Dead " + (str)(dead) + "\n"
    allMut = "all "+ (str)(allMutations) + "\n"
    first = "sg "+ (str)(sg_count) + "\n"
    second = "igi " + (str)(igi_count) + "\n"
    third = "aa " + (str)(aa_count) + "\n"
    fourth = "it " + (str)(it_count) +"\n"
    fifth = "a " + (str)(a_count) + "\n"
    #fifth = "dr " + (str)(dr_count) + "\n"
    seventh = "gu " + (str)(gu_count) + "\n"
    eigth = "ai " + (str)(ai_count) + "\n"
    uniquePhenos = "Unique phenotypes: " + str(set(mutInts)) + "\n"

    summary = total + cancer + aliveString + deadString + reg + allMut + first + second + third + fourth + fifth + seventh + eigth + uniquePhenos
   # 
    output = numpy.save(open(outfilename, 'w'), summary);



         

            