import os
import csv
import numpy
import matplotlib
import matplotlib.pylab as pylab
import matplotlib.pyplot
import pdb
import datetime
#i, j, state, quis, nec, apop, agg1, agg2, agg3, al, mut, sg, igi, aa, it, a, gu, ai, gp = numpy.loadtxt(infilename, dtype=int, unpack=True)

#Order of phenotypes in binary lowest to highest bit
#Has been updated on April 19
one = "glyco phenotype"
two = "avoids immunity"
three = "genome unstable"
four = "angiogenesis"
five = "ignores telomere"
six = "avoids apoptosis"
seven = "ignores growth inhibition"
eight = "self grows"
nine = "dead" 
ten = "mut"
eleven = "alive"
twelve = "agg3"
thirteen = "agg2"
fourteen = "agg1"
fifteen = "apop"
sixteen = "nec"
seventeen = "quis"

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1
    
def make_r(x):
    x = int(x)
    if (x < 128): #Cell is dead so colour black #April 19 had to swap this.. now the 1 is alive.. so 128 or greater is alive
        return 0
    if (x == 128): 
        return 100
    if (x*4 < 255):
        return x*4
    if (x*3 < 255):
        return x*3
    return x*2    
    
def make_g(x):
    x = int(x) 
    if (x < 128): #Cell is dead so colour black
        return 0
    if (x == 128): #cell is healthy and alive 
        return 149
    else:
        g = x*3
        if (g > 255):
            g = g/4
            return g
        return g
        
def make_b(x):
    x = int(x)
    if (x < 128): #Cell is dead so colour black
        return 0
    if (x == 128): #cell is healthy and alive 
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
    #Updated for glyco pheno April 19th
    if (get_bit(mutBin, 8) == 1):
        pheno = pheno + nine + ","
    if (pheno == "Mutations: dead"):
        return "Dead"
    if (pheno == "Mutations: "):
        return "Normal"
    else:
        return pheno
        
def decode_phenos(phenoNumbers):
    output = "Phenotypes by number: \n"
    for x in phenoNumbers:
        output = output + "Pheno number: " + str(x) + " is " + convert_to_hallmarks(x) + " and occurs " + str(phenoNumbers[x]) + " times \n"
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
fileTime = str(now.month) + "_" + str(now.day) + "_" + str(now.year) + "_" + str(now.hour)+ "_" + str(now.minute)
newLocation = "C:\\Users\\Jenna\\Documents\\Visual Studio 2012\\Projects\\HallSim2\\Output\\" + fileTime
#newLocation = "C:\\Users\\Jenna\\Documents\\Visual Studio 2012\\Projects\\HallSim2\\Output\\" + str(now.month) + "_" + str(now.day) + "_" + str(now.year) + "_" + str(now.hour)+ "_" + str(now.minute)
os.chdir("C:\Users\Jenna\Documents\Visual Studio 2012\Projects\HallSim2\Output\Working2") #Working3    #Changes to the director with the input files
os.makedirs(newLocation) #Make the directory for the output
newLocation = newLocation + "\\"
folder = "C:\Users\Jenna\Documents\Visual Studio 2012\Projects\HallSim2\Output\Working2"

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
    outfilename = newLocation + 'summary_' + filename;
    #picturefilename = 'tumour_' + run + '.jpg'

    
    #Load in data
    i, j, state, quis, nec, apop, agg1, agg2, agg3, al, mut, sg, igi, aa, it, a, gu, ai, gp = numpy.loadtxt(infilename, dtype=int, unpack=True)
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
    numQuis = 0
    numApop = 0
    numNec = 0
    imdead = 0
    numGly = 0
    numAgg1 = 0
    numAgg2 = 0
    numAgg3 = 0
    mutInts = []
    phenotypes = {}
    phenoHalls = {}
    regAlive = 0
    
    #State can be:
    #	ALIVE, DEAD, QUIS, APOP, NEC, IM_DEAD, AGG1, AGG2, AGG3, GLY
    #     0     1      2    3     4      5      6      7     8    9
    
    #Need to add new arrays to keep track of 
    # Glyco cells
    # Quiescent cells
    # And then different dead cells 
    # April 19th.. haven't done this yet
    # Also watch out a few lines down.. the mut value includes the state which is no longer just 0 or 1... might have to do something about that
    
    for x in range(0, total):
        # get the binary value of the mutations
        #Al is the new state... either 0 or 1... 1 is alive, 0 is dead
        mutValue =  "0b" + str(al[x]) + str(sg[x]) + str(igi[x]) + str(aa[x]) + str(it[x]) + str(a[x]) + str(gu[x]) + str(ai[x])
        #i, j, state, quis, nec, apop, agg1, agg2, agg3, al, mut, sg, igi, aa, it, a, gu, ai, gp
        #mutValue = "0b" + str(quis[x]) + str(nec[x]) + str(al[x]) + str(sg[x]) + str(igi[x]) + str(aa[x]) + str(it[x]) + str(a[x]) + str(gu[x]) + str(ai[x])
        #Convert into integer
        mutInt = int(mutValue, 2)
        #Keep track of how many of each phenotype are present
        if (mutInt in phenotypes):
            phenotypes[mutInt]+=1
        else:
            phenotypes[mutInt] = 1
            
        #Keep track of what each phenotype represents
        if (mutInt not in phenoHalls):
            phenoHalls[mutInt] = convert_to_hallmarks(mutValue)
            
        #mutInt is a unique integer for every phenotype...
        #Keep track of each cells phenotype
        mutInts.append(str(mutInt))

        #Tally various combinations and cell types
        #April 19
        if state[x] == 2:
            numQuis += 1
        if state[x] == 3:
            numApop += 1
        if state[x] == 4:
            numNec += 1
        if state[x] == 5:
            imdead += 1
        if state[x] == 9:
            numGly += 1
        if state[x] == 6:
            numAgg1 += 1
        if state[x] == 7:
            numAgg2 += 1
        if state[x] == 8:
            numAgg3 += 1
        if a[x] == 1 and igi[x] == 1 and state[x] == 0:
            a_igi_Cellsi.append(i[x])
            a_igi_Cellsj.append(j[x])
        if a[x] == 1 and aa[x] == 1 and state[x] == 0:
            a_aa_Cellsi.append(i[x])
            a_aa_Cellsj.append(j[x])
        if sg[x] == 1 and state[x] == 0:
            sg_count += 1
        #if sg[x] == 1 and igi[x] == 1 and aa[x] ==1 and it[x] == 1 and gu[x] == 1 and state[x] == 0:
         #   allMutations += 1
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
            #April 19
        if mut[x] == 1 and (state[x] == 0 or state[x] == 6 or state[x] == 7 or state[x] == 8 or state[x] == 9):
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
        if state[x] == 0 and mut[x] == 0:
            regAlive += 1
    

    import re        
    currentStep = re.sub("Cells_[a-z]_01_2_1_", "", filename)
    
    #Cells_a_01_2_1_4000_it_2


    #Attempting to order stuff
    import operator
    #This creates a list of phenotype tuples that sorted smallest to end.. so the last 10 are the top 10 phenotypes
    #Print these out into their own file
    sorted_phenos = sorted(phenotypes.iteritems(), key=operator.itemgetter(1)) 
    #sorted phenos = [(phenoInt, numberOfOccurances)]
    last = len(sorted_phenos)
    if (last > 10):
        biggestPheno1 = str(sorted_phenos[last-1][0]) + " " + str(sorted_phenos[last-1][1]) + "\n"
        biggestPheno2 = str(sorted_phenos[last-2][0]) + " " + str(sorted_phenos[last-1][1]) + "\n"
        biggestPheno3 = str(sorted_phenos[last-3][0]) + " " + str(sorted_phenos[last-1][1]) + "\n"
        biggestPheno4 = str(sorted_phenos[last-4][0]) + " " + str(sorted_phenos[last-1][1]) + "\n"
        biggestPheno5 = str(sorted_phenos[last-5][0]) + " " + str(sorted_phenos[last-1][1]) + "\n"
        biggestPheno6 = str(sorted_phenos[last-6][0]) + " " + str(sorted_phenos[last-1][1]) + "\n"
        biggestPheno7 = str(sorted_phenos[last-7][0]) + " " + str(sorted_phenos[last-1][1]) + "\n"
        biggestPheno8 = str(sorted_phenos[last-8][0]) + " " + str(sorted_phenos[last-1][1]) + "\n"
        biggestPheno9 = str(sorted_phenos[last-9][0]) + " " + str(sorted_phenos[last-1][1]) + "\n"
        biggestPheno10 = str(sorted_phenos[last-10][0]) + " " + str(sorted_phenos[last-1][1]) + "\n"
        top10Phenos = biggestPheno1 + biggestPheno2 + biggestPheno3 + biggestPheno4 + biggestPheno5 + biggestPheno6 + biggestPheno7 + biggestPheno8 + biggestPheno9 + biggestPheno10
        outfilenameTop = newLocation + 'Top10_' + filename;
        numpy.save(open(outfilenameTop, 'w'), top10Phenos);


       
    currentStep = re.sub("_it_\d+.txt", "", currentStep)
    currentStep = (int)(currentStep)/2000
 
    import Image as im 
    import ImageFilter as imf 
    frame = im.new('RGBA', (300,300), "white") 
    c = frame.load() 
    for m in range(0,pylab.size(mutInts)): 
        #Convert the mutation into a colour
        r = make_r(mutInts[m])
        g = make_g(mutInts[m])
        b = make_b(mutInts[m])
        if (type(r) != int):
            r = 150
            print "error: returned r was not an int.. subbing in 150"
        if (type(g) != int):
            g = 150
            print "error: returned g was not an int.. subbing in 150"
        if (type(b) != int):
            b = 150
            print "error: returned b was not an int.. subbing in 150"
        c[int(i[m]),300-int(j[m])] = (r, g, b, 255)

    print frame.size 
    frame = frame.resize([int(3*s) for s in frame.size]) 
    frame = frame.filter(imf.SMOOTH_MORE) 
    picturefilename = newLocation + 'picture_' + filename[:-4] + '.jpg'
    #Save the picture 
    frame.save(picturefilename)
       
            

#   # 
    regCells = total - cancer
#   #         
    total = "Total " + (str)(total) + "\n"
    cancer = "Cancer " + (str)(cancer) + "\n"
    reg = "Noncancer " + (str)(regCells) + "\n"
    aliveString = "Alive " + (str)(alive) + "\n"
    regAliveString = "HealthyAlive " + (str)(regAlive) + "\n"
    deadString = "Random Dead " + (str)(dead) + "\n"
    apopString = "Apoptosis " + (str)(apop) + "\n"
    quisString = "Quis " + (str)(quis) + "\n"
    necString = "Necrotic " + (str)(nec) + "\n"
    imString = "Immunity death " + (str)(imdead) + "\n"
    glyString = "Glyco phenotype " + (str)(gly) + "\n"
    agg1String = "Agg1 " + (str)(agg1) + "\n"
    agg2String = "Agg2 " + (str)(agg2) + "\n"
    agg3String = "Agg3 " + (str)(agg3) + "\n"

    unique = "Unique phenotypes: "+ str(set(mutInts)) + "\n"
    breakdown = decode_phenos(phenotypes)
    first = "sg "+ (str)(sg_count) + "\n"
    second = "igi " + (str)(igi_count) + "\n"
    third = "aa " + (str)(aa_count) + "\n"
    fourth = "it " + (str)(it_count) +"\n"
    fifth = "a " + (str)(a_count) + "\n"
    #fifth = "dr " + (str)(dr_count) + "\n"
    seventh = "gu " + (str)(gu_count) + "\n"
    eigth = "ai " + (str)(ai_count) + "\n"
    
    #Write out all of the data we are interested in
    summary = total + cancer + apopString + quisString + necString + imString + glyString + agg1String + agg2String + agg3String + aliveString + regAliveString + deadString + reg + unique + breakdown + "Individual: " + first + second + third + fourth + fifth + seventh + eigth
   


    output = numpy.save(open(outfilename, 'w'), summary);


    


#return fileTime? For putting the parameters in the right folder?
         

            