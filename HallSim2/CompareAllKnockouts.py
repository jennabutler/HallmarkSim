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

numReps = 37
allCancer = [19, 41, 66, 91, 118, 146, 169, 199, 214, 244, 260, 291, 304, 337, 361,
377, 411, 420, 450, 478, 475, 500, 483, 493, 492, 529, 614, 748, 910, 1075,
1246, 1422, 1597, 1777, 1954, 2133]
combo_b = [8, 22, 36, 53, 74, 94, 111, 134, 150, 176, 194, 223, 236, 265, 291,
301, 332, 345, 356, 376, 362, 376, 353, 337, 338, 419, 561, 732, 917, 1114,
1312, 1510, 1711, 1912, 2113, 2311]
combo_c = [18, 39, 59, 78, 100, 114, 129, 151, 163, 182, 197, 211, 219, 229, 246,
259, 271, 280, 298, 300, 304, 325, 333, 325, 355, 437, 534, 667, 791, 938, 1048, 1177, 1309, 1471, 1619, 1781]
combo_d = [12,  31,  50,   71,   95, 121,  144,  176,  197,  233,  257,  292,  311,  351,  382,
404,  448,  460,  486,  514,  499,  516,  483,  447,  385,  449,  587,  766,  956, 1145,
1331, 1516, 1706, 1894, 2083, 2271]
combo_v = [16,36,60,82,105,130,151,181,206,237,266,297,320,347,376
,405,437,457,493,507,509,531,495,466,441,459,524,654,768,920
,1057,1220,1371,1538,1709,1877]
combo_s = [8,22,37,55,71,92,111,139,156,180,193,220,233,255,276
,291,315,321,335,354,353,363,341,314,307,349,433,542,669,824
,994,1159,1333,1507,1681,1852]
combo_p = [14,30,46,61,78,94,111,129,138,157,168,189,201,219,240
,240,262,270,283,299,295,303,289,275,269,355,492,656,823,988
,1156,1331,1508,1682,1863,2039]
combo_o = [16,40,64,88,112,134,153,175,189,210,221,242,254,269,287
,295,314,323,332,342,333,336,327,311,268,325,433,568,717,871
,1025,1185,1349,1518,1683,1845]
combo_m = [10,27,49,71,93,114,133,154,164,183,194,217,228,247,260
,269,290,298,307,321,310,323,306,310,309,340,418,536,665,797
,944,1089,1236,1382,1531,1679]
combo_i = [11,26,43,64,86,109,129,159,181,213,238,271,290,322,347
,369,403,415,436,459,437,450,418,366,379,476,618,785,952,1126
,1303,1482,1666,1850,2037,2225]
combo_e = [14,36,64,92,120,151,177,205,226,244,264,291,302,323,339,357,380,391,397,408,414,402,352,273,83,47,15,9,0,0,0,0,0,0,0,0]
combo_f = [15,35,58,79,100,123,141,168,180,208,213,244,254,274,292,297,329,330,337,352,337,332,279,148,44,22,0,0,0,0,0,0,0,0,0,0]
combo_g = [4,11,21,35,49,65,81,100,116,137,146,168,185,199,217,224,247,259
,260,273,261,252,215,153,47,26,26,11,0,0,0,0,0,0,0,0]
combo_h = [16,34,54,75,96,117,130,149,160,178,187,207,208,230,240
,243,261,262,279,288,286,290,282,271,292,381,522,673,833,996
,1156,1324,1489,1655,1824,1992]
combo_j  = [17,41,64,89,119,145,168,202,220,256,265,304,308,341,363,372,410,392
,415,412,386,319,219,32,0,0,0,0,0,0,0,0,0,0,0,0]
combo_k = [14,34,60,89,121,150,175,207,230,263,282,324,338,380,412
,433,471,474,510,535,518,523,489,433,451,539,697,877,1061,1249
,1442,1639,1844,2043,2245,2453]
combo_l = [11,26,43,64,86,109,129,159,181,213,238,271,290,322,347
,369,403,415,436,459,437,450,418,366,379,476,618,785,952,1126
,1303,1482,1666,1850,2037,2225]
combo_n = [18,38,61,83,105,127,143,174,185,214,223,246,252,271,282,283,309,300
,301,316,286,282,240,127,32,21,0,0,0,0,0,0,0,0,0,0]
combo_q = [15,36,59,84,109,133,156,186,203,234,248,284,290,326,353,356,394,387
,409,435,417,423,373,251,94,77,0,0,0,0,0,0,0,0,0,0]
combo_r = [16,34,56,80,108,138,169,208,232,269,292,331,347,383,418
,429,472,499,503,537,524,525,486,437,345,351,398,482,560,655
,760,877,1006,1130,1254,1380]
combo_u = [22,46,73,101,128,155,180,207,230,252,271,297,314,345,354,376,412,410
,427,448,428,416,377,223,149,107,74,56,0,0,0,0,0,0,0,0]
combo_t = [29,63,103,139,185,221,249,295,314,358,373,426,443,476,514,532,567,599
,605,616,603,606,527,407,218,139,6,0,0,0,0,0,0,0,0,0]


runs = range(1, numReps)
runs = numpy.array(runs)
print runs
import pdb
pylab.clf()
pylab.show()
pylab.draw()
allCancer = numpy.array(allCancer)
combo_b = numpy.array(combo_b)
combo_c = numpy.array(combo_c)
combo_d = numpy.array(combo_d)
ca = matplotlib.pyplot.scatter(runs, allCancer, color="red")
c1 = matplotlib.pyplot.scatter(runs, combo_b, color="BlueViolet")
c2 = matplotlib.pyplot.scatter(runs, combo_c, color="Coral")
c3 = matplotlib.pyplot.scatter(runs, combo_d, color="DarkCyan")
c4 = matplotlib.pyplot.scatter(runs, combo_e, color="DarkMagenta")
c5 = matplotlib.pyplot.scatter(runs, combo_f, color="GoldenRod")
c6 = matplotlib.pyplot.scatter(runs, combo_g, color="DimGray")
c7 = matplotlib.pyplot.scatter(runs, combo_h, color="MediumSeaGreen")
c8 = matplotlib.pyplot.scatter(runs, combo_i, color="MediumVioletRed")
c9 = matplotlib.pyplot.scatter(runs, combo_j, color="Yellow")
c10 = matplotlib.pyplot.scatter(runs, combo_k, color="Tomato")
c11 = matplotlib.pyplot.scatter(runs, combo_l, color="DarkOliveGreen")
c12 = matplotlib.pyplot.scatter(runs, combo_m, color="DarkOrange")
c13 = matplotlib.pyplot.scatter(runs, combo_n, color="DarkSlateBlue")
c14 = matplotlib.pyplot.scatter(runs, combo_o, color="DarkSlateGray")
c15 = matplotlib.pyplot.scatter(runs, combo_p, color="DeepPink")
c16 = matplotlib.pyplot.scatter(runs, combo_q, color="Khaki")
c17 = matplotlib.pyplot.scatter(runs, combo_r, color="LawnGreen")
c18 = matplotlib.pyplot.scatter(runs, combo_s, color="Maroon")
c19 = matplotlib.pyplot.scatter(runs, combo_t, color="MediumTurquoise")
c20 = matplotlib.pyplot.scatter(runs, combo_u, color="Navy")
c21 = matplotlib.pyplot.scatter(runs, combo_v, color="Peru")
from pylab import *

pylab.legend((ca, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15, c16, c17, c18, c19, c20, c21), 
('All Hallmarks', 'SG & IGI', 'SG & AA', 'SG & IT', 'SG & A', 'SG & GU', 'SG & AI', 'IGI & AA', 'IGI & IT', 'IGI & A', 'IGI & GU', 'IGI & AI',
'AA & IT', 'AA & A', 'AA & GU', 'AA & AI', 'IT & A', 'IT & GU', 'IT & AI', 'A & GU', 'A & AI', 'GU & AI'), scatterpoints=1,loc='upper left', fontsize=11, ncol=2)
pylab.xlabel("Simulation counter (step size = 400)", fontsize=11)
pylab.ylabel("Cancer cell count", fontsize=11)
pylab.ylim(-100,2500)
pylab.xlim(0, 40)
pictureFileName = "..\\Stats\\" + "ComparePhenoGrowth.png"
#pylab.savefig(pictureFileName, dpi=150)
