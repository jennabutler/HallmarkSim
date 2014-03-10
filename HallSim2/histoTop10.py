#top 10
import numpy
import matplotlib
import matplotlib.pylab as pylab
import matplotlib.pyplot
import pdb
from collections import Counter


phenos = [128, 20, 0, 144, 4, 16, 160, 136, 192, 52, 128, 20, 0, 4, 16, 144, 130, 136, 132, 22, 
128, 160, 4, 0, 32, 36, 132, 136, 164, 130, 128, 22, 4, 0, 144, 160, 54, 130, 178, 132, 
128, 4, 0, 136, 132, 68, 196, 130, 192, 8, 128, 4, 0, 20, 22, 132, 144, 192, 130, 2, 
128, 4, 0, 132, 20, 136, 144, 192, 64, 130, 128, 4, 0, 144, 132, 28, 192, 20, 16, 136, 
128, 6, 4, 134, 0, 130, 160, 132, 192, 2,  128, 4, 0, 132, 68, 160, 192, 36, 64, 
128, 4, 0, 136, 192, 8, 160, 12, 36, 128, 4, 0, 22, 20, 144, 86, 132, 82, 160,
128, 4, 0, 132, 20, 192, 144, 160, 68, 64, 128, 4, 0, 132, 160, 144, 136, 192, 68, 20]



from collections import Counter
import numpy as np
import matplotlib.pyplot as plt
from operator import itemgetter


c = Counter(phenos).items()
c.sort(key=itemgetter(1))

font = {'family' : 'sanserif',
        'color'  : 'black',
        'weight' : 'normal',
        'size'   : 12,
        }

font2 = {'family' : 'sansserif',
        'color'  : 'black',
        'weight' : 'normal',
        'size'   : 16,
        }
labels, values = zip(*c)

pylab.show()
pylab.draw()

indexes = np.arange(0, 2*len(labels), 2)
width = 2
plt.bar(indexes, values, width=2, color="blueviolet")
plt.xlabel("Phenotype identifier", fontdict=font)
plt.ylabel("Number of occurances in top 10 \n phenotypes for cancerous tumours", fontdict=font)
#plt.title("Number of occurances for different phenotypes \n in top 10 subclones of a tumour", fontdict=font2)
plt.xticks(indexes + width * 0.5, labels, rotation='vertical')
pictureFileName2 = "..\\Stats\\"  + "Phenos.png"
pylab.savefig(pictureFileName2, dpi=150)


#fig.set_size_inches(18.5,10.5)
#plt.savefig('test2png.png',dpi=100)

