#! /usr/bin/env python

# imports of external packages to use in our code
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

#setting matplotlib ticksize
matplotlib.rc('xtick', labelsize=14) 
matplotlib.rc('ytick', labelsize=14) 

#set matplotlib global font size
matplotlib.rcParams['font.size']=14

# import our Random class from python/Random.py file
sys.path.append(".")
from python.MySort import MySort


# helper functions


def PoisProb(x, rate):
    return (rate**x) * np.exp(-rate) / math.factorial(x)

def PlotHypotheses(array0, array1, title, alpha):

    N0 = len(array0)
    N1 = len(array1)

    hmin = min(array0[0], array1[0])
    hmax = max(array0[N0-1], array1[N1-1])

    fig = plt.figure(figsize=[12,7])
    ax = plt.axes()
    plt.hist(array0, N0 + 1, density=True, facecolor='b', alpha=0.5, label='$P(\lambda | H0)$')
    plt.hist(array1, N1 + 1, density=True, facecolor='g', alpha=0.5, label='$P(\lambda | H1)$')
    
    ax.set_yscale('log')

    ax.set_xlabel('$\lambda = \log [ \mathcal{L}(H1) / \mathcal{L}(H0) ]$')
    ax.set_ylabel('Probability')

    plt.legend()

    return fig, ax

# main function for our CookieAnalysis Python code
if __name__ == "__main__":
   
    haveInput = [False, False]

    InputFile = [None, None]

    alpha = 0.

    for i in range(1,len(sys.argv)):
        if sys.argv[i] == '-h' or sys.argv[i] == '--help':
            continue

        # seeing if we have input files
        if sys.argv[i] == '-H0':
            InputFile[0] = sys.argv[i+1]
            haveInput[0] = True

        if sys.argv[i] == '-H1':
            InputFile[1] = sys.argv[i+1]
            haveInput[1] = True
    
    if '-h' in sys.argv or '--help' in sys.argv or not np.all(haveInput):
        print ("Usage: %s [options] -H0 [input file for H0] -H1 [input file for H1]" % sys.argv[0])
        print ("  options:")
        print ("   --help(-h)          print options")
        print ("   -alpha [number]     significance of test")
        print
        sys.exit(1)
    
    # reading in data from files

    Nmeas = 0
    rate = []
    counts = []
    need_rate = True
    
    for h in range(2):
        
        need_rate = True
        this_hyp = []
        
        with open(InputFile[h]) as ifile:
            for line in ifile:
                if need_rate:
                    need_rate = False
                    rate.append(float(line))
                    continue
            
                lineVals = line.split()
                Nmeas = len(lineVals)
                
                this_exp = []
                for m in range(Nmeas):
                    this_exp.append(float(lineVals[m]))
                this_hyp.append(this_exp)

        counts.append(this_hyp)


    LLR = []

    for h in range(2):
        
        this_hyp = []

        Nexp = len(counts[h])

        for e in range(Nexp):
            Nmeas = len(counts[h][e])

            LogLikeRatio = 0.

            for m in range(Nmeas):

                LogLikeRatio += np.log( PoisProb( counts[h][e][m], rate[1] ) ) # LLR for H1

                LogLikeRatio -= np.log( PoisProb( counts[h][e][m], rate[0] ) ) # LLR for H0

            this_hyp.append(LogLikeRatio)

        LLR.append(this_hyp)
    
    Sorter = MySort()

    LLR[0] = Sorter.DefaultSort(LLR[0])
    LLR[1] = Sorter.DefaultSort(LLR[1])

    plot_title = "{} measurements / experiment with rates {:.2f}, {:.2f} counts / sec".format(Nmeas, rate[0], rate[1])

    fig, ax = PlotHypotheses(LLR[0], LLR[1], plot_title, alpha)


    fig.savefig('Hypotheses.pdf')

    fig.show()

