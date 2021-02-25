#! /usr/bin/env python

# imports of external packages to use in our code
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scipy.special as special

#setting matplotlib ticksize
matplotlib.rc('xtick', labelsize=14) 
matplotlib.rc('ytick', labelsize=14) 

#set matplotlib global font size
matplotlib.rcParams['font.size']=14

# import our Random class from python/Random.py file
sys.path.append(".")
from python.MySort import MySort


# helper functions

# PoisProb approximates the Poisson probability
# using the Stirling approximation
def PoisProb(x, rate):

    logP = x * np.log(rate) - rate - special.gammaln(x + 1)

    P = np.exp(logP)
    
    return P

# given two sets of data and a significance level, plots the histograms
# with significance
def PlotHypotheses(array0, array1, title, alpha):

    N0 = len(array0)
    N1 = len(array1)

    hmin = min(array0[0], array1[0])
    hmax = max(array0[N0-1], array1[N1-1])

    fig = plt.figure(figsize=[12,7])
    ax = plt.axes()
    n1, b1, p1 = plt.hist(array0, N0 + 1, density=True, facecolor='b', alpha=0.5, label='$P(\lambda | H0)$')
    n2, b2, p2 = plt.hist(array1, N1 + 1, density=True, facecolor='g', alpha=0.5, label='$P(\lambda | H1)$')
    
    ax.set_yscale('log')

    ax.set_xlabel('$\lambda = \log [ \mathcal{L}(H1) / \mathcal{L}(H0) ]$')
    ax.set_ylabel('Probability')

    plt.legend()
    
    
    lambda_crit = 0
    beta = 0
    
    print(alpha)
    print(1/N0)
    # TODO figure out this block TODO
    if(alpha > 1/N0):
        lambda_crit = array0[min(int((1-alpha)*N0), N0-1)]
        beta = 0
        
        below_alpha = []
        above_alpha = []

        for i in range(N0):
            if array0[i] > lambda_crit:
                above_alpha.append(array0[i])
        for i in range(N1):
            if array1[i] < lambda_crit:
                beta +=1

                below_alpha.append(array1[i])
        beta /= N1
        
        below_alpha = np.array(below_alpha)
        above_alpha = np.array(above_alpha)

        w1 = below_alpha / N0
        w2 = above_alpha / N1

        #plt.hist(below_alpha, bins=b1, weights=w1, alpha=0.5, color='darkblue')
        #plt.hist(above_alpha, bins=b2, weights=w2, alpha=0.5, color='lime')
        
        plt.axvline(lambda_crit)
        plt.text(lambda_crit, ax.get_ylim()[1] * 0.8, 'alpha = {:.3f}'.format(alpha))
    fig.savefig('Hypotheses.pdf')
    fig.show()

    return

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

        if sys.argv[i] == '-alpha':
            alpha = float(sys.argv[i + 1])
    
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
    
    # loop over all hypotheses (only 2)
    for h in range(2):
        
        need_rate = True
        this_hyp = []
        
        with open(InputFile[h]) as ifile:
            
            # parse each line
            for line in ifile:
                
                # first line is the rate parameter
                if need_rate:
                    need_rate = False
                    rate.append(float(line))
                    continue
            
                # each line is a different experiment
                lineVals = line.split()
                Nmeas = len(lineVals)
                
                this_exp = []
                
                # need to go through all measurements to convert them from string to float
                for m in range(Nmeas):
                    this_exp.append(float(lineVals[m]))
                this_hyp.append(this_exp)

        counts.append(this_hyp)


    LLR = []
    
    # loop over all hypotheses
    for h in range(2):
        
        this_hyp = []

        Nexp = len(counts[h])

        # loop over all experiments
        for e in range(Nexp):
            Nmeas = len(counts[h][e])

            LogLikeRatio = 0.

            # loop over all measurements to calculate the LLR
            for m in range(Nmeas):
    
                # LLR is a sum; one contributes positive, other negative
                LogLikeRatio += np.log( PoisProb( counts[h][e][m], rate[1] ) ) # LLR for H1

                LogLikeRatio -= np.log( PoisProb( counts[h][e][m], rate[0] ) ) # LLR for H0

            this_hyp.append(LogLikeRatio)

        LLR.append(this_hyp)
    
    # sort the data
    Sorter = MySort()

    LLR[0] = Sorter.DefaultSort(LLR[0])
    LLR[1] = Sorter.DefaultSort(LLR[1])

    plot_title = "{} measurements / experiment with rates {:.2f}, {:.2f} counts / sec".format(Nmeas, rate[0], rate[1])

    # plot the histogram
    PlotHypotheses(LLR[0], LLR[1], plot_title, alpha)

