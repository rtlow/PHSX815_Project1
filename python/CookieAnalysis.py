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

# main function for our CookieAnalysis Python code
if __name__ == "__main__":
   
    haveInput = False

    for i in range(1,len(sys.argv)):
        if sys.argv[i] == '-h' or sys.argv[i] == '--help':
            continue

        InputFile = sys.argv[i]
        haveInput = True
    
    if '-h' in sys.argv or '--help' in sys.argv or not haveInput:
        print ("Usage: %s [options] [input file]" % sys.argv[0])
        print ("  options:")
        print ("   --help(-h)          print options")
        print
        sys.exit(1)
    
    Nmeas = 1
    times = []
    times_avg = []

    need_rate = True
    
    with open(InputFile) as ifile:
        for line in ifile:
            if need_rate:
                need_rate = False
                rate = float(line)
                continue
            
            lineVals = line.split()
            Nmeas = len(lineVals)
            t_avg = 0
            for v in lineVals:
                t_avg += float(v)
                times.append(float(v))

            t_avg /= Nmeas
            times_avg.append(t_avg)

    Sorter = MySort()

    times = Sorter.DefaultSort(times)
    times_avg = Sorter.DefaultSort(times_avg)
    # try some other methods! see how long they take
    # times_avg = Sorter.BubbleSort(times_avg)
    # times_avg = Sorter.InsertionSort(times_avg)
    # times_avg = Sorter.QuickSort(times_avg)

    # ADD YOUR CODE TO PLOT times AND times_avg HERE
    
    N = len(times_avg) 

    n_bins = 100

    weights = np.ones_like(times_avg)/N
    
    title = "{} measurements / experiment with rate {:.2f} cookies / day".format(Nmeas, rate)
    
    median = times_avg[int(0.5 * N)]

    plt.figure(figsize=[12,7])
    ax = plt.axes()

    plt.hist(times_avg, bins=n_bins, weights=weights, color='springgreen', alpha=0.5, ec='black')
    
    ax.axvline(median)
    plt.text(median-0.005, 0.001, 'median = {:.2f}'.format(median), rotation=90, color='k')
    
    sigma = 1

    while sigma <= 3:
        P = math.erf(sigma/np.sqrt(2.))/2.

        low = times_avg[int((0.5-P)*N)]
        high = times_avg[int((0.5+P)*N)]

        ax.axvline(low)
        plt.text(low-0.005, 0.001, '-{} $\sigma$ = {:.3f}'.format(sigma, low), rotation=90, color='k')
        ax.axvline(high)
        plt.text(high-0.005, 0.001, '+{} $\sigma$ = {:.3f}'.format(sigma, high), rotation=90, color='k')

        sigma += 1

    plt.yscale('log')
    
    plt.title(title)
    
    xaxis = "Average time between missing cookies [days]"
    
    plt.xlabel(xaxis)

    plt.ylabel('Probability')
    
    plt.grid(True, axis='y')
    
    plt.savefig('time_hist.pdf')

    plt.show()
    
    
