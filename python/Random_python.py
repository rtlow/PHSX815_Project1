#! /usr/bin/env python

# imports of external packages to use in our code
import sys
import numpy as np
import matplotlib.pyplot as plt
from Random import *

# main function for this Python code
if __name__ == "__main__":
    # if the user includes the flag -h or --help print the options
    if '-h' in sys.argv or '--help' in sys.argv:
        print ("Usage: %s [-seed number]" % sys.argv[0])
        print
        sys.exit(1)

    # default seed
    seed = 5555

    fname = None

    # read the user-provided seed from the command line (if there)
    if '-seed' in sys.argv:
        p = sys.argv.index('-seed')
        seed = sys.argv[p+1]

    # read the user-provided filename (if there)
    if '-file' in sys.argv:
        p = sys.argv.index('-file')
        fname = sys.argv[p+1]


    # class instance of our Random class using seed
    random = Random(seed)

    # create some random data
    N = 10000

    # an array of random numbers using our Random class
    myx = []
    for i in range(0,N):
        myx.append(random.Poisson(rate=1000))

    # create histogram of our data
    n, bins, patches = plt.hist(myx, 50, density=True, facecolor='g', alpha=0.75)

    # plot formating options
    plt.xlabel('x')
    plt.ylabel('Probability')
    plt.title('Normally distributed random number')
    plt.grid(True)

    # save the figure if called for
    if fname is not None:
        plt.savefig(fname)

    # show figure (program only ends once closed )
    plt.show()
