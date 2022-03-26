#!/usr/bin/python3
import matplotlib.pyplot as plt
import numpy as np
import sys
import getopt
import researchdata1258


# check if we run this as a main program
if __name__ == "__main__":
    subj = 1
    startsec = 2

    helptext = 'usage: {} -p participant -s startsec -f file -h'.format(sys.argv[0])

    try:
        # Gather the arguments
        all_args = sys.argv[1:]
        opts, arg = getopt.getopt(all_args, 'p:s:f:')
        # Iterate over the options and values
        for opt, arg_val in opts:
            if '-p' in opt:
                subj = int(arg_val)
            elif '-s' in opt:
                startsec = int(arg_val)
            elif '-h' in opt:
                raise getopt.GetoptError()
            else:
                raise getopt.GetoptError()
    except getopt.GetoptError:
        print (helptext)
        sys.exit(2)

    ep = researchdata1258.Evoked_potentials(subj,startsec)
    t,avg = ep.get_averaged_ep()
    plt.figure("P300")
    plt.plot(t,avg)
    plt.xlabel("t/ms")
    plt.ylabel("P300/uV")
    plt.show()
