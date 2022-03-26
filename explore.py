#!/usr/bin/python3
import matplotlib.pyplot as plt
import numpy as np
import sys
import getopt
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from scipy import signal
import researchdata1258

# check if we run this as a main program
if __name__ == "__main__":
    subj = 1
    startsec = 60
    endsec = False
    task = researchdata1258.Tasks.TASKS[0]
    usePlotly = True

    helptext = 'usage: {} -p participant -s startsec -e endsec -f noiseredfile.tsv -t task -m -h'.format(sys.argv[0])

    try:
        # Gather the arguments
        all_args = sys.argv[1:]
        opts, arg = getopt.getopt(all_args, 'p:s:e:f:t:m')
        # Iterate over the options and values
        for opt, arg_val in opts:
            if '-p' in opt:
                subj = int(arg_val)
            elif '-s' in opt:
                startsec = int(arg_val)
            elif '-e' in opt:
                endsec = int(arg_val)
            elif '-t' in opt:
                task = arg_val
            elif '-m' in opt:
                usePlotly = False
            elif '-h' in opt:
                raise getopt.GetoptError()
            else:
                raise getopt.GetoptError()
    except getopt.GetoptError:
        print (helptext)
        print ("Option -m switches over to matplotlib. Default is plotly.")
        sys.exit(2)

    data = researchdata1258.Tasks(subj,task,startsec,endsec)

    plt.plot(data.t,data.ch1)
    plt.show()
