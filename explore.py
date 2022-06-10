#!/usr/bin/python3
import matplotlib.pyplot as plt
import numpy as np
import sys
import getopt
from scipy import signal
import researchdata1258

# check if we run this as a main program
if __name__ == "__main__":
    subj = 1
    startsec = 60
    endsec = False
    task = researchdata1258.Tasks.TASKS[0]

    helptext = 'usage: {} -p participant -s startsec -e endsec -f noiseredfile.tsv -t task -m -h'.format(sys.argv[0])
    helptext = helptext + "\nOption -m switches over to matplotlib. Default is plotly."

    try:
        # Gather the arguments
        all_args = sys.argv[1:]
        opts, arg = getopt.getopt(all_args, 'p:s:e:f:t:mh')
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
            elif '-h' in opt:
                raise getopt.GetoptError(helptext)
            else:
                raise getopt.GetoptError(helptext)
    except getopt.GetoptError as err:
        print (err)
        sys.exit(2)

    data = researchdata1258.Tasks(subj,task,startsec,endsec)

    plt.figure("Timedomain: {} {}".format(subj,task))
    plt.xlabel("time/sec")
    plt.ylabel("EEG/uV")
    plt.plot(data.t,data.ch1 * 1E6)
    plt.figure("Spectrum: {} {}".format(subj,task))
    plt.xlabel("Frequency/Hz")
    plt.ylabel("Power/V^2")
    fx, fy = signal.periodogram(data.ch1,data.Fs,scaling="spectrum",nfft=data.Fs)
    plt.xlim([0,100])
    plt.plot(fx,fy)
    plt.show()
