#!/usr/bin/python3
import matplotlib.pyplot as plt
import numpy as np
import sys
import getopt
from scipy import signal
import researchdata1258

subj = 1
startsec = 60
endsec = False
task = False

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

fig, axs = plt.subplots(len(researchdata1258.Tasks.TASKS), 2,sharex='col')

def plotTask(r,tn):
    data = researchdata1258.Tasks(subj,tn,startsec,endsec)
    fx, fy = signal.periodogram(data.ch1,data.Fs,scaling="spectrum",nfft=(data.Fs*10))
    axs[r,0].set_xlabel("time/sec")
    axs[r,0].set_ylabel("{}\nEEG/uV".format(tn))
    axs[r,0].set_ylim([-300,300])
    axs[r,0].set_xlim([0,60])
    axs[r,0].plot(data.t,data.ch1 * 1E6)
    axs[r,1].set_xlabel("F/Hz")
    axs[r,1].set_ylabel("P/1E10 V^2")
    axs[r,1].set_xlim([0,100])
    axs[r,1].set_ylim([1E-14,1E-10])
    axs[r,1].semilogy(fx,fy)
    

for t in range(len(researchdata1258.Tasks.TASKS)):
    tn = researchdata1258.Tasks.TASKS[t]
    plotTask(t,tn)
    
plt.show()
