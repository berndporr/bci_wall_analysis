#!/usr/bin/python3
import matplotlib.pyplot as plt
import numpy as np
import sys
import getopt
from scipy import signal
import researchdata1258

subj = 20
startsec = 60
endsec = False
task = False

helptext = 'usage: {} -p participant -s startsec -f noiseredfile.tsv -t task -m -h'.format(sys.argv[0])
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
            startsec = float(arg_val)
        elif '-t' in opt:
            task = arg_val
        elif '-h' in opt:
            raise getopt.GetoptError(helptext)
        else:
            raise getopt.GetoptError(helptext)
except getopt.GetoptError as err:
    print (err)
    sys.exit(2)

def plotTask(r,tn,axs):
    data = researchdata1258.Tasks(subj,tn,startsec)
    fx, fy = signal.periodogram(data.ch1,data.Fs,scaling="spectrum",nfft=int(data.Fs*10))
    axs[r,0].set_xlabel("time/sec")
    axs[r,0].set_ylabel("{}\nEEG/uV".format(tn))
    if researchdata1258.Tasks.TASKS[0] in tn:
        axs[r,0].set_ylim([-1000,1000])
    else:
        axs[r,0].set_ylim([-200,200])
    axs[r,0].set_xlim([0,60])
    axs[r,0].plot(data.t,data.ch1 * 1E6)
    axs[r,1].set_xlabel("F/Hz")
    axs[r,1].set_ylabel("P/1E10 V^2")
    axs[r,1].set_xlim([0,100])
    axs[r,1].set_ylim([1E-15,1E-10])
    axs[r,1].semilogy(fx,fy)

if task:
    fig, axs = plt.subplots(1, 2, squeeze=False)
    plotTask(0,task,axs)
else:
    fig, axs = plt.subplots(len(researchdata1258.Tasks.TASKS), 2,sharex='col')
    for t in range(len(researchdata1258.Tasks.TASKS)):
        tn = researchdata1258.Tasks.TASKS[t]
        plotTask(t,tn,axs)

plt.show()
