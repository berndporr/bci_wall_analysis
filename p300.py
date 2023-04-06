#!/usr/bin/python3
import matplotlib.pyplot as plt
import numpy as np
import sys
import getopt
import researchdata1258
from scipy import signal

subj = 20
startsec = 60

helptext = 'usage: {} -p participant -s startsec -f file -h'.format(sys.argv[0])

try:
    # Gather the arguments
    all_args = sys.argv[1:]
    opts, arg = getopt.getopt(all_args, 'p:s:h')
    # Iterate over the options and values
    for opt, arg_val in opts:
        if '-p' in opt:
            subj = int(arg_val)
        elif '-s' in opt:
            startsec = int(arg_val)
        elif '-h' in opt:
            raise getopt.GetoptError(helptext)
        else:
            raise getopt.GetoptError(helptext)
except getopt.GetoptError as err:
    print (err)
    sys.exit(2)

ep = researchdata1258.Evoked_potentials(subj)

fig, axs = plt.subplots(1, 2)
axs[0].set_xlabel("time/sec")
axs[0].set_ylabel("EEG/uV")
axs[0].set_ylim([-200,200])
#axs[0].set_xlim([0,60])
axs[0].plot(ep.t[ep.initial_samples_to_ignore:-ep.final_samples_to_ignore],
            ep.eeg[ep.initial_samples_to_ignore:-ep.final_samples_to_ignore] * 1E6, label="EEG")
axs[0].plot(ep.oddball_samples/ep.Fs,np.ones(len(ep.oddball_samples))*100,"|",
            label="oddballl events")
axs[0].legend()

t,avg = ep.get_averaged_ep()
axs[1].plot(t,avg * 1E6)
axs[1].set_xlabel("t/ms")
axs[1].set_ylabel("P300/uV")

plt.show()
