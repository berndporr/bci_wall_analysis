#!/usr/bin/python3
import matplotlib.pyplot as plt
import numpy as np
import sys
import getopt
from scipy import signal
import researchdata1258
import noise_wall
from snr import SNR

subj = 20
startsec = 60
endsec = False
task = 'read'
p300daqlatency = 100 # ms

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
        elif '-t' in opt:
            task = arg_val
        elif '-h' in opt:
            raise getopt.GetoptError(helptext)
        else:
            raise getopt.GetoptError(helptext)
except getopt.GetoptError as err:
    print (err)
    sys.exit(2)

def plotTask(r,axs,minF,maxF):
    print("\nTask:",task)
    data = researchdata1258.Tasks(subj,task,startsec=startsec,band_low=minF,band_high=maxF)
    noisewall = noise_wall.NoiseWall(subj,task,startsec=startsec,minF=minF,maxF=maxF)
    noisewall.calcNoiseWall()
    snr = SNR(subj,task,startsec=startsec,minF=minF,maxF=maxF)
    snr.calcSNR()
    print("rho = {}, Wall = {}dB, SNR = {}dB.".format(noisewall.rho,noisewall.SNRwall,snr.snrvalue))
    tn = "no filter"
    if minF and maxF:
        if (minF < 0) and (maxF < 0):
            tn = "d/dt"
        else:
            tn = "{}-{}Hz".format(minF,maxF)
    fx, fy = signal.periodogram(data.ch1,data.Fs,scaling="spectrum",nfft=int(data.Fs*10))
    axs[r,0].set_xlabel("time/sec")
    axs[r,0].set_ylabel("{}\nEEG/uV".format(tn))
    noiseMin = np.ones(len(data.ch1))*np.sqrt(noisewall.noiseVarMin)
    noiseMax = np.ones(len(data.ch1))*np.sqrt(noisewall.noiseVarMax)
    axs[r,0].set_xlim([0,60])
    axs[r,0].plot(data.t,data.ch1 * 1E6, color='blue')
    axs[r,0].plot(data.t,noiseMin * 1E6, label="\u03C3 min", linestyle='--', linewidth=1, color='orange')
    axs[r,0].plot(data.t,noiseMax * 1E6, label="\u03C3 max", linestyle='-', linewidth=1, color='red')
    axs[r,0].plot(data.t,-noiseMin * 1E6, linestyle='--', linewidth=1, color='orange')
    axs[r,0].plot(data.t,-noiseMax * 1E6, linestyle='-', linewidth=1, color='red')
    axs[r,0].legend()
    axs[r,1].set_xlabel("F/Hz")
    axs[r,1].set_ylabel("P/1E10 V^2")
    axs[r,1].set_xlim([0,100])
    axs[r,1].set_ylim([1E-15,1E-10])
    axs[r,1].semilogy(fx,fy,color='blue')
    mindB = -20
    axs[r,2].bar(["SNR_wall","SNR"],[noisewall.SNRwall+20,snr.snrvalue-mindB],bottom=mindB)
    axs[r,2].set_ylim([mindB,20])
    axs[r,2].set_ylabel("dB")

franges = [
    [False,False],
    [0.1,3],
    [8,18],
    [8,12],
    [-1,-1]
]

fig, axs = plt.subplots(len(franges),3,sharex='col',gridspec_kw={'width_ratios': [4, 1, 1]})
fig.suptitle("Task: "+task)
for fr in range(len(franges)):
    plotTask(fr,axs,franges[fr][0],franges[fr][1])

ep = researchdata1258.Evoked_potentials(subj)

fig, axs = plt.subplots(1, 2)
fig.suptitle("Task: "+task)
axs[0].set_xlabel("time/sec")
axs[0].set_ylabel("EEG/uV")
axs[0].set_ylim([-200,200])
#axs[0].set_xlim([0,60])
axs[0].plot(ep.t[ep.initial_samples_to_ignore:-ep.final_samples_to_ignore],
            ep.eeg[ep.initial_samples_to_ignore:-ep.final_samples_to_ignore] * 1E6, label="EEG", color='blue')
axs[0].plot(ep.oddball_samples/ep.Fs,np.ones(len(ep.oddball_samples))*100,"|",
            label="oddball events", color='orange')
axs[0].legend()

t,avg = ep.get_averaged_ep()
axs[1].plot(t-p300daqlatency,avg * 1E6, color='blue')
axs[1].set_xlabel("t/ms")
axs[1].set_ylabel("P300/uV")
snr = SNR(subj,task,startsec=startsec)
snr.calcSNR()
signalPwr = snr.calcSignalPower()
signalAmpl = np.ones(len(t))*np.sqrt(signalPwr)
axs[1].plot(t,signalAmpl * 1E6, color='blue')

plt.show()
