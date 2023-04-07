#!/usr/bin/python3
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal
import math as math
from scipy.interpolate import interp1d
import scipy.stats as stats
import researchdata1258
import bci_wall
import snr
import sys
import getopt


subjectsOK = [1,3,4,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]

startsec = 30

def doStats(EEGsignal_min_f=False,EEGsignal_max_f=False):
    if EEGsignal_min_f or EEGsignal_max_f:
        print("Stats with min_f={}, max_f={}:".format(EEGsignal_min_f,EEGsignal_max_f))
    else:
        print("Stats with the full frequency range:")
    wall_mean=[]
    wall_stddev=[]
    snr_mean=[]
    snr_stddev=[]
    pval = []
    pval_for_significance = 0.05
    for e in researchdata1258.Tasks.TASKS:
        print("\n"+e)
        wall_tmp = []
        snr_tmp = []
        for subj in subjectsOK:
            print(e,subj,":")
            noiseWall = bci_wall.NoiseWall(subj,e,startsec=startsec,minF=EEGsignal_min_f,maxF=EEGsignal_max_f)
            noiseWall.calcNoiseWall()
            wall_tmp.append(noiseWall.SNRwall)
            s = snr.SNR(subj,e,startsec=startsec,minF=EEGsignal_min_f,maxF=EEGsignal_max_f)
            s.calcSNR()
            snr_tmp.append(s.snrvalue)
            print("Wall = {}, SNR = {}.".format(noiseWall.SNRwall,s.snrvalue))
        wall_mean.append(np.mean(wall_tmp))
        wall_stddev.append(np.std(wall_tmp))
        snr_mean.append(np.mean(snr_tmp))
        snr_stddev.append(np.std(snr_tmp))
        # test how likely it is that the average SNR has hit the SNR wall
        t, p = stats.ttest_rel(snr_tmp, wall_tmp)
        if t < 0:
            p = 1
        print("Experiment {} has p={}, t={}".format(e,p,t))
        pval.append(p)
        

    index = np.arange(len(researchdata1258.Tasks.TASKS))
    height = 0.35
    fig, ax = plt.subplots()
    baseline = 20
    xleft = np.ones(len(researchdata1258.Tasks.TASKS)) * -baseline
    wall_mean_shift = [x+baseline for x in wall_mean]
    snr_mean_shift = [x+baseline for x in snr_mean]
    rects_wall = ax.barh(index+height*1.1,wall_mean_shift,height,left=xleft,align='edge',color='b',xerr=wall_stddev)
    rects_snr = ax.barh(index,snr_mean_shift,height,color='y',left=xleft,align='edge',xerr=snr_stddev)
    ax.set_xlabel('dB')
    ax.set_title('SNR vs SNR wall, BP:{:.1f}-{:.1f} Hz'.format(EEGsignal_min_f,
                                                               EEGsignal_max_f));
    ax.set_yticks(index + height / 2)
    ax.set_yticklabels(researchdata1258.Tasks.TASKS)
    ax.set_xlim([-20,20])
    ax.legend((rects_wall, rects_snr), ('Wall', 'SNR'))
    for i in range(len(researchdata1258.Tasks.TASKS)):
        s = " (p={:.3f})".format(pval[i])
        if (pval[i] < 0.05):
            s = s + "*"
        xpos = max([wall_mean[i],snr_mean[i]]) + 1
        ax.text(xpos, i + .25, s, color='blue', fontweight='bold')

    

helptext = 'usage: {} -m -n -e -d -a lowF -b highF'.format(sys.argv[0])
helptext = helptext + "\n  -m: 8-18 Hz"
helptext = helptext + "\n  -n: 8-12 Hz"
helptext = helptext + "\n  -e: 0.1-3 Hz"
helptext = helptext + "\n  -d: 1st order highpass"

minF = False
maxF = False

try:
    # Gather the arguments
    all_args = sys.argv[1:]
    opts, arg = getopt.getopt(all_args, 'a:b:mnedh')
    # Iterate over the options and values
    for opt, arg_val in opts:
        if '-m' in opt:
            minF = 8
            maxF = 18
        elif '-n' in opt:
            minF = 8
            maxF = 12
        elif '-e' in opt:
            minF = 0.1
            maxF = 3
        elif '-d' in opt:
            minF = -1
            maxF = -1
        elif '-a' in opt:
            minF = float(arg_val)
        elif '-b' in opt:
            maxF = float(arg_val)
        elif '-h' in opt:
            raise getopt.GetoptError(helptext)
        else:
            raise getopt.GetoptError(helptext)
except getopt.GetoptError as err:
    print (err)
    sys.exit(2)

doStats(minF,maxF)
plt.show()
