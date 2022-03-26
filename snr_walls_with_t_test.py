import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal
import math as math
from scipy.interpolate import interp1d
import scipy.stats as stats
import researchdata1258
import noisewall
import snr

def doStats(EEGsignal_min_f,EEGsignal_max_f):
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
        for subj in range(1,20):
            noiseWall = noisewall.NoiseWall(subj,e)
            noiseWall.calcNoiseWall()
            wall_tmp.append(noiseWall.SNRwall)
            s = snr.SNR(subj,e)
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
        

    index = np.arange(len(experiments))
    height = 0.35
    fig, ax = plt.subplots()
    baseline = 20
    xleft = np.ones(len(experiments)) * -baseline
    wall_mean_shift = [x+baseline for x in wall_mean]
    snr_mean_shift = [x+baseline for x in snr_mean]
    rects_wall = ax.barh(index+height*1.1,wall_mean_shift,height,left=xleft,align='edge',color='b',xerr=wall_stddev)
    rects_snr = ax.barh(index,snr_mean_shift,height,color='y',left=xleft,align='edge',xerr=snr_stddev)
    ax.set_xlabel('dB')
    ax.set_title('SNR vs SNR wall, EEG:{:.1f}-{:.1f} Hz, BP:{:.1f}-{:.1f} Hz, NR = {:.1f}'.format(EEGsignal_min_f,
                                                                                                  EEGsignal_max_f,
                                                                                                  filter_low_f,
                                                                                                  filter_high_f,
                                                                                                  noiseReduction))
    ax.set_yticks(index + height / 2)
    ax.set_yticklabels(experiments)
    ax.set_xlim([-20,20])
    ax.legend((rects_wall, rects_snr), ('Wall', 'SNR'))
    for i in range(len(experiments)):
        s = " (p={:.3f})".format(pval[i])
        if (pval[i] < 0.05):
            s = s + "*"
        xpos = max([wall_mean[i],snr_mean[i]]) + 1
        ax.text(xpos, i + .25, s, color='blue', fontweight='bold')

doStats(4,35)

doStats(8,18)

plt.show()
