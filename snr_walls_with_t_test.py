from NoiseWall import NoiseWall
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal
import math as math
from scipy.interpolate import interp1d
import scipy.stats as stats

def doStats(low_f,high_f):
    experiments = ["lie_relax","blink","eyescrunching","raisingeyebrows","jaw","readingsitting","readinglieing","flow","sudoku","wordsearch","templerun"]
    wall_mean=[]
    wall_stddev=[]
    snr_mean=[]
    snr_stddev=[]
    pval = []
    pval_for_significance = 0.05
    for e in experiments:
        print("\n"+e)
        wall_tmp = []
        snr_tmp = []
        for subj in range(2,28):
            noiseWall = NoiseWall(subj,e)
            noiseWall.noiseReduction = 6
            if noiseWall.dataok:
                noiseWall.filterData(low_f,high_f)
                try:
                    noiseWall.doAllCalcs()
                    print("Noise wall of subj %s = %f dB" % (subj,noiseWall.getSNRwall()))
                    print("SNR of subj %s = %f dB" % (subj,noiseWall.getSNR()))
                    wall_tmp.append(noiseWall.getSNRwall())
                    snr_tmp.append(noiseWall.getSNR())
                except NoiseWall.NoiseWallException as err:
                    print("subj"+str(subj)+"="+str(err))
        wall_mean.append(np.mean(wall_tmp))
        wall_stddev.append(np.std(wall_tmp))
        snr_mean.append(np.mean(snr_tmp))
        snr_stddev.append(np.std(snr_tmp))
        # test how likely it is that the average SNR has hit the SNR wall
        t, p = stats.ttest_rel(snr_tmp, wall_tmp)
        if t < 0:
            p = 1
        print("Experiment %s has p=%f, t=%f" % (e,p,t))
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
    ax.set_title('SNR vs SNR wall, p = average SNR has hit the SNR wall')
    ax.set_yticks(index + height / 2)
    ax.set_yticklabels(experiments)
    ax.set_xlim([-20,20])
    ax.legend((rects_wall, rects_snr), ('Wall', 'SNR'))
    for i in range(len(experiments)):
        s = " (p=%0.03f)" % pval[i]
        xpos = max([wall_mean[i],snr_mean[i]]) + 1
        ax.text(xpos, i + .25, s, color='blue', fontweight='bold')
    plt.show()




doStats(4,35)
#doStats(8,13)
#doStats(4,50)
