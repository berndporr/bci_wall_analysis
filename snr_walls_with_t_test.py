# -*- coding: utf-8 -*-
"""
Created on Sun Jul 29 21:45:30 2018

(C) 2018 Wanting Huang <172258368@qq.com>
(C) 2018 Bernd Porr <bernd.porr@glasgow.ac.uk>

GNU GENERAL PUBLIC LICENSE
Version 3, 29 June 2007
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal
import math as math

import scipy.stats as stats

# Directory of the dataset (http://researchdata.gla.ac.uk/676/):
global dataset676dir
dataset676dir = "../dataset_676"

class NoiseWall:

    ## generate EEG power from a paralysed persion in a certain
    ## frequency band
    def generateParalysedEEGVariance(self,band_low = 0,band_high = 0):
        ## from paper
        if (band_high == 0):
            band_high = self.fs/2
        self.pureEEGVar = 0
        for f in np.arange(band_low,band_high,1.0):
            p = -11.5 - f*0.02
            if ((f>7) and (f<13)):
                # alpha peak
                p = p + 1
            power_density = 10**p
            self.pureEEGVar = self.pureEEGVar + power_density
        
        #print("Paralysed EEG variance is %e" % self.pureEEGVar)

    ## Loads the data from the database
    def loadDataFromFile(self,subdir):

        self.data=np.loadtxt(subdir+"emgeeg.dat")
        self.zero_data=np.loadtxt(subdir+"zero_time_data.dat")
        self.zero_video=np.loadtxt(subdir+"zero_time_video.dat")
        self.artefact=np.loadtxt(subdir+"artefact.dat")
        self.relaxed=np.loadtxt(subdir+"dataok.dat")

        self.t=self.data[:,0]          #timestamp or sample # (sampling rate fs=1kHz)
        self.eeg=self.data[:,1]        #eeg
        self.emg=self.data[:,2]        #emg
        self.trigger=self.data[:,3]    #switch 

        #AMPLIFER GAIN IS 500, SAMPLING RATE IS 1kHz
        self.eeg=self.eeg/500
        self.emg=self.emg/500
        self.fs=1000
        self.T=1/self.fs
        self.t=np.arange(0,self.T*len(self.eeg),self.T)

    ## Filters out known powerline interference
    def filterData(self,band_low=0,band_high=0):
        # smooth it at 100Hz cutoff
        bLP,aLP = signal.butter(4,100/self.fs*2)
        self.eeg = signal.lfilter(bLP,aLP,self.eeg);

        ## highpass at 1Hz and 50Hz notch
        bfilt50hz,afilt50hz = signal.butter(2,[49/self.fs*2,51/self.fs*2],'stop')
        bhp,ahp = signal.butter(4,0.5/self.fs*2,'high')
        self.eeg = signal.lfilter(bhp,ahp,signal.lfilter(bfilt50hz,afilt50hz,self.eeg));
        #emg_clean = signal.lfilter(bhp,ahp,signal.lfilter(bfilt50hz,afilt50hz,emg));
        #trigger_clean = signal.lfilter(bhp,ahp,signal.lfilter(bfilt50hz,afilt50hz,trigger));

        ## strange 80Hz interference
        bfilt80hz,afilt80hz = signal.butter(2,[78/self.fs*2,82/self.fs*2],'stop')
        self.eeg = signal.lfilter(bfilt80hz,afilt80hz,self.eeg);
        #emg_clean = signal.lfilter(bfilt80hz,afilt80hz,emg_clean);
        #trigger_clean = signal.lfilter(bfilt80hz,afilt80hz,trigger_clean);

        ## strange 25 Hz interference
        bfilt25hz,afilt25hz = signal.butter(2,[24/self.fs*2,26/self.fs*2],'stop')
        self.eeg = signal.lfilter(bfilt25hz,afilt25hz,self.eeg);
        #emg_clean = signal.lfilter(bfilt25hz,afilt25hz,emg_clean);
        #trigger_clean = signal.lfilter(bfilt25hz,afilt25hz,trigger_clean);

        ## do we just look at a specific band?
        if (band_high > 0) and (band_low > 0) and (band_low < band_high):
            bfilt1hz,afilt10hz = signal.butter(2,[1/self.fs*2,10/self.fs*2],'bandpass')
            self.eeg = signal.lfilter(bfilt1hz,afilt10hz,self.eeg)
        #FILTER COMPLETE

        self.generateParalysedEEGVariance(band_low,band_high)

    #noise variance without an artefact / activity
    def calcNoiseVarRelaxed(self):
        dt=self.zero_data-self.zero_video
        t1=int(self.fs*(self.relaxed[0]+dt))
        t2=int(self.fs*(self.relaxed[1]+dt))
        yRelaxed=self.eeg[t1:t2]
        relaxedPower=np.var(yRelaxed) - self.pureEEGVar
        if (relaxedPower < 0):
            return 0
        return relaxedPower
    
    #noise variance with artefacts / activity
    def calcNoiseVarArtefact(self):
        #ARTEFACTS beginning / stop
        tbeginVideo=self.artefact[:,0]
        tendVideo=self.artefact[:,1]
        dt=self.zero_data-self.zero_video
        tbegin=tbeginVideo+dt
        tend=tendVideo+dt
        artefactPowerList=[]

        for i in range(len(tbegin)):
            t1=tbegin[i]
            t2=tend[i]
            t1=int(self.fs*t1)
            t2=int(self.fs*t2)
            signalWithArtefact=self.eeg[t1:t2]
            artefactPower = np.var(signalWithArtefact) - self.pureEEGVar
            artefactPowerList.append(artefactPower)

        averageArtefactPower = np.mean(artefactPowerList)
        if (averageArtefactPower < 0):
            return 0
        return averageArtefactPower

    def calcRho(self):
        noiseVarArtefact = self.calcNoiseVarArtefact()
        noiseVarRelaxed = self.calcNoiseVarRelaxed()
        rho = math.sqrt(noiseVarArtefact / noiseVarRelaxed)
        return rho
   
    def calcNoiseWall(self):
        p = self.calcRho()
        SNRwall = 10 * math.log10(p - 1/p)
        return SNRwall
    
    def calcSNR(self):
        noiseVariance = self.calcNoiseVarArtefact() / self.calcRho()
        SNR= self.pureEEGVar / noiseVariance
        SNR= 10 * math.log10(SNR)
        return SNR   #nominal noise variance of signal
    
    
def calcNoiseWall(subj,experiment,band_low=0,band_high=0):
    
    subj = "%02d" % subj
    criterion=np.loadtxt(dataset676dir+"/experiment_data/subj"+subj+"/"+"all_exp_ok.dat", dtype=bytes).astype(str)
    #print(criterion)
    if criterion == 'False':
        return 0
    else:
        noiseWall = NoiseWall()
        noiseWall.loadDataFromFile(dataset676dir+"/experiment_data/subj"+subj+"/"+experiment+"/")
        noiseWall.filterData(band_low,band_high)
        #noiseWall.plotData()
        noiseWall.calcNoiseVarRelaxed()
        noiseWall.calcNoiseVarArtefact()
        return (noiseWall.calcNoiseWall(),noiseWall.calcSNR())



    
experiments = ["lie_relax","blink","eyescrunching","raisingeyebrows","jaw","readingsitting","readinglieing","flow","sudoku","wordsearch","templerun"]
Num=0
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
        wall = calcNoiseWall(subj,e)

        if not (wall == 0):
            print("subject%02d: Noise wall of %s = %f dB" % (subj,e,wall[0]))
            print("subject%02d: SNR of %s = %f dB" % (subj,e,wall[1]))
            wall_tmp.append(wall[0])
            snr_tmp.append(wall[1])

    wall_mean.append(np.mean(wall_tmp))
    wall_stddev.append(np.std(wall_tmp))
    snr_mean.append(np.mean(snr_tmp))
    snr_stddev.append(np.std(snr_tmp))
    t, p = stats.ttest_rel(snr_tmp, wall_tmp)
    # one sided: p is half
    p = p / 2
    # we reject anything where the SNR is less than the Wall
    if (t<0):
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
ax.set_title('SNR vs SNR wall')
ax.set_yticks(index + height / 2)
ax.set_yticklabels(experiments)
ax.set_xlim([-20,20])
ax.legend((rects_wall, rects_snr), ('Wall', 'SNR'))
for i in range(len(experiments)):
    s = ""     
    if (pval[i] < pval_for_significance):
        s = s + "* p=%0.03f" % pval[i]
    else:
        s = s + " (p=%0.03f)" % pval[i]
    xpos = max([wall_mean[i],snr_mean[i]]) + 1
    ax.text(xpos, i + .25, s, color='blue', fontweight='bold')
plt.show()
