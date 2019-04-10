# -*- coding: utf-8 -*-
"""
(C) 2018 Wanting Huang <172258368@qq.com>
(C) 2018-2019 Bernd Porr <bernd.porr@glasgow.ac.uk>

GNU GENERAL PUBLIC LICENSE
Version 3, 29 June 2007

This script calculates the noise walls
of EEG recordings during different tasks.

Thus, the script calculates if it possible
at all to detect a change in EEG at all
buried in EMG noise.

The dataset used is: http://researchdata.gla.ac.uk/676
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal
import math as math
from scipy.interpolate import interp1d
import scipy.stats as stats

# Directory of the dataset (http://researchdata.gla.ac.uk/676/):
global dataset676dir
dataset676dir = "../dataset_676"

# Calculates the noise wall from one subject doing one experiment
class NoiseWall:

    # Exception arising from the calculations
    # These arise from sanity checks and point
    # to possibly wrong/bad recordings so that they
    # can be excluded.
    class NoiseWallException(Exception):
        pass

    DATA_INVALID = "Data invalid"
    MIN_VAR_NEG = "Min variance less than pure EEG var"
    MAX_VAR_NEG = "Max variance less than pure EEG var"
    MIN_VAR_LARGER_THAN_MAX_VAR = "Min variance larger than max variance"
    MIN_VAR_UNUSUALLY_HIGH = "Unusually high Min variance"

    # Constructor with subject number and experiment
    # Loads the data if a complete dataset exists for this experiment
    def __init__(self,subj,experiment):
        self.noiseReduction = 1
        s = "%02d" % subj
        d = np.loadtxt(dataset676dir+"/experiment_data/subj"+s+"/"+"all_exp_ok.dat", dtype=bytes).astype(str)
        self.dataok = d in ['True','true','ok','OK']
        self.subdir = dataset676dir+"/experiment_data/subj"+s+"/"+experiment+"/"
        if self.dataok:
            self.loadDataFromFile(self.subdir)

    def calcParalysedEEGVariance(self,filename,band_low,band_high):
        a = np.loadtxt(filename)
        f = a[:,0]
        p = a[:,1]
        psd = interp1d(f, p, kind='cubic')
        bandpower = 0
        for f2 in np.arange(band_low,band_high,1.0):
            bandpower = bandpower + 10**psd(f2)
        return bandpower

    def generateParalysedEEGVariance(self,band_low,band_high):
        self.pureEEGVar = 0
        self.pureEEGVar = self.pureEEGVar + self.calcParalysedEEGVariance("sub1a.dat",band_low,band_high)
        self.pureEEGVar = self.pureEEGVar + self.calcParalysedEEGVariance("sub1b.dat",band_low,band_high)
        self.pureEEGVar = self.pureEEGVar + self.calcParalysedEEGVariance("sub1c.dat",band_low,band_high)
        self.pureEEGVar = self.pureEEGVar + self.calcParalysedEEGVariance("sub2a.dat",band_low,band_high)
        self.pureEEGVar = self.pureEEGVar + self.calcParalysedEEGVariance("sub2b.dat",band_low,band_high)
        self.pureEEGVar = self.pureEEGVar + self.calcParalysedEEGVariance("sub2c.dat",band_low,band_high)
        self.pureEEGVar = self.pureEEGVar / 6.0

    # Loads the data from the database
    # this is a private function and there is no need to call it
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
        corrfactor = 2 # adjusting the average EEG from the exeriment to 10E-13 V^2/Hz from E.M.Whitham
        self.eeg=self.eeg/500*corrfactor
        self.fs=1000
        self.T=1/self.fs
        self.t=np.arange(0,self.T*len(self.eeg),self.T)

    # Filters out known powerline interference and limits the EEG
    # to the band of interest.
    def filterData(self,band_low=0,band_high=0):
        # smooth it at 50Hz cutoff
        bLP,aLP = signal.butter(4,50/self.fs*2)
        self.eeg = signal.lfilter(bLP,aLP,self.eeg);

        ## highpass at 1Hz and 50Hz notch
        bfilt50hz,afilt50hz = signal.butter(2,[49/self.fs*2,51/self.fs*2],'stop')
        bhp,ahp = signal.butter(4,0.5/self.fs*2,'high')
        self.eeg = signal.lfilter(bhp,ahp,signal.lfilter(bfilt50hz,afilt50hz,self.eeg));
        #emg_clean = signal.lfilter(bhp,ahp,signal.lfilter(bfilt50hz,afilt50hz,emg));
        #trigger_clean = signal.lfilter(bhp,ahp,signal.lfilter(bfilt50hz,afilt50hz,trigger));

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


    def getMinNoiseVarEEGChunk(self):
        dt=self.zero_data-self.zero_video
        t1=int(self.fs*(self.relaxed[0]+dt))
        t2=int(self.fs*(self.relaxed[1]+dt))
        return self.eeg[t1:t2]
        
    # Minimal noise variance taken from a stretch before the experiment
    # starts.
    def calcNoiseVarMin(self):
        yMin=self.getMinNoiseVarEEGChunk()
        yMinVar = np.var(yMin)
        self.noiseVarMin= (yMinVar - self.pureEEGVar) / self.noiseReduction
        if (self.noiseVarMin < 0):
            raise self.NoiseWallException(self.MIN_VAR_NEG)
        if (yMinVar**0.5) > 50E-6:
            raise self.NoiseWallException(self.MIN_VAR_UNUSUALLY_HIGH)
        
    
    # Maximum noise variance taken from a range of section where the
    # variance is highest = where an artefact happens.
    # The median of all variances is returned to elimiate outliers.
    def calcNoiseVarMax(self):
        #ARTEFACTS beginning / stop
        tbeginVideo=self.artefact[:,0]
        tendVideo=self.artefact[:,1]
        dt=self.zero_data-self.zero_video
        tbegin=tbeginVideo+dt
        tend=tendVideo+dt
        maxVarList=[]

        for i in range(len(tbegin)):
            t1=tbegin[i]
            t2=tend[i]
            t1=int(self.fs*t1)
            t2=int(self.fs*t2)
            signalWithArtefact=self.eeg[t1:t2]
            artefactVar = (np.var(signalWithArtefact) - self.pureEEGVar)  / self.noiseReduction
            maxVarList.append(artefactVar)

        self.noiseVarMax = np.median(maxVarList)
        if (self.noiseVarMax < 0):
            raise self.NoiseWallException(self.MAX_VAR_NEG)

    # Calculates the noise uncertainty as the ratio between the
    # min variance and the max variance.
    def calcRho(self):
        self.rho = np.sqrt( self.noiseVarMax / self.noiseVarMin )

    # Calculates the noise wall in decibel
    def calcNoiseWall(self):
        if (self.rho < 1):
            raise self.NoiseWallException(self.MIN_VAR_LARGER_THAN_MAX_VAR)
        self.SNRwall = 10 * math.log10(self.rho - 1/self.rho)

    # Calculates the SNR in decibel
    def calcSNR(self):
        noiseVariance = self.noiseVarMin * self.rho
        self.SNR= 10 * math.log10( self.pureEEGVar / noiseVariance )

    # Do all calculations in one go
    def doAllCalcs(self):
        if not self.dataok:
            raise self.NoiseWallException(self.DATA_INVALID)
        self.calcNoiseVarMin()
        self.calcNoiseVarMax()
        self.calcRho()
        self.calcNoiseWall()
        self.calcSNR()

    # Returns the SNR wall. Only if the average SNR is higher EEG can
    # be detected at all.
    def getSNRwall(self):
        return self.SNRwall

    # Returns the average SNR of the signal.
    def getSNR(self):
        return self.SNR

