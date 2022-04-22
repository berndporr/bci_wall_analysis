#!/usr/bin/python3
import matplotlib.pyplot as plt
import numpy as np
import sys
import p300
import scipy.signal as signal
import getopt
import os
import researchdata1258

VEPstartTime = 0.35 # sec
VEPendTime = 0.45 # sec

class NoiseWall:
    class NoiseWallException(Exception):
        pass
    
    def __init__(self,subj,task,startsec=False,minF=False,maxF=False):
        self.subj = subj
        self.task = task
        self.startsec = startsec
        self.minF = minF
        self.maxF = maxF
    
    # Calculates the noise uncertainty as the ratio between the
    # min variance and the max variance.
    def calcRho(self):
        task = researchdata1258.Tasks(self.subj,self.task,band_low=self.minF,band_high=self.maxF)
        y = task.ch1
        winSize = 1000
        self.noiseVarMin = 1E10
        self.noiseVarMax = -1E10
        for i in range(len(y)-winSize-1):
            v = np.var(y[i:i+winSize])
            if (v > self.noiseVarMax) and (v > 0):
                self.noiseVarMax = v
            if (v < self.noiseVarMin) and (v > 0):
                self.noiseVarMin = v
        print(self.noiseVarMin,self.noiseVarMax)
        self.rho = np.sqrt( self.noiseVarMax / self.noiseVarMin )

    # Calculates the noise wall in decibel
    def calcNoiseWall(self):
        self.calcRho()
        if (self.rho < 1):
            raise self.NoiseWallException(self.MIN_VAR_LARGER_THAN_MAX_VAR)
        self.SNRwall = 10 * np.log10(self.rho - 1/self.rho)




    
# check if we run this as a main program
if __name__ == "__main__":
    subj = 1
    task = researchdata1258.Tasks.TASKS[0]
    minF = False
    maxF = False

    helptext = 'usage: {} -p participant -s startsec -t task -h'.format(sys.argv[0])

    try:
        # Gather the arguments
        all_args = sys.argv[1:]
        opts, arg = getopt.getopt(all_args, 'p:s:t:a:b:h')
        # Iterate over the options and values
        for opt, arg_val in opts:
            if '-p' in opt:
                subj = int(arg_val)
            elif '-s' in opt:
                startsec = int(arg_val)
            elif '-a' in opt:
                minF = int(arg_val)
            elif '-b' in opt:
                maxF = int(arg_val)
            elif '-t' in opt:
                task = arg_val
            elif '-h' in opt:
                raise getopt.GetoptError(helptext)
            else:
                raise getopt.GetoptError(helptext)
    except getopt.GetoptError as err:
        print (err)
        sys.exit(2)

    noisewall = NoiseWall(subj,task,minF=minF,maxF=maxF)
    noisewall.calcNoiseWall()
    
    print("Subject: {}, Task: {}, Wall= {}dB".format(subj,task,noisewall.SNRwall))
