#!/usr/bin/python3
import matplotlib.pyplot as plt
import numpy as np
import sys
import p300
import scipy.signal as signal
import getopt
import os
import researchdata1258

subjectsOK = [1,3,4,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]

SNRbandMin = 5 # Hz
SNRbandMax = 100 # Hz

VEPstartTime = 0.3 # sec
VEPendTime = 0.5 # sec

class SNR:
    def __init__(self,subj,task,startsec=False):
        self.subj = subj
        self.task = task
        self.startsec = startsec
    
    def calcNoisePower(self):
        task = researchdata1258(self.subj,self.task)
        y = task.ch1
        return np.var(y)

    def calcSNR(self):
        NoisePwr = self.calcNoisePower()
        vep = p300.calcVEP(self.subj,self.startsec)
        SignalPwr = np.median(vep[int(self.fs*VEPstartTime):int(self.fs*VEPendTime)]**2)
        print("Signal Power:",SignalPwr)
        print("NoisePwr:",NoisePwr)
        snr = np.log10(SignalPwr/NoisePwr)*10
        return snr

def calcAllSNRimprovemements(startsec,
                             noisefolder,
                             fs,
                             filtered_filename):
    beforeArray = np.array([])
    afterArray = np.array([])
    for subj in subjectsOK:
        print("Subject",subj)
        snr = SNR(subj=subj,startsec=startsec,fs=fs,folder=noisefolder,noisered_filename=filtered_filename)
        snrdnf, wdnf = snr.calcSNRdnf()
        snrinner, winner = snr.calcSNRinner()
        impr = snrdnf-snrinner
        print("SNR improvement: {} - {} = {}".format(snrinner,snrdnf,impr))
        beforeArray = np.append(beforeArray,snrinner)
        afterArray = np.append(afterArray,snrdnf)
    imprArray = afterArray - beforeArray
    snrdiff_av = np.mean(imprArray)
    snrdiff_sd = np.std(imprArray)
    return beforeArray,afterArray,snrdiff_av,snrdiff_sd


# check if we run this as a main program
if __name__ == "__main__":
    subj = 1

    helptext = 'usage: {} -p participant -s startsec -t task -h'.format(sys.argv[0])

    try:
        # Gather the arguments
        all_args = sys.argv[1:]
        opts, arg = getopt.getopt(all_args, 'p:s:f:t:')
        # Iterate over the options and values
        for opt, arg_val in opts:
            if '-p' in opt:
                subj = int(arg_val)
            elif '-s' in opt:
                startsec = int(arg_val)
            elif '-t' in opt:
                taskfolder = arg_val
            elif '-h' in opt:
                raise getopt.GetoptError()
            else:
                raise getopt.GetoptError()
    except getopt.GetoptError:
        print (helptext)
        sys.exit(2)
