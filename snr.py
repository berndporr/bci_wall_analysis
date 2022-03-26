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

class SNR:
    def __init__(self,subj,task,startsec=False):
        self.subj = subj
        self.task = task
        self.startsec = startsec
    
    def calcNoisePower(self):
        task = researchdata1258.Tasks(self.subj,self.task)
        y = task.ch1
        return np.var(y)

    def calcEPPower(self):
        ep = researchdata1258.Evoked_potentials(self.subj)
        t,p300 = ep.get_averaged_ep()
        p300peak = p300[int(ep.Fs*VEPstartTime):int(ep.Fs*VEPendTime)]
        return np.median(p300peak**2)        

    def calcSNR(self):
        NoisePwr = self.calcNoisePower()
        SignalPwr = self.calcEPPower()
        print("Signal Power:",SignalPwr)
        print("NoisePwr:",NoisePwr)
        self.snrvalue = np.log10(SignalPwr/NoisePwr)*10

    
# check if we run this as a main program
if __name__ == "__main__":
    subj = 1
    task = researchdata1258.Tasks.TASKS[0]

    helptext = 'usage: {} -p participant -s startsec -t task -h'.format(sys.argv[0])

    try:
        # Gather the arguments
        all_args = sys.argv[1:]
        opts, arg = getopt.getopt(all_args, 'p:s:t:')
        # Iterate over the options and values
        for opt, arg_val in opts:
            if '-p' in opt:
                subj = int(arg_val)
            elif '-s' in opt:
                startsec = int(arg_val)
            elif '-t' in opt:
                task = arg_val
            elif '-h' in opt:
                raise getopt.GetoptError()
            else:
                raise getopt.GetoptError()
    except getopt.GetoptError:
        print (helptext)
        sys.exit(2)

    snr = SNR(subj,task)
    snr.calcSNR()
    print("Subject: {}, Task: {}, SNR= {}dB".format(subj,task,snr.snrvalue))
