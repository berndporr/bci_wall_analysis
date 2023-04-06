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

desyncPowerReductionPercent = 40 # percent

class SNR:
    def __init__(self,subj,task,startsec=False,minF=False,maxF=False):
        self.subj = subj
        self.task = task
        self.startsec = startsec
        self.minF = minF
        self.maxF = maxF

    def calcNoisePower(self):
        task = researchdata1258.Tasks(self.subj,self.task,band_low=self.minF,band_high=self.maxF)
        y = task.ch1
        return np.var(y)

    def calcP300power(self):
        ep = researchdata1258.Evoked_potentials(self.subj)
        t,p300 = ep.get_averaged_ep()
        p300peak = p300[int(ep.Fs*VEPstartTime):int(ep.Fs*VEPendTime)]
        idx = np.argmax(p300peak)
        return p300peak[idx]**2        

    def calcSignalPower(self):
        if self.minF and self.maxF:
            p300full = self.calcP300power()
            p300reduced = p300full * (100-desyncPowerReductionPercent)/100
            pSignal = p300full - p300reduced
            print("Conscious {}% reduction (or boost) of power (desync): {} - {} = {}"
                  .format(desyncPowerReductionPercent,p300full,p300reduced,pSignal))
            return pSignal
        else:
            print("Power from P300")
            return self.calcP300power()

    def calcSNR(self,band_low=False,band_high=False):
        NoisePwr = self.calcNoisePower()
        SignalPwr = self.calcSignalPower()
        print("Signal Power = {} V^2, Amplitude = {} uV".format(SignalPwr,(SignalPwr**0.5)*1E6))
        print("Noise Power = {}V^2, Noise amplitude = {}V".format(NoisePwr,np.sqrt(NoisePwr)))
        self.snrvalue = np.log10(SignalPwr/NoisePwr)*10

    
# check if we run this as a main program
if __name__ == "__main__":
    subj = 20
    task = researchdata1258.Tasks.TASKS[0]
    a = False
    b = False

    helptext = 'usage: {} -p participant -s startsec -t task -a [min Hz] -b [max Hz] -h'.format(sys.argv[0])

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
                a = float(arg_val)
            elif '-b' in opt:
                b = float(arg_val)
            elif '-t' in opt:
                task = arg_val
            elif '-h' in opt:
                raise getopt.GetoptError(helptext)
            else:
                raise getopt.GetoptError(helptext)
    except getopt.GetoptError as err:
        print (err)
        sys.exit(2)

    snr = SNR(subj,task,minF=a,maxF=b)
    snr.calcSNR()
    print("Subject: {}, Task: {}, SNR= {} dB".format(subj,task,snr.snrvalue))
