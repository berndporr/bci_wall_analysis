import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal
import math as math
from scipy.interpolate import interp1d
import scipy.stats as stats

class ParalysedEEG:
    def __init__(self,f_min,f_max):
        self.f_min = f_min
        self.f_max = f_max

    def readParalysedEEGVarianceFromWhithamEtAl(self,filename):
        a = np.loadtxt(filename)
        f = a[:,0]
        p = a[:,1]
        psd = interp1d(f, p, kind='cubic')
        self.bandpower = 0
        for f2 in np.arange(int(self.f_min),int(self.f_max)):
            self.bandpower = self.bandpower + ( 10**psd(f2) )
        return self.bandpower

    def getPureEEGVar(self):
        self.pureEEGVar = 0
        self.pureEEGVar = self.pureEEGVar + self.readParalysedEEGVarianceFromWhithamEtAl("sub1a.dat")
        self.pureEEGVar = self.pureEEGVar + self.readParalysedEEGVarianceFromWhithamEtAl("sub1b.dat")
        self.pureEEGVar = self.pureEEGVar + self.readParalysedEEGVarianceFromWhithamEtAl("sub1c.dat")
        self.pureEEGVar = self.pureEEGVar + self.readParalysedEEGVarianceFromWhithamEtAl("sub2a.dat")
        self.pureEEGVar = self.pureEEGVar + self.readParalysedEEGVarianceFromWhithamEtAl("sub2b.dat")
        self.pureEEGVar = self.pureEEGVar + self.readParalysedEEGVarianceFromWhithamEtAl("sub2c.dat")
        self.pureEEGVar = self.pureEEGVar / 6.0
        return self.pureEEGVar
