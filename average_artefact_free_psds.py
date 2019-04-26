from NoiseWall import NoiseWall
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal

experiments = ["lie_relax","blink","eyescrunching","raisingeyebrows","jaw","readinglieing"]

fs = 1000
favgeeg = np.zeros(fs)
n = 0

for e in experiments:
    print("\n"+e)
    for subj in range(15,28):
        print(subj)
        noiseWall = NoiseWall(subj,e)
        if noiseWall.dataok:
            eeg = noiseWall.getMinNoiseVarEEGChunk()
            feeg = np.fft.fft(eeg) / len(eeg)
            feeg = np.abs(feeg)
            feeg = feeg*feeg
            feeg = np.log10(feeg)
            feeg = signal.resample(feeg,fs)
            favgeeg = favgeeg + feeg
            n = n + 1
            plt.plot(feeg)
            plt.xlim([0,100])
            plt.ylim([-16,-10])

plt.figure(2)
plt.plot(favgeeg / n)
plt.xlim([0,100])
plt.ylim([-16,-10])


plt.show()
