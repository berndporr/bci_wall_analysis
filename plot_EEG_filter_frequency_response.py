from NoiseWall import NoiseWall
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal

band_low = 5
band_high = 20

noiseWall = NoiseWall(20,"jaw")
noiseWall.filterData(band_low,band_high)
plt.plot(noiseWall.eegFilterFrequencyResponse**2)
plt.xlabel("Hz")
plt.ylabel("gain (linear)")
plt.title("EEG filter response")

plt.figure()
noiseWall = NoiseWall(20,"jaw")
noiseWall.filterData(1,95)
plt.plot(noiseWall.eegFilterFrequencyResponse**2)
plt.xlabel("Hz")
plt.ylabel("gain (linear)")
plt.title("EEG filter response")


plt.figure()
noiseWall = NoiseWall(20,"jaw")
noiseWall.filterData(-1)
plt.plot(noiseWall.eegFilterFrequencyResponse**2)
plt.xlabel("Hz")
plt.ylabel("gain (linear)")
plt.title("EEG filter response")


plt.show()
