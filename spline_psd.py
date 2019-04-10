import sys
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import numpy as np
a = np.loadtxt(sys.argv[1])
x = a[:,0]
y = a[:,1]
f = interp1d(x, y, kind='cubic')


psd = interp1d(x, y, kind='cubic')
bandpower = 0
for f2 in np.arange(1,90,1.0):
    bandpower = bandpower + 10**psd(f2)
print("Average EEG voltage is (uV):",round((bandpower**0.5)*1E6))

x2 = np.linspace(1,90,100)
plt.plot(x2,f(x2))
plt.ylim([-16,-10])
plt.xlim([0,100])
plt.show()
