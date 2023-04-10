# BCI-walls

These scripts calculate the SNr and the SNR-walls

## Dataset

Download the dataset from https://researchdata.gla.ac.uk/1258/
and unpack the zip outwidth of this repository so that the directory
structure looks like this:
```
--+--deepNeuronalFilter/eeg_filter
  |
  +--gla_researchdata_1258
```

## Pre-requisites

The following python modules are required:
 - matplotlib
 - numpy
 - plotly
 - scipy

## Commandline parameters

All python commandline tools have the same parameters:
 - -p participant number between 1 and 20 (default is 20)
 - -s The second at which the data analysis starts (default is 60sec)
 - -t The task which can be: "jawclench", "read", "colour", "wordsearch", "sudoku", "phoneApp", "lyingEC", "lyingEO"
 - -a [minF] -b [maxF] the optional frequency band when detecting bandlimited events (default: full spectrum)
 - -h Help

All parameters are optional and have default values set.

## explore.py

Plots the data in the timedomain and frequency domain.

## snr.py

Calculates the signal to noise ratio from one subject doing one task.

## p300.py

Calculates the signal power of the P300 peak which is used in the SNR calculations.

## noise_wall.py

Calculates the SNR wall.

## bci_wall_one_subj_analysis.py

Calculates the SNR-walls and SNRs for one subject (default = 20).

## bci_walls_with_t_test.py

Calculates all SNR and noise wall values for all task over all subjects and performs t-tests if it is significantly possible to detect conscious EEG changes.

Additional pre-set parameters for the frequency bands:
 -  -m: 8-18 Hz (used for BCI)
 -  -n: 8-12 Hz (alpha frequency range used for BCI)
 -  -d: derivative (1st order highpass used for BCI)
 -  -f: 0.1-fs/2 (full range)
 -  -e: 0.1-3 Hz (eyeblink frequency range)

# Credits

Bernd Porr
Lucía Muñoz Bohollo
