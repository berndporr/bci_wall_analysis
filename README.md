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
 - -p participant number between 1 and 20 (default 1)
 - -s The second at which the data analysis starts (default is zero)
 - -t The task which can be: "jawclench", "read", "colour", "wordsearch", "sudoku", "phoneApp", "lyingEC", "lyingEO"
 - -a [minF] -b [maxF] the optional frequency band when detecting bandlimited events (default: full spectrum)
 - -h Help

All parameters are optional and have default values set.

## explore.py

Plots the data in the timedomain and frequency domain

## snr.py

Calculates the signal to noise ratio from one subject doing one task

## p300.py

Calculates the signal power of the P300 peak which is used in the SNR calculations.

## noisewall.py

Calculates the noise wall from the data.

## snr_walls_with_t_test.py

Calculates all SNR and noise wall values for all task over all subjects.

Additional pre-set parameters for the frequency bands:
 -  -w: 4..35 Hz
 -  -n: 8..18 Hz
 -  -a: 8..12 Hz
 -  -l: 1..20 Hz
Default is full range (engery detector)

# Credits

Bernd Porr
Lucía Muñoz Bohollo
