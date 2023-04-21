# BCI-walls

These scripts calculate SNR and the SNR-walls of EEG

[![DOI](https://zenodo.org/badge/146623395.svg)](https://zenodo.org/badge/latestdoi/146623395)

## Dataset

Download the dataset from https://researchdata.gla.ac.uk/1258/
and unpack the zip outwidth of this repository so that the directory
structure looks like this:
```
--+--bci_wall_analysis
  |
  +--gla_researchdata_1258
```

## Pre-requisites

The following python modules are required:
 - matplotlib
 - numpy
 - scipy

## Commandline parameters

All python commandline tools have the same parameters:
 - -p participant number between 1 and 20 (default is 20)
 - -s The second at which the data analysis starts (default is 60sec)
 - -t The task which can be: "jawclench", "read", "colour", "wordsearch", "sudoku", "phoneApp", "lyingEC", "lyingEO"
 - -h Help

All parameters are optional and have default values set.

## `explore.py`

Plots the data in the timedomain and frequency domain of one subject and for all tasks (or one task if specified).

## `snr.py`

Calculates the signal to noise ratio from one subject doing one task. It can be used both as a module and main program.

## `p300.py`

Plots the signal power of the P300 peak which is used in the SNR calculations.

## `noise_wall.py`

Calculates the SNR-wall. It can be used both as a module and main program.

## `bci_wall_one_subj_analysis.py`

Calculates the SNR-walls and SNRs for one subject.

## `bci_walls_with_t_test.py`

Calculates all SNR and SNR-wall values for all task over all subjects and
performs t-tests if it is significantly possible to detect conscious EEG changes.

Pre-set parameters for the frequency bands:
 -  -f: full range
 -  -m: 8-18 Hz (used for BCI for motor imagination)
 -  -n: 8-12 Hz (alpha frequency range used for BCI)
 -  -d: derivative (1st order highpass used for BCI)
 -  -e: 0.1-3 Hz (eyeblink frequency range)

Alternatively, set the frequency range yourself:
 - -a [minF] -b [maxF] frequency band when detecting bandlimited events (default: full spectrum)

## `do_all_stats.sh`

Runs five `bci_walls_with_t_test.py` scripts for all pre-set parameters in separate
processes.

# Credits

Bernd Porr

Lucía Muñoz Bohollo
