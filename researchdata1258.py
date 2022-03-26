import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as signal

class Tasks:
    TASKS = ["jawclench", "read", "colour", "wordsearch", "sudoku", "phoneApp", "lyingEC", "lyingEO"]
    Fs = 500
    
    def __init__(self,_participant,_task,startsec=False,endsec=False,filterData=True):
        '''
        _participant is the participant number
        _task is one of the tasks from the TASKS array
        filterData if set to true the EEG data is filtered:
        Highpass filter at 1Hz to remove Dc.
        Notch filter betwwen 48 and 52Hz to remove 50Hz noise.
        Bandstop filter betwwen 148 and 152Hz to remove 150Hz interference.
        
        Butterworth filter used. Show flat frequency response in the passband.
        '''
        self.participant = _participant
        self.task = _task
        fullpath = "../gla_researchdata_1258/EEG_recordings/participant{:03d}/{}.tsv".format(self.participant,self.task)
        self.data = np.loadtxt(fullpath)
        
        self.ch1 = self.data[:,7]
        self.ch2 = self.data[:,8]
        self.t = np.linspace(0,len(self.ch1)/self.Fs,len(self.ch1))
        if not filterData:
            return

        # Remove DC
        bHigh,aHigh = signal.butter(4,1/self.Fs*2,'high')
        self.ch1 = signal.lfilter(bHigh,aHigh,self.ch1);
        self.ch2 = signal.lfilter(bHigh,aHigh,self.ch2);
        
        # Remove 50Hz noise
        b50,a50 = signal.butter(4,[48/self.Fs*2,52/self.Fs*2],'stop')
        self.ch1 = signal.lfilter(b50,a50,self.ch1);
        self.ch2 = signal.lfilter(b50,a50,self.ch2);
        
        # Remove 150Hz interference
        b150,a150 = signal.butter(4,[148/self.Fs*2,152/self.Fs*2],'stop')
        self.ch1 = signal.lfilter(b150,a150,self.ch1);
        self.ch2 = signal.lfilter(b150,a150,self.ch2);

        self.startsec = startsec
        if startsec:
            a = startsec * self.Fs
        else:
            a = 10 * self.Fs
        if endsec:
            b = endsec * self.Fs
        else:
            b = -1

        self.ch1 = self.ch1[a:b]
        self.ch2 = self.ch2[a:b]


class Evoked_potentials:
    Fs = 250
    HPfc = 0.5
        
    def __init__(self,_participant,startsec=False,do_filter_data = True):
        """
        Loads the P300 or VEP of one Participant.
        _participant is the integer number of the participant
        _ep should be taken from self.EPs and is either "VEP" or "P300"
        The following instance variables are then available:
        self.t is an array of timestamps of the length of the recording
        self.eeg is the eeg data. If do_filter_data is True it's filtered.
        self.oddball_flags contains either 0 or 1 and is one if an evoked potential has been triggered
        self.oddball_samples contains the sample numbers where a an oddball has been triggered
        self.initial_samples_to_ignore is the number of samples to ignore because of filtering
        """
        self.participant = _participant
        fullpath = "../gla_researchdata_1258/EEG_recordings/participant{:03d}/rawp300.tsv".format(self.participant)
        self.data = np.loadtxt(fullpath)
        self.startsec = startsec
        if startsec:
            a = startsec * self.Fs
        else:
            a = 2 * self.Fs
        self.eeg = self.data[a:,0]
        self.oddball_flags = self.data[a:,2]
        self.oddball_samples = np.argwhere(self.oddball_flags > 0.5)
        self.initial_samples_to_ignore = 0
        self.t = np.linspace(0,len(self.data)/self.Fs,len(self.data))
        if do_filter_data:
            self.__filter_data()
        
    def __filter_data(self):
        # Remove DC
        bHigh,aHigh = signal.butter(2,0.25/self.Fs*2,'high')
        self.eeg = signal.lfilter(bHigh,aHigh,self.eeg);
        # for VEP
        self.initial_samples_to_ignore = int(self.Fs / self.HPfc) * 3

        # Remove 50Hz noise
        b50,a50 = signal.butter(4,[48/self.Fs*2,52/self.Fs*2],'stop')
        self.eeg = signal.lfilter(b50,a50,self.eeg);

        # Remove 150Hz interference
        b100,a100 = signal.butter(4,[98/self.Fs*2,102/self.Fs*2],'stop')
        self.eeg = signal.lfilter(b100,a100,self.eeg);

        
    def get_averaged_ep(self):
        """
        Calculates the evoked potential usign the oddball samples and
        averages over them.
        Returns: timestamps,vep
        """
        self.navg = int(self.Fs)
        self.avg = np.zeros(self.navg)
        
        n = 0
        for [ob] in self.oddball_samples:
            if ((ob+self.navg) < len(self.eeg)) and (ob > self.initial_samples_to_ignore):
                self.avg = self.avg + self.eeg[int(ob):int(ob+self.navg)]
                n = n + 1
                
        avg = self.avg / n
        
        time = np.linspace(0,self.navg/self.Fs,self.navg) * 1000
        return time,avg

    def get_stimulus_stats(self):
        dt = np.array([])
        for i in range(len(self.oddball_samples)-1):
            dt = np.append(dt,(self.oddball_samples[i+1] - self.oddball_samples[i])/self.Fs)
        return np.mean(dt),np.std(dt),int(np.round(min(dt))),int(np.round(max(dt)))

