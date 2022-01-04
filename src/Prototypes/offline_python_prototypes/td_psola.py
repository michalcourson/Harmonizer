#modified from https://github.com/sannawag/TD-PSOLA
import numpy as np
from numpy.fft import fft, ifft
import matplotlib.pyplot as plt
from scipy.io import wavfile

rollover = np.zeros(1024)
frame_offset = 0
peak_offset = 0
last_P1 = -1
period_memory = []
period_voiced = []
period_graph = []
periods = [40,40,40]
NOISE_GATE_THRESHOLD = 500

def signalGen(freq, Fs, dur):
    wav = np.zeros(Fs*dur)
    t = np.linspace(0, dur, Fs*dur)
    for i in range(t.size):
        wav[i] = np.cos(2*np.pi*freq*t[i])
    return wav

def shift_pitch(signal, fs, f_ratio):
    peaks = find_peaks(signal, fs)
    new_signal = psola(signal, peaks, f_ratio)
    return new_signal


def find_peaks(signal, fs, max_hz=950, min_hz=75, analysis_win_ms=40, max_change=1.5, min_change=0.666):
    global peak_offset
    global period_memory
    global periods
    global NOISE_GATE_THRESHOLD
    N = len(signal)
    min_period = fs // max_hz
    max_period = fs // min_hz

    # compute pitch periodicity
    sequence = int(1024) 
    periods.pop(0)
    periods.append(compute_periods_per_sequence(signal[2048:], sequence, min_period, max_period)[0])
    

    # simple hack to avoid octave error: assume that the pitch should not vary much, restrict range -- could ax if needed
    mean_period = np.mean(periods)
    max_period = int(mean_period * 1.1)
    min_period = int(mean_period * 0.9)
    periods[2] = compute_periods_per_sequence(signal[2048:], sequence, min_period, max_period)[0]

    #append period to memory
    period_memory.append(periods[2])

    #segment out unvoiced segments
    if len(period_memory) > 5:
        mean_period = np.mean(period_memory[-5:])
        for i, p in enumerate(periods):
            if abs(p - mean_period) > .1*mean_period or np.max(signal[1024*i:1024*(i+1)]) < NOISE_GATE_THRESHOLD:
                periods[i] = 0

    #append voiced/unvoiced period to array
    period_voiced.append(periods[2])

    #correct unvoiced
    for i, p in enumerate(periods):
        if p == 0:
            periods[i] = period_voiced[np.max(np.nonzero(period_voiced))]

    #moving average
    periods[2] = (periods[1] + periods[2]) >>1

    #append to graph list for debugging
    period_graph.append(periods[2])
    
    #get peaks
    peaks = [np.argmax(signal[:int(periods[0]*1.1)])]
    while True:
        prev = peaks[-1]
        idx = prev // sequence  # current autocorrelation analysis window
        if prev + int(periods[idx] * max_change) >= N:
            break
        # find maximum near expected location
        peaks.append(prev + int(periods[idx] * min_change) +
                np.argmax(signal[prev + int(periods[idx] * min_change): prev + int(periods[idx] * max_change)]))

    #correct peaks for phase consistency
    if peak_offset != 0:
        for peak in peaks:
            if peak > 1024:
                peaks += (peak_offset - peak)
                break
    for peak in peaks:
        if peak > 2048:
            peak_offset = peak - 1024
            break
    return np.array(peaks)



#maybe replace with pYIN (https://code.soundsoftware.ac.uk/projects/pyin)
def compute_periods_per_sequence(signal, sequence, min_period, max_period):
    offset = 0  # current sample offset
    periods = []  # period length of each analysis sequence

    while offset < len(signal):
        fourier = fft(signal[offset:offset+sequence])
        fourier[0] = 0  # remove DC component
        autoc = ifft(fourier * np.conj(fourier)).real
        autoc_peak = min_period + np.argmax(autoc[min_period: max_period])
        periods.append(autoc_peak)
        offset += sequence
    return periods


def psola(signal, peaks, f_ratio):
    global rollover
    global frame_offset
    global last_P1
    N = len(signal)

    #new peaks the easy way
    new_peaks = []
    for i in range(len(peaks)):
        if peaks[i] > 1024:
            new_period = (peaks[i+1] - peaks[i]) / f_ratio
            break
    new_peaks = np.arange(1024 + frame_offset, 3072 + frame_offset, new_period)
    new_peaks = np.asarray(new_peaks, np.int16)


    # PSOLA
    new_signal = np.zeros(N)
    new_signal[1024:1024+len(rollover)] += rollover
    for j in range(len(new_peaks)):
        if new_peaks[j] < 1024: #could just get rid of those some values
            continue
        if new_peaks[j] > 2048:
            rollover = np.copy(new_signal[2048:])

        # find the corresponding old peak index
        i = np.argmin(np.abs(peaks - new_peaks[j]))

        # get the window length
        if last_P1 != -1:
            if f_ratio > 1:
                P1 = last_P1 if j == 0 else new_peaks[j] - new_peaks[j-1]
            else:
                P1 = last_P1 if j == 0 else peaks[i] - peaks[i-1]
        else:
            if f_ratio > 1:
                P1 = new_peaks[j+1] - new_peaks[j] if j == 0 else new_peaks[j] - new_peaks[j-1]
            else:
                P1 = peaks[i+1] - peaks[i] if i == 0 else peaks[i] - peaks[i-1]
        
        #adjust for edge cases
        if peaks[i] + P1 >= 3072:
            P1 = 3072 - peaks[i]
        if new_peaks[j] + P1 >= 3072:
            P1 = 3072 - new_peaks[j]
        
            
        #window and add
        if peaks[i] - P1 < 0:
            sig_to_window = signal[peaks[i+1] - P1:peaks[i+1] + P1]
        else:
            sig_to_window = signal[peaks[i] - P1:peaks[i] + P1]
        window_len = len(sig_to_window)
        window = np.hanning(window_len)
        

        if new_peaks[j] - P1 < 0:
            new_signal[0: new_peaks[j] + P1] += (window * sig_to_window)[int(np.ceil(window_len/2)) - new_peaks[j]:]
        elif new_peaks[j] + P1 >= 3072:
            new_signal[new_peaks[j] - P1:3072] += (window * sig_to_window)[0:len(new_signal[new_peaks[j] - P1:3072])]
        else:
            new_signal[new_peaks[j] - P1: new_peaks[j] + P1] += window * sig_to_window
        
        #get the frame to return, as well as set values for phase consistency
        if(new_peaks[j]) > 2048 and j != 0:
            ret_sig = new_signal[1024: 2048]
            frame_offset = new_peaks[j] - 2048
            last_P1 = P1
            break
        
    return ret_sig


if __name__=="__main__":
    # Load audio
    
    fs = 44100

    #TODO: MODIFY THE SETUP CODE ACCORDING TO COMMENTS BELOW.
    
    #IF YOU WANT TO USE A GENERATED SIGNAL, UNCOMMENT THE TWO LINES BELOW
    #freq = 200                               #TODO: SET FREQUENCY OF SIGNAL GENERATOR
    #orig_signal = signalGen(freq, fs, 2)    #TODO: IF YOU WANT TO TEST A SINUSOID, USE THIS. ARG 1 IS NOTE FREQUENCY, ARG 2 IS SAMPLE FREQUENCY (44.1kHz), ARG 3 IS SIGNAL LENGTH IN SECONDS
    
    #IF YOU WANT TO USE AN INPUT WAV FILE, UNCOMMENT THESE LINES
    #fs, orig_signal = wavfile.read("<INSERT FILE NAME HERE>")   #TODO: OTHERWISE, SUBSTITUTE THE FILENAME WITH SOME OTHER FILE
    
    #Convert stereo signals to mono
    if orig_signal.shape[1] == 2:
        orig_signal = (orig_signal[...,0] + orig_signal[...,1])/2
    
    #Capture the length of the signal
    N = len(orig_signal)


    #TODO: CHANGE THE VALUE TO SHIFT TO. THIS IS THE PITCH SHIFT AMOUNT AS A RATIO
    f_ratio = 2 ** (-12 / 12)
        

    # Shift pitch
    orig_signal = orig_signal[0:1024*int(len(orig_signal)*1024)]
    np.insert(orig_signal, 0, np.zeros(1024))
    #orig_signal.prepend(zeros(1024))
    #new_signal = shift_pitch(orig_signal, fs, f_ratio)
    new_signal = np.empty(len(orig_signal))
    for i in np.arange(0,len(orig_signal) - 3072, 1024):
        new_signal[i:i+1024] = shift_pitch(orig_signal[i:i+3072],fs, f_ratio)

    # Plot
    plt.style.use('ggplot')
    plt.plot(orig_signal[:-1])
    plt.show()
    plt.plot(new_signal[:-1])
    x = np.arange(0,len(new_signal), 1024)
    plt.scatter(x, new_signal[x], color='blue')
    plt.figure()
    plt.plot(period_graph)
    plt.show()
    # Write to disk
    new_signal16 = np.asarray( new_signal, np.int16)
    wavfile.write("sax2_transposed_{:01.2f}.wav".format(f_ratio), fs, new_signal16)
