	
'''
The program accepts either a .wav file as input or generates an ideal sin wave at a specific frequency 
to be used as input.  The pitch shift ratio (denoted by f_ratio), synthesis hop size and analysis hop size 
are all defined at the beginning of the process.  The program starts in the analysis stage by calling the stft() 
function in order to acquire RFFT data from the entire signal; a hanning window is applied to the input, and 
each RFFT is separated by the analysis hopsize.  

After the RFFTs have been taken, the function doit() is called and the synthesis stage begins.  This function 
computes the unwrapped-phase difference between respective frequencies bins in neighboring rows of the previously 
calculated STFT matrix.  The unwrapped-phase difference is then multiplied by the pitch shift ratio and added to 
the previous synthesized result’s phase to get the newly synthesized phase for the current iteration.  This newly 
synthesized phase, along with the magnitude data from the current analysis window, is used to compute the newly 
synthesized RFFT result.  These synthesized RFFT results get stored within the 2-dimensional synth[] array. 

When the synth[] array has been completely filled, the script then calls the istft() function in order to get the 
synthesized result into the time domain.  The istft() function takes the inverse RFFT on each row of the 2-dimensional 
synth[] array and stores them in a  1-dimensional array to be outputted.  Consecutive inverse RFFT results get 
overlapped and added by the synthesis hopsize amount.  The resulting output is the time stretched version of the input. 

The final step is to resample the time stretched output in order to perform the desired pitch shift.  This is performed 
with the sf.write() function.  The filename is the first argument, the time stretched output is the second argument, 
and sampling frequency of the output is the third argument.  The sampling frequency must be the initial sampling frequency 
of the input multiplied by the pitch shift ratio in order for the output signal to be the correct duration and pitch. 

'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import wavfile
import scipy.signal.windows as window
import scipy.io 
from scipy.fft import fft, ifft, fftshift, ifftshift
import cmath
import soundfile as sf


#generates a sine wav at the frequency freq and duration dur
def signalGen(freq, Fs, dur):
    wav = np.zeros(Fs*dur)
    t = np.linspace(0, dur, Fs*dur)
    for i in range(t.size):
        wav[i] = np.cos(2*np.pi*freq*t[i])
    return wav

#takes an stft of the signal x, win is the window applied before transforming
#over is the hopsize
def stft(x, win, over):
    numBins = int((len(x)/over) - ((len(win))/over))
    dfts = np.zeros((numBins, int(len(win)/2 + 1)), dtype = complex)
    raw = np.zeros((numBins, len(win)))

    for q in range(numBins):
        newfirstw = x[over*q:len(win)+over*q] * win
        raw[q] = newfirstw
        dfts[q] = np.fft.rfft(newfirstw)

    return dfts, raw


#takes inverse stft in the same manner as the stft function
def istft(X, win, over):
    numBins = (X.shape)[0]
    print(numBins)
    reconstructed = np.zeros(numBins*over)
    for q in range(numBins):
        if(over*q + len(win) > numBins*over):
            continue
        result = np.fft.irfft(X[q])
        result = result.real
        result = result * win
        reconstructed[over*q:over*q + len(win)] += result

    return reconstructed

#function to plot a spectrogram from stft data
def plotSpectrogram(freq, x, whichSignal, win, over, Fs):
    yaxis = (Fs/len(win))*(np.arange(len(win)/2))
    specxaxis = np.arange(int((len(x)/over) - ((len(win) - over)/over)))
    numBins = int((len(x)/over) - ((len(win) - over)/over))
    newPart = np.zeros((numBins, len(win)))

    for q in range(numBins):
        newPart[q] = 10*np.log10(fftshift(abs(freq[q])))
        
    print(newPart.shape)
    plt.pcolormesh(specxaxis*(over/Fs)*1000, yaxis/1000, (newPart[:,int(len(win)/2):int(len(win))].T), shading = 'gouraud', cmap='jet')
    if whichSignal == 0:
        plt.title("Spectrogram of STFT of input_signal Signal Using Rectangular Window")
    elif whichSignal == 1: 
        plt.title("Spectrogram of STFT of Anote Signal Using Rectangular Window")
    plt.xlabel("Time (ms)")
    plt.ylabel("Frequency (kHz)")
    plt.show()


#calculates the angle x wrapped in the range -pi to pi
def princargx(x):
    return (x + np.pi)%(2*np.pi) - np.pi

#calculates unwrapped phase difference between two stft frames
def unwrappedPhase(ph1, ph2, hopsize, Fs, k, num_bins):
    h = hopsize/Fs
    fk = (Fs*(k))/(num_bins)
    return (h*fk*2*np.pi) + princargx(ph2 - ph1 - (h*fk*2*np.pi))

#returns a phase shifted number
def phaseShift(num, angle):
    return num * complex(np.cos(angle), np.sin(angle))

#calculates new phases for every frequency bin in each stft frame
def doit(stft, analysis_hopsize, synthesis_hopsize, Fs, num_bins, f_ratio):
    synth = np.copy(stft)
    for i in range(1,stft.shape[0]):
        for q in range(stft.shape[1]):
            phase_diff = unwrappedPhase((np.angle(stft[i-1,q])), (np.angle(stft[i, q])), analysis_hopsize, Fs, q, num_bins)
            synth[i,q] = phaseShift(np.abs(synth[i,q]), np.angle(synth[i-1,q]) + phase_diff*f_ratio) 
    return synth



Fs = 44100

#TODO: MODIFY THE SETUP CODE ACCORDING TO COMMENTS BELOW.

#IF YOU WANT TO USE A GENERATED SIGNAL, UNCOMMENT THE TWO LINES BELOW
#freq = 200                               #TODO: SET FREQUENCY OF SIGNAL GENERATOR
#input_signal = signalGen(freq, Fs, 2)    #TODO: IF YOU WANT TO TEST A SINUSOID, USE THIS. ARG 1 IS NOTE FREQUENCY, ARG 2 IS SAMPLE FREQUENCY (44.1kHz), ARG 3 IS SIGNAL LENGTH IN SECONDS

#IF YOU WANT TO USE AN INPUT WAV FILE, UNCOMMENT THESE LINES
#Fs, input_signal = wavfile.read("<INSERT FILE NAME HERE>")   #TODO: OTHERWISE, SUBSTITUTE THE FILENAME WITH SOME OTHER FILE

#Convert stereo signals to mono
if input_signal.shape[1] == 2:
    input_signal = (input_signal[...,0] + input_signal[...,1])/2

#TODO: YOU CAN PLAY WITH DIFFERENT VALUES OF WINDOW SIZE. WE SET IT TO 2048 HERE, BUT TO 1024 IN THE C++ IMPLEMENTATIONS.
window_size = 2048

#Perform test shifts for -12 to 13 semitones from base frequency of note.
semitones = list(range(-12,13))

for s in semitones:
    f_ratio = 2 ** (s / 12)
    synthesis_hopsize = int(window_size/4)
    analysis_hopsize = int(synthesis_hopsize/f_ratio)
    rectwind = np.ones(window_size)
    hannwind = window.hann(window_size)
    hop = int((hannwind.size)/2)
    stftinput_signal, raw = stft(input_signal, hannwind, analysis_hopsize)


    synth = doit(stftinput_signal, analysis_hopsize,synthesis_hopsize , Fs, len(hannwind), f_ratio)
    output = istft(synth, hannwind, synthesis_hopsize)
    sf.write("shifted_{}.wav".format(s+12), output, int(Fs*f_ratio))

