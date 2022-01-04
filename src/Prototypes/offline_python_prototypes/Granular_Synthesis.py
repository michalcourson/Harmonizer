# -*- coding: utf-8 -*-
"""
Created on Sat Oct 23 15:24:40 2021

GRANULAR SYNTHESIS in python

@purpose
The focus is on real-time pitch shifting implementation using python
Real time implying low-latency and with a live signal as opposed to recorded
The algorithm will be converted to C++ to port/send to an  microcontroller board
The board is to be used bare metal (no operating system or anything)
At the moment the  teensy is  going to be used to interface the MIDI with the code

"""

import struct
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import wavfile
import wave
import sys

###############################################
#               steps
###############################################
#1)resample original signal
#    * figure out length of singal
#    *figure out number of samples
#    *OR samples per time frame
#2)buffer
#3)window (tapered) /corssfade apply
##################################################

        #Functions
        
#This function is simple fractional resampling as a linear interpolation between the current and next sample
def get_subsample(x, t):
    n = int(t)
    a = 1.0 - (t - n)
    try:
        return a * x[n] + (1 - a) * x[n + 1] 
    except IndexError:
        try:
            return a * x[n]
        except IndexError:
            return 0
        
#This is running through the original signal and sampling it at a desired rate realtive to the wanted pitch 
def resample(x, factor):
    # length of the output signal after resampling
    n_out = int(np.floor(len(x) / factor))
    y = np.zeros(n_out)
    for n in range(0, n_out):
        y[n] = get_subsample(x, float(n) * factor)
    return y

#This is converting miliseconds into a number of samples, which is merely a formula
def milliseconds2samples(ms, sampleRate):
    ms2Samp = int(float(sampleRate) * float(ms) / 1000.0)
    return ms2Samp
        

#This is the tapering window to be used for the grains for smooth overlap
def win_taper(N, overlap):
    R = int(N * overlap / 2)
    r = np.arange(0, R) / float(R)
    win = np.r_[r, np.ones(N - 2*R), r[::-1]]
    stride = N - R - 1
    return win, stride

#This is the overlap between grainsusing a tapering envelope/window
def GS_pshift(x, factor, grain_size, overlap=0.5):
    N = len(x)
    y = np.zeros(N)
    # size of input buffer given target ouptut grain size and resampling factor
    input_chunk_size = int(grain_size * factor + 0.5)
    win, stride = win_taper(grain_size, overlap)
    for n in range(0, len(x) - max(input_chunk_size, grain_size), stride):
        w = resample(x[n:n+input_chunk_size], factor)
        y[n:n+grain_size] += w * win
    return y

#This creates a cosine signal at a fixed frequency
def signalGen(freq, Fs, dur):
    wav = np.zeros(Fs*dur)
    t = np.linspace(0, dur, Fs*dur)
    for i in range(t.size):
        wav[i] = np.cos(2*np.pi*freq*t[i])
    return wav

                        #Main Code
 #TODO: MODIFY THE SETUP CODE ACCORDING TO COMMENTS BELOW.
    
#IF YOU WANT TO USE A GENERATED SIGNAL, UNCOMMENT THE TWO LINES BELOW
freq = 200                               #TODO: SET FREQUENCY OF SIGNAL GENERATOR
sampleRate= 48000
dataArray = signalGen(freq,sampleRate,2)    #TODO: IF YOU WANT TO TEST A SINUSOID, USE THIS. ARG 1 IS NOTE FREQUENCY, ARG 2 IS SAMPLE FREQUENCY (44.1kHz), ARG 3 IS SIGNAL LENGTH IN SECONDS
plotType =0     #TODO: THIS IS A FLAG TO DO THE CORRESPONDING TYPE OF PLOTTING
    
#IF YOU WANT TO USE AN INPUT WAV FILE, UNCOMMENT THE FOLLOWING LINES
#wavFileName ="<INSERT FILE NAME HERE>"
#sampleRate, dataArray = wavfile.read(wavFileName)   #TODO: OTHERWISE, SUBSTITUTE THE FILENAME WITH SOME OTHER FILE
#plotType =1    #TODO: THIS IS A FLAG TO DO THE CORRESPONDING TYPE OF PLOTTING
    
#Convert stereo signals to mono
if dataArray.shape[0] == 2:
    dataArray = (dataArray[...,0] + dataArray[...,1])/2


#  You can change these values to see the different effects.
grainInms = 120
pitch = 1.7   # 1=original pitch, <1 = lower in pitch, >1 = higher in pitch
overlapGrainWindow =.5 #25%

#convert to float and put between range of [-1, 1], max of range of int (-32768 to 32767)
range1Neg1 = dataArray /32767.0

grain_size = milliseconds2samples(grainInms, sampleRate)

pitchShiftedDataArray = GS_pshift(range1Neg1, pitch, milliseconds2samples(grainInms, sampleRate), overlapGrainWindow)

pitchShiftedDataArray = (pitchShiftedDataArray *2**15)

wavfile.write("PitchShiftedOutput.wav", sampleRate, pitchShiftedDataArray.astype(np.int16))
                        #End Of Main Code


                #Visualing and Analyzing the audio signals

#######################signwave plotting
if plotType==0 :
    #Input Signal
    plt.plot( dataArray[0:])
    plt.xlim(0,10000)
    # display the plot
    plt.show()

    #Output Signal
    wavFileName ="PitchShiftedOutput.wav"
    sampleRate, dataArray = wavfile.read(wavFileName)
    plt.plot( dataArray[0:])
    plt.xlim(0,10000)
    # display the plot
    plt.show()
    
else:
    spf = wave.open("//engin-labs.m.storage.umich.edu/allanahm/windat.v2/Documents/OSR_us_000_0010_8k.wav", "r")
    
    # Extract Raw Audio from Wav File
    signal = spf.readframes(-1)
    signal = np.frombuffer(signal, dtype='int16')
    fs = spf.getframerate()
    
    # If Stereo
    if spf.getnchannels() == 2:
        print("Just mono files")
        sys.exit(0)
    
    
    Time = np.linspace(0, len(signal) / fs, num=len(signal))
    
    
    # read audio samples
    input_data = wavfile.read(wavFileName)
    audio = input_data[1]
    
    
    # plot the first 1024 samples of original signal
    plt.plot(audio[0:])
    # label the axes
    plt.ylabel("Amplitude")
    plt.xlabel("Time")
    # set the title  
    plt.title("Audio Signal of Wav File") 
    # display the plot
    plt.show()
    
    
    wav_file = wave.open(wavFileName)
    totalFrames=wav_file.getframerate()
    data = wav_file.readframes(totalFrames)
    wav_file.close()
    data = struct.unpack('{n}h'.format(n=totalFrames), data)
    data = np.array(data)
    data_fft = np.fft.fft(data)
    # This will give us the frequency we want
    
    frequencies = np.abs(data_fft)
    print("The frequency is {} Hz".format(np.argmax(frequencies)))
    
    
    #Original Signal Analysis and Plots
    plt.plot(data[:])
    plt.title("Original Audio Wave")
    plt.ylabel("Amplitude")
    plt.show()
    
    plt.plot(frequencies)
    plt.title("Frequencies found")
    plt.xlim(0,3000)
    plt.ylabel("Quantity")
    plt.xlabel("Frequency")
    plt.show()
    
    
    #Output File Analysis and Plots (same methods as above)
    wav_file = wave.open("PitchShiftedOutput.wav", 'r')
    totalFrames=wav_file.getframerate()
    data = wav_file.readframes(totalFrames)
    wav_file.close()
    data = struct.unpack('{n}h'.format(n=totalFrames), data)
    data = np.array(data)
    data_fft = np.fft.fft(data)
    # This will give us the frequency we want
    
    frequencies = np.abs(data_fft)
    print("The frequency is {} Hz".format(np.argmax(frequencies)))
    
    plt.plot(data[:])
    plt.title("Output Audio Wave")
    plt.ylabel("Amplitude")
    plt.show()
    
    plt.plot(frequencies)
    plt.title("Output Frequencies Found")
    plt.ylabel("Quantity")
    plt.xlabel("Frequency")
    plt.xlim(0,3000)
    plt.show()
