/*
    This file contains implementations of our functions to process the input audio
    by pitch shifting it. The amount of shifting is dictated by the variables "period"
    and "new_period" corresponding to the input pitch and the target pitch. In it's current 
    form, this object will process in the signal using an algorithm called TD-PSOLA 
    unless the target pitch is around an octave or below the input pitch, where it will be
    processed by an algorithm called phase vocoding. There is a crossfade region where
    the outputs of both algorithm are mixed together to make the transition from one 
    algorithm to the other smooth. 
*/
#include <Arduino.h>
#include "effect_pitch_shift.h"
#include <cmath>
#include <complex>

using namespace std;



static void copy_buffer(void *destination, const void *source) {
    const uint16_t *src = ( const uint16_t * )source;
    uint16_t *dst = (  uint16_t * )destination;
    for (int i=0; i < AUDIO_BLOCK_SAMPLES; i++)  *dst++ = (*src++);
}

inline void AudioPitchShift::transmit_push(audio_block_t * block){
    transmit_queue[(transmit_head + transmit_num) & 15] = block;
    ++transmit_num;
}

void AudioPitchShift::update(void){

    audio_block_t *block;
    block = receiveReadOnly();
    if (!block) return;
    if ( !enabled ) {
        release( block );
        return;
    }

    //add new samples to the rolling buffer of input audio
    copy_buffer( AudioBuffer, AudioBuffer+AUDIO_BLOCK_SAMPLES);
    copy_buffer( AudioBuffer+AUDIO_BLOCK_SAMPLES, AudioBuffer+(2*AUDIO_BLOCK_SAMPLES));
    copy_buffer( AudioBuffer+(2*AUDIO_BLOCK_SAMPLES), block->data);
    release(block);
    
    //only process signal if there is a target pitch
    if (new_period != 0.0) process( );


    if (transmit_num != 0){
        audio_block_t * transmit_block = transmit_queue[transmit_head];
        transmit_head = (transmit_head + 1) & 15;
        --transmit_num;
        //fade in and out if this object is starting or stopping playback - reduces poping
        if(fade_in){
            double gain = 0;
            for(int i = 0; i < 128; ++i){
                transmit_block->data[i] *= gain;
                gain += 1/128.0;
            }
            fade_in = false;
        }
        if(fade_out){
            double gain = 1;
            for(int i = 0; i < 128; ++i){
                transmit_block->data[i+896] *= gain;
                gain -= 1/128.0;
            }
            fade_out = false;
        }
        transmit(transmit_block);
        release(transmit_block);
    }else{
        //if there is no processed audio, output a block of 0's
        audio_block_t * transmit_block = allocate();
        for(int i = 0; i < AUDIO_BLOCK_SAMPLES; ++i){
            transmit_block->data[i] = 0;
        }
        transmit(transmit_block);
        release(transmit_block);
    }
}

void AudioPitchShift::process(void){

    audio_block_t *transmit_block = allocate();
    double ratio = period/new_period;
    //since it is fast, run td_psola every time
    td_psola();
    if (ratio < .7){
        //run phase vocoding if needed, and phase vocoding and td_psola as needed
        phasevocode();
        double td_gain = 0; //TODO WAS 0

        if(ratio > .5){
            td_gain = (ratio - .5)*5;
        }

        for(int i = 0; i < AUDIO_BLOCK_SAMPLES; ++i){
            transmit_block->data[i] = td_gain * TD_OutputBuffer[i+AUDIO_BLOCK_SAMPLES] + (1-td_gain)*PV_OutputBuffer[i];
        }
    }else{
        copy_buffer(transmit_block->data,  TD_OutputBuffer+(AUDIO_BLOCK_SAMPLES));
    }
    transmit_push(transmit_block);

    //rotate rolling output buffers
    copy_buffer(TD_OutputBuffer+AUDIO_BLOCK_SAMPLES,  TD_OutputBuffer+(2*AUDIO_BLOCK_SAMPLES));
    for(int j = 2048; j < 3072; ++j){
        TD_OutputBuffer[j] = 0;
    }
    copy_buffer(PV_OutputBuffer,  PV_OutputBuffer+AUDIO_BLOCK_SAMPLES);
    for(int j = 1024; j < 2048; ++j){
        PV_OutputBuffer[j] = 0;
    }
}

void AudioPitchShift::clear_peaks(){
    for(int i = 0; i < 32; ++i){
        peaks[i] = 0;
        new_peaks[i] = 0;
    }
    num_peaks = 0;
    num_new_peaks = 0;
}

//find the closest input peak to a given output peak
int AudioPitchShift::find_closest_peak(uint16_t val){
    int upper_indx = 0;
    for(int i = 0; i < 16; ++i){
        if(peaks[upper_indx] > val) break;
        ++upper_indx;
    }
    if(upper_indx == 0){
        return 0;
    }
    int diff1 = peaks[upper_indx] - val;
    int diff2 = val - peaks[upper_indx - 1];
    if (diff1 < 0) diff1 *= -1;
    if (diff2 < 0) diff2 *= -1;
    if (diff1 < diff2){
        return upper_indx;
    }else{
        return upper_indx - 1;
    }
}

void AudioPitchShift::td_psola(){
    const uint16_t N = 3072;
    clear_peaks();
    
    //set pitch marking for both the input and the output
    int peak_num = 0;
    while(true){
        if (orig_peak_offset + peak_num*period >= N) break;
        uint16_t num = uint16_t(orig_peak_offset+ peak_num*period);
        if (num_peaks >= 32) break;
        peaks[num_peaks++] = num;
        ++peak_num;
    }
    for(int i =0; i < num_peaks; ++i){
        if(peaks[i] > 2048){
            orig_peak_offset = peaks[i] - 1024;
            break;
        }
    }

    peak_num = 1;
    while(true){
        if(new_peak_offset + peak_num*new_period >= N) break;
        uint16_t num = uint16_t(new_peak_offset + peak_num*new_period);
        if(num_new_peaks >= 32) break;
        new_peaks[num_new_peaks++] = num;
        ++peak_num;
    }
    for(int i = 0; i < num_new_peaks; ++i){
        if(new_peaks[i] > 2048){
            new_peak_offset = new_peaks[i] - 1024;
            break;
        }
    }
    
    for (int i = 0; i < num_new_peaks; ++i){
        if (new_peaks[i] < 1024) continue;
        
        int peak_indx = find_closest_peak(new_peaks[i]);

        //set an appropriate window lenght
        int window_length;
        if(period > new_period){
            window_length = int(new_period);
        }else{
            window_length = int(period);
        }
        if (last_width != -1 && i == 0){
            window_length = last_width;
        }

        //check for possible array out of bounds, adjust window lenght as needed
        if(peaks[peak_indx] + window_length >= 3072){
            window_length = 3072 - peaks[peak_indx];
        }
        if(new_peaks[i] + window_length >= 3072){
            window_length = 3072 - new_peaks[i];
        }
        if(peaks[peak_indx] - window_length < 0){
            window_length = peaks[peak_indx];
        }
        if(new_peaks[i] - window_length < 0){
            window_length = new_peaks[i];
        }

        //use a linear window centered on the input pitch marking, and copy it over to the output pitch marking
        double multiplier = 0;
        double step = 1.0/(double)window_length;
        for (int j = -window_length; j < 0; j++) {
            TD_OutputBuffer[new_peaks[i] + j] += int16_t(multiplier * AudioBuffer[peaks[peak_indx] + j]);
            multiplier += step;
        }
        for (int j = 0; j < window_length; j++) {
            TD_OutputBuffer[new_peaks[i] + j] += int16_t(multiplier * AudioBuffer[peaks[peak_indx] + j]);
            multiplier -= step;
        }

        //if there is a full block of output, return
        if(new_peaks[i] > 2048 && i !=0){
            last_width = window_length;
            return;
        }

    }
    
}

//wrapper for fft function
inline void AudioPitchShift::transform(float32_t*FFT_array, float32_t *time_array){
    arm_rfft_fast_f32(&fft_instance, time_array, FFT_array, 0);
}

//wrapper for ifft function
inline void AudioPitchShift::inverseTransform(float32_t*FFT_array, float32_t *time_array){
    arm_rfft_fast_f32(&fft_instance, FFT_array, time_array, 1);
}

//resamples the timestretched buffer to fit the block size
void AudioPitchShift::resample(double ratio, int samp){
    int n_out = int(samp/ratio);
    double t = 0;
    double a = 0; 
    int n = 0;
    for(int q = 0; q < n_out ; q++){
        t = (double(q)) * ratio;
        n = int(t);
        a = 1.0 - (t - n); 
        PV_OutputBuffer[q] += int16_t((a * TimeStretchedBuffer[n]) + ((1.0 - a)*(TimeStretchedBuffer[n + 1]))); 
    }
    for(int q = 0; q < 8184; q++){
        TimeStretchedBuffer[q] = 0;
    }
}

//creates hanning window of length FFT_LENGTH
void AudioPitchShift::hanning(){
    for(double i = 0; i < (FFT_LENGTH/2); i++){
        hann[int(i)] = 1 - cos(PI*(i/FFT_LENGTH));
    }
    for(double j = -(FFT_LENGTH/2); j < 0; j++){
        hann[(FFT_LENGTH/2) + (int(j) + (FFT_LENGTH/2))] = 1 - cos(PI*((j + 1)/FFT_LENGTH));
    }
}

//computes a phase wrapped between -pi and pi
double AudioPitchShift::princargx(double val){
    double z = 0;
    double v = val + PI; 
    double two_pi = 2*PI;
    z = v - two_pi * floor(v/two_pi);
    return (z - PI); 
}

void AudioPitchShift::phasevocode(){
    double f_ratio = period/new_period;

    //compute initialization values
    if(first_run){
        anlHop = 128;
        Fs = 44100; 
        synthHop = anlHop*f_ratio; 
        arm_rfft_fast_init_f32(&fft_instance, FFT_LENGTH);
    }
    
    int num_iters = 1024/anlHop;
    //synthesis hopsize is dictated by the ratio of target pitch to input pitch
    synthHop = int(anlHop*f_ratio);
    for(int i = 0; i < num_iters; i++){

        //copy time domain samples from the audio buffer
        for(int q = 0; q < FFT_LENGTH; ++q){
            if(i==0 && first_run){
                Time_Init[q] = AudioBuffer[((i)*anlHop) + q + 1024]* hann[q];
            }
            Time_Next[q] = AudioBuffer[((i+1)*anlHop) + q + 1024]* hann[q];
        }

        //transform time domain samples to frequency domain
        if(i==0 && first_run){
            transform(FFT_Init, Time_Init);
            for (int f = 0; f < FFT_LENGTH; ++f){
                FFT_Synth_Init[f] = FFT_Init[f];
            }
        }
        transform(FFT_Next, Time_Next);

        double h = anlHop/(double)Fs; 

        //calculate phase rotation for each frequency bin
        for(int j = 0; j < FFT_LENGTH/2; j++){
            complex<float32_t> init(FFT_Init[2*j], FFT_Init[2*j + 1]);
            complex<float32_t> next(FFT_Next[2*j], FFT_Next[2*j + 1]);
            complex<float32_t> synth(FFT_Synth_Init[2*j], FFT_Synth_Init[2*j + 1]);
            phase_init = arg(init);
            phase_next = arg(next);
            synthPhase_init = arg(synth);
            double fk = (Fs*(j))/(double)FFT_LENGTH; 
            double prin = princargx(phase_next - phase_init - (h*fk*2*PI));
            phasediff = (h*fk*2*PI) + prin;

            mag_next = abs(next);

            FFT_Synth_Next[2*j] = mag_next * cos(synthPhase_init + phasediff*f_ratio);
            FFT_Synth_Next[2*j + 1] = mag_next * sin(synthPhase_init + phasediff*f_ratio);
        }

        for (int k = 0; k < FFT_LENGTH; k++){
            FFT_Synth_Init[k] = FFT_Synth_Next[k];
        }

        //tranform frequency domain synthesized output to time domain
        inverseTransform(FFT_Synth_Next, Time_Synth); 

        //copy time domain samples to time stretched buffer
        for (int j = 0; j < FFT_LENGTH; j++){
            TimeStretchedBuffer[synthHop*i + j] += int16_t(Time_Synth[j]*hann[j]);
        }
        
        for (int j = 0; j < FFT_LENGTH; j++){
            FFT_Synth_Next[j] = 0;
            FFT_Init[j] = FFT_Next[j];
        }
    }
    
    //resample time stretched buffer to final output buffer
    int samps = synthHop*(num_iters - 1) + FFT_LENGTH;
    resample(f_ratio, samps);
    first_run = false;
} 

void AudioPitchShift::begin( ) {
    __disable_irq( );
    first_run      = true;
    enabled        = true;
    hanning();
    __enable_irq( );
}

//sets input period
void AudioPitchShift::set_period(double p){
    __disable_irq( );
    period = p;
    __enable_irq( );
}

//sets output period as well as fade_in and fade_out flags
void AudioPitchShift::set_new_period(double p){
    __disable_irq( );
    if(new_period == 0 && p != 0){
        fade_in = true;
        fade_out = false;
        first_run = true;
    }
    if(p == 0){
        fade_out = true;
    }
    new_period = p;
    __enable_irq( );
}