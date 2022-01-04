#include "effect_pitch_shift.h"
#include "FftRealPair.hpp"
#include <iostream>
#include <vector>
#include <cmath>

//static audio
#define PI 3.14159265358

using namespace std;

static void copy_buffer(void *destination, const void *source) {
    const int16_t *src = ( const int16_t * )source;
    int16_t *dst = (  int16_t * )destination;
    for (int i=0; i < AUDIO_BLOCK_SAMPLES; i++)  *dst++ = (*src++);
}

void AudioPitchShift::clear_peaks(){
    for(int i = 0; i < 32; ++i){
        peaks[i] = 0;
        new_peaks[i] = 0;
    }
    num_peaks = 0;
    num_new_peaks = 0;
}

inline void AudioPitchShift::transmit_push(audio_block_t * block){
    transmit_queue[(transmit_head + transmit_num) & 15] = block;
    ++transmit_num;
}

audio_block_t* AudioPitchShift::update(audio_block_t* block){
    if (!block) return nullptr;
    
    if ( !enabled ) {
        return nullptr;
    }
    blocklist1[state++] = block;
    

    if ( state >= AUDIO_SHIFT_BLOCKS/2 ) {
        for ( int i = 0; i < AUDIO_SHIFT_BLOCKS/2; i++ ) copy_buffer( AudioBuffer+( i * 0x80 ), AudioBuffer+( (i+8) * 0x80 ));
        for ( int i = 0; i < AUDIO_SHIFT_BLOCKS/2; i++ ) copy_buffer( AudioBuffer+( (i+8) * 0x80 ), blocklist1[i]->data );
        if ( !first_run ) process( );
        first_run = false;
        state = 0;
    }
    if (transmit_num != 0){
        audio_block_t * transmit_block = transmit_queue[transmit_head];
        transmit_head = (transmit_head + 1) & 15;
        --transmit_num;
        return transmit_block;
    }else{
        return nullptr;
    }
}

void AudioPitchShift::process(void){
    if(method == 0){
        td_psola();
    }
    else if (method == 1){
        phasevocode();
    }
    
}

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

double AudioPitchShift::phase(double &r, double &i){
    if (r >= 0){
        if (i == 0){
            return 0;
        }
        else{
            return atan(i/r);
        }
    }
    else if (i < 0){
        return atan(i/r) - PI;
    }
    else{
        return atan(i/r) + PI; 
    }
    
}

void AudioPitchShift::resample(double ratio, int samp){
    int n_out = int(samp/ratio);
    double t = 0;
    double a = 0; 
    int n = 0;
    for(int q = 0; q < n_out ; q++){
        t = (double(q)) * ratio;
        n = int(t);
        a = 1.0 - (t - n); 
        OutputBuffer[q + 512] += int16_t((a * PHVCOutputBuffer[n]) + ((1.0 - a)*(PHVCOutputBuffer[n + 1]))); 
    }
    for(int q = 0; q < 8192; q++){
        PHVCOutputBuffer[q] = 0;
    }
}


double AudioPitchShift::mag(double &r, double &i){
    return sqrt(pow(r,2) + pow(i,2));
}

double AudioPitchShift::princargx(double val){
    double z = 0;
    double v = val + PI; 
    if ( v < 0 ){
        v = -1 * v;
        z = fmod(v, 2*PI); 
        z = 2*PI - z;
    } 
    else{    
        z = fmod(v, 2*PI);
    }
    return (z - PI); 
}

void AudioPitchShift::td_psola(){
    const uint16_t N = 2048;
    

    clear_peaks();
    int peak_num = 0;
    while(true){	
        if (orig_peak_offset + peak_num*period >= N) break;	
        uint16_t num = uint16_t(orig_peak_offset+ peak_num*period);	
        if (num_peaks >= 32) break;	
        peaks[num_peaks++] = num;	
        ++peak_num;	
    }	
    for(int i =0; i < num_peaks; ++i){	
        if(peaks[i] > 1536){	
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
        if(new_peaks[i] > 1536){	
            new_peak_offset = new_peaks[i] - 1024;	
            break;	
        }	
    }

    for (int i = 0; i < num_new_peaks; ++i){
        if (new_peaks[i] < 512) continue;
        
        int peak_indx = find_closest_peak(new_peaks[i]);

        int window_length;
        if(period > new_period){
            window_length = int(new_period);
        }else{
            window_length = int(period);
        }
        if (last_width != -1 && i == 0){
            window_length = last_width;
        }

        if(peaks[peak_indx] + window_length >= 2048){
            window_length = 2048 - peaks[peak_indx];
        }
        if(new_peaks[i] + window_length >= 2048){
            window_length = 2048 - new_peaks[i];
        }
        if(peaks[peak_indx] - window_length < 0){
            window_length = peaks[peak_indx];
        }
        if(new_peaks[i] - window_length < 0){
            window_length = new_peaks[i];
        }
        double multiplier = 0;
        double step = 1.0/(double)window_length;
        for (int j = -window_length; j < 0; j++) {
            OutputBuffer[new_peaks[i] + j] += int16_t(multiplier * AudioBuffer[peaks[peak_indx] + j]);
            multiplier += step;
        }
        for (int j = 0; j < window_length; j++) {
            OutputBuffer[new_peaks[i] + j] += int16_t(multiplier * AudioBuffer[peaks[peak_indx] + j]);
            multiplier -= step;
        }

        if(new_peaks[i] > 1536 && i !=0){
            last_width = window_length;
            for(int j = 0; j < 8; j++){
                audio_block_t * transmit_block = new audio_block_t;
                copy_buffer(transmit_block->data,  OutputBuffer+( (j+4) * 0x80 ));
                transmit_push(transmit_block);
            }
            for(int j = 0; j < 8; ++j){
                copy_buffer(OutputBuffer+( (j) * 0x80 ),  OutputBuffer+( (j+8) * 0x80 ));
            }
            for(int j = 1024; j < 2048; ++j){
                OutputBuffer[j] = 0;
            }
            return;
        }

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
    for(int q = 0; q < FFT_LENGTH; q++){ //printing to terminal to verify window created correctly 
        std::cout << "Hann element " << q << " : " << hann[q] << std::endl;
    }
    
}

/* 
When phasevocode() is first called (when process_step == 0), the hanning window, along with the 
analysis hopsize, synthesis hopsize, sampling rate and FFT length are initialized. Four arrays 
(realInit, imagInit, realNext and imagNext) are used for storing FFT results from consective 
windows of the input signal.  The FFT is taken on both consecutive windows using the Fft::transform()
function. 
The synthesized FFT results are stored in four arrays (synthRealInit, synthImagInit, synthRealNext, 
synthImagNext).  Within the synthesis loop (starts on line 281), phase and magnitude calculations 
are performed on the FFT arrays in order to populate the synthesized results in synthRealNext and
synthImagNext.  Results from these arrays are then copied into synthRealInit and synthImagInit before
the Fft::inverseTransform() function is called (this function modifies the arrays passed in by reference)
so the next iteration has the information needed to perform the necessary phase calculations. 
After the inverse transform has been taken on the synthesized results, overlap and add is performed 
on the output buffer using the synthesis hopsize (to define the overlap amount).
Resampling is then performed on the output buffer using the resample() function.  This step of the 
code is what actually pitch shifts the output and makes the time duration of the output equivalent to 
the time duration of the input. 
*/


void AudioPitchShift::phasevocode(){

    if(process_step == 0){
        hanning();
        anlHop = 64;
        Fs = 44100; 
        FFT_length = 256; 
        numBins = int(numSamp/anlHop) - int(FFT_length/anlHop) + 1; 
        synthHop = anlHop*(period/new_period); 
        outputSize = synthHop*numBins; 
    }
    
    for(int i = 0; i < (AUDIO_BLOCK_SAMPLES*8)/analHop; i++){
        for(int q = 0; q < FFT_LENGTH; ++q){
            realInit[q] = AudioBuffer[(i*analHop) + q]; //ADD WINDOW by adding '* hann[q]'
            imagInit[q] = 0;
            realNext[q] = AudioBuffer[((i+1)*analHop) + q]; //ADD WINDOW by adding '* hann[q]'
            imagNext[q] = 0; 
        }

        Fft::transform(realInit, imagInit);
        Fft::transform(realNext, imagNext);

        if ((i == 0) & (process_step == 0)){
            for (int f = 0; f < FFT_LENGTH; ++f){
                synthRealInit[f] = realInit[f];
                synthImagInit[f] = imagInit[f];
            }
        }
        for(int j = 0; j < FFT_length; j++){
            if(j <= FFT_length/2){
                phase_init = phase(realInit[j], imagInit[j]);
                phase_next = phase(realNext[j], imagNext[j]);
                synthPhase_init = phase(synthRealInit[j], synthImagInit[j]);
                phasediff = ((anlHop/Fs)*((Fs*j)/FFT_length)*2*PI) + princargx(phase_next - phase_init - ((anlHop/Fs)*((Fs*j)/FFT_length)*2*PI));
                mag_next = mag(realNext[j], imagNext[j]);
                synthRealNext[j] = mag_next * cos(synthPhase_init + phasediff*(period/new_period));
                synthImagNext[j] = mag_next * sin(synthPhase_init + phasediff*(period/new_period)); 
            }
            else if(j > FFT_length/2){
                synthRealNext[j] = synthRealNext[(FFT_LENGTH/2) - (j - ((FFT_LENGTH)/2))];
                synthImagNext[j] = -1*(synthImagNext[(FFT_LENGTH/2) - (j - ((FFT_LENGTH)/2))]);
            }
        }

        if( (i == 0) & (process_step == 0)){
            Fft::inverseTransform(synthRealInit, synthImagInit); 
        }
        for (int j = 0; j < FFT_length; j++){
            synthRealInit[j] = synthRealNext[j];
            synthImagInit[j] = synthImagNext[j];
        }

        Fft::inverseTransform(synthRealNext, synthImagNext);

        for (int j = 0; j < FFT_LENGTH; j++){
            PHVCOutputBuffer[synthHop*i + j] = PHVCOutputBuffer[synthHop*i + j] + (synthRealNext[j]);
        }
        
        for (int j = 0; j < FFT_LENGTH; j++){
            synthRealNext[j] = 0;
            synthImagNext[j] = 0; 
            realInit[j] = 0;
            imagInit[j] = 0;
            realNext[j] = 0;
            imagNext[j] = 0;
        }
        
    }
    
    resample((period/new_period), synthHop*15 + 256);
    for(int j = 0; j < 8; j++){
        audio_block_t * transmit_block = new audio_block_t;
        copy_buffer(transmit_block->data,  OutputBuffer+( (j+4) * 0x80 ));
        transmit_push(transmit_block);
    }
    for(int j = 0; j < 8; ++j){
        copy_buffer(OutputBuffer+( (j) * 0x80 ),  OutputBuffer+( (j+8) * 0x80 ));
    }
    for(int j = 1024; j < 2048; ++j){
        OutputBuffer[j] = 0;
    }
    
    process_step++;

} 

void AudioPitchShift::print_audio_buff(){
    std::cout << "size of AudioBuffer:  " << sizeof(AudioBuffer) << std::endl; 
}

void AudioPitchShift::begin( ) {
    first_run      = true;
    enabled        = true;
    state          = 0;
}


void AudioPitchShift::set_period(double p){
    period = p;
}

void AudioPitchShift::set_new_period(double p){
    new_period = p;
}

void AudioPitchShift::set_method(int mthd){
    if (mthd == 1){
        method = 0;
    }
    else if (mthd == 2){
        method = 1;
    }
}

void Fft::transform(double *real, double *imag) {
    transformRadix2(real, imag); 
}

void Fft::inverseTransform(double *real, double *imag) {
	transform(imag, real);
    for (int i = 0; i < FFT_LENGTH; i++){
        real[i] = (real[i])/FFT_LENGTH;
        imag[i] = (imag[i])/FFT_LENGTH; 
    }
}

static size_t reverseBits(size_t val, int width) {
	size_t result = 0;
	for (int i = 0; i < width; i++, val >>= 1)
		result = (result << 1) | (val & 1U);
	return result;
}

void Fft::transformRadix2(double *real, double *imag) {
	// Length variables
    
	size_t n = FFT_LENGTH;
	int levels = 0;  // Compute levels = floor(log2(n))
	for (size_t temp = n; temp > 1U; temp >>= 1)
		levels++;
	if (static_cast<size_t>(1U) << levels != n)
		throw std::domain_error("Length is not a power of 2");
	
	// Trigonometric tables
    double cosTable[n/2];
    double sinTable[n/2];
	for (size_t i = 0; i < n / 2; i++) {
		cosTable[i] = std::cos(2 * M_PI * i / n);
		sinTable[i] = std::sin(2 * M_PI * i / n);
	}
	
	// Bit-reversed addressing permutation
	for (size_t i = 0; i < n; i++) {
		size_t j = reverseBits(i, levels);
		if (j > i) {
			std::swap(real[i], real[j]);
			std::swap(imag[i], imag[j]);
		}
	}
	
	// Cooley-Tukey decimation-in-time radix-2 FFT
	for (size_t size = 2; size <= n; size *= 2) {
		size_t halfsize = size / 2;
		size_t tablestep = n / size;
		for (size_t i = 0; i < n; i += size) {
			for (size_t j = i, k = 0; j < i + halfsize; j++, k += tablestep) {
				size_t l = j + halfsize;
				double tpre =  real[l] * cosTable[k] + imag[l] * sinTable[k];
				double tpim = -real[l] * sinTable[k] + imag[l] * cosTable[k];
				real[l] = real[j] - tpre;
				imag[l] = imag[j] - tpim;
				real[j] += tpre;
				imag[j] += tpim;
			}
		}
		if (size == n)  // Prevent overflow in 'size *= 2'
			break;
	}
}
