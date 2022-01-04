#ifndef AudioPitchShift_h_
#define AudioPitchShift_h_

#include <cstdint>
#include <cmath>
#include <string>
#include "FftRealPair.hpp"


#define AUDIO_SHIFT_BLOCKS 16
#define AUDIO_BLOCK_SAMPLES 128
#define FFT_LENGTH 256  
#define FS 44100 

typedef struct audio_block_struct {
	uint8_t  ref_count;
	uint8_t  reserved1;
	uint16_t memory_pool_index;
	int16_t  data[AUDIO_BLOCK_SAMPLES];
} audio_block_t;

class AudioPitchShift {
public:
    AudioPitchShift(void) : enabled(false){}

    void begin(void);

    audio_block_t* update(audio_block_t*);

    void set_period(double);
    
    void set_new_period(double);
    
    void set_method(int);
    
    int numSamp = 0;

private:

    void process(void);

    void td_psola(void);
    
    //function which creates hanning window of size FFTLength
    void hanning(void); 
    
    //phasevocode function
    void phasevocode(void); 
    
    void print_audio_buff(void); 
    
    int find_closest_peak(uint16_t val);
    
    double phase(double &r, double &i);
    
    void resample(double ratio, int samp);
    
    double mag(double &r, double &i);
    
    double princargx(double val);

    void clear_peaks(void);

    inline void transmit_push(audio_block_t*);
    
    void transform(double &real, double &imag); 
    
    void inverseTransform(double &real, double &imag); 
    
    void transformRadix2(double &real, double &imag); 

    int16_t new_peak_offset = 0;
    int16_t orig_peak_offset = 0;
    int16_t last_width = -1;
    uint8_t state;
    bool first_run, enabled, method; 
    double period = 0;
    double new_period = 0;
    double hann[FFT_LENGTH];  //hanning window of size FFTLength 
    uint16_t peaks[32];	
    uint16_t new_peaks[32];	
    uint16_t num_peaks = 0;	
    uint16_t num_new_peaks = 0;
    uint16_t transmit_num = 0;
    uint16_t transmit_head = 0;
    uint16_t process_step = 0; //for debugging purposes, will get rid of later
    uint16_t analHop = 64; //used to index Audio buffer
    uint16_t synthHop = 0;
    double anlHop = 0; //needed for floating point calcs
    double Fs = 0; 
    double FFT_length = 0; 
    double phase_init = 0;
    double phase_next = 0;
    double synthPhase_init = 0;
    double phasediff = 0;
    double mag_next = 0;
    int numBins = 0; 
    int outputSize = 0;
    int16_t  AudioBuffer[AUDIO_SHIFT_BLOCKS*128] __attribute__ ( ( aligned ( 4 ) ) );
    int16_t  OutputBuffer[AUDIO_SHIFT_BLOCKS*128] __attribute__ ( ( aligned ( 4 ) ) );
    int16_t PHVCOutputBuffer[AUDIO_SHIFT_BLOCKS*128*4] __attribute__ ( ( aligned (4) ) ); 
    double realInit[FFT_LENGTH];  //Vector for the real part of FFT 
    double imagInit[FFT_LENGTH];  //Vector for imag part of FFT 
    double realNext[FFT_LENGTH];  
    double imagNext[FFT_LENGTH];  
    double synthRealInit[FFT_LENGTH]; 
    double synthImagInit[FFT_LENGTH]; 
    double synthRealNext[FFT_LENGTH];
    double synthImagNext[FFT_LENGTH]; 
    audio_block_t *transmit_queue[16];
    audio_block_t *blocklist1[AUDIO_SHIFT_BLOCKS];
    audio_block_t *inputQueueArray[1];
};
#endif
