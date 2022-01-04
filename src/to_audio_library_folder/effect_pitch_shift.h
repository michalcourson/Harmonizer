#ifndef AudioPitchShift_h_
#define AudioPitchShift_h_

#include <cstdint>
#include "Audio.h"
#include "AudioStream.h"

#define AUDIO_SHIFT_BLOCKS 3
//#define AUDIO_BLOCK_SAMPLES 128
#define FFT_LENGTH 512  //Daniel added this
#define FSS 44100 //Daniel added this

class AudioPitchShift : public AudioStream {
public:
    AudioPitchShift(void) : AudioStream(1, inputQueueArray), enabled(false){}
    void begin(void);
    virtual void update(void);
    void set_period(double);
    void set_new_period(double);
    int numSamp = 0;
private:
    void process(void);
    void clear_peaks(void);
    int find_closest_peak(uint16_t val);
    void td_psola(void);
    void hanning(void); 
    void resample(double ratio, int samp);
    double princargx(double val);
    void phasevocode(void); 

    inline void transmit_push(audio_block_t*);
    inline void transform(float32_t*FFT_array, float32_t *time_array);
    inline void inverseTransform(float32_t*FFT_array, float32_t *time_array);
    
    int16_t new_peak_offset = 1024;
    int16_t orig_peak_offset = 1024;
    int16_t last_width = -1;
    int16_t algorithm = 0;
    bool first_run, enabled;
    bool fade_in, fade_out;
    double period = 0;
    double new_period = 0;
    double hann[FFT_LENGTH];  
    uint16_t peaks[32];
    uint16_t new_peaks[32];
    uint16_t num_peaks = 0;
    uint16_t num_new_peaks = 0;
    uint16_t transmit_num = 0;
    uint16_t transmit_head = 0;
    uint16_t anlHop = 512; 
    uint16_t synthHop; 
    double Fs = 0;  
    double phase_init = 0;
    double phase_next = 0;
    double synthPhase_init = 0;
    double phasediff = 0;
    double mag_next = 0;
    int16_t  AudioBuffer[3072] __attribute__ ( ( aligned ( 4 ) ) );
    int16_t  TD_OutputBuffer[3072] __attribute__ ( ( aligned ( 4 ) ) );
    int16_t  PV_OutputBuffer[3072] __attribute__ ( ( aligned ( 4 ) ) );
    int16_t TimeStretchedBuffer[2046*4] __attribute__ ( ( aligned (4) ) ); 
    float32_t Time_Init[FFT_LENGTH];
    float32_t Time_Next[FFT_LENGTH];
    float32_t Time_Synth[FFT_LENGTH];
    float32_t FFT_Init[FFT_LENGTH];
    float32_t FFT_Next[FFT_LENGTH];
    float32_t FFT_Synth_Next[FFT_LENGTH];
    float32_t FFT_Synth_Init[FFT_LENGTH];
    arm_rfft_fast_instance_f32 fft_instance;
    audio_block_t *transmit_queue[16];
    audio_block_t *blocklist1[AUDIO_SHIFT_BLOCKS];
    audio_block_t *inputQueueArray[1];
};
#endif
