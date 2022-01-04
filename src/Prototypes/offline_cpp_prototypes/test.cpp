#include "effect_pitch_shift.h"
#include "AudioFile.h"
#include <iostream> 
#include <fstream>
using namespace std;

int main(){
    
    string filename, outfile, alg; 
    double periodInit, periodNew;
    
    cout << "Hello user!  Please input the name of the .wav file you wish to be processed.  Make sure to include the file extension!" << endl; 
    cin >> filename; 
    cout << "Please type in the name of the output file you would like to save to. Make sure to include the file extension!" << endl; 
    cin >> outfile; 
    cout << "Please type in '1' if you wish to process the file with TDPSOLA, or type in '2' to process the file with Phase-Vocoding" << endl; 
    cin >> alg; 
    while ( (alg != "1") & (alg != "2")){
        cout << "Your input was incorrect! Please type in '1' if you wish to process the file with TDPSOLA, or type in '2' to process the file with Phase-Vocoding" << endl; 
        cin >> alg;
    }
    cout << "Please input the base frequency you would like to shift the input from"  << endl; 
    cin >> periodInit;
    cout << "Please input the target frequency you would like to shift to" << endl; 
    cin >> periodNew; 
    
    
    AudioFile<double> audioFile;
    AudioFile<double> audioFileNumeroDos; 
    AudioFile<double> audioFileNumeroTres; 
    AudioPitchShift pitch;
    pitch.begin();
    if(alg == "1"){ pitch.set_method(1); }
    else if (alg == "2"){ pitch.set_method(2); }
    pitch.set_period(44100/periodInit); 
    pitch.set_new_period(44100/periodNew); 
    
    

    audioFile.load(filename); 
    int channel = 0;
    int numSamples = audioFile.getNumSamplesPerChannel();
    pitch.numSamp = numSamples; 

    AudioFile<double>::AudioBuffer out_buff;
    AudioFile<double>::AudioBuffer out_buff_phvcd;
    AudioFile<double>::AudioBuffer out_buff_phvcd_resamp; 
    out_buff.resize(1);
    out_buff[0].resize(numSamples);
    int out_indx = 0;
    out_buff_phvcd.resize(1);
    out_buff_phvcd_resamp.resize(1); 
    
    audio_block_t* output_block;
    for(int i = 0; i < numSamples; i += 128){
        audio_block_t *block = new audio_block_t;
        for(int j = 0; j < 128; ++j){
            block->data[j] = int16_t(audioFile.samples[channel][i+j] * 0x10000);
        }
        output_block = pitch.update(block);
        if(output_block){
            for(int j = 0; j < 128; ++j){
                out_buff[0][out_indx++] = (output_block->data[j]) / (double)0x10000;
            }
        }
    }
    bool ok = audioFile.setAudioBuffer(out_buff);
    if(ok){
        audioFile.save(outfile);
    }
    cout << "Number of samples in female_scale.wav: " <<  numSamples << endl; 
}
