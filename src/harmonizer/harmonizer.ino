/*
	This is the driver Teensy sketch for our harmonizer. In this file, 
	four pitch shifting objects are created and mixed together with
	an unprocessed version of the input signal. This sketch also handles
	MIDI messages from an external host, as well as handling messages 
	where the teensy is acting as it's own host. The MIDI notes are then
	delegated to the pitch shifting objects by means of a queue structure.
*/
#include <USBHost_t36.h>
#include <Audio.h>
#include <Wire.h>
#include <SD.h>
#include <SPI.h>
#include <SerialFlash.h>


#define NOTE_LIMIT 4
  
const int myInput = AUDIO_INPUT_LINEIN;

struct MidiNote {
	byte channel, note, velocity;
};

bool operator==(const MidiNote& LHS, const MidiNote RHS) {
	return LHS.channel == RHS.channel && LHS.note == RHS.note;
}

//queue structure that handles note delegation between the pitch shifting objects
class NoteQueue {
private:
	static const size_t max_elts = NOTE_LIMIT;
	size_t num_elts = 0;
	size_t start_location = 0;

	MidiNote note_array[max_elts];

public:
	NoteQueue() {
		for (size_t i = 0; i < max_elts; ++i) note_array[i] = { 0,0,0 };
	}

	//Pop the first item from the note queue
	void pop_front() {
		if (!empty()) {
			at_helper(0) = { 0,0,0 };
			start_location = (start_location + 1) % max_elts;
			--num_elts;
		}
	}

	void push_back(const MidiNote& note) {
		if (!contains(note)) {
			if (num_elts == max_elts) pop_front();
			at_helper(num_elts) = note;
			++num_elts;
		}
	}

	const MidiNote& at(const size_t index) {
		return at_helper(index);
	}

	bool remove(const MidiNote& note) {
		//Find the element. If it isn't in the queue, return false.
		bool found = false;
		size_t found_index = 0;

		//Go through and find. If found, shift everything left 1
		for (size_t i = 0; i < size(); ++i) {
			if (at_helper(i) == note) {
				found = true;
				num_elts--;
			}
			if (found) at_helper(i) = at_helper(i + 1);
		}

		//Zero out the last element just because
		if (found) at_helper(size()) = { 0,0,0 };

		return found;
	}

	MidiNote& operator[](const size_t index) {
		return at_helper(index);
	}

	size_t size() const {
		return num_elts;
	}

	bool empty() {
		return num_elts == 0;
	}

	bool contains(const MidiNote& note) {
		for (size_t i = 0; i < size(); ++i) {
			if (at(i) == note)  return true;
		}
		return false;
	}

	void clear() {
		while (!empty()) pop_front();
	}

private:
	MidiNote& at_helper(const size_t index) {
		return note_array[(start_location + index) % max_elts];
	}

};

USBHost myusb;
USBHub hub1(myusb);
MIDIDevice midi1(myusb);
NoteQueue MIDINotes;

AudioInputI2S         audioInput;         // audio shield: mic or line-in
AudioOutputI2S        audioOutput;        // audio shield: headphones & line-out
AudioAnalyzeNoteFrequency analysis;

AudioPitchShift pitch1;
AudioPitchShift pitch2;
AudioPitchShift pitch3;
AudioPitchShift pitch4;

AudioMixer4 pitchmixer;
AudioMixer4 mastermixer;

//send input to a pitch estimation object
AudioConnection c1(audioInput, 0, analysis, 0);

//send input to four different pitch shifting objects
AudioConnection c3(audioInput, 0, pitch1, 0);
AudioConnection d3(audioInput, 0, pitch2, 0);
AudioConnection d4(audioInput, 0, pitch3, 0);
AudioConnection d5(audioInput, 0, pitch4, 0);

//route the pitch shifters to a mixer
AudioConnection c4(pitch1, 0, pitchmixer, 0);
AudioConnection c5(pitch2, 0, pitchmixer, 1);
AudioConnection c6(pitch3, 0, pitchmixer, 2);
AudioConnection c7(pitch4, 0, pitchmixer, 3);

//route the pitch mixer and original input into a second mixer
AudioConnection c8(audioInput, 0, mastermixer, 0);
AudioConnection c9(pitchmixer, 0, mastermixer, 1);

//send the master mixer to the output
AudioConnection c10(mastermixer, 0, audioOutput, 0);
AudioConnection c11(mastermixer, 0, audioOutput, 1);

AudioControlSGTL5000 audioShield;

void setup() {
  Serial.begin(115200);
  delay(300);
  pinMode(LED_BUILTIN, OUTPUT);

  // allocate memory for the audio library
  AudioMemory(128);
  audioShield.enable();
  audioShield.volume(0.8 );
  
  //initialize mixer gains
  pitchmixer.gain(0,.75);
  pitchmixer.gain(1,.75);
  pitchmixer.gain(2,.75);
  pitchmixer.gain(3,.75);

  mastermixer.gain(0,.5);
  mastermixer.gain(1, 1);
  mastermixer.gain(2, 0);
  mastermixer.gain(3, 0);

  myusb.begin();
  usbMIDI.begin();

  //Host Mode MIDI handlers
  midi1.setHandleNoteOff(OnNoteOff);
  midi1.setHandleNoteOn(OnNoteOn);
  midi1.setHandleControlChange(OnControlChange);

  //Device Mode MIDI handlers
  usbMIDI.setHandleNoteOff(OnNoteOff);
  usbMIDI.setHandleNoteOn(OnNoteOn);
  usbMIDI.setHandleControlChange(OnControlChange);

  // Initialize the processing objects
  analysis.begin(.15);
  pitch1.begin();
  pitch2.begin();
  pitch3.begin();
  pitch4.begin();
  pitch1.set_new_period(0);
  pitch2.set_new_period(0);
  pitch3.set_new_period(0);
  pitch4.set_new_period(0);
  pitch1.set_period(300);
  pitch2.set_period(300);
  pitch3.set_period(300);
  pitch4.set_period(300);
}

enum class MIDItype { NONE, HOST_MODE, DEVICE_MODE };
MIDItype miditype = MIDItype::NONE;

void loop()
{	
	//set input pitch estimations when they are available
	if(analysis.available()){
		float freq = analysis.read();
    Serial.println(freq);
		float period = 44100/freq;
		pitch1.set_period(period);
		pitch2.set_period(period);
		pitch3.set_period(period);
		pitch4.set_period(period);
	}

	//Host Mode: if MIDI keyboard is plugged in, use midi1
	if (midi1) {
		//Did we switch to this mode? If so clear the queue, turn on the LED for confirmation
		if (miditype != MIDItype::HOST_MODE) {
			miditype = MIDItype::HOST_MODE;
			MIDINotes.clear();
			digitalWrite(LED_BUILTIN, HIGH);
			Serial.println("Booting into Host Mode");
		}

		myusb.Task();
		midi1.read();
	}

	//Device Mode: if no MIDI keyboard plugged in, use usbMIDI
	else {
		//Did we switch to this mode? If so clear the queue, turn off the LED for confirmation
		if (miditype != MIDItype::DEVICE_MODE) {
			miditype = MIDItype::DEVICE_MODE;
			MIDINotes.clear();
			digitalWrite(LED_BUILTIN, LOW);
			Serial.println("Booting into Device Mode");
		}

		usbMIDI.read();
	}

}

void OnNoteOn(byte channel, byte note, byte velocity)
{
	MidiNote midi_note = { channel, note, velocity };
	MIDINotes.push_back(midi_note);
	updateTargets();
}

void OnNoteOff(byte channel, byte note, byte velocity)
{
	MidiNote midi_note = { channel, note, velocity };
	MIDINotes.remove(midi_note);
	updateTargets();
}

void OnControlChange(byte channel, byte control, byte value)
{
	//TODO: How can we use control changes to mess with params?
	Serial.print("Control Change, ch=");
	Serial.print(channel);
	Serial.print(", control=");
	Serial.print(control);
	Serial.print(", value=");
	Serial.print(value);
	Serial.println();
}

void updateTargets() {
	//set pitch targets based on MIDI queue
	pitch1.set_new_period(noteToPeriod(MIDINotes[0].note));
	pitch2.set_new_period(noteToPeriod(MIDINotes[1].note));
	pitch3.set_new_period(noteToPeriod(MIDINotes[2].note));
	pitch4.set_new_period(noteToPeriod(MIDINotes[3].note));

	Serial.println();
}


float noteToPeriod(int note) {
	if (note == 0) return 0.0;

	float a = 440; //frequency of A
	float f =  (a / 32) * pow(2, ((note - 9) / 12.0));

	return 44100. / f;
}
