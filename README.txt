Team members:
	Audrey Cooke
	audeo@umich.edu
 
	Michal Courson 
	mcourson@umich.edu
 
	Daniel Kohler
	danikohl@umich.edu
 
	Allanah Matthews
	allanahm@umich.edu
 
Zipfile Contents:
COPYRIGHT.TXT
	--Document which highlights the necessary copyright provisions for the Harmonizer project

2021-Team-Harmonizer-final.pdf
	--Final report for the Harmonizer project.  

2021-Team-Harmonizer-slides.pdf
	--Final presentation slide-deck to be presented in class. 
 
2021-Team-Harmonizer-poster.pdf
	--Final project poster to be presented at the Fall 2021 CoE Design Expo

vidoes 
	--Folder which includes real time demonstration videos of the Harmonizer for our final presentation slides

src 
	--Folder which contains the source code for our project
		--Includes code written for the offline python and offline C++ prototypes of our algorithms.  
		Each offline folder includes its own USAGE.TXT to outline how to use the code.

		--The OVERVIEW.TXT file provides an overview of each embedded code file needed for the Teensy.
		It also gives a basic overview of the software operation of the system. 

		--The USAGE.TXT file contains user-interface information which details the designated inputs 
		and outputs for the Harmonizer on the Teensy itself. 

		--The INSTALL.TXT file contains instructions for software installation on the Teensy.  It also 
		details the designated folder locations for specific files in order for the Harmonizer to successfully run in real time. 

		--The to_audio_library_folder houses the files which need to be stored in the Audio folder for the Teensy.

		--The to_core_teensy4_folder houses the files which need to be stored in the CoreTeensy4 folder for the Teensy. 

		--The harmonizer folder includes the .ino file which needs to be compiled on the Teensy in order for the system to operate in real time. 
data 
	--The audio subfolder includes results of pitch shifting a tritonal signal (at a note frequency of C3) using the TD-PSOLA algorithm, 
	the Phase-Vocoding algorithm, and both of these algorithms together with crossfading implemented at the switching threshold.  
	The results from each of these tests can be compared to Target_Output.wav

	--The figures subfolder includes images and plots which are used in the final presentation and final report. 

design 
	--This folder houses subfolders which are designated for files associated with the PCB design, the package design and images of the assembly process.
  
	--The pdf for bill of materials is also included in this folder (not included in its own subfolder). 
