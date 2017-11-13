November 12, 2017

This repository contains two folders: (1) DER_C++ (source code for DER) and (2) DER_Book_Examples (utility files to generate figures in MATLAB)

(1) DER_C++:
The source code for DER can be found here: http://www.cs.columbia.edu/cg/elastic_coiling/
We used the above code to simulate a few simple problems. If you want to use DER, please go to the above URL and cite the appropriate paper (M. K. Jawed et al, Proc. Natl. Acad. Sci. 2014).

Our examples are intended to accompany the following book: 
"A Primer on the Kinematics of Discrete Elastic Rods" by M. K. Jawed, A. Novelia, and O. O'Reilly

This code works only with Linux or OS X. It does not work on Windows (you may try to use cygwin or something similar).
Unfortunately, this code-base is not actively maintained. It requires the following libraries. We have noted the version number of these libraries that are known to work.
* Boost 1.47 http://www.boost.org/users/history/version_1_47_0.html
* TCLAP 1.2.1 http://sourceforge.net/projects/tclap/files/tclap-1.2.1.tar.gz/download
* LIBPNG 1.0.60 http://sourceforge.net/projects/libpng/files/libpng10/older-releases/1.0.60/
* ZLIB 1.2.7 http://sourceforge.net/projects/libpng/files/zlib/1.2.7/
* Compiler: gcc-4.4 and g++-4.4

You may choose to edit your "bashrc" (on Linux) or "bash_profile" (on OS X) file to add the following lines; make sure to change these lines to reflect to file paths on your computer:
export CC=/usr/bin/gcc-4.4
export CXX=/usr/bin/g++-4.4
export CPP=/usr/bin/g++-4.4
export TCLAP_INC_DIR=/usr/local/include/tclap

How to compile and run?
Create a directory called "build" (or give your favorite name) under DER_Book_Codes/DER_C++.
"cd" to that directory in the terminal or command line. Then follow the following steps:

(a) Run "$ cmake .."
(b) Run "$ ccmake .." and change the necessary values, e.g. set EIGEN_INCLUDE_PATH to "../BASim/src/Eigen2".
Use "CMAKE_BUILD_TYPE" as "Release" (instead of "Debug") for faster performance.
Usually, you would set "USE_LAPACK" to "ON" and "USE_MKL" to "OFF".
Then, as mentioned at the bottom of the "ccmake" window, hit "c" to configure (may need to do it twice) and "g" to generate & exit.
You may ignore warnings.
(c) Run "$ make" and wait a minute or two.

If you get errors (ignore the warnings), please make sure that all the libraries are installed and linked properly.
Unfortunately, we are unable to provide assistance with installation.

A primer on code-base: the "Problems" are written in "DER_C++/Apps/BASimulator/Problems/".
Each problem requires one ".hh" file and one ".cc" file.
The functions associated with computing the elastic/inertial/gravity forces are located in BASim/src/Physics/ElasticRods/

Run this command to get a list of problems in this code:
$ Apps/BASimulator/BASimulator -p

Run this command to get help:
$ Apps/BASimulator/BASimulator -h

To run the specific problems in the book, please run commands provided near the end of this file.
Let us go over some basics before running the problems:
* Note that a simulation window will open once you execute these commands.
* A text file will be created by the simulator in the "build" folder that contains the necessary data
* Hit "Spacebar" to start the simulation. Hit "Spacebar" again if you want to pause.
* Hit "1" to turn on/off the material frame, hit "2" to turn on/off the reference frame. By default, material frame may be ON.
* Use your mouse and "Control/Command" keys to change perspective. Especially, press the left mouse button and move around the cursor in the simulator.
* Hit "q" to quit.
* Right click on the simulator to see more options.
* Each problem requires an "options" file and is specified by the "-f" flag (as shown in the commands below).
* Sample options files are provided in the "assets" folder. Feel free to play around by changing Young's Modulus, length, number of vertices, etc.
* Set "end-time" to "-1" in the options file if you do not want the simulator to exit on its own.

(i) To run the cantilever beam example in Chapter 1:
$ Apps/BASimulator/BASimulator -r 1 -f ../assets/ch1_cantilever.txt &
// This will create "ch1_helixTerminalMomentData.txt" file

(ii) To run the rod under terminal moment example in Chapter 1:
$ Apps/BASimulator/BASimulator -r 2 -f ../assets/ch1_helixTerminalMomentData.txt &
// This will create "ch1_helixTerminalMomentData.txt" file

(iii) To run the helix uncoiling example in Chapter 5:
$ Apps/BASimulator/BASimulator -r 3 -f ../assets/ch5_helixUncoiling.txt &
// This will create "ch5_helixUncoilingData.txt" file

(iv) To run the vibrating rod example in Chapter 8:
$ Apps/BASimulator/BASimulator -r 4 -f ../assets/ch8_vibratingCable.txt &
// This will create "ch8_vibratingCableData.txt" and "ch8_vibratingCableEnergyData.txt" files

(2) DER_Book_Examples:
This folder contains a few subfolders containing each example in the book. Each subfolder contains a 00README_First.txt with instructions.

Khalid Jawed
khalidjm@seas.ucla.edu
