# solar-eclipse-8.5.1
solar-eclipse-8.5.1 is built and installed using the install_solar.sh script located in the main directory.  
It's recommended that prior to running the script the options should be viewed by running: install_solar.sh --help
Prior to running the script please ensure you have the Intel Math Kernel library installed.  A required input of the
script is the location of the lib and include folders for MKL.  By default MKL usually installs itself in /opt/
which means the path that would be entered into the script is /opt/intel/mkl .  It's recommended that the g++ and gcc compilers
be version 7 or less (so long as c99 and c++11 are available).  Building solar on Mac computers will most likely require third
party g++,gcc, and gfortran compilers downloaded through homebrew. By default the name of the startup script to use solar is named solar.
This can be changed using the --script_name option.
