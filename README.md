# NovacPPP
This is the code for the Novac Post Processing Program (aka NovacPPP or NovacP3), 
which is a program used to batch-process spectra after they have been collected in order to improve
the spectral evaluation and/or the estimated fluxes. Compared to the NovacProgram is the NovacPPP allowed
to take some  more time evaluating the spectra in order to achieve a better analysis, it has more and more advanced
options for performing the evaluations and it is possible to run the program on all collected Novac data in one
run.

The program has during 2018 been modified to be platform independent and has been tested to build on both:
* Windows 10, using Visual Studio 2017, 2019 and 2022 (community edition)
* Ubuntu 16.04, and Ubuntu 18.04 using gcc and make
The build system used to achieve this is [CMake](https://cmake.org/) and scripts are included to build the program in
Windows using visual Studio and in Ubuntu using gcc and make. 

In order to make the program platform independent was the user interface removed entirely. The program is now console only.

The program depends on [Poco libraries](https://pocoproject.org/) for performing common networking and filesystem tasks 
and Poco must therefore be installed on the system prior to building.

The main tool is the NovacPPP executable, which is found in the PPPExe folder. The NovacPPP tool is configured by a set of 
xml-files which must reside in a sub directory named _configuration_ in the same directory as the executable. 
Examples of such configurations can be found in the directories _Debug/configuration_ and _Release/configuration_ in the
source tree. 

## Setting up

1. Clone the repository from your own fork

2. The SpectralEvaluation is setup as a git sub-module.
Before you can start building the program you need to setup the sub-module.
This is done by opening a command prompt and typing the two commands 'git submodule init' and 'git submodule update'

3. Copy the pre-built binaries of Poco to ../poco/install/ (this is where CMake will search for them)

4. Run the desired build script to create Visual Studio solution or a Unix make file.
Build_UnixMake.sh builds the entire solution (and creates executables) for Linux.
Build_VisualStudio20XX.cmd creates a Visual Studio solution which can then be opened in Visual Studio to build an executable.

