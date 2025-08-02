
Scilab/Scicos Toolbox for Adaptive Optics  
===========================================


Note
=====

A legacy code from https://sciao.sourceforge.net/

This repository is obsolete and is not maintained anymore. 

scitao is a toolbox based on Scilab/Scicos environment for the simulation 
of wave optics, especially for the simulation of adaptive optics.

Installation
============

To Install this toolbox: 
you must have scilab (V4) installed in your system.

We Suppose here that <PATH> stands for  the path of the directory
containing this README file.

- On Unix/Linux systems
	
1. The scitao toolbox need some third party Libraries to support,
	you should first download and install them on you system. All the 
	required three third party libraries are:

	(1) Cfitsio v2.401 : 
		http://heasarc.gsfc.nasa.gov/docs/software/fitsio/fitsio.html
		
	(2) Fastest Fourier Transform in the West v3.0
		http://www.fftw.org/
		
		note: need compile two versions of this library: the (default) 
		double version, and a float version.
		
	(3) LAPACK v 3.0+ (May 2000 release):
		http://www.netlib.org/lapack/
		
	note: there are some functions of LAPACK library functions in standard Scilab
	but it not complete, we need the complete LAPACK library to surport scitao toolbox.

2. Then you can install scitao toolbox follow the standard linux install sequences (./configure && make && make install). See the file 	INSTALL for further instructions on how to build and install it.
	
	Users can also execute the following instruction within Scilab:	
	exec <PATH>/loader.sce  before using the toolbox, he can also put it 
	in his .scilab startup file for automatic loading.

- On Windows systems

	We have prepare a setup binary file, you need only click it and 
	all needed file will be installed on the proper directive. If you like to 
	compile the source code (need c/c++ compiler),you can execute file 
	builder.sce to build it.

	Users can also execute the following instruction within Scilab:
	exec <PATH>\loader.sce. before using the toolbox, he can also
	put it in his .scilab startup file for automatic loading.
     
		Note: to load scitao toolbox, you should copy all the 
      *.dll files in the directory of src/arroyo/win-util into the 
      SCI/bin directory. 


Contents
========

- README.txt				: This file
- AUTHORS.txt				: The information of the AUTHORS
- COPYING.txt				: GPL license
- INSTALL.txt				: installation instructions (for Linux systems)

- loader.sce				: Installation script
- builder.sce				: Script for building library
- win-builder.sce			: a scilab script for building under windows system.

- src					: Directory of c/c++ routines

	- arroyo				: Arroyo c++ class library      
	- lightPipes			: lightPipes function library     
	- wrapper				: some c wrapper function library for Arroyo library     
	- interf				: optical interface function     
	- scicos				: computational function for optical scicos diagrams
     
- macros					: Directory of optics and scicos-block interface function

     - *.sci				: Source versions
     - *.bin				: Precompiled binary versions (generated)
     - names				: Table of functions (generated)
     - lib				: Scilab library binary save (generated)
     
- man					: Directory for help.

     - optics				: help files for optical Scilab functions 
     - optics_scicos			: help files for optical Scicos-blocks 
     - builder.sce			: Script for building help
     - loader.sce			: Script for load help 
     
- demos					: demos directory.

- examples				: Contain some simple examples for learn and use.
