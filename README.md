## Small package for Spectral-Band / Multi-Wavelength pyrometry

  This package is based on part of the code that was written for laboratory setup for measuring spectral emissivity. The details of the setup itself are described in an article that (I hope) will be published in the near future. 
  
  This repository codes can be used to process experimentally measured thermal emission spectra of real surfaces, that is, those for which the spectral emissivity depends on the wavelength. In particular, it allows one to calculate the surface temperature of a real object with an unknown emissivity from its thermal radiation intensity (for correct determination of intensity, the measured spectrum should be corrected according to the spectral sensitivity of the receiving system, thats why we need a laboratory setup!). Multiwavelength pyrometry approach is realized in the general module `BandPyrometry.jl`. 
The package also includes `Planck.jl` module for calculating the blackbody thermal emission intensity (aka Planck function) and its derivatives. 
And also, two small modules called `JDXreader.jl` and `Pyrometers.jl`. The first one is for reading files in the JCAMP-DX=4.24 format. This format is very quite common in spectroscopy and has a large number of variations, so my module does not pretend to be a fully developed reader for this format, but it worked well for FTIR spectrometer which was available for me. The second one is for partial radiation pyrometers. Several common pyrometers types together with their working spectral ranges are brought together. 
  
  In `/test` folder there are two [Pluto](https://plutojl.org/) notebooks, which can be used as examples of package usage.


  Full documentation is available at https://manarom.github.io/BandPyrometry.jl/

# Installation 

1) Download [Julia](https://julialang.org/downloads)

#### For usage

2a) Clone this repository to your local machine 

3a) use `include("BandPyrometry_folder\src\BandPyrometry.jl)` in REPL or other module's to bring BandPyrometry module to the global scope

4a) use `using .BandPyrometry` to bring the module and its content to the corresponding namespace

#### For development

2c) Clone this repository to `username/.julia/dev/`.

3c) Enter the package manager in REPL by pressing `]`  then add the package by typing `dev BandPyrometry`


This is my first repository on GitHub, so some problems are quite possible, but I hope that these codes will be useful to someone.
