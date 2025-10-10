## Small package for Spectral-Band / Multi-Wavelength pyrometry

  This package is based on part of the code that was written for laboratory setup for measuring spectral emissivity. The details of the setup itself are described in the [article](http://dx.doi.org/10.1007/s00340-024-08331-9), the preprint is also [available](http://dx.doi.org/10.21203/rs.3.rs-4766080/v1)  
  
This repository contains code for processing experimentally measured thermal emission spectra of real surfaces, where the spectral emissivity varies with wavelength. Specifically, it enables the calculation of the surface temperature of an object with unknown emissivity based on its thermal radiation intensity. To ensure accurate intensity determination, the measured spectrum must be corrected according to the spectral sensitivity of the detection system, which requires a laboratory setup. The multiwavelength pyrometry method is implemented in `BandPyrometry.jl` module.

Previously, this package also included  `Planck.jl`, `JDXreader.jl` and `Pyrometers.jl`  modules for calculating the blackbody thermal emission, reading spectroscopic format JCAMP-DX and  virtual partial radiation pyrometers calculator. These packages are now moved to independent repositories [PlanckFunctions.jl](https://github.com/Manarom/PlanckFunctions.jl.git), [JCAMPDXir.jl](https://github.com/Manarom/JCAMPDXir.jl.git) and [RadiationPyrometers.jl](https://github.com/Manarom/RadiationPyrometers.jl.git) respectively.
  
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

