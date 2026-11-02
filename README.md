## Small package for Spectral-Band / Multi-Wavelength pyrometry

This package is based on part of the code written for a laboratory setup designed to measure spectral emissivity. The details of the setup are described in this [article](http://dx.doi.org/10.1007/s00340-024-08331-9), and the preprint is also [available](http://dx.doi.org/10.21203/rs.3.rs-4766080/v1).

> **Warning:** The development of this package has been suspended in favor of the [Pyrometers.jl](https://github.com/Manarom/Pyrometers.jl) package, which now includes various types of convenient pyrometers and a multi-wavelength pyrometer.

This repository contains code for processing experimentally measured thermal emission spectra of real surfaces where the spectral emissivity varies with wavelength. Specifically, it enables the calculation of the surface temperature of an object with unknown emissivity based on its thermal radiation intensity. To ensure accurate intensity determination, the measured spectrum must be corrected according to the spectral sensitivity of the detection system, which requires a laboratory setup calibration. The multiwavelength pyrometry method is implemented in the `BandPyrometry.jl` module.

Previously, this package also included the `Planck.jl`, `JDXreader.jl`, and `Pyrometers.jl` modules for calculating blackbody thermal emission, reading the spectroscopic JCAMP-DX format, and calculating virtual partial radiation pyrometers. These modules have now been moved to independent repositories: [PlanckFunctions.jl](https://github.com/Manarom/PlanckFunctions.jl.git), [JCAMPDXir.jl](https://github.com/Manarom/JCAMPDXir.jl.git), and [Pyrometers.jl](https://github.com/Manarom/Pyrometers.jl.git) respectively.

This package also depends on the unregistered package [`ScaledPolynomials.jl`](https://github.com/Manarom/ScaledPolynomials.jl.git), which provides a custom implementation of Bernstein basis polynomials (among others) useful for the constrained optimization of spectral emissivity.

The `.//notebooks` folder contains two [Pluto](https://plutojl.org/) notebooks that can be used as examples of package usage.

Full documentation is available at https://manarom.github.io/BandPyrometry.jl/

# Installation

1) Download [Julia](https://julialang.org/downloads)

#### For usage

2a) Using the Julia package manager:

```julia
import Pkg
# ScaledPolynomials.jl is an unregistered package for polynomials
Pkg.add(url = "https://github.com/Manarom/ScaledPolynomials.jl.git") 
# Link to the package itself
Pkg.add(url = "https://github.com/Manarom/BandPyrometry.jl.git") 
```

#### For development

2c) Clone this repository to `username/.julia/dev/`.

3c) Enter the package manager in the REPL by pressing `]` and then add the package by typing `dev BandPyrometry`.