# Multiwavelength pyrometry package

## More detailed description
  This package is based on part of the code that was written to work with a laboratory installation for measuring spectral emissivity. The details of the laboratory setup itself are described in an article that I hope will be published in the near future. This program can be used to process experimentally measured thermal emission spectra of real surfaces, that is, those for which the spectral emissivity depends on the wavelength. In particular, it allows one to calculate the surface temperature of a real object with an unknown emissivity from the spectrum of its thermal radiation intensity (for a correct determination of the intensity, the measured spectrum must be corrected for the spectral sensitivity of the receiving system using a reference emitter, thats why we need the laboratory setup!), this approach is implemented in the main module 'BandPyrometry.jl'. 
In addition to BandPyrometry.jl module, the package also includes a 'Planck.jl' module for calculating the blackbody thermal emission intensity viz Planck function and its derivatives. And also a small module for reading files in the Dcamp Zhdeixreader format, which is often used in spectroscopy, and a small module with a set of basic types of partial radiation pyrometers Pyrometers zh.l.
In addition to the codes themselves, the test folder also contains two pluto laptops that can be used as examples

This is my first repository on GitHub, so some problems are quite possible, but I hope that these programs will be useful to someone.

Documentation is available at https://manarom.github.io/BandPyrometry/
