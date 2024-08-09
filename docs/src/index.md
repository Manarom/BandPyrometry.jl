
# BandPyrometry.jl

## General description

```@raw html
<p>Detailed description of the underlying mathematics is available at:<a href="assets/supplementary_v_0_0_1.pdf">Download PDF</a>.</p
```

- Planck.jl module contains several functions to evaluate: 
    - Blackbody spectral intensity (Planck function)
    - The first and the second derivatives of Planck function with respect to the wavelength and temperature 
    - Planck function fast integration over the specified wavelength region  

- BandPyrometry.jl module provides:
    - EmPoint type formulates the least-square problem to optimize the temperature of blackbody 
    - BandPyrometryPoint type formulates the least-square problem to optimize the temperature and spectral emissivity approximation to the real surface (with emissivity dependent on wavelength) formulates the least-square problem to optimize the temperature of blackbody 
    - Methods to solve the least-square optimization problem using zero,first and second order methods using optimizers provided by Optim.jl package 

- JDXreader.jl module provides:
    - Function to read files in spectroscopic format (JCAMP-DX=4.24 version) 

- Pyrometers.jl module
    - Brings Pyrometer type (partial radiation pyrometer)

## Contact

To contact me, please do it through the [GitHub repository](https://github.com/Manarom/BandPyrometry).

## License

Copyright (c) 2024 Roman Mironov

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.

