# Functions and type to fit the temperature of blackbody and real surface thermal emission

This Module introduces two main types: 
    - EmPoint for blackbody thermal emission 
    - BandPyrometryPoint for real surface (with spectrally dependent emissivity) 
This types can be used to solve the least-square optimization problem to fit the temperature 
and spectral emissivity 

```@autodocs
    Modules = [BandPyrometry]
    Order   = [:function, :type]
```