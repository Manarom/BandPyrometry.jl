### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 1621271d-075d-4034-a023-6d5ee1b3bee9
try 
	using ScaledPolynomials
	using Plots,StaticArrays,Interpolations,  PlutoUI, LinearAlgebra, PrettyTables, Optimization,OptimizationOptimJL
	using BandPyrometry 
	import PlanckFunctions as Planck
	using DelimitedFiles
catch notfound 
	import Pkg 
	Pkg.activate(@__DIR__())
	Pkg.instantiate()
	using Plots , StaticArrays, Interpolations,  PlutoUI, LinearAlgebra, PrettyTables, Optimization,OptimizationOptimJL
	using BandPyrometry 
	import PlanckFunctions as Planck
	using DelimitedFiles
	using LaTeXStrings
end

# ╔═╡ 30743a02-c643-4bdc-837e-b97299f9520a
md"""
#  Pluto notebook for testing the `BandPyrometry.jl` module  

##### Package **`BandPyrometry.jl`** contains methods to obtain the real surface's temperature from its thermal emission spectrum using multiwavelenght or band-pyrometry method. 
____________________
### Installation

To run this notebook, you need:
1) Install `julia` language itself from its official [download page](https://julialang.org/downloads) 
2) Install [Pluto](https://plutojl.org/) notebook from `julia` REPL by entering the following commands line-by-line:
```julia
import Pkg
Pkg.add("Pluto")
using Pluto
Pluto.run()
```
The last line will launch the Pluto starting page in your default browser 
3) Copy the entire GitHub [repository](https://github.com/Manarom/BandPyrometry.jl.git) to your local folder
4) Open this notebook file `BandPyrometry_test_git.jl` in `Pluto` by providing the full path to the *"Open a notebook"* text field on `Pluto`'s starting page. As far as `Pluto` has its own package manager, it will automatically install all necessary dependancies, which are marked in `using` cell of this file . 

"""

# ╔═╡ c167b60a-ae39-46e1-89b2-bef27b205490
plot_common_args = (grid = true, gridlinewidth=3, gridstyle = :dot,minorgrid=true, box = :on, linewidth = 3)

# ╔═╡ 171409eb-22b5-4bc5-a8e2-eac0932a24f3
PlutoUI.TableOfContents(indent=true, depth=4, aside=true)

# ╔═╡ d442014a-20e6-4be4-ac7f-f13de329dec5
md"""
## Part I. Introduction


#### I.I. The `blackbody` vs `real surface` thermal emission 
_______________________

All heated bodies emit thermal radiation. According to Planck's law, the spectrum of ideal emitter (called the *blackbody*) is governed solaly by its temperature. The blackbody spectral intensity can be calculated as follows \


``I_{blackbody}(\lambda , T) =  \frac{C_1}{\lambda ^5} \cdot \frac{1}{e^{\frac{C_2}{\lambda T } } - 1}``, \


where ``C_1`` = $(Planck.C₁), ``W \cdot μm/m² \cdot sr`` and ``C_2`` = $(Planck.C₂), ``μm \cdot K`` , ``\lambda`` - wavelength in ``\mu m``, ``T`` - temperature in Kelvins

A real surface thermal emission intensity is lower than the one of the blackbody. The fraction of blackbody thermal radiation intensity emitted by a real surface is characterized by directional spectral emissivity ``\epsilon (\lambda, T,\vec{\Omega})``:

``I_{real\ \ surface}(\lambda , T, \vec{\Omega} ) = \epsilon (\lambda, T,\vec{\Omega}) \cdot  I_{blackbody}(\lambda , T)``, \


here  ``\vec{\Omega}`` stays for direction. \


It is interesting that, unlike the blackbody, the real surface thermal emission (in general) depends  on the direction of radiation. Therefore, the most general characteristic for thermal radiation of a real surface is the `directional spectral emissivity`.

In [Planck.jl](https://manarom.github.io/BandPyrometry.jl/planck/) module there are several functions to calculate the blackbody thermal emission spectra (and various derivatives, integrals etc.).
The following figure show the impact of spectral emissivity on the real surface thermal emission intensity.
"""

# ╔═╡ 27b3c586-9eb0-4a51-b9ca-a9c0379fccdf
@bind  λ_BB PlutoUI.combine() do Child
	md"""
	Blackbody spectral range, ``\mu m`` : \
	``\lambda_{left}`` = $(
		Child(Slider(0.1:0.1:50,default=0.1,show_value = true))
	)   -- 
	 $(
		Child(Slider(0.1:0.1:50,default=25.0,show_value = true))
	)  ``\lambda_{right}`` 
	"""
end

# ╔═╡ f22d22b6-5d98-4cc4-998f-a53e92809618
@bind  T_BB PlutoUI.combine() do Child
	md"""
	Blackbody temperatures: \
	T₁ = $(
		Child(Slider(-273.0:1:2500,default=300,show_value = true))
	) ``^o C`` \
	T₂ = $(
		Child(Slider(-273.0:1:2500,default=900,show_value = true))
	)  ``^o C`` 
	"""
end

# ╔═╡ 7cc110e9-7655-4dfc-b1e0-ab3905866425
@bind  scales_BB PlutoUI.combine() do Child
	md"""
	Scales : \
	xscale = $(
		Child(Select([:identity,:ln,:log10]))
	)   yscale  
	 $(
		Child(Select([:identity,:ln,:log10]))
	)   
	"""
end

# ╔═╡ 8a066ee5-80e9-462f-9a61-15851468aa63
md"""
#### I.II. Partial radiation pyrometry
_______________________

As far as the `blackbody` thermal radiation energy strongly depends on temperature, this quantity can be used to measure the temperature of a real surface. This is the general idea of partial radiation pyrometry: **measure intensity to get the temperature**. 

Figure shows various typical partial radiation pyrometers available in `RadiationPyrometers.jl` together with their spectral tange:

![GitHub image](https://github.com/Manarom/RadiationPyrometers.jl/blob/main/notebooks/Pyrometers.png?raw=true)
"""

# ╔═╡ 9b08b767-7e8f-4483-9f2f-226022ce10e4
md"""
	The legend of the figure above shows the `measured` temperature for emissivity set to one for several common pyrometers types.  Spectral regions are shown with different colours. This figure was generated using [RadiationPyrometers.jl](https://github.com/Manarom/RadiationPyrometers.jl), a small package, which provides several function to work with `virtual` partial radiation pyrometers. More details on this picture are available in  `RadiationPyrometers.jl` package description.
	"""

# ╔═╡ 824a6af1-3f70-4de0-8a94-6c63663a546a
md"""
	#### I.III. Multiwavelength pyrometry
	_______________________

	Unlike classical partial radiation pyrometry, which requires setting a constant emissivity in some relatively narrow wavelength region, the **multiwavelength pyrometry**  in theory, allows one to obtain the temperature of a surface without knowing the emissivity. More about multiwavelength pyrometry can be found e.g. in [Multi-spectral pyrometry—a review](https://iopscience.iop.org/article/10.1088/1361-6501/aa7b4b). 

	In order to achieve this goal, the **multiwavelength pyrometry** assumes that the dependence of emissivity on wavelength in some spectral range can be described by some relatively small number (``N``) of parameters. Hence, if you have measured thermal radiation intensity at relatively large number (``M``) of wavelengths, and if ``M>N+1``, you can formulate the optimization problem in space of ``N+1`` optimization variables viz ``\vec{x}= \begin{bmatrix} a_0 , \dots,  a_{N-1} , T\end{bmatrix}``, here `` \begin{bmatrix} a_0 , \dots,  a_{N-1}  \end{bmatrix}`` are ``N``  emissivity approximation coefficients, and the ``(N+1)``'th optimization variable ``T`` is the temperature.	The **multiwavelength pyrometry** optimization problem has several features that can be utilized in order to obtain a computationally-effective algorithm:

	* First, the emissivity approximation is a linear problem, this means that emissivity approximation coefficients can be taken independent of both the optimization variables and the independent variables (wavelength)
	* Second, a highly non-linear term (viz Planck function) depends on only one of the optimizaiton variables
	* And third, the target function is the product of linear and non-linear optimization problems

	Mathematical consequencies of these features are described in this repository supplementary materials [download pdf](https://manarom.github.io/BandPyrometry.jl/assets/supplementary_v_0_0_1.pdf) in more details.

	In `BandPyrometry.jl` , the real surface thermal emission is approximated as a product of Planck function (ideal surface thermal emission) and linear (with respect to the optimization variables e.g. polynomial coefficients) approximation of spectral emissivity.	

	In `BandPyrometry.jl` the spectral emissivity is approximated as a linear combination of basis functions:
	"""

# ╔═╡ 5cc20c03-6c6e-4425-b974-242f69fe29be
L"""
\epsilon(λ)=\sum_{n=0}^{N-1}a_n \cdot \phi_n(λ)
	"""

# ╔═╡ ba2d0691-aeef-4d0a-808a-0c01a8b49e12
md"""
	where ``\phi_n`` is the basis function column vector, e.g. for standard basis it is:
	"""

# ╔═╡ cbaf05db-3265-4e2d-8ea4-759442f798e0
L"""
\mathbf{ \phi_n(λ)  = \begin{bmatrix} λ₁ⁿ  \\ \vdots \\ λ_{M}ⁿ \end{bmatrix}
}
"""

# ╔═╡ c2b82e9d-cae8-465c-bca0-160599e06102
md"""
	Full emission spectrum of a real surface is calculated as:
	"""

# ╔═╡ 478c8bb8-aa3c-4afe-b787-4924525858fa
L"""
	\mathbf{ 
	I_{real\:surface} (\lambda,T) =I_{blackbody}(\lambda,T) \cdot  \epsilon(λ)=I_{blackbody}(\lambda,T) \cdot \sum_{n=0}^{N-1}a_n \cdot \phi_n(λ)
	} 
		"""


# ╔═╡ 0831fddd-d2db-4ada-b293-00848d3673fb
md"Now, the optimization problem can be formulated:"

# ╔═╡ 4accbaec-4e08-43d3-8f36-2217b9394e86
L"""
	\begin{gather}

	\vec{x}^*=argmin\{F(\vec{x})\} \\

	F(\vec{x})=\sum_{i=1}^{M}[y_i - I_{blackbody}(\lambda_i,T) \cdot (\sum_{n=0}^{N-1}a_n \cdot \phi_n(λ_i))]^2 \\

	\vec{x} = [\vec{a},T]^t
	\end{gather} 

	"""

# ╔═╡ c47637c6-b243-4d8b-8234-40c68608939c
md"where ``\vec{a}`` is the column vector of emissivity approximation (  ``[]^t`` stays for transposition), ``\vec{x}^*`` is the local minimum"

# ╔═╡ d55793f5-45ff-4968-b2e1-8b846de8b91f
md"""

	To solve this optimization problem `BandPyrometry.jl` package provides functions to evaluate the discrepancy function ``F``, ``\nabla F`` and ``\nabla^2F`` which are needed to solve the optimization problem using zero, first or second order optimization algorithms. It also has special type to work with the emissivity linear approximation.
	"""

# ╔═╡ da6d5178-c083-4baa-91c3-e62f735bf808
md"""
#### I.IV. Emissivity approximation functions
_______________________

To approximate the emissivity `BandPyrometry.jl`  uses several polynomial bases:

* Standard basis (from [`Polynomials.jl`](https://juliamath.github.io/Polynomials.jl/stable/) package)
* Chebyshev basis (from [`Polynomials.jl`](https://juliamath.github.io/Polynomials.jl/stable/) package)
* Legendre polynomials (from [`LegendrePolynoials.jl`](https://jishnub.github.io/LegendrePolynomials.jl/stable/) package)
* Trigonometric basis: ``\phi_n(λ) = sin(\pi \cdot n \lambda)``  for odd n, and  ``\phi_n(λ) = cos(\pi \cdot n \lambda)`` for even n
* Bernstein basis of degree ``D``: ``\phi_{k}^{n}(\lambda) = (\begin{matrix} n \\ k \end{matrix}) (\frac {\lambda - a}{b - a})^k(\frac {b-\lambda}{b-a})^{n-k}`` - Bernstein basis function for ``\lambda \in [\lambda_a...\lambda_b]`` , ``k \in [0...n]`` , ``a = -1``, ``b = 1``

All basis vectors are stored in a structure called `VanderMatrix`, this type: 

1. Stores basis vectors for selected polynomial type and degreee  - columns of matrix ``V``: ``V= \begin{bmatrix} \vec{\phi_1} , \dots,  \vec{\phi_n} \end{bmatrix}``

2. The resulting emissivity can be calculated as a product of `VanderMatrix` and the vector of emissivity approximation polynomial coefficients vector: ``\vec{\epsilon}=V\cdot\vec{a}``

Independent variable vector ``\lambda`` is stored in normalized (in order to fit withing the range of -1...1) form, `VanderMatrix` also stores all data needed to return the original  ``\lambda`` vector. Matrix ``V`` can  be used to fit (in a least-square sense) by solving the overdetermined system of equations: ``\vec{a}= V^{\dagger} \cdot \vec{\epsilon}``, where  ``[]^{\dagger}`` means pseudo-inverse. The package implements this in  `polyfit` function (show docs) $(@bind show_polyfit_docs CheckBox(default=false)). 
"""

# ╔═╡ ae970ce5-fc5d-40d8-9f08-5dce0f99a509
show_polyfit_docs ? @doc(BandPyrometry.polyfit) : nothing

# ╔═╡ 529b07d7-e622-4816-8de9-e31581ea96a6
@bind  λ_fit_vand PlutoUI.combine() do Child
	md"""
	Emissivity fitting spectral range, ``\mu m`` : \
	``\lambda_{left}`` = $(
		Child(Slider(0.1:0.1:16,default=4.0,show_value = true))
	)   -- 
	 $(
		Child(Slider(0.1:0.1:16,default=8.0,show_value = true))
	)  ``\lambda_{right}`` 
	"""
end

# ╔═╡ 5880b468-64ec-4efa-9f6a-c39ecb090440
md"Select polynomial basis type $(@bind em_approx_poly_type Select(collect(keys(BandPyrometry.ScaledPolynomials.SUPPORTED_POLYNOMIAL_TYPES)),default = :stand))"

# ╔═╡ 988274c4-ad6a-42df-aead-5f40e0000998
md"Set polynomial degree : $(@bind poly_fit_degree Select(0:10,default=3)) (the polynomial degree = numer of basis functions - 1, thus, zero order polynomial is an all-units vector)"

# ╔═╡ 111122e9-1260-4f00-aba6-0fdf6cf04c4e
begin 
	N1 = 25
	λ_fit_vec = collect(range(λ_fit_vand...,length=N1)) 
	rt_emissivity_data = readdlm("real_surface_emissivity.txt") # loading file, actually this file is two-column data with no headers, but this is ok
	rt_emissivity_interpolation = linear_interpolation(rt_emissivity_data[:,1],rt_emissivity_data[:,2],extrapolation_bc = Interpolations.Flat())
	
	interpolated_values =rt_emissivity_interpolation(λ_fit_vec) # interpolating data at λ points
	PolyType = BandPyrometry.ScaledPolynomials.SUPPORTED_POLYNOMIAL_TYPES[em_approx_poly_type]{poly_fit_degree + 1 , Float64}
	Vander_test = BandPyrometry.ScaledPolynomials.VanderMatrix(SVector{N1}(λ_fit_vec),PolyType()) # creating new matrix 
	(a_fit_check,fitted_value,goodness_of_fit) = BandPyrometry.ScaledPolynomials.polyfit(Vander_test,λ_fit_vec,interpolated_values)# fitting polynomial coefficients
	plot(rt_emissivity_data[:,1],rt_emissivity_data[:,2],label = "ϵ real")
	scatter!(λ_fit_vec,fitted_value, label="ϵ fitted"; plot_common_args...)
	xlabel!("Wavelength, μm")
	ylabel!("Emissivity")
end

# ╔═╡ 4d6337aa-cfc7-4154-a395-5aa53e23d01a
begin 
	T₁ = T_BB[1] + Planck.Tₖ # converting to Kelvins
	T₂ = T_BB[2] + Planck.Tₖ # converting to Kelvins
	λ_bb = collect(range(λ_BB...,length=500));
	ibb1 = Planck.ibb.(λ_bb,T₁ );ibb2 =  Planck.ibb.(λ_bb,T₂);
	p_bb = plot(λ_bb,ibb1,label="blackbody T=$(T₁),K", xscale=scales_BB[1],yscale =scales_BB[2],fillrange=0, fillalpha=0.3)
	plot!(λ_bb,ibb2,label="blackbody T=$(T₂),K", xscale=scales_BB[1],yscale =scales_BB[2],fillrange=0, fillalpha=0.3)
	plot!(λ_bb,ibb1.*rt_emissivity_interpolation(λ_bb),label="real surface T=$(T₁),K", xscale=scales_BB[1],yscale =scales_BB[2],fillrange=0, fillalpha=0.2; plot_common_args...)	
	plot!(λ_bb,ibb2.*rt_emissivity_interpolation(λ_bb),label="real surface T=$(T₂),K", xscale=scales_BB[1],yscale =scales_BB[2],fillrange=0, fillalpha=0.2)
	xlabel!("Wavelength, μm")
	ylabel!("Spectral intensity, W/m²⋅sr⋅μm")
end

# ╔═╡ 03d76e64-ebf4-432b-b9be-d4cb26275f55
begin 
	e_real_BB = rt_emissivity_interpolation(λ_bb) # this spectral emissivity was measured up to 18 μm, thus for higher wavelengths it uses flat extrapolation
	plot(λ_bb, e_real_BB, label=nothing, linewidth=3.0; plot_common_args...)
	title!("Real (measured) surface spectral emissivity")
	xlabel!("Wavelength, μm");ylabel!("Spectral emissivity (ϵ)")
end

# ╔═╡ efbc4a5f-8853-47b1-8842-c77b060d2de7
pretty_table(HTML,hcat(["a$(i)" for i in 0:1:poly_fit_degree],a_fit_check), top_left_string="Table of the coefficients of emissivity linear approximation in the band from $(λ_fit_vand[1]) to $(λ_fit_vand[2]) μm using $(em_approx_poly_type) bases type,  the goodness of fit = $(goodness_of_fit)"  , column_labels = ["coeff","val"])

# ╔═╡ 38300c92-e5e6-4d5e-a394-aa1a47cfd757
md"""
	## Part II. `BandPyrometry.jl` testing

	In this notebook, the "measured" thermal emission spectrum is calculated as the same model  as the `BandPyrometryPoint.jl` internal representation, further this spectrum is fitted using optimization tools provided by **`Optimization.jl`** package, some random noise can be added (`optionally`). 

	#### II.I. "Measured" emissivity generation
	"""

# ╔═╡ 8a54d856-2298-4225-81ac-23bf66f35136
md"Measured spectrum fitting region:"

# ╔═╡ d14c4a0d-b641-457c-b1ed-491ff047a468
N = 50 # number of spectral points

# ╔═╡ 03a7a670-3ddd-4205-b49a-0488b55886ab
	BENCHMARK_DATA = [1.4	0.177627029	0.746041814;
			1.554268293	0.308975548	1.123202266;
			1.715243902	0.476864722	1.503503265;
			1.869512195	0.673742272	1.810258761;
			2.030487805	0.965854253	2.043896955;
			2.184756098	1.298184072	2.630983628;
			2.345731707	1.55399753	2.719482727;
			2.5	1.708743729	2.720488281];

# ╔═╡ 278c2da0-b396-493e-bdcd-fc2a539780b6
md"Select the polynomial type = $(@bind poly_type confirm(Select(collect(keys(BandPyrometry.ScaledPolynomials.SUPPORTED_POLYNOMIAL_TYPES)),default = :bernsteinsym)))"

# ╔═╡ d90da478-40b3-4677-82b9-fbfd79a72ef0
md"""Set the "measured" spectrum emissivity approximation coefficients:"""

# ╔═╡ 22eab67b-868a-44fd-9d40-69234f1ecb43
@bind  a_real PlutoUI.combine() do Child
	md"""
	a0 = $(
		Child(Slider(-1:0.01:1,default=0.6,show_value = true))
	) \
	a1 = $(
		Child(Slider(-1:0.01:1,default=0.8,show_value = true))
	)\
	a2 = $(
		Child(Slider(-1:0.01:1,default=0.8,show_value = true))
	)\
	a3 = $(
		Child(Slider(-1:0.01:1,default=0.8,show_value = true))
	)\
	a4 = $(
		Child(Slider(-1:0.01:1,default=0.8,show_value = true))
	)\
	a5 = $(
		Child(Slider(-1:0.01:1,default=0.8,show_value = true))
	)\
	"""
end

# ╔═╡ 8ef42759-fb53-41af-904e-8916924415fa
md"Set polynomial degree : $(@bind real_poly_degree confirm(Select(0:6,default=3))) (the polynomial degree= numer of basis functions-1, thus zero order polynomial is constant)"

# ╔═╡ a4dc0b5e-0ee7-4c22-8deb-eafdeb672207
md"""
The following two figures show:

1)basis vectors for selected polynomial (columns of `VanderMatrix.v`): ``V``

2)the resulting emissivity, calculated as a product of `VanderMatrix` and the vector of emissivity approximation polynomial: ``\vec{\epsilon}= \begin{bmatrix} \vec{\phi_1} , \dots,  \vec{\phi_n} \end{bmatrix}\cdot\vec{a} = V\cdot\vec{a}``

"""

# ╔═╡ 80dacea5-ea46-479f-b0d8-da9a9cf41aa2
md"""
#### II.II. Generating the "measured" spectral intensity
"""

# ╔═╡ 321172de-dcf6-4f51-ab31-edb37b8c3983
md"""
Real values of the optimization variables for "measured" thermal emission spectrum:  """

# ╔═╡ 0232ff3c-daa9-4e86-a6c7-582d66a16cb9
md" Use real emissivity $(@bind is_real_emissivity CheckBox(default = false))"

# ╔═╡ a7861092-024a-4011-820f-79835473a281
md" Use benchmark data $(@bind is_benchmark_data CheckBox(default = false))"

# ╔═╡ 2ef835b8-f62b-4254-9337-e7aa4d44c584
md"""
#### II.III. Starting optimization variables vector:
"""

# ╔═╡ 7f56174f-03f2-4996-a158-efcfc9ce9979
md"Set polynomial degree : $(@bind poly_degree Select(0:11,default=8)) (the polynomial degree= numer of basis functions-1, thus zero order polynomial is constant)"

# ╔═╡ 09f16e07-0f21-4bcc-a6e8-cdba3097d57b
@bind  a_starting PlutoUI.combine() do Child
	md"""
	a0 = $(
		Child(Slider(-1:0.01:1,default=0.3,show_value = true))
	) \
	a1 = $(
		Child(Slider(-1:0.01:1,default=0.1,show_value = true))
	)\
	a2 = $(
		Child(Slider(-1:0.01:1,default=0.2,show_value = true))
	)\
	a3 = $(
		Child(Slider(-1:0.01:1,default=0.2,show_value = true))
	)\
	a4 = $(
		Child(Slider(-1:0.01:1,default=0.2,show_value = true))
	)\
	a5 = $(
		Child(Slider(-1:0.01:1,default=0.2,show_value = true))
	)\
	a6 = $(
		Child(Slider(-1:0.01:1,default=0.1,show_value = true))
	)\
	a7 = $(
		Child(Slider(-1:0.01:1,default=0.2,show_value = true))
	)\
	a8 = $(
		Child(Slider(-1:0.01:1,default=0.2,show_value = true))
	)\
	a9 = $(
		Child(Slider(-1:0.01:1,default=0.2,show_value = true))
	)\
	a10 = $(
		Child(Slider(-1:0.01:1,default=0.2,show_value = true))
	)\
	
	"""
end

# ╔═╡ acc77ec3-6afb-4e76-ad2f-ec137555d1bc
md"""
Starting temperature , K = $(@bind T_starting	Slider(range(273,2573,length=1000), show_value=true,default=600))
"""

# ╔═╡ 808f8601-90fc-4cb1-b455-46cfc11b8fc7
begin 
	x_starting = MVector{poly_degree + 2}(zeros(poly_degree + 2))
	x_starting[end] = T_starting
	x_starting[1:end-1] .= a_starting[1:poly_degree+1]
	x_starting
end

# ╔═╡ 9016369d-bc8b-4907-bdb2-f4e81c444d30
md"""
#### II.IV. Solving the optimization problem
"""

# ╔═╡ 6c38418d-9c4d-4bb6-a386-dd0bde9af8d9
md"""
"Measured" spectrum temperature, K = $(@bind T_to_fit Slider(range(273,2573,length=1000), show_value=true,default=1000))
"""

# ╔═╡ 36c655f8-141b-4472-96d7-01c8fd1f1515
x_real_data = [a_real[1:real_poly_degree+1]...,T_to_fit];

# ╔═╡ c4df3ad4-9f1a-414c-b02c-8161f007ccd5
L"""%$(x_real_data)"""

# ╔═╡ 01cd1f0d-15b8-474b-a05a-eec840c54fff
md"""
	add noise to the "measured" spectrum = $(@bind noise_amplitude  Slider(0.0:0.0001:0.5,default = 0.0,show_value = true))
	"""

# ╔═╡ 6f17606c-e52e-4913-87f6-56190d209308
@bind  lam_region confirm(PlutoUI.combine() do Child
	md"""
	λₘᵢₙ ,μm = $(
		Child(Slider(0.01:0.1:30,default=8.0,show_value = true))
	) \
	λₘₐₓ ,μm = $(
		Child(Slider(0.0:0.1:30,default=13.0,show_value = true))
	) 
	"""
end)

# ╔═╡ 3b4308dc-da18-43ab-9f50-2c8bc05acb78
begin
	if !is_benchmark_data
		λ = SVector{N}(collect(range(lam_region[1],lam_region[2],length=N)))

		RealPolyType = BandPyrometry.ScaledPolynomials.SUPPORTED_POLYNOMIAL_TYPES[poly_type]{real_poly_degree + 1,Float64}
		
		VVV= BandPyrometry.ScaledPolynomials.VanderMatrix(λ,RealPolyType())
		Ib = Planck.ibb.(λ,T_to_fit) # bb spectrum
		ϵ_data = if is_real_emissivity 
				rt_emissivity_interpolation(Vector(λ))
			else
				VVV*x_real_data[1:end-1]# "measured" spectrum emissivity 
			end
		I_data = Ib.*ϵ_data .+ noise_amplitude*Ib.*(0.5 .-rand(length(λ)))# adding noise to the "measured" spectrum
	else
		
	end
	
end;

# ╔═╡ 4f2db18c-6f48-4c41-9a53-470042decf5f
begin 
	p_em = plot(λ, VVV*[i for i in a_real[1:real_poly_degree + 1]],label=nothing,linewidth=6; plot_common_args...)
	title!(raw"""Generated "real" surface spectral emissivity""")
    xlabel!("Wavelength, μm")
	ylabel!("ϵ")
	p_em
end

# ╔═╡ d17a5db7-0125-48b5-9016-76da6d72c673
begin
	
	
	p_mode = plot(VVV.xi,VVV.v[:,1], label="n= "*string(0),linewidth=3)
	for (i,V) in enumerate(eachcol(VVV.v[:,2:end]))
		plot!(VVV.xi,V, label="n= "*string(i),linewidth=3;plot_common_args...)
	end
	title!("Emissivity fitting modes for $(poly_type) -type polynomial")
	xlabel!("Normalized wavelength")
	ylabel!("Normalized amplitude")
	p_mode
end

# ╔═╡ 1c60da73-a088-45ac-a08d-885bb1d98dcc
begin 
	CheckType = BandPyrometry.ScaledPolynomials.SUPPORTED_POLYNOMIAL_TYPES[poly_type]{poly_degree + 1,Float64}
		
	VV_check = BandPyrometry.ScaledPolynomials.VanderMatrix(λ,CheckType())
	Ib_check = Planck.ibb.(λ,T_starting) # bb spectrum
	ϵ_data_check = VV_check*x_starting[1:end-1]# "measured" spectrum emissivity 
	I_data_check = Ib_check.*ϵ_data_check
	plot(λ,I_data,label="Real")
	plot!(λ,I_data_check,label = "check",title = "Spectral intensity",xlabel="Wavelength,μm",ylabel = "I"; plot_common_args...)
end

# ╔═╡ 08cad9ec-ce1c-4c2e-9758-58a1988cd1d3
begin 
	plot(λ,ϵ_data_check,title= "Spectral emissivity",xlabel = "wavelength, μm", ylabel = "ϵ", label = "Starting emissivity")
	plot!(λ,ϵ_data,label = "real emissivity";plot_common_args...)
end	

# ╔═╡ e9b7178b-f321-477f-bc63-060b9c96f4a0
# BandPyrometry.fitting_error(band_fit)

# ╔═╡ c3fc6fbf-5358-4d68-acf1-40fbc82c13dc
md"""Choose the optimization method $(@bind optim_type confirm(Select(["NelderMead","BFGS","GradientDescent", "LBFGS","NewtonTrustRegion","ParticleSwarm","Newton","IPNewton" ],default = "LBFGS")))"""

# ╔═╡ 50517cca-c1c9-4ecc-b6c7-e0549ea1e14a
@bind constraint_type confirm(Select(["constraint","unconstraint" ],default="unconstraint"))

# ╔═╡ 0ca01e99-bf48-4e24-b937-c3863eb03f50
is_constraint = constraint_type=="constraint";

# ╔═╡ c01d5d8c-06f7-42b2-b845-e3ec0e82078e
md"""
Spectral emissivity range:


ϵₗ = $(@bind eps_lower confirm(Slider(0:1e-3:1,show_value = true,default = 0.1)))

ϵᵣ = $(@bind eps_upper confirm(Slider(0:1e-3:1,show_value = true,default = 0.99)))
"""

# ╔═╡ c9253144-b8c6-4a24-98d4-3faa4175d96a
md"""
Constraint temperature: $(@bind is_temp_constr CheckBox(false))

Tₗ = $(@bind t_lower confirm(Slider(300:1.0:5000,show_value = true,default = 500)))

Tᵣ = $(@bind t_upper confirm(Slider(300:1.0:5000,show_value = true,default = 1500)))
"""

# ╔═╡ 15025715-1ff3-4e7c-8136-cc442b77528a
md"""
Selected optimization method is $(optim_type)
the actual optimizer is
$(BandPyrometry.optimizer_switch(optim_type,is_box_constraint=is_constraint,is_lagrange_constraint = false)) 
"""
# this function returns the methods which was actually used (some methods from optim_type list does not support constraints)

# ╔═╡ 6118cebc-d18e-42cc-ac1b-aa1ad95cd719
begin 
	
	band_fit =BandPyrometry.BandPyrometryPoint(I_data,λ,x_starting,polynomial_type=poly_type)# creating BandPyrometryPoint from the data to be fitted, wavelength region and starting optimization variables vector x_starting
	e_range = (eps_lower,eps_upper)
	t_range = is_temp_constr ? (t_lower,t_upper) : nothing
	@time out = BandPyrometry.fit_T!(band_fit,optimizer_name=optim_type,
										  is_lagrange_constraint = false,
										 is_box_constraint = is_constraint,
										 emissivity_range = e_range,
										 temperature_range = t_range) # this function runs the optimizaiton is_box_constraint=is_constraint,
end # starting values

# ╔═╡ c3e9e5e5-0d1b-4329-a5b6-020e0e7321b8
band_fit.hessian[end]

# ╔═╡ 950a61f6-fe83-47a3-b65f-fd4b70ff12ba
plot(band_fit;plot_common_args...)

# ╔═╡ a73dc68d-7e41-4748-b0bb-458d6ef73305
T_fitted = BandPyrometry.temperature(band_fit)

# ╔═╡ c947d126-092b-4b92-ab56-373d5a387908
begin 
	p3 = plot(λ, I_data, lw=2,label="Data T=  $(T_to_fit) ");
	plot!(p3,λ, band_fit.Ic, lw=2,label="Fitting T= $(T_fitted)", ls=:dot; plot_common_args...)
	title!("Initial data vs fitting results")
	xlabel!("Wavelength, μm") 
	ylabel!("Emission intensity")
end

# ╔═╡ 095cbaa2-9dfd-45d9-81a2-c075f7ba4838
(lb1,ub1) = BandPyrometry.evaluate_box_constraints(band_fit, e_range,t_range)

# ╔═╡ 73c85c85-2fc4-48ac-b5f4-4c35edcf935c
begin 
	a_vect_fitted = out[2]
	a_vect_real = x_real_data[1:end-1]
end;

# ╔═╡ ebc0b228-01c8-42e4-9388-8194be1dd669
md"""
	Goodness of fit:

	"""

# ╔═╡ 5e46ca94-1ca2-4af1-88c3-54c7e86d1aef
if length(x_real_data) == length(a_vect_fitted) + 1
	L"""
	\begin{bmatrix} T_{fitted} = %$(T_fitted) & T_{real} = %$(T_to_fit) & \Delta T /T = %$(abs(T_fitted-T_to_fit)/T_to_fit) \\ |\vec{a}_{fitted}| = %$(norm(a_vect_fitted)) & |\vec{a}_{real}| = %$(norm(x_real_data[1:end-1])) & |\vec{a}_{fitted} -\vec{a}_{real}|= %$(norm(a_vect_fitted - a_vect_real))\end{bmatrix}

"""
end

# ╔═╡ e9b7169f-0d75-44bb-8505-ccfde1bdce98
e_fitted_spectrum = band_fit.ϵ;

# ╔═╡ fead76f5-c0b5-4fcb-a39e-c97ded034653
begin
	p1 = plot(λ,ϵ_data,label="Real emissivity "); #a=  $(map(i->@sprintf("%.2f",i),x_data[1:end-1]) )
	scatter!(p1,λ, e_fitted_spectrum,label="Fitted emissivity") #, a=  $(map(i->@sprintf("%.2f",i), a_vect_fitted ))
	if is_constraint 
		plot!(p1,λ,eps_lower*λ./λ,label = "lower bound",linestile = :dash,linewidth = 3)
		plot!(p1,λ,eps_upper*λ./λ,label = "upper bound",linestile = :dash,linewidth = 3;plot_common_args...)
	end
	ylims!(0, 1.2)
	xlabel!("Wavelength, μm") 
	ylabel!(L"\epsilon")
end

# ╔═╡ 3f1e1a4c-48bf-44fa-a146-020dde04d2ff
md"Final discrepancy value: $(BandPyrometry.disc(out[4].u, band_fit ))"

# ╔═╡ 027f22c1-5309-49d4-86ac-53791c2e0eaf
begin	
	if is_benchmark_data
		ni = 1
		r = 1:ni:8
		N_bench = length(r)
		#NPOINTS = 30
		SV = SVector{N,Float64}
		#lam_bench = SV(data[r,1])
		lam_init = BENCHMARK_DATA[r,1]
		lam_bench = SV(range(extrema(lam_init)...,N)) #linear_interpolation
		I_interp = linear_interpolation(lam_init,BENCHMARK_DATA[r,2])(Vector(lam_bench))
		I_meas_bench = SV(1e3*(I_interp))
		Planck.ibb.(lam_bench,1050.0)


		p_bench = plot(lam_bench,I_meas_bench, label = "init", markershape=:diamond)
		T_out = Float64[]
		e_bench = Matrix{Float64}(undef,N,7)
		for n in 2:N_bench
		    x_start = zeros(n)
		    x_start[1] = 0.5
		    x_start[end] = 1000.0
		    bm_bench = BandPyrometry.BandPyrometryPoint(I_meas_bench,lam_bench, SVector{n,Float64}(x_start))
		    out_bench = BandPyrometry.fit_T!(bm_bench,
		            optimizer_name="ParticleSwarm",
		            is_constraint=true)
		    push!(T_out,out_bench[1])
		    e_bench[:,n-1] = bm_bench.ϵ
		    @show out_bench[4]
		    plot!(p_bench,lam_bench,bm_bench.Ic,label = string(n))
		end
		@show T_out
		p=plot(e_bench)
		p_bench
	end
end

# ╔═╡ Cell order:
# ╟─30743a02-c643-4bdc-837e-b97299f9520a
# ╠═1621271d-075d-4034-a023-6d5ee1b3bee9
# ╠═c167b60a-ae39-46e1-89b2-bef27b205490
# ╟─171409eb-22b5-4bc5-a8e2-eac0932a24f3
# ╟─d442014a-20e6-4be4-ac7f-f13de329dec5
# ╟─27b3c586-9eb0-4a51-b9ca-a9c0379fccdf
# ╟─f22d22b6-5d98-4cc4-998f-a53e92809618
# ╟─7cc110e9-7655-4dfc-b1e0-ab3905866425
# ╠═4d6337aa-cfc7-4154-a395-5aa53e23d01a
# ╠═03d76e64-ebf4-432b-b9be-d4cb26275f55
# ╠═8a066ee5-80e9-462f-9a61-15851468aa63
# ╟─9b08b767-7e8f-4483-9f2f-226022ce10e4
# ╟─824a6af1-3f70-4de0-8a94-6c63663a546a
# ╟─5cc20c03-6c6e-4425-b974-242f69fe29be
# ╟─ba2d0691-aeef-4d0a-808a-0c01a8b49e12
# ╟─cbaf05db-3265-4e2d-8ea4-759442f798e0
# ╟─c2b82e9d-cae8-465c-bca0-160599e06102
# ╟─478c8bb8-aa3c-4afe-b787-4924525858fa
# ╟─0831fddd-d2db-4ada-b293-00848d3673fb
# ╟─4accbaec-4e08-43d3-8f36-2217b9394e86
# ╟─c47637c6-b243-4d8b-8234-40c68608939c
# ╟─d55793f5-45ff-4968-b2e1-8b846de8b91f
# ╟─da6d5178-c083-4baa-91c3-e62f735bf808
# ╟─ae970ce5-fc5d-40d8-9f08-5dce0f99a509
# ╟─529b07d7-e622-4816-8de9-e31581ea96a6
# ╟─5880b468-64ec-4efa-9f6a-c39ecb090440
# ╟─988274c4-ad6a-42df-aead-5f40e0000998
# ╠═111122e9-1260-4f00-aba6-0fdf6cf04c4e
# ╠═efbc4a5f-8853-47b1-8842-c77b060d2de7
# ╟─38300c92-e5e6-4d5e-a394-aa1a47cfd757
# ╟─8a54d856-2298-4225-81ac-23bf66f35136
# ╟─d14c4a0d-b641-457c-b1ed-491ff047a468
# ╟─03a7a670-3ddd-4205-b49a-0488b55886ab
# ╟─278c2da0-b396-493e-bdcd-fc2a539780b6
# ╟─d90da478-40b3-4677-82b9-fbfd79a72ef0
# ╟─22eab67b-868a-44fd-9d40-69234f1ecb43
# ╟─8ef42759-fb53-41af-904e-8916924415fa
# ╟─a4dc0b5e-0ee7-4c22-8deb-eafdeb672207
# ╟─4f2db18c-6f48-4c41-9a53-470042decf5f
# ╟─d17a5db7-0125-48b5-9016-76da6d72c673
# ╟─80dacea5-ea46-479f-b0d8-da9a9cf41aa2
# ╟─36c655f8-141b-4472-96d7-01c8fd1f1515
# ╟─321172de-dcf6-4f51-ab31-edb37b8c3983
# ╟─c4df3ad4-9f1a-414c-b02c-8161f007ccd5
# ╟─0232ff3c-daa9-4e86-a6c7-582d66a16cb9
# ╟─a7861092-024a-4011-820f-79835473a281
# ╟─3b4308dc-da18-43ab-9f50-2c8bc05acb78
# ╟─2ef835b8-f62b-4254-9337-e7aa4d44c584
# ╠═08cad9ec-ce1c-4c2e-9758-58a1988cd1d3
# ╠═7f56174f-03f2-4996-a158-efcfc9ce9979
# ╟─09f16e07-0f21-4bcc-a6e8-cdba3097d57b
# ╠═1c60da73-a088-45ac-a08d-885bb1d98dcc
# ╟─acc77ec3-6afb-4e76-ad2f-ec137555d1bc
# ╟─808f8601-90fc-4cb1-b455-46cfc11b8fc7
# ╟─9016369d-bc8b-4907-bdb2-f4e81c444d30
# ╟─6c38418d-9c4d-4bb6-a386-dd0bde9af8d9
# ╟─01cd1f0d-15b8-474b-a05a-eec840c54fff
# ╠═c947d126-092b-4b92-ab56-373d5a387908
# ╟─6f17606c-e52e-4913-87f6-56190d209308
# ╟─c3e9e5e5-0d1b-4329-a5b6-020e0e7321b8
# ╠═fead76f5-c0b5-4fcb-a39e-c97ded034653
# ╠═950a61f6-fe83-47a3-b65f-fd4b70ff12ba
# ╟─e9b7178b-f321-477f-bc63-060b9c96f4a0
# ╟─c3fc6fbf-5358-4d68-acf1-40fbc82c13dc
# ╟─50517cca-c1c9-4ecc-b6c7-e0549ea1e14a
# ╠═0ca01e99-bf48-4e24-b937-c3863eb03f50
# ╟─c01d5d8c-06f7-42b2-b845-e3ec0e82078e
# ╟─c9253144-b8c6-4a24-98d4-3faa4175d96a
# ╟─15025715-1ff3-4e7c-8136-cc442b77528a
# ╟─6118cebc-d18e-42cc-ac1b-aa1ad95cd719
# ╠═a73dc68d-7e41-4748-b0bb-458d6ef73305
# ╠═095cbaa2-9dfd-45d9-81a2-c075f7ba4838
# ╟─73c85c85-2fc4-48ac-b5f4-4c35edcf935c
# ╠═ebc0b228-01c8-42e4-9388-8194be1dd669
# ╠═5e46ca94-1ca2-4af1-88c3-54c7e86d1aef
# ╟─e9b7169f-0d75-44bb-8505-ccfde1bdce98
# ╠═3f1e1a4c-48bf-44fa-a146-020dde04d2ff
# ╠═027f22c1-5309-49d4-86ac-53791c2e0eaf
