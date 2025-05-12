### A Pluto.jl notebook ###
# v0.20.8

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

# ╔═╡ fc6d97d3-3775-42d2-a991-09d875b0dd38
begin# we need the project structrue inside the notebook
	notebook_path = @__DIR__()# this notebook local path
	project_path = abspath(joinpath(notebook_path,".."))#project's path
	import Pkg
	Pkg.activate(project_path)
	sources_path = joinpath(project_path,"src")# it is supposed that sources are in separate folder \project_folder\src\
end;

# ╔═╡ c6f4b5c0-e007-4d4f-8b79-b3380cc5e10c
using StaticArrays, 
	PlutoUI, Plots, 
	Optimization,
	OptimizationOptimJL,BenchmarkTools, 
	LinearAlgebra, Interpolations,Printf,
	MKL,LinearAlgebra, LegendrePolynomials, 
	Polynomials, ForwardDiff,Revise, LaTeXStrings, PrettyTables

# ╔═╡ 16039ccb-7cb1-4f44-be4a-03fdb0d08873
using Main.BandPyrometry #this line returns not defined error on the first Pluto run (probably, because of the Pluto running all "using"'s before the cells) just re-run this cell manually

# ╔═╡ 500f0cb0-bf70-11ee-19a8-35930e422cab
md"""
# Test of Planck.jl, EmPoint type from BandPyrometry.jl and JDXreader.jl
"""

# ╔═╡ 40e2112d-58cc-4650-8c6b-53b23a9fe470
md"""
##### Contents

* Planck functions
* Optimization of BB temperature with and without the EmPoint type 
* JDXreader test

"""

# ╔═╡ 77e9ee5b-b7f5-4216-8745-ca0286ab9024
includet(joinpath("../src","BandPyrometry.jl")) # the package file is includet using Revise methods, thus the code of the main file can be modified

# ╔═╡ 27ca1398-96ab-4d24-9c6d-63e02da49a94
import .BandPyrometry.Planck as PL

# ╔═╡ 6f5a8d02-d519-47c4-bdcd-7e52bf92098c
import Main.BandPyrometry as BP

# ╔═╡ 9e047346-b39e-4179-8f8a-bb1188941606
md"""
### Testing main functionality of Planck module functions

 Plotting blackbody spectrum and its derivarives
"""

# ╔═╡ 7098ee9a-1dc9-43fb-b554-35d87174fc64
md" Set wavelength region: "

# ╔═╡ f1664b05-3dfb-4aa9-8b81-1e3e1ed488c8
md"wavelength region left boundary,  λₗ = $(@bind lam_left confirm(Slider(0.1:0.1:60,default=0.1,show_value = true)))"

# ╔═╡ 323ed168-3ebb-4f9a-bdd9-2db3dace879c
md"wavelength region right boundary λᵣ = $(@bind lam_right confirm(Slider(0.1:0.1:60,default=lam_left+10,show_value = true)))"

# ╔═╡ 6beb18fe-4e7c-4851-ab94-b5e6336bb8e6
@bind  T_values PlutoUI.combine() do Child
	md"""
	T₁ = $(
		Child(Slider(-200.0:10:2500,default=300,show_value = true))
	) 
	T₂ = $(
		Child(Slider(-200.0:10:2500,default=900,show_value = true))
	)
	T₃ = $(
		Child(Slider(-200.0:10:2500,default=1500,show_value = true))
	)
	"""
end

# ╔═╡ 536278e5-8cb7-4787-9018-b862f0adfc07
T = [+(t,PL.Tₖ) for t in T_values]; # convert all temperatures to Kelvin

# ╔═╡ baf332fc-a805-4981-b569-70324a879a99
lam = range(lam_left,lam_right,length=100); # full region (as far as EmPoint uses StaticArrays, the length of lam range should not be too high)

# ╔═╡ b3d56970-6fe3-4834-89df-c931c07269c9
pretty_table(HTML,hcat(T,PL.∫ibbₗ.(T,λₗ = lam_left,λᵣ=lam_right)), header = ["Temperature", "Fractional power"], top_left_str="Table of fraction of total power within spectral range $(lam_left) ... $(lam_right)")

# ╔═╡ 3de53f57-ca6f-433c-b3b7-5a59250895e9
@bind p_func Select([("blackbody intensity","I", PL.ibb),("first derivative with λ","dI/dλ", PL.∇ₗibb), ("first derivative with T","dI/dT", PL.∇ₜibb),("second derivative with T","d²I/dT²", PL.∇²ₜibb),("second derivative  with λ","d²I/dλ²", PL.∇²ₗibb) ])

# ╔═╡ 04261201-2729-4e94-8244-d30ccc22a19d
begin 
	p=plot(lam,p_func[3].(lam,T[1]),
		title=p_func[1],
		legend=true,
		xlabel = "λ [μm]", 
		ylabel =p_func[2],
		label =@sprintf("T=%5.1f K",T[1]))
	for i in 2:3
		plot!(lam,p_func[3].(lam,T[i]),
		title=p_func[1],
		legend=true,
		label =@sprintf("T=%5.1f K",T[i]))
	end
	p
end

# ╔═╡ 3a16c8c8-17d9-4c36-9bc1-99a081a32c33
md"""
### BB temperature fitting using Planck module and EmPoint type from BandPyrometry module
_________

Least square discrepancy:\
``r(T)=\sum_{i}(I_{bb}(\lambda_i, T) -I_{bb}(\lambda_i, T_{real}))^2 ``

``I_{bb}(\lambda, T)``  is the blackbody spectrum
"""

# ╔═╡ 05c7a161-4346-4464-b41f-66ca7b3e149d
md" Blackbody real temperature ``T_{real}``  to be  fitted (Celsius)  $(@bind Ttr confirm(Slider(0.0:10:2500,default=500,show_value = true)))"

# ╔═╡ a0b68a76-8eb0-43af-927a-2ffc1abcea98
# "Measured" BB spectrum for temperature Ttr
Imeas = PL.ibb.(lam,Ttr);

# ╔═╡ a9214e4c-5467-4cca-81b3-059590b121b7
em_point = BP.EmPoint(Imeas, collect(lam));# emissivity point 

# ╔═╡ ae5738ce-9c29-42f3-8b9c-a465c8bf2952
# residual function constructor	
residual_func(T,l)=0.5*norm( Imeas .- PL.ibb.(l,T))^2;

# ╔═╡ c52a950c-c15e-4188-a5ce-1ba3100d43b9
# gradient function calculated directly using Planck module functions
function grad_fun!(G,T,l) 
	amat = repeat(l,1,3)
	PL.a₁₂₃!(amat,l,T[end])
	Ib = copy(l)
	PL.ibb!(Ib,l,amat)
	∇I = copy(l)
	PL.∇ₜibb!(∇I,T[end],amat,Ib)
	G[1] = -dot(∇I,(Imeas -Ib))
end;

# ╔═╡ f3013bb6-eab2-4256-b244-ac3cd272cc42
# Hessian direct implementation	
function hess_fun!(H,T,l) 
	(i1,i2,i3) = PL.Dₜibb(l,T)# returns planck,first and second derivatives
	H[1] = -dot(i3,Imeas .-i1) + dot(i2,i2);
end;

# ╔═╡ 9bc6c540-53c8-4a4d-9194-629d9e0277ca
md" Set the value of ``T_{try}`` to calculate ``r(T_{try}),\nabla_T r(T_{try}),\nabla_T^2 r(T_{try}) `` manually (Celsius) $@bind Ttry Slider(0.0:0.1:2500,default=500.0,show_value = true)"

# ╔═╡ fad07c97-14ee-4df8-9914-72309fa55d24
begin
	TtryK = Ttry+PL.Tₖ
	discr_values = [0.0,0.0,"-"]
	grad_values  = [0.0,0.0,0.0]
	hess_values = [0.0,0.0,0.0]
	lam_vec = collect(lam)
	discr_values[1] = residual_func(TtryK,lam_vec) # direct discrepancy calculation
	discr_values[2] =BP.disc([TtryK],em_point) # calculating the discrepancy inside  the EmPoint

	#calculating first order derivatives
	grad_fun!(view(grad_values,1),TtryK,lam_vec) # direct calculation of the gradient 
	G = [0.0]
	BP.grad!(G,[TtryK],em_point) # EmPoint calculation of the gradient
	grad_values[2] = G[]
	auto_fo_der(t) =  ForwardDiff.derivative(x->residual_func(x,lam_vec),t)
	grad_values[3] = auto_fo_der(TtryK)
	# calculating second order derivatives
	hess_values[3] = ForwardDiff.derivative(auto_fo_der, TtryK) #using autodiff
	hess_fun!(view(hess_values,1),[TtryK],lam_vec)
	BP.hess!(G,[TtryK],em_point) # EmPoint calculation of the second order derivative
	hess_values[2] =G[]
	pretty_table(HTML,hcat(["Direct","EmPoint", "AutoDiff"],discr_values,grad_values,hess_values),header = ["calculated using:","discrepancy", "gradient", "hessian"],top_left_str= "Table of discrepancy and its derivatives values, calculated for T=$(TtryK),K" )
end

# ╔═╡ 439c48f9-69e7-46bc-8c42-f53ad5456772
begin 
	benchmark_measures = Dict{String,NamedTuple{(:t_direct,:t_emPoint),Tuple{Float64, Float64}}}();# dictionary stores results of benchmark tests
	bm_res = Base.RefValue{BenchmarkTools.Trial}();# stores current benchmark test results
	t1= Ref(0.0)
	t2 = Ref(0.0)
end;# rerun this cel 

# ╔═╡ e7d5617b-1da9-4b96-bd7a-0490e8e4ceca
md"""
#### Try zero order method fitting = $(@bind is_zero_order_run  confirm(CheckBox()))
	"""

# ╔═╡ c7caa7f6-fae1-4e07-aae8-a11312d74f09
md"""Zero order optimizer $(@bind zero_order_method Select([("NelderMead",NelderMead),("ParticleSwarm",ParticleSwarm)]))
"""

# ╔═╡ e04ad4d7-c26a-4753-baa8-921121c27f08
if is_zero_order_run
	begin 
		#optimization problem in terms of primitive implementation
		prob_zo= OptimizationProblem(OptimizationFunction(residual_func,Optimization.AutoForwardDiff()), [235.0],collect(lam),lb=[20.0],ub=[2500.0]); 
		sol_zo = solve(prob_zo,  zero_order_method[2]());#solution
		#optimization problem solution using BandPyrometry.EmPoint object 
		em_point_zo = BP.EmPoint(copy(Imeas), collect(lam));# emissivity point 
		em_prob_zo= OptimizationProblem(OptimizationFunction(BP.disc,grad=BP.grad!), [235.0],em_point_zo);# optimization problem with Empoint implementation
		em_sol_zo = solve(em_prob_zo,  zero_order_method[2]()) # solution using 
	end
end;

# ╔═╡ dd0b1fbb-c56f-4ac8-8ad8-57d6add2ddcc
is_zero_order_run ? md"Tres direct = $( sol_zo.u[])" : nothing

# ╔═╡ ccdce2ef-cbc9-4339-ab8b-7eaab5479c97
is_zero_order_run ? md"Tres EmPoint =  $( em_sol_zo.u[])" : nothing

# ╔═╡ 88dd837e-bb77-4ec4-87d8-1134b0f741ce
md"show zero order output = $(@bind is_zero_order_out  CheckBox())"

# ╔═╡ ede76282-d27a-41a2-ac92-101c26129346
if is_zero_order_run && is_zero_order_out
	sol_zo.original
end

# ╔═╡ 6cc50dda-bb5c-4138-a520-973e8dae3df0
if is_zero_order_run && is_zero_order_out
	em_sol_zo.original
end

# ╔═╡ 8c13923d-612b-44d5-8046-bf4aa4bc175e
if is_zero_order_run&&is_zero_order_out
	begin
		plot(lam,Imeas, label="Measured spectrum Treal=$(Ttr)")
		plot!(lam, PL.ibb(lam, sol_zo.u[]),label="Direct solution T=$(sol_zo.u[])")
		scatter!(lam, PL.ibb(lam, em_sol_zo.u[]),label="EmPoint solution T=$(em_sol_zo.u[])")
		xlabel!("Wavelength, μm")
		ylabel!("Blackbody intensity, a.u.")
	end
end

# ╔═╡ 1bbc82eb-16e3-47a2-814a-0a69bd803c75
md"run zero order benchmark = $(@bind is_zero_order_bench  CheckBox())"

# ╔═╡ 7a4e3653-ee72-404c-a849-0e4c389b65f5
if is_zero_order_run && is_zero_order_bench
	begin 
		bm_res[] = @benchmark solve(prob_zo,  zero_order_method[2]())
		t1[]=mean(bm_res[].times)
		bm_res[]
	end
end

# ╔═╡ 8b70cb0a-ce22-4b11-bb75-4240b2b2b36b
if is_zero_order_run && is_zero_order_bench
	begin 
		bm_res[] = @benchmark solve(em_prob_zo,  zero_order_method[2]())
		t2[]=mean(bm_res[].times)
		bm_res[]
	end
end

# ╔═╡ f80f2317-c3b1-4d06-b030-fcbb565ddd82
if is_zero_order_run && is_zero_order_bench
	benchmark_measures[zero_order_method[1]]=(t_direct=t1[],t_emPoint=t2[])
end;

# ╔═╡ 59af65c0-1520-4a8e-9f51-69c1446a723b
md"""
	#### Try first order method fitting = $(confirm(@bind is_first_order_run  CheckBox()))
	"""

# ╔═╡ e498d53f-ade6-4085-9b3f-5468e8e99721
@bind first_order_method Select([("GradientDescent",GradientDescent),("BFGS",BFGS),("LBFGS",LBFGS)])

# ╔═╡ 2ea9247a-e667-447d-8ed1-6100f76876a5
if is_first_order_run
	begin 
		prob_fo= OptimizationProblem(OptimizationFunction(residual_func,grad=grad_fun!), [235.0],collect(lam)); 
		sol_fo = solve(prob_fo,  first_order_method[2]());
		# optimization problem with Empoint implementation
		em_point_fo = BP.EmPoint(copy(Imeas), collect(lam));# emissivity point 
		em_prob_fo= OptimizationProblem(OptimizationFunction(BP.disc,
			grad=BP.grad!), [835.0],em_point_fo);
		em_sol_fo = solve(em_prob_fo,  first_order_method[2]()); # solution using 
	end;
end

# ╔═╡ 6974f06e-e0a9-48aa-ac4b-3c2b55b13734
is_first_order_run ? md"Tres direct implementation= $( sol_fo.u[])" : nothing

# ╔═╡ ab8d29d4-e221-4223-b5a8-dda4e2a32c0f
is_first_order_run ? md"Tres EmPoint= $( em_sol_fo.u[])" : nothing

# ╔═╡ 608f923e-c327-4252-a245-49ada28cc874
md"show first order output = $(@bind is_first_order_out  CheckBox())"

# ╔═╡ e691fe83-64e9-4b85-a32b-93dc0c669376
is_first_order_run&&is_first_order_out ? sol_fo.original : nothing

# ╔═╡ 6e8ec563-8797-4a40-89bb-c52a3d7802b8
is_first_order_run&&is_first_order_out ? em_sol_fo.original : nothing

# ╔═╡ 5730698e-7376-4723-829d-f2cc602d072a
if is_first_order_run&&is_first_order_out
	begin
		p1 = plot(lam,Imeas, label="Measured spectrum Treal=$(Ttr)")
		title!("First order method $(first_order_method[1])")
		plot!(lam, PL.ibb(lam, sol_fo.u[]),label="Direct solution T=$(sol_fo.u[])")
		scatter!(lam, PL.ibb(lam, em_sol_fo.u[]),label="EmPoint solution T=$(em_sol_fo.u[])")
		xlabel!("Wavelength, μm")
		ylabel!("Blackbody intensity, a.u.")
	end
end

# ╔═╡ c34c44bf-1485-4f7d-8715-16fadc7b79bc
md"run first order benchmark = $(@bind is_first_order_bench  CheckBox())"

# ╔═╡ 21137cef-18d9-42cb-90d7-f565a1b243eb
if is_first_order_run&&is_first_order_bench
	begin 
		bm_res[] = @benchmark solve(prob_fo,  first_order_method[2]())
		t1[]=mean(bm_res[].times)
		bm_res[]
	end
end

# ╔═╡ d76e765f-ff56-4c7e-b61e-0b49e88f558a
if is_first_order_run&&is_first_order_bench
	begin 
		bm_res[] = @benchmark solve(em_prob_fo,  first_order_method[2]())
		t2[]=mean(bm_res[].times)
		bm_res[]
	end
end

# ╔═╡ 7b32578c-de60-47fa-ae3f-628bb1c887ba
if is_first_order_run&&is_first_order_bench
	benchmark_measures[first_order_method[1]]=(t_direct=t1[],t_emPoint=t2[])
end;

# ╔═╡ 95d1107f-07c5-4c31-904a-1471cb51ad33
md"""
	#### Try second order method fitting = $(@bind is_second_order_run  confirm(CheckBox()))
		"""

# ╔═╡ 631a1822-5bf0-4713-b048-a17b93bfa66e
@bind second_order_method Select([("Newton",Newton), ("NewtonTrustRegion",NewtonTrustRegion)])

# ╔═╡ ae12d403-1149-426e-ab13-5afc5e5615d9
if is_second_order_run
	begin 
		prob_so= OptimizationProblem(OptimizationFunction(residual_func,
			grad=grad_fun!,hess=hess_fun!), [235.0],collect(lam)); #optimization problem in terms of primitive implementation ,lb=[20.0],ub=[2500.0]
		sol_so = solve(prob_so,  second_order_method[2]());
		# optimization problem with Empoint implementation
		em_point_so = BP.EmPoint(copy(Imeas), collect(lam));# emissivity point 
		em_prob_so= OptimizationProblem(OptimizationFunction(BP.disc,
			grad=BP.grad!, hess=BP.hess!), [835.0],em_point_so);
		em_sol_so = solve(em_prob_so,  second_order_method[2]()) # solution using 
	end;
end

# ╔═╡ c2071134-3355-4479-a2d7-55d2141bb7ff
is_second_order_run ? md"Tres direct implementation= $( sol_so.u[])" : nothing

# ╔═╡ d0afa1de-101b-4eea-b4c3-75d2b4bc27d9
is_second_order_run ? md"Tres EmPoint= $( em_sol_so.u[])" : nothing

# ╔═╡ 70f59af9-f609-4b9b-b23e-606bc3b2efa2
md"show second order output = $(@bind is_second_order_out  CheckBox())"

# ╔═╡ 51988df7-5be5-478f-8cda-b4999b32ff6f
is_second_order_run&&is_second_order_out ? sol_so.original : nothing

# ╔═╡ 207929ea-af4a-4d4c-b2d2-59a0dc927ba6
is_second_order_run&&is_second_order_out ? em_sol_so.original : nothing

# ╔═╡ e9ea4596-7d7b-4e87-87d2-88fb40a1b6ad
if is_second_order_run&&is_second_order_out
	begin
		p2 = plot(lam,Imeas, label="Measured spectrum Treal=$(Ttr)")
		title!("First order method $(second_order_method[1])")
		plot!(lam, PL.ibb(lam, sol_so.u[]),label="Direct solution T=$(sol_so.u[])")
		scatter!(lam, PL.ibb(lam, em_sol_so.u[]),label="EmPoint solution T=$(em_sol_so.u[])")
		xlabel!("Wavelength, μm")
		ylabel!("Blackbody intensity, a.u.")
	end
end

# ╔═╡ 52894f69-76ab-4045-b70c-67bea0043279
md"run second order benchmark = $(@bind is_second_order_bench  CheckBox())"

# ╔═╡ 70476796-ae99-490b-a77f-7b609b9e5b5f
if is_second_order_bench
	begin 
		bm_res[] = @benchmark solve(prob_so,  second_order_method[2]())
		t1[]=mean(bm_res[].times)
		bm_res[]
	end
end

# ╔═╡ 7e8ef8b5-630b-4150-aa46-d58537906103
if is_second_order_run&&is_second_order_bench
	begin 
		bm_res[] = @benchmark solve(em_prob_so,  second_order_method[2]())
		t2[]=mean(bm_res[].times)
		bm_res[]
	end
end

# ╔═╡ 3cd039d3-3926-4889-a743-a8b908bf1796
if is_second_order_run&&is_second_order_bench
	benchmark_measures[second_order_method[1]]=(t_direct=t1[],t_emPoint=t2[])
end;

# ╔═╡ 8cfa738e-05cc-4d86-b40c-86442d14b4b1
if is_zero_order_bench || is_first_order_bench || is_second_order_bench
	t_dir_vec = [ b[2].t_direct/1e6 for b in benchmark_measures]
	t_em_vec = [ b[2].t_emPoint/1e6 for b in benchmark_measures]
	optims_vec = collect(keys(benchmark_measures))
	pretty_table(HTML,hcat(optims_vec,t_dir_vec,t_em_vec),header = ["Optimizer", "time Direct, ms", "time EmPoint, ms"])
end

# ╔═╡ c200b8df-3580-4cb5-9da4-bf5e2113bccb
md"Select benchmarks to compare"

# ╔═╡ f9d71607-f558-4f90-b5a0-b4c445c97f2e
@bind plot_data_selector  confirm(MultiSelect(collect(keys(benchmark_measures))))

# ╔═╡ a5565a17-6572-4946-8147-7a9c7ca203f2
md"""
	Load jdx file? $(@bind is_jdx_loaded confirm(CheckBox()))
	"""

# ╔═╡ 92b7b55f-b50a-48cf-8d6e-0a95ae7b2465
md"This figure show measured spectrum loaded using JDXreader"

# ╔═╡ 3f4d5d0e-8f35-4f57-8b0f-2aae19b2f7be
if is_jdx_loaded
begin
	filedir = joinpath(pwd(),"JCAM_test_file.txt")
	if isfile(filedir)
		xydata = BandPyrometry.JDXreader.read_jdx_file(filedir);
		pj = plot(xydata.x,xydata.y,label = nothing)
		xlabel!(xydata.headers["XUNITS"])
		ylabel!(xydata.headers["YUNITS"])
		title!(xydata.headers["DATATYPE"])
		pj
	end
end;
end

# ╔═╡ 21d63be3-b945-452b-91d8-63efb3568b95
function plot_bencnhmark(benchmark_measures)
	performed_tests = collect(keys(benchmark_measures))
	N =  2*length(performed_tests)
	if N==0
		return nothing
	end
	r_c = RGB(1,0,0)
	b_c = RGB(0.2, 0.2, 1)
	bars_coords =collect(1.0:N)
	bars_coords[1:2:end] .=bars_coords[1:2:end] .+0.1
	bars_coords[2:2:end] .=bars_coords[2:2:end] .-0.1
	bars_values = copy(bars_coords)
	bars_values[1:2:end] .= [x[2].t_direct/1e6 for x in benchmark_measures]
	bars_values[2:2:end] .= [x[2].t_emPoint/1e6 for x in benchmark_measures]
	colrs = [isodd(i) ? b_c : r_c  for i in 1:N]
	ppp = plot(bar(bars_coords, bars_values, 
		label = "",
		alpha = 0.5, 
		bar_width = 0.8,fillcolor = colrs),		
		#background_color = :transparent,
		background_color = :white,
		#tickfontcolor= RGB(1, 1, 1),
		tickfontsize=10,
		xticks=(1:2:(2*length(performed_tests)),performed_tests),
		xrotation = -10
	)
	ylims!(ppp,0.9*minimum(bars_values),1.1*maximum(bars_values))
	title!("""Performance of direct impelementation
				vs EmPoint type""")
	temp = [-10 -20; -30 -40 ]
	plot!(temp, color=[b_c r_c], label=["Direct" "EmPoint"],linewidth=25)
	ylabel!("Soluion time, ms")
	return ppp
end;

# ╔═╡ 6f7497ac-d157-4ed5-8ed8-2a80f267efad
begin 
	selected_benchs = Dict(zip(plot_data_selector,[benchmark_measures[k] for k in plot_data_selector])) 
	plot_bencnhmark(selected_benchs)
end

# ╔═╡ Cell order:
# ╟─500f0cb0-bf70-11ee-19a8-35930e422cab
# ╟─40e2112d-58cc-4650-8c6b-53b23a9fe470
# ╠═fc6d97d3-3775-42d2-a991-09d875b0dd38
# ╠═c6f4b5c0-e007-4d4f-8b79-b3380cc5e10c
# ╟─77e9ee5b-b7f5-4216-8745-ca0286ab9024
# ╟─16039ccb-7cb1-4f44-be4a-03fdb0d08873
# ╠═27ca1398-96ab-4d24-9c6d-63e02da49a94
# ╠═6f5a8d02-d519-47c4-bdcd-7e52bf92098c
# ╟─9e047346-b39e-4179-8f8a-bb1188941606
# ╟─7098ee9a-1dc9-43fb-b554-35d87174fc64
# ╟─f1664b05-3dfb-4aa9-8b81-1e3e1ed488c8
# ╟─323ed168-3ebb-4f9a-bdd9-2db3dace879c
# ╟─6beb18fe-4e7c-4851-ab94-b5e6336bb8e6
# ╟─536278e5-8cb7-4787-9018-b862f0adfc07
# ╟─baf332fc-a805-4981-b569-70324a879a99
# ╟─b3d56970-6fe3-4834-89df-c931c07269c9
# ╟─3de53f57-ca6f-433c-b3b7-5a59250895e9
# ╟─04261201-2729-4e94-8244-d30ccc22a19d
# ╟─3a16c8c8-17d9-4c36-9bc1-99a081a32c33
# ╟─05c7a161-4346-4464-b41f-66ca7b3e149d
# ╟─a0b68a76-8eb0-43af-927a-2ffc1abcea98
# ╟─a9214e4c-5467-4cca-81b3-059590b121b7
# ╟─ae5738ce-9c29-42f3-8b9c-a465c8bf2952
# ╟─c52a950c-c15e-4188-a5ce-1ba3100d43b9
# ╟─f3013bb6-eab2-4256-b244-ac3cd272cc42
# ╟─9bc6c540-53c8-4a4d-9194-629d9e0277ca
# ╟─fad07c97-14ee-4df8-9914-72309fa55d24
# ╠═439c48f9-69e7-46bc-8c42-f53ad5456772
# ╟─e7d5617b-1da9-4b96-bd7a-0490e8e4ceca
# ╟─c7caa7f6-fae1-4e07-aae8-a11312d74f09
# ╟─e04ad4d7-c26a-4753-baa8-921121c27f08
# ╟─dd0b1fbb-c56f-4ac8-8ad8-57d6add2ddcc
# ╟─ccdce2ef-cbc9-4339-ab8b-7eaab5479c97
# ╟─88dd837e-bb77-4ec4-87d8-1134b0f741ce
# ╟─ede76282-d27a-41a2-ac92-101c26129346
# ╟─6cc50dda-bb5c-4138-a520-973e8dae3df0
# ╟─8c13923d-612b-44d5-8046-bf4aa4bc175e
# ╟─1bbc82eb-16e3-47a2-814a-0a69bd803c75
# ╟─7a4e3653-ee72-404c-a849-0e4c389b65f5
# ╟─8b70cb0a-ce22-4b11-bb75-4240b2b2b36b
# ╟─f80f2317-c3b1-4d06-b030-fcbb565ddd82
# ╟─59af65c0-1520-4a8e-9f51-69c1446a723b
# ╟─e498d53f-ade6-4085-9b3f-5468e8e99721
# ╟─2ea9247a-e667-447d-8ed1-6100f76876a5
# ╟─6974f06e-e0a9-48aa-ac4b-3c2b55b13734
# ╟─ab8d29d4-e221-4223-b5a8-dda4e2a32c0f
# ╟─608f923e-c327-4252-a245-49ada28cc874
# ╟─e691fe83-64e9-4b85-a32b-93dc0c669376
# ╟─6e8ec563-8797-4a40-89bb-c52a3d7802b8
# ╟─5730698e-7376-4723-829d-f2cc602d072a
# ╟─c34c44bf-1485-4f7d-8715-16fadc7b79bc
# ╟─21137cef-18d9-42cb-90d7-f565a1b243eb
# ╟─d76e765f-ff56-4c7e-b61e-0b49e88f558a
# ╟─7b32578c-de60-47fa-ae3f-628bb1c887ba
# ╟─95d1107f-07c5-4c31-904a-1471cb51ad33
# ╟─631a1822-5bf0-4713-b048-a17b93bfa66e
# ╟─ae12d403-1149-426e-ab13-5afc5e5615d9
# ╟─c2071134-3355-4479-a2d7-55d2141bb7ff
# ╟─d0afa1de-101b-4eea-b4c3-75d2b4bc27d9
# ╟─70f59af9-f609-4b9b-b23e-606bc3b2efa2
# ╟─51988df7-5be5-478f-8cda-b4999b32ff6f
# ╟─207929ea-af4a-4d4c-b2d2-59a0dc927ba6
# ╟─e9ea4596-7d7b-4e87-87d2-88fb40a1b6ad
# ╟─52894f69-76ab-4045-b70c-67bea0043279
# ╟─70476796-ae99-490b-a77f-7b609b9e5b5f
# ╟─7e8ef8b5-630b-4150-aa46-d58537906103
# ╟─3cd039d3-3926-4889-a743-a8b908bf1796
# ╟─8cfa738e-05cc-4d86-b40c-86442d14b4b1
# ╟─c200b8df-3580-4cb5-9da4-bf5e2113bccb
# ╠═f9d71607-f558-4f90-b5a0-b4c445c97f2e
# ╟─6f7497ac-d157-4ed5-8ed8-2a80f267efad
# ╟─a5565a17-6572-4946-8147-7a9c7ca203f2
# ╟─92b7b55f-b50a-48cf-8d6e-0a95ae7b2465
# ╠═3f4d5d0e-8f35-4f57-8b0f-2aae19b2f7be
# ╟─21d63be3-b945-452b-91d8-63efb3568b95
