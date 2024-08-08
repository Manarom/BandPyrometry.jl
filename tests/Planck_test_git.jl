### A Pluto.jl notebook ###
# v0.19.43

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ fc6d97d3-3775-42d2-a991-09d875b0dd38
using StaticArrays, 
PlutoUI, Plots, 
Optimization,
OptimizationOptimJL,BenchmarkTools, 
LinearAlgebra, Interpolations,Printf,
MKL,LinearAlgebra, LegendrePolynomials, 
Polynomials, ForwardDiff, LaTeXStrings

# ╔═╡ 500f0cb0-bf70-11ee-19a8-35930e422cab
md"""
# Test of Planck.jl, EmPoint type from BandPyrometry.jl and JDXreader.jl
"""

# ╔═╡ 40e2112d-58cc-4650-8c6b-53b23a9fe470
md"""
##### Contents

* Planck functions
* Optimization of BB temperature with and without EmPoint type 
* JDXreader test

	"""
# ╔═╡ aba2d996-2cea-4484-81d7-b33d56d525e4
# ENTER HERE PATH TO THE FOLDER WITH BandPyrometry.jl, Planck.jl and JDXreader.jl 
cur_dir = Ref("working_dir")

# ╔═╡ 3930e657-b34d-4d4b-b7a2-596118f4ce74
is_files_imported_flag = Ref(false)

# ╔═╡ 77e9ee5b-b7f5-4216-8745-ca0286ab9024
needed_files = ["BandPyrometry.jl", "JDXreader.jl"];

# ╔═╡ 83663678-22c4-4b88-8513-3487b12f9736
files_needed = Dict(zip(needed_files,["",""]));

# ╔═╡ 0abbeb1a-062d-43fb-b622-91ad63c2719a
# trying to find necessary files
if !is_files_imported_flag[]
	begin 
		for (k,val) ∈ files_needed
				global cur_dir,files_needed
				files_needed[k]= joinpath(cur_dir[],k)
		end
	
		cd(cur_dir[])
		for k in needed_files
			include(files_needed[k])
		end
		is_files_imported_flag[]=true
	end;
end

# ╔═╡ 6ddac362-fc91-4881-8f29-8a835870b781
#names(Main)

# ╔═╡ 27ca1398-96ab-4d24-9c6d-63e02da49a94
import Main.Planck as PL

# ╔═╡ 6f5a8d02-d519-47c4-bdcd-7e52bf92098c
import Main.BandPyrometry as BP

# ╔═╡ 9e047346-b39e-4179-8f8a-bb1188941606
md"""
## Plotting blackbody spectrum and its derivarives
"""

# ╔═╡ 7098ee9a-1dc9-43fb-b554-35d87174fc64
md" Set wavelength region "

# ╔═╡ f1664b05-3dfb-4aa9-8b81-1e3e1ed488c8
md"wavelength region left boundary,  λₗ = $(confirm(@bind lam_left Slider(0.1:0.1:60,default=0.1,show_value = true)))"

# ╔═╡ 323ed168-3ebb-4f9a-bdd9-2db3dace879c
md"wavelength region right boundary λᵣ = $(confirm(@bind lam_right Slider(0.1:0.1:60,default=lam_left+10,show_value = true)))"

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
T = map((X)->+(X,PL.Tₖ),collect(T_values) ); # convert all temperatures to Kelvin

# ╔═╡ baf332fc-a805-4981-b569-70324a879a99
lam = range(lam_left,lam_right,length=100); # full region (as far as EmPoint uses StaticArrays, the length of lam range should not be too high)

# ╔═╡ 3de53f57-ca6f-433c-b3b7-5a59250895e9
@bind p_func Select([("blackbody intensity","I", PL.ibb),("blackbody intensity first derivative by wavelength","dI/dλ", PL.∇ₗibb), ("blackbody intensity first derivative by  temperature","dI/dT", PL.∇ₜibb),("blackbody intensity second derivative by temperature","d²I/dT²", PL.∇²ₜibb),("blackbody intensity second derivative  by wavelength","d²I/dλ²", PL.∇²ₗibb) ])

# ╔═╡ 0713ba68-9c5d-464a-ad38-ea05d5239511
md" fraction of total power within spectral range $(lam_left) ... $(lam_right):"

# ╔═╡ d8924882-d84c-417b-a4ec-9be6a32d992a
@. "for T="*string(T)*" : "*string(PL.∫ibbₗ.(T,λₗ = lam_left,λᵣ=lam_right))

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
Least square discrepancy:\
r(T)= ∑(I(λᵢ,T)-I(λᵢ,Tᵣₑₐₗ))²\
I(λᵢ,Tᵣₑₐₗ) - BB spectrum
"""

# ╔═╡ 05c7a161-4346-4464-b41f-66ca7b3e149d
md" BB temperature to  fit (Celsius)  $(confirm(@bind Ttrue Slider(0.0:10:2500,default=500,show_value = true)))"

# ╔═╡ f1968bb1-dc06-47a9-939f-4b76215c9bfe
Ttr=Ttrue+PL.Tₖ; # converting Celsius to Kelvins

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
	G[1] = -dot(∇I,(Imeas -Ib));
end;

# ╔═╡ 9bc6c540-53c8-4a4d-9194-629d9e0277ca
md" Ttry $@bind Ttry Slider(0.0:0.1:2500,default=500.0,show_value = true)"

# ╔═╡ 0f93061c-5c62-427d-b1fb-0f9b8066a786
" discrepancy value direct calculation $(residual_func(Ttry+PL.Tₖ,lam))"

# ╔═╡ ec4494d0-61e4-4afa-a135-5a6a3e71ba84
" discrepancy value em_point $( BP.disc([Ttry+PL.Tₖ],em_point) ) "

# ╔═╡ eb4e98ae-8891-4d54-9581-ffd54e925e17
begin 
	G = [0.0]
	grad_fun!(G,Ttry+PL.Tₖ,collect(lam)) 
	md" direct gradient calculation $(G[end])"
end

# ╔═╡ 28429534-2221-44ce-8e75-0c6146e137e6
function grad_fun2!(G,T,l) 
	(i1,i2) = PL.Dₜibb(l,T); #∇ₜibb!(y,λ::Vector,T::Vector)
	#G[1] = dot(2*PL.∇ₜibb.(l,T),PL.ibb.(l,T) - PL.ibb.(l,Ttr));
	G[1] = dot(i2,(i1 - PL.ibb.(l,Ttr)))
end;

# ╔═╡ 66263f0a-3921-478e-813b-29c33f3d557e
begin 
	BP.grad!(G,Ttry+PL.Tₖ,em_point)
	md" gradient calculation using EmPoint $(G[end])"
	
end

# ╔═╡ f3879d0d-45b0-4820-8153-19859d235ef1
begin 
	grad_fun2!(G,[Ttry+PL.Tₖ],collect(lam)) 
	md" direct gradient calculation 2nd version $(G[end])"
end

# ╔═╡ 2ba39927-635d-4054-ad28-cff311199263
begin 
	ForwardDiff.derivative(x->residual_func(x,collect(lam)),Ttry+PL.Tₖ)
end

# ╔═╡ 439c48f9-69e7-46bc-8c42-f53ad5456772
begin 
	benchmark_measures = Dict{String,NamedTuple{(:t_direct,:t_emPoint),Tuple{Float64, Float64}}}();# dictionary stores results of benchmark tests
	bm_res = Base.RefValue{BenchmarkTools.Trial}();# stores current benchmark test results
	t1= Ref(0.0)
	t2 = Ref(0.0)
end;

# ╔═╡ e7d5617b-1da9-4b96-bd7a-0490e8e4ceca
md"""
#### Try zero order method fitting = $(confirm(@bind is_zero_order_run  CheckBox()))
	"""

# ╔═╡ c7caa7f6-fae1-4e07-aae8-a11312d74f09
@bind zero_order_method Select([("NelderMead",NelderMead),("ParticleSwarm",ParticleSwarm)])

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
		p0 = plot(lam,Imeas, label="Measured spectrum Treal=$(Ttr)")
		plot!(lam, PL.ibb(lam, sol_zo.u[]),label="Direct solution T=$(sol_zo.u[])")
		plot!(lam, PL.ibb(lam, em_sol_zo.u[]),label="EmPoint solution T=$(em_sol_zo.u[])")
		p0
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
end

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
		plot!(lam, PL.ibb(lam, em_sol_fo.u[]),label="EmPoint solution T=$(em_sol_fo.u[])")
		p1
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
end

# ╔═╡ 95d1107f-07c5-4c31-904a-1471cb51ad33
md"""
	#### Try second order method fitting = $(confirm(@bind is_second_order_run  CheckBox()))
		"""

# ╔═╡ 631a1822-5bf0-4713-b048-a17b93bfa66e
@bind second_order_method Select([("Newton",Newton), ("NewtonTrustRegion",NewtonTrustRegion)])

# ╔═╡ 773f5479-a7de-4a81-b513-0545913b37bd
# Hessian direct implementation	
function hess_fun!(H,T,l) 
	(i1,i2,i3) = PL.Dₜibb(l,T);
	i1 -=PL.ibb.(l,Ttr);
	H[1] = dot(2*i3,i1) + dot(2*i2,i2); 
end

# ╔═╡ ae12d403-1149-426e-ab13-5afc5e5615d9
if is_second_order_run
	begin 
		prob_so= OptimizationProblem(OptimizationFunction(residual_func,grad=grad_fun!,hess=hess_fun!), [235.0],collect(lam)); #optimization problem in terms of primitive implementation ,lb=[20.0],ub=[2500.0]
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
		plot!(lam, PL.ibb(lam, em_sol_so.u[]),label="EmPoint solution T=$(em_sol_so.u[])")
		p2
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
end

# ╔═╡ c758524b-04f3-4699-8594-cb65fe30f7a4
#L"\begin{bmatrix} Ttr=%$(Ttr) & %$(em_point.Tib[]) \\ %$(2) & %$(3) \end{bmatrix}"

# ╔═╡ 8cfa738e-05cc-4d86-b40c-86442d14b4b1
if is_zero_order_bench || is_first_order_bench || is_second_order_bench
	benchmark_measures
end

# ╔═╡ a5565a17-6572-4946-8147-7a9c7ca203f2
md"""
	Load jdx file? $(@bind is_jdx_loaded CheckBox())
	"""

# ╔═╡ 3f4d5d0e-8f35-4f57-8b0f-2aae19b2f7be
if is_jdx_loaded
begin 
	filedir = open_dialog("select file");
	if isfile(filedir)
		xydata = Main.JDXreader.read_jdx_file(filedir);
		pj = plot(xydata.x,xydata.y)
		xlabel!(xydata.headers["XUNITS"])
		ylabel!(xydata.headers["YUNITS"])
		title!(xydata.headers["DATATYPE"])
		pj
	end
end;
end

# ╔═╡ 21d63be3-b945-452b-91d8-63efb3568b95
function plot_becnhmark()
	performed_tests = collect(keys(benchmark_measures))
	if isempty(benchmark_measures)
		return nothing
	end
	direct_bars_values = [x[2].t_direct/1e6 for x in benchmark_measures]
	em_bars_values = [x[2].t_emPoint/1e6 for x in benchmark_measures]
	direct_bars_coords=1:3:3*length(performed_tests)
	em_bars_coords=2:3:3*length(performed_tests)
	p = plot(bar(direct_bars_coords, direct_bars_values, 
		label = "",
		alpha = 0.5, 
		bar_width = 1),
		bar(em_bars_coords, em_bars_values, 
		fillcolor=RGB(1,0,0),
		alpha = 0.5, 
		label = "",
		bar_width = 1),
		background_color = RGB(0.2, 0.2, 0.2),
		tickfontcolor= RGB(1, 1, 1),
		tickfontsize=15,
		ticks=performed_tests
		)
	ylims!(p,0.9minimum(em_bars_values),1.1maximum(direct_bars_values))
	return p
end

# ╔═╡ 725efcfe-6514-42d9-a92a-5b8cdbf8caa9
ppp = plot_becnhmark()

# ╔═╡ Cell order:
# ╟─500f0cb0-bf70-11ee-19a8-35930e422cab
# ╠═40e2112d-58cc-4650-8c6b-53b23a9fe470
# ╠═fc6d97d3-3775-42d2-a991-09d875b0dd38
# ╠═d8d114ae-2b48-4332-a656-7d5d6c2eebc2
# ╠═f28df1b7-8054-40f5-9254-67e4d15de96a
# ╠═139370dc-04aa-4687-97ee-8bd9c956d967
# ╠═aba2d996-2cea-4484-81d7-b33d56d525e4
# ╠═3930e657-b34d-4d4b-b7a2-596118f4ce74
# ╠═77e9ee5b-b7f5-4216-8745-ca0286ab9024
# ╠═83663678-22c4-4b88-8513-3487b12f9736
# ╠═0abbeb1a-062d-43fb-b622-91ad63c2719a
# ╟─6ddac362-fc91-4881-8f29-8a835870b781
# ╟─27ca1398-96ab-4d24-9c6d-63e02da49a94
# ╟─6f5a8d02-d519-47c4-bdcd-7e52bf92098c
# ╟─9e047346-b39e-4179-8f8a-bb1188941606
# ╟─7098ee9a-1dc9-43fb-b554-35d87174fc64
# ╟─f1664b05-3dfb-4aa9-8b81-1e3e1ed488c8
# ╟─323ed168-3ebb-4f9a-bdd9-2db3dace879c
# ╟─6beb18fe-4e7c-4851-ab94-b5e6336bb8e6
# ╟─536278e5-8cb7-4787-9018-b862f0adfc07
# ╟─baf332fc-a805-4981-b569-70324a879a99
# ╟─3de53f57-ca6f-433c-b3b7-5a59250895e9
# ╟─0713ba68-9c5d-464a-ad38-ea05d5239511
# ╟─d8924882-d84c-417b-a4ec-9be6a32d992a
# ╟─04261201-2729-4e94-8244-d30ccc22a19d
# ╟─3a16c8c8-17d9-4c36-9bc1-99a081a32c33
# ╟─05c7a161-4346-4464-b41f-66ca7b3e149d
# ╟─f1968bb1-dc06-47a9-939f-4b76215c9bfe
# ╟─a0b68a76-8eb0-43af-927a-2ffc1abcea98
# ╟─a9214e4c-5467-4cca-81b3-059590b121b7
# ╟─ae5738ce-9c29-42f3-8b9c-a465c8bf2952
# ╠═c52a950c-c15e-4188-a5ce-1ba3100d43b9
# ╠═9bc6c540-53c8-4a4d-9194-629d9e0277ca
# ╠═0f93061c-5c62-427d-b1fb-0f9b8066a786
# ╠═ec4494d0-61e4-4afa-a135-5a6a3e71ba84
# ╟─eb4e98ae-8891-4d54-9581-ffd54e925e17
# ╠═28429534-2221-44ce-8e75-0c6146e137e6
# ╠═66263f0a-3921-478e-813b-29c33f3d557e
# ╠═f3879d0d-45b0-4820-8153-19859d235ef1
# ╠═2ba39927-635d-4054-ad28-cff311199263
# ╠═439c48f9-69e7-46bc-8c42-f53ad5456772
# ╠═e7d5617b-1da9-4b96-bd7a-0490e8e4ceca
# ╠═c7caa7f6-fae1-4e07-aae8-a11312d74f09
# ╠═e04ad4d7-c26a-4753-baa8-921121c27f08
# ╠═dd0b1fbb-c56f-4ac8-8ad8-57d6add2ddcc
# ╠═ccdce2ef-cbc9-4339-ab8b-7eaab5479c97
# ╟─88dd837e-bb77-4ec4-87d8-1134b0f741ce
# ╠═ede76282-d27a-41a2-ac92-101c26129346
# ╠═6cc50dda-bb5c-4138-a520-973e8dae3df0
# ╠═8c13923d-612b-44d5-8046-bf4aa4bc175e
# ╠═1bbc82eb-16e3-47a2-814a-0a69bd803c75
# ╠═7a4e3653-ee72-404c-a849-0e4c389b65f5
# ╠═8b70cb0a-ce22-4b11-bb75-4240b2b2b36b
# ╠═f80f2317-c3b1-4d06-b030-fcbb565ddd82
# ╠═59af65c0-1520-4a8e-9f51-69c1446a723b
# ╠═e498d53f-ade6-4085-9b3f-5468e8e99721
# ╠═2ea9247a-e667-447d-8ed1-6100f76876a5
# ╠═6974f06e-e0a9-48aa-ac4b-3c2b55b13734
# ╠═ab8d29d4-e221-4223-b5a8-dda4e2a32c0f
# ╠═608f923e-c327-4252-a245-49ada28cc874
# ╠═e691fe83-64e9-4b85-a32b-93dc0c669376
# ╠═6e8ec563-8797-4a40-89bb-c52a3d7802b8
# ╠═5730698e-7376-4723-829d-f2cc602d072a
# ╠═c34c44bf-1485-4f7d-8715-16fadc7b79bc
# ╠═21137cef-18d9-42cb-90d7-f565a1b243eb
# ╠═d76e765f-ff56-4c7e-b61e-0b49e88f558a
# ╠═7b32578c-de60-47fa-ae3f-628bb1c887ba
# ╟─95d1107f-07c5-4c31-904a-1471cb51ad33
# ╠═631a1822-5bf0-4713-b048-a17b93bfa66e
# ╠═773f5479-a7de-4a81-b513-0545913b37bd
# ╠═ae12d403-1149-426e-ab13-5afc5e5615d9
# ╟─c2071134-3355-4479-a2d7-55d2141bb7ff
# ╟─d0afa1de-101b-4eea-b4c3-75d2b4bc27d9
# ╟─70f59af9-f609-4b9b-b23e-606bc3b2efa2
# ╠═51988df7-5be5-478f-8cda-b4999b32ff6f
# ╟─207929ea-af4a-4d4c-b2d2-59a0dc927ba6
# ╟─e9ea4596-7d7b-4e87-87d2-88fb40a1b6ad
# ╟─52894f69-76ab-4045-b70c-67bea0043279
# ╠═70476796-ae99-490b-a77f-7b609b9e5b5f
# ╟─7e8ef8b5-630b-4150-aa46-d58537906103
# ╠═3cd039d3-3926-4889-a743-a8b908bf1796
# ╠═c758524b-04f3-4699-8594-cb65fe30f7a4
# ╠═8cfa738e-05cc-4d86-b40c-86442d14b4b1
# ╠═a5565a17-6572-4946-8147-7a9c7ca203f2
# ╠═3f4d5d0e-8f35-4f57-8b0f-2aae19b2f7be
# ╠═725efcfe-6514-42d9-a92a-5b8cdbf8caa9
# ╠═21d63be3-b945-452b-91d8-63efb3568b95
