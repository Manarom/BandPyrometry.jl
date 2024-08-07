using Revise,BenchmarkTools,Plots,LinearAlgebra, InteractiveUtils
include("Planck.jl")
module TestPlanck
 import ..Planck
 using ..BenchmarkTools, ..LinearAlgebra, ..InteractiveUtils
 function test_ibb()
    lam = collect(range(0.1,30,length=10_000))
    T = rand(273.0:50.0:3073.0,20) # random temperatures vector
    t_cur_ar = T[rand(1:1:20,1)]
    i_1 = similar(lam)
    t_cur = t_cur_ar[1]
    amat = Matrix{Float64}(undef, length(lam),3)
    Planck.a₁₂₃!(amat,lam,t_cur)
    ibb_ethalon = Planck.ibb.(lam,t_cur) 
# TESTING PLANCK FUNCTION
   # ibb(λ::Vector{Float64}, T::Vector{Float64})
    display("ibb( $(typeof(lam)) , $(typeof(t_cur)) )")
    display(@code_warntype Planck.ibb(lam,t_cur))
    display(@benchmark Planck.ibb($lam,$t_cur))
    display(@show norm(Planck.ibb(lam,t_cur) - ibb_ethalon))
   #  ibb(λ::Vector{Float64}, T) (with broadcasting)
    display("ibb.( $(typeof(lam)) , $(typeof(t_cur)) )")
    display(Meta.@lower Planck.ibb.(lam,t_cur))
    display(@benchmark Planck.ibb.($lam,$t_cur))
    display(@show norm(Planck.ibb(lam,t_cur) - ibb_ethalon))
   # ibb(λ::Vector{Float64}, T::Base.RefValue{Float64})
    t_cur_ref = Ref(t_cur)
    display("ibb.( $(typeof(lam)) , $(typeof(t_cur_ref)) )")
    display(Meta.@lower Planck.ibb.(lam,t_cur_ref))
    display(@benchmark Planck.ibb.($lam,$t_cur_ref))

    #display("ibb(λ::Vector{Float64}; amat)")
    #display(Meta.@lower Planck.ibb(lam;amat=amat))
    #display(@benchmark Planck.ibb($lam;amat=$amat)) 
    
 end
end