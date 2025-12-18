using BandPyrometry
using Test,StaticArrays,ForwardDiff#,LinearAlgebra,Statistics
import PlanckFunctions as PL
import BandPyrometry as BP
#include(joinpath(@__DIR__,"tests data","TestingData.jl")) # TestingData.benchmark_data
#=
temperature # Kelvins
lower # lower wavelength limit μm
upper # upper wavelength limit μm
total_radiance # W/(m²⋅sr)  
peak_wavelength # μm
peak_spectral_radiance # W/m²/sr/μm
band_radiance # W/m²/sr total radiance within the band lower-
spectrum # 
=#
N = 30
lam = collect(range(1.0,5.0,30))
T_real_test = [1453.57, 653.567, 2376.1]
λ = SVector{N,Float64}(lam)
λ_norm = BP.normalize_x(λ)[1] # normalized lambda vector to check gradients 
@testset "BandPyrometry.jl" begin
    # testing EmPoint and derivatives, compare to numerical derivarives evaluation
    for T_real in T_real_test
        I_measured = PL.ibb.(λ,T_real)
        em = BP.EmPoint(I_measured,λ)
        f_em(T) = begin
            0.5*sum(x->x^2,I_measured .- PL.ibb.(λ,T[]))
        end
        @test f_em([1000.0]) ≈ BP.disc(1000.0,em) # discrepancy
        g = [0.0];
        BP.grad!(g,[1000.0],em);
        @test ForwardDiff.gradient(f_em,[1000.0])[] ≈ g[] #gradient
        BP.hess!(g,[1000.0],em)
        @test ForwardDiff.hessian(f_em,[1000.0])[] ≈ g[] #hessian
        @test em() ≈ T_real # testint temperature fitting
    end


    
    em_coeff_real = [0.5, 0.01,4.2e-3]
    x_vect_test = [1,0.01,0.02,1000] # values of derivarives evaluation
    for T_real in T_real_test
        f_bp(x) = begin 
            0.5*sum(x->x^2,@. I_measured - PL.ibb(λ,x[end])*(x[1] + x[2]*λ_norm + x[3]*λ_norm^2))
        end
        e =@. (em_coeff_real[1] + em_coeff_real[2]*λ_norm + em_coeff_real[3]*λ_norm^2)
        I_measured =@. PL.ibb(λ,T_real)*e
        
        # function value test 
        bp = BandPyrometry.BandPyrometryPoint(I_measured,λ,MVector{4}([0.9,0.0,0.0,678.8]))
        @test BandPyrometry.disc(x_vect_test,bp) ≈ f_bp(x_vect_test)[]

        # gradient evaluation test
        g_num = ForwardDiff.gradient(f_bp,x_vect_test)
        g = similar(g_num)
        BandPyrometry.grad!(g,x_vect_test,bp)
        @test all(g .≈ g_num)

        # hessian evaluation test 
        h_num = ForwardDiff.hessian(f_bp,x_vect_test)
        h = similar(h_num)
        BandPyrometry.hess!(h,x_vect_test,bp)
        @test all(h .≈ h_num)
        @test T_real ≈ bp() # fitting
        @test all(bp.x_em_vec .≈ [em_coeff_real...,T_real] )
    end
end
