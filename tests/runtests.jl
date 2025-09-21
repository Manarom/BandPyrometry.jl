using BandPyrometry
using PlanckFunctions
using Test #,LinearAlgebra,Statistics
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

@testset "BandPyrometry.jl" begin
    # Write your tests here.
    for point in values(TestingData.benchmark_data) #iterating over data for various temperatures
        @show T = point.temperature
        λₗ = point.lower #lower wavelength
        λᵣ = point.upper #upper wavelength
        sp = point.spectrum 
        # total radiance test
        @test PlanckFunctions.power(T) ≈ point.total_radiance rtol=1e-4 # (possible due to some difference in Stefan constant there is a discrepancy)
        #peak wavelength test
        @test PlanckFunctions.λₘ(T) ≈ point.peak_wavelength rtol=1e-6 
        #simple Planck function
        @test PlanckFunctions.ibb(PlanckFunctions.λₘ(T),T) ≈ point.peak_spectral_radiance rtol=1e-6
        #testing the in-band radiance
        @test band_power(T,λₗ =λₗ ,λᵣ =λᵣ ) ≈ point.band_radiance rtol=1e-4
        # testing least-square difference between calculated spectra
        points_number = size(sp,1)
        n = sqrt(sum(i->i^2, ibb.(sp[:,1],T) .-sp[:,2]))/points_number
        @test n ≈ 0 atol=1e-2
    end
end
