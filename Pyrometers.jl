#include("BandPyrometry.jl")
using LinearAlgebra, #
MKL, # using MKL turns default matrix multiplication library from openBLAS to mkl  
Optimization,
OptimizationOptimJL

include("Planck.jl")
module Pyrometers
    using   ..Optimization,
            ..OptimizationOptimJL,
            ..LinearAlgebra
    import  ..Planck
    const pyrometers_types = Dict(
                    "P"=> [2.0, 2.6],
                    "M"=> [3.4], 
                    "D"=>[3.9],
                    "L"=> [4.6],
                    "E"=>[4.8, 5.2],
                    "F"=>[7.9],
                    "K"=>[8.0, 9.0],
                    "B"=> [9.1,14.0]
                    );

    struct Pyrometer # this type supports methods for radiative pyrometers
        type::String
        λ::Vector{Float64}
        ϵ::Base.RefValue{Float64}
        Pyrometer(type::String) = begin
            haskey(pyrometers_types,type) ?  new(type,
                                                pyrometers_types[type],
                                                Ref(0.0)) : error("There is no such pyrometer type")
        end
    end
    function Base.isless(p1::Pyrometer,p2::Pyrometer) # is used to sort the vector of pyrometers
        return all(p1.λ.<p2.λ)
    end
    function wavelength_number()
        return mapreduce(x->length(x),+,pyrometers_types)
    end
    function full_wavelength_range()
        #sz = mapreduce(x->length(x),+,pyrometers_types)
        λ = Vector{Float64}()
        pyr_names = Vector{String}()
        for l in pyrometers_types
            append!(λ,l[2])
            push!(pyr_names,l[1])
            if length(l[2])>1
                push!(pyr_names,l[1])
            end
        end
        inds = sortperm(λ)
        pyr_names.=pyr_names[inds]
        λ.=λ[inds]
        return λ,pyr_names
    end
    function produce_pyrometers()
        pyr_vec = Vector{Pyrometer}()
        for l in pyrometers_types
            push!(pyr_vec,Pyrometer(l[1]))
        end
        sort!(pyr_vec) # sorting according to the wavelength increase
        return pyr_vec
    end
    function set_emissivity(p::Pyrometer,em_value::Float64)
        if !(0.0<em_value<=1.0)
            em_value = 1.0
        end
        p.ϵ[]=em_value
    end
    function fit_ϵ(p::Vector{Pyrometer},Treal::Float64,Tmeasured::Vector{Float64})
        return length(p)==length(Tmeasured) ? fit_ϵ.(p,Tmeasured,Treal) : error("Vectors must be of the same size")
    end
    function fit_ϵ_wavelength(p::Vector{Pyrometer},Treal::Float64,Tmeasured::Vector{Float64})
        if length(p)!=length(Tmeasured)
            error("Vectors must be of the same size")
        end
        e_out = Vector{Float64}()
        for (pr,t) in zip(p,Tmeasured)
            append!(e_out,fit_ϵ_wavelength(pr,t,Treal))
        end    
        return e_out
    end
    function fit_ϵ(p::Pyrometer,Tmeasured::Float64,Treal::Float64)
        # fits the emssivity of the pyrometer
        # Tmeasured - is the temperature measured by the pyrometer 
        # Treal  - is the real temperature of the surface
        # the real temperature of the surface should be higher due to the emissivity
        if Treal<Tmeasured
            (Treal,Tmeasured)=(Tmeasured,Treal)
        end
        # tup = (p,Tmeasured,Treal)
        if length(p.λ)==1 # for pyrometers with fixed wavelength
            fun = OptimizationFunction((ϵ,tup)->norm(Planck.ibb(tup[1].λ[1],tup[2]) - ϵ[1]*Planck.ibb(tup[1].λ[1],tup[3])),AutoForwardDiff())
        else # pyrometer with limited band
            fun = OptimizationFunction((ϵ,tup)->norm(Planck.band_power(tup[2],λₗ=tup[1].λ[1],λᵣ=tup[1].λ[2]) .- ϵ[1]*Planck.band_power(tup[3],λₗ=tup[1].λ[1],λᵣ=tup[1].λ[2])),AutoForwardDiff())
        end
        prob = OptimizationProblem(
            fun, #fun
           [0.1],#starting vectors
           (p,Tmeasured,Treal),
           lb=[0.00],
           ub=[1.0])
        res = solve(prob,NelderMead())
        set_emissivity(p,res.u[1])
        return p.ϵ[]
    end
    function fit_ϵ_wavelength(p::Pyrometer,Tmeasured::Float64,Treal::Float64) # this is the same as fit_ϵ with the exception that 
        # this function returns a vector of values, if pyrometer is single wavelength it returns one -element array
        e_out = similar(p.λ)
        e_out .=fit_ϵ(p,Tmeasured,Treal)
        return e_out
    end
    function switch_the_type(λ::Float64)
        if ( (λ>=2)&&(λ <=2.6)) return "P"#;//0 1 - 0
            elseif isapprox(λ,3.4,atol=0.05) return "M"#;//2 - 1
            elseif isapprox(λ,3.9,atol=0.05) return "D"#;//3 - 2
            elseif isapprox(λ,4.6,atol=0.05) return "L"#;//4 - 3
            elseif ( (λ >=4.8)&& (λ<=5.2) ) return "E"#;// 5 6 - 4
            elseif isapprox( λ,7.9,atol=0.05) return "F"#;//7 - 5
            elseif ( (λ >=8)&& (λ<=9)) return "K"#;//8 9 - 6
            elseif ( (λ >=9.1)&&(λ<=14) ) return "B"#;//10 11 - 7
            else return "";
        end
    end
end