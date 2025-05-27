
module Pyrometers
    using   Optimization,
            OptimizationOptimJL,
            LinearAlgebra,
            StaticArrays 
    import  PlanckFunctions as Planck
    """
    Default pyrometers types 

    """
    const pyrometers_types = Dict(
                    "P"=> [2.0, 2.6],
                    "M"=> [3.4], 
                    "D"=>[3.9],
                    "L"=> [4.6],
                    "E"=>[4.8, 5.2],
                    "F"=>[7.9],
                    "K"=>[8.0, 9.0],
                    "B"=> [9.1,14.0]
    )

    struct Pyrometer{N} # this type supports methods for radiative pyrometers
        type::String
        λ::SVector{N,Float64}
        ϵ::Base.RefValue{Float64}
        """
    Pyrometer(type::String)

    Pyrometer object Constructor, the type of pyrometer can chosen from the pyrometers_types dictionary
    Input:
        type - pyrometer type, must be member of pyrometers_types 
"""
        Pyrometer(type::String) = begin
            haskey(pyrometers_types,type) ?  new(type,
                                                pyrometers_types[type],
                                                Ref(1.0)) : error("Unknown pyrometer type")
        end
        Pyrometer(;type::String,λ::Vector{Float64},ϵ::Float64) = begin
            new(type,)
        end
    end
    """
    measure(p::Pyrometer,i::Float64)

    Calculates the "measured" temperature from "mesaured" intensity by fitting the Planck function.
    The intensity units should be consistent with PlanckFunctions.ibb(λ,T) function,
    it should be in [W/m²⋅sr⋅μm]
    Input:
        p - pyrometer object
        i - measured intensity in [W/m²⋅sr⋅μm]
            returns temperature "measured" by this pyrometer  
        (optional)
        T_starting  - starting temperature value
"""
    function measure(p::Pyrometer,i::Float64;T_starting::Float64=600.0)
        tup = (p,i)
        if length(p.λ)==1 # single wavelength pyrometer
            fun = OptimizationFunction((t,tup)-> norm(tup[1].ϵ[]*Planck.ibb(tup[1].λ[1],t[]) - tup[2]))
        else# narrow band pyrometer
            fun = OptimizationFunction((t,tup)-> norm(tup[1].ϵ[]*Planck.band_power(t[],λₗ=tup[1].λ[1],λᵣ=tup[1].λ[2]) - tup[2]))
        end
        prob = OptimizationProblem(
                fun, #fun
                [T_starting],#starting temperature
                tup)
        return solve(prob,NelderMead()).u[]   
    end
    """
    Base.isless(p1::Pyrometer,p2::Pyrometer)

    Vector of Pyrometer objects can be sorted using isless
"""
function Base.isless(p1::Pyrometer,p2::Pyrometer) # is used to sort the vector of pyrometers
        return all(p1.λ.<p2.λ)
    end
    """
    wavelength_number()

    Returns the length of wavelengths vector covered by the pyrometers_types
"""
function wavelength_number()
        return mapreduce(x->length(x),+,pyrometers_types)
    end
    """
    full_wavelength_range()

    Creates the wavelengths vector covered by default pyrometers see [`pyrometers_types`](@ref) 
"""
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
    """
    produce_pyrometers()

    Creates the vector of all supported pyrometers 
"""
function produce_pyrometers()
        pyr_vec = Vector{Pyrometer}()
        for l in pyrometers_types
            push!(pyr_vec,Pyrometer(l[1]))
        end
        sort!(pyr_vec) # sorting according to the wavelength increase
        return pyr_vec
    end
    """
    set_emissivity(p::Pyrometer,em_value::Float64)

    Setter for spectral emissivity
"""
function set_emissivity(p::Pyrometer,em_value::Float64)
        if !(0.0<em_value<=1.0)
            em_value = 1.0
        end
        p.ϵ[]=em_value
    end

    """
    fit_ϵ(p::Vector{Pyrometer},Treal::Float64,Tmeasured::Vector{Float64})

    Optimizes the emissivity of the pyrometer to make measured by the pyrometer temperature fit
    fit the real temperature
    
    Input:
        p - pyrometer objects vector , [Nx0]
        Treal - real temperature of the surface, Kelvins
        Tmeasured - temperatures measured by the pyrometer, Kelvins, [Nx0]
"""
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
    """
    fit_ϵ(p::Pyrometer,Tmeasured::Float64,Treal::Float64)

    Optimizes the emissivity of the pyrometer to make measured by the pyrometer temperature fit
    fit the real temperature
    
    Input:
        p - pyrometer object
        Treal - real temperature of the surface, Kelvins
        Tmeasured - temperature measured by the pyrometer, Kelvins
"""
function fit_ϵ(p::Pyrometer,Tmeasured::Float64,Treal::Float64;optimizer = NelderMead())
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
        res = solve(prob,optimizer)
        set_emissivity(p,res.u[1])
        return p.ϵ[]
    end
    """
    fit_ϵ_wavelength(p::Pyrometer,Tmeasured::Float64,Treal::Float64)

    Fits emissivity and returns it as a vector of the same size as pyrometer's wavelength region
    Some pyrometers has 2-wavelength, other work on a single wavelength, for two-wavelength pyrometers
    e_out will be a two-element vector
    Input:
        p - pyrometer object
        Tmeasured - temperature measured by the pyrometer, Kelvins
        Treal - real temeprature of the surface, Kelvins
     
"""
function fit_ϵ_wavelength(p::Pyrometer,Tmeasured::Float64,Treal::Float64) # this is the same as fit_ϵ with the exception that 
        # this function returns a vector of values, if pyrometer is single wavelength it returns one -element array
        e_out = similar(p.λ)
        e_out .=fit_ϵ(p,Tmeasured,Treal)
        return e_out
    end
    """
    switch_the_type(λ::Float64)

    Returns the type (which can be used as a key of pyrometers_types dict) depending on the wavelengh
    Input:
        λ - wavelength in μm
"""
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