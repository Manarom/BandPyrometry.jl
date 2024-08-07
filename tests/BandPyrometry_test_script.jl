# BandPyrometry test Script

using Plots, Revise, DelimitedFiles, Interpolations,DataInterpolations
cur_folder = @__DIR__
upper_folder = join(split(cur_folder,"\\")[1:end-1],"\\")
if isfile( joinpath(upper_folder,"Planck.jl"))
    #includet(upper_folder*"Planck.jl")
    includet(joinpath(upper_folder,"BandPyrometry.jl"))
end
module i1
    using ..Planck,..Plots,..Interpolations
    import ..BandPyrometry
    λ_full = collect(range(2,18,length=200)) # full spectral range
    e_real = 0.5*(λ_full.^0) .+ 1e-2*λ_full .+ 1e-3*λ_full.^2 
    p1 = plot(λ_full, e_real) 
    display(p1 )
    Treal = 1973.4
    Imeasured_full = ibb.(λ_full,Treal).*e_real
    p2 = plot(λ_full,Imeasured_full)
    display(p2)
    Iint = linear_interpolation(λ_full,Imeasured_full)
    polynomial_type = "leg"#"chebT","trig", "leg"
    initial_x = [0.5,0,1767.3]
    run1(;optimizer="Default",
        is_constraint::Bool=false,
        λ_left::Float64=0.0,
        λ_right::Float64=Inf,x_init::Union{Vector{Float64},Nothing}=initial_x) = begin
            flag  = (λ_full.<λ_right) .&   (λ_full.>λ_left) 
            n = sum(flag)
            if n==0
                return nothing
            end
            λ = Vector{Float64}(undef,n)
            λ .= λ_full[flag]
            bp1 = BandPyrometry.BandPyrometryPoint(
                Imeasured_full[flag],
                λ,
                x_init,
                polynomial_type=polynomial_type
            )
            fitting_out = BandPyrometry.fit_T(bp1,optimizer_name=optimizer,is_constraint=is_constraint)
        return (λ,Imeasured_full, fitting_out...,bp1) 
    end
end
i1.run1(optimizer="NelderMead",
is_constraint=false,
λ_left=5.0,
λ_right=9.0,x_init=[0.5,574.5])