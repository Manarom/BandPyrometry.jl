
# THIS SHOULD BE EVALUATED IN THE MAIN MODULE
include("Planck.jl") # Brings Planck module 
#
module BandPyrometry
    using LinearAlgebra, #
    MKL, # using MKL turns default LinearAlgebra from library from openBLAS to mkl  
    Optimization,
    OptimizationOptimJL, 
    Interpolations,
    StaticArrays, #,
    Plots,
    ..Planck # ATTENTION Planck module should be in the scope!!
    import Polynomials,LegendrePolynomials
    const optim_dic = Base.ImmutableDict("NelderMead"=>NelderMead,
                            "Newton"=>Newton,
                            "BFGS"=>BFGS,
                            "GradientDescent"=>GradientDescent,
                            "NewtonTrustRegion"=>NewtonTrustRegion,
                            "ParticleSwarm"=>ParticleSwarm,
                            "Default"=>NelderMead,
                            "LBFGS"=>LBFGS,
                            "IPNewton"=>IPNewton) # list of supported optimizers
    const support_constraint_optimizers = ["NelderMead", 
                                            "LBFGS",
                                            "IPNewton",
                                            "ParticleSwarm"]
    const support_lagrange_constraints = ["IPNewton"]
    """
    optimizer_switch(name::String,is_constraint::Bool)

    Switch the appropriate optimizer
"""
function optimizer_switch(name::String;is_constraint::Bool=false,
            is_lagrange_constraint::Bool=false)
        if is_lagrange_constraint
            name = filter(x -> ==(name,x),support_lagrange_constraints)
            default_ = optim_dic[support_lagrange_constraints[1]]
            return length(name)>=1 ? get(optim_dic,name[1],default_ ) :
                                                                default_            
        elseif is_constraint
            name = filter(x -> ==(name,x),support_constraint_optimizers)
            return length(name)>=1 ? get(optim_dic,name[1],optim_dic["Default"]) :
            optim_dic["Default"]
        else
            return get(optim_dic,name,optim_dic["Default"])
        end

end
    include("BandPyrometryTypes.jl") # Brings types and functions for working with types
    ## BAND PYROMETRY POINT METHODS
    """
        Evaluates lower box boundary
    """
    function lower_box_constraints(bp::BandPyrometryPoint) 
        # methods calculates box constraints 
        # of the fesible region dumb version
        lower_bond_vector = copy(bp.x)
        a = @view lower_bond_vector[1:end-1]
        fill!(a,-1.0)
        lower_bond_vector[end] =200 # 200 Kelvins limit for the temperature
        return lower_bond_vector
    end
    """

        Evaluates upper box boundary

    """    
    function upper_box_constraints(bp::BandPyrometryPoint) # methods calculates box constraints 
        # of the fesible region dumb version
        upper_bond_vector = copy(bp.x)
        #max_wavelength = max(bp.e_p.λ)
        max_ind = argmax(bp.e_p.λ)
        a = @view upper_bond_vector[1:end-1]
        vm = @view bp.vandermonde.v[max_ind,:]
        a .= 1.0 ./vm
        upper_bond_vector[end] =5000.0
        return upper_bond_vector
    end
    """
        Evaluates maximum emissivity in the whole wavelength range 
        This function is used in the constraints

    """
    function em_cons!(constraint_value::AbstractArray,
                            x::AbstractVector, 
                            bp::BandPyrometryPoint)
        # evaluate the constraints on emissivity (it should not be greater than one in a whole spectra range)
        feval!(bp,x)  
        constraint_value.=extrema(bp.ϵ) # (minimum,maximum) values of the emissivity 
        return nothing
        #   in a whole spectrum range
    end
    """
        Fills emissivity for the current BandPyrometry point
    """
    function emissivity!(bp::BandPyrometryPoint,x::AbstractVector)
        a = @view x[1:end-1] #emissivity approximation variables
        bp.ϵ .= bp.vandermonde*a
        return bp.ϵ
    end
    """
        Fills emissivity and emission spectra 
    """
    function feval!(bp::BandPyrometryPoint,x::AbstractVector)
        # evaluates residual vector
        #a = @view x[1:end-1] #emissivity approximation variables
        feval!(bp.e_p,x[end]) # refreshes planck function values
        if x!=bp.x_em_vec
            if bp.is_has_Iₛᵤᵣ
                bp.Ic .= (bp.e_p.Ib .- bp.e_p.Iₛᵤᵣ).*emissivity!(bp,x) # I=(Ibb-Isur)*ϵ
            else
                bp.Ic .= bp.e_p.Ib.*emissivity!(bp,x) # I=Ibb*ϵ
            end
            bp.x_em_vec.=x
        end
        return bp.Ic
    end
    """
        Fills emissivity, emission spectra and evaluates residual vector
    """
    function residual!(bp::BandPyrometryPoint,x::AbstractVector)
        feval!(bp,x)   # feval! calculates function value only if current x is not the same as 
        bp.r .=bp.e_p.I_measured .- bp.Ic # measured data - calculated 
        bp.e_p.r[] = 0.5*norm(bp.r)^2 # discrepancy value
        return bp.r # returns residual vector
    end
    
    """
    disc(x::AbstractVector,bp::BandPyrometryPoint)

    Fills discrepancy value, bo.e_p strores the residual function norm
"""
function  disc(x::AbstractVector,bp::BandPyrometryPoint)
        residual!(bp,x)
        return bp.e_p.r[]# returns current value of discrepancy
    end

    function jacobian!(x::AbstractVector,bp::BandPyrometryPoint) # evaluates Planck function
        ∇!(x[end],bp.e_p) # refresh Planck function first derivative
        if x!=bp.x_jac_vec
            J1 = @view bp.jacobian[:,1:end-1] # Jacobian without temperature derivatives
            J2 = @view bp.jacobian[:,end] # Last column of the jacobian 
            #a  = @view (x,1,end-1)
            J1 .= bp.e_p.Ib.*bp.vandermonde
            J2 .= bp.e_p.∇I.*emissivity!(bp,x)
            bp.x_jac_vec.=x
        end
    end   

    function grad!(g::AbstractVector,x::AbstractVector,bp::BandPyrometryPoint)
        residual!(bp,x)
        jacobian!(x,bp) # calculated Jₘ
        g .= -transpose(bp.jacobian)*bp.r # calculates gradient ∇f = -Jₘᵀ*r
        return nothing
    end
    function hess_approx!(ha, x::AbstractVector,bp::BandPyrometryPoint)
        # calculates approximate hessian which is Hₐ = Jᵀ*J (J - Jacobian)
        if x!=bp.x_hess_approx
            jacobian!(x,bp)
            bp.hessian_approx .= transpose(bp.jacobian)*bp.jacobian 
            # this matrix is always symmetric positive definite
        end
        ha .= bp.hessian_approx
        return nothing
    end
    function hess!(h,x::AbstractVector,bp::BandPyrometryPoint)
        if x!=bp.x_hess_vec
            hess_approx!(h,x,bp) # filling approximate hessian Jᵀ*J
            ∇²!(x,bp.e_p) # refreshes second derivative of the Planck function
            # H = Ha - Hm, Ha is approximate Hessian
            # Hm_vec = Vᵀ*I'ᴰ*r - vector Hm,
            # V - Vandermonde matrix, I'ᴰ - first 
            # derivative diagonal matrix,
            # r - residual vector
            last_hess_col = @view bp.h[1:end-1,end] 
            # view of last column of hessian without 
            # initial formula: Hm_vec = Vᵀ*I'ᴰ*r  => transpose(V)*diagm(I')*r 
            # A*diagm(b) <=> A.*transpose(b) <=> transpose(Aᵀ.*b) 
            # Hm_vec = (V.*I')ᵀ*r
            last_hess_col .-= transpose(bp.vandermonde.*bp.e_p.∇I)*bp.r
            bp.h[end,1:end-1] .=last_hess_col # the sample
            # only right-down corner of hessian contains the second derivative
            # hm = rᵀ*(∇²Ibb)ᴰ*V*a
            bp.h[end,end] .-= dot(bp.r.*bp.e_p.∇²I,bp.ϵ) # dot product
        end
        h.=bp.h
        return nothing
    end

    # EMISSION POINT METHODS
    lower_box_constraints(::EmPoint) = [0.0] # emissivity point has fixed boundaries
    upper_box_constraints(::EmPoint) =[5000.0]
    feval!(e::EmPoint,T::AbstractArray) = feval!(e,T[end])
    function feval!(e::EmPoint,t::Float64) # fills planck spectrum
        if !isapprox(t,e.Tib[]) # if current temperature is the same as the last recorded, 
            #a₁₂₃!(e_obj.amat,e_obj.λ,t) # filling amat
            Planck.a₁₂₃!(e.amat,e.λ,t) #fills amatrix
            Planck.ibb!(e.Ib, e.λ, e.amat) 
            e.Tib[] = t # save the current temperature
        end
        return e.Ib
    end

    residual!(e::EmPoint,T::AbstractArray) = residual!(e::EmPoint,T[end])
    function residual!(e::EmPoint,t::Float64)
        feval!(e,t)
        if !isapprox(t,e.Tri[]) # if current temperature is the same as the last recorded, 
            e.ri .= e.I_measured .- e.Ib
            e.r[] =0.5* norm(e.ri)^2 # discrepancy value
            e.Tri[]=t
        end
        return e.ri # returns residual vector
    end
    
    function  disc(T,e::EmPoint)
        residual!(e,T)
        return e.r[] # returns current value of discrepancy
    end

    ∇!(T::AbstractVector,e::EmPoint) = ∇!(T[end],e)
    function ∇!(t::Float64,e::EmPoint) # evaluates Planck function
        feval!(e,t)
        if !isapprox(t,e.T∇ib[]) # current temperature is not equal to the temperature of gradient calculation
            Planck.∇ₜibb!(e.∇I,t, e.amat,e.Ib)
            e.T∇ib[] = t # refresh gradient calculation temperature
        end
        return e.∇I
    end

    function grad!(g::AbstractVector,T ,e::EmPoint)
        ∇!(T,e)
        residual!(e,T)
        if !isapprox(T[end],e.Tgrad[])
            g[1] = -dot(e.ri,e.∇I) # filling gradient vector
            e.Tgrad[] = T[1]
        end
        return nothing
    end


    ∇²!(T::AbstractVector,e::EmPoint)=∇²!(T[end],e)
    function ∇²!(t::Float64,e::EmPoint)
        ∇!(t,e)
        if !isapprox(t,e.T∇²ib[])
           Planck.∇²ₜibb!(e.∇²I,t,e.amat,e.∇I) 
           e.T∇²ib[] = t # ref value
        end
        return e.∇²I
    end
    function hess!(h,T,e::EmPoint) # calculates hessian of a simple Planck function fitting
        t = T[end]
        ∇²!(t,e)
        if !isapprox(t,e.Thess[])
            Base.@inbounds h[1]= (dot(e.∇I,e.∇I) - dot(e.ri,e.∇²I)) # fill hessian vector from current value
            e.Thess[] = T[1]
        end
        return nothing
    end
    # OPTIMIZATION TOOLS
    function fit_T(point::Union{EmPoint,BandPyrometryPoint};
        optimizer_name="Default",is_constraint::Bool=false,
        is_lagrange_constraint::Bool=false)
        optimizer = optimizer_switch(optimizer_name,
                                    is_constraint = is_constraint,
                                    is_lagrange_constraint = is_lagrange_constraint)
        if is_lagrange_constraint
            fun=OptimizationFunction(disc,grad=grad!,hess=hess!,cons=em_cons!)
        else
            fun=OptimizationFunction(disc,grad=grad!,hess=hess!) 
        end
        # by default all derivatives are supported
        if point isa EmPoint
            starting_vector = MVector{1}([235.0])
        else
            starting_vector = copy(point.x);
        end
        if is_lagrange_constraint
            probl= OptimizationProblem(fun, 
                            starting_vector,
                            point, 
                            lcons = [0.0,0.0], 
                            ucons = [1.0,1.0])         
        elseif is_constraint
            probl= OptimizationProblem(fun, 
                            starting_vector,
                            point, 
                            lb=lower_box_constraints(point),
                            ub=upper_box_constraints(point))
        else
            probl= OptimizationProblem(fun, 
                                starting_vector,
                                point)           
        end
        results = solve(probl,optimizer())
        return (results.u, results,optimizer)
    end
    # PLOTTING TOOLS
    function internal_fields(st,prop::Vector{Symbol})
        data=st
        for s ∈ prop
            data = Base.getproperty(data,s)
            if data isa AbstractVector
                return data
            end
        end
        
    end
    function bp_plot(plot_handle,
        data_struct,
        x::Union{Vector{Symbol}, Symbol},
        y::Union{Vector{Symbol}, Symbol}) 
        if x isa Vector
            x_data = internal_fields(data_struct,x)
        else
            x_data = Base.getproperty(data_struct,x)
        end
        if y isa Vector
            y_data = internal_fields(data_struct,y)    
        else
            y_data = Base.getproperty(data_struct,y)
        end

        if isnothing(plot_handle)
            plot_handle = plot(plot_handle,x_data,y_data);
        else
            plot!(plot_handle,x_data,y_data)
        end
        display(plot_handle)
        return (plot_handle,x_data,y_data)
    end
end

