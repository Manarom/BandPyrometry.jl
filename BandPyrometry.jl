#cur_dir = @__DIR__
#cd(cur_dir)
#import Pkg
#Pkg.activate(cur_dir)]


# THIS SHOULD BE EVALUATED IN THE MAIN MODULE
include("Planck.jl") # Brings Planck module 
#
module BandPyrometry
    using LinearAlgebra, #
    MKL, # using MKL turns default LinearAlgebra from library from openBLAS to mkl  
    Optimization,
    OptimizationOptimJL, 
    DataInterpolations,
    StaticArrays, #,
    Plots,
    ..Planck # ATTENTION Planck module should be in the scope!!
    #Makie,
    #GLMakie
    const optim_dic = Base.ImmutableDict("NelderMead"=>NelderMead,
                            "Newton"=>Newton,
                            "BFGS"=>BFGS,
                            "GradientDescent"=>GradientDescent,
                            "NewtonTrustRegion"=>NewtonTrustRegion,
                            "ParticleSwarm"=>ParticleSwarm, # Particle swarm methos supports only box_type constraints
                            "Default"=>NelderMead,
                            "LBFGS"=>LBFGS,
                            "IPNewton"=>IPNewton) # list of supported optimizers

    const support_constraint_optimizers = ["NelderMead", "LBFGS","IPNewton","ParticleSwarm"]
    
    #include("BandPyrometryTypes.jl") # Brings types
"""
EmPoint type stores data about thermal emission spectrum and its first and second derivatives
it also stores "Measurements " vector which further can be fitted? it also provides the constructor
EmPoint(I_measured,λ) -  I_measured is a measured spectrum
                      -  λ - wavelength vector (in μm)

"""

struct EmPoint{VectType,AmatType} 

    I_measured :: VectType#MVector{N,Float64}# data to fit
    λ:: VectType #MVector{N,Float64}  # wavelength vector (it is constant during the optimization)
    Ib::VectType#MVector{N,Float64} # Planck function values vector ????
    ri::VectType#MVector{N,Float64} # discrepancy vector
    r::Base.RefValue{Float64} # discrepancy value
    ∇I::VectType#MVector{N,Float64} # first derivative value
    ∇²I::VectType#MVector{N,Float64} # second derivative vector
    amat::AmatType#MMatrix{N,3,Float64,L} # intermediate private data 
    #                           temparatures of :
    Tib::Base.RefValue{Float64} # Intensity evaluation
    Tri::Base.RefValue{Float64} # Residual vector evaluation
    T∇ib::Base.RefValue{Float64} # Planck derivative evaluation
    Tgrad::Base.RefValue{Float64} # Gradient of emission discrapancy function evaluation
    T∇²ib::Base.RefValue{Float64} # Planck function second derivative evaluation
    Thess::Base.RefValue{Float64} # Discrapancy function evaluation

    points_number::Int64 # number of points in wavelength
    EmPoint(I_measured::AbstractVector,λ::AbstractVector) = begin 
       #new(yin,similar(yin),similar(yin)) # fills fi vector with the same values
       points_number = length(λ)
       L = points_number*3
       VectType = MVector{points_number,Float64}
       MatType = MMatrix{points_number,3,Float64,3*points_number} 
       @assert length(I_measured)==points_number
       new{VectType,MatType}(MVector{points_number}(I_measured), # measured value
            VectType(λ),# wavelength
            VectType(undef),#Ib::AbstractVector# Planck function values vector
            VectType(undef),#ri::AbstractVector  # discrepancy vector
            Ref(maxintfloat(Float64)),#r::Base.RefValue{Float64}# discrepancy value
            VectType(undef),#∇I::AbstractVector # first derivative value
            VectType(undef),#∇²I::AbstractVector # second derivative vector
            triplicate_columns(λ,MatType),#amat::Matrix{Float64} # intermediate private data 
            Ref(Float64(0)),#Ti::Float64 # this field stores the temperature of 
            Ref(Float64(0)), # Tri
            Ref(Float64(0)),# T∇ibb
            Ref(Float64(0)),#Tgrad
            Ref(Float64(0)),# T∇²ibb
            Ref(Float64(0)),#Tsec
            points_number#points_number::Int64 # number of points in wavelength            
       ) # calling the constructor

    end
end

function triplicate_columns(a::AbstractVector,T)
    return T(repeat(a,1,3))
end
"""
BandPyrometryPoint type stores data about thermal emission spectrum of a real body with 
emissivity polynomial approximation, and  its first and second derivatives
it also stores "Measurements " vector which further can be fitted? it also provides the constructor
BandPyrometryPoint(I_measured,λ,initial_x,polynomial_type) where
                    -  I_measured is a measured spectrum
                    -  λ - wavelength vector (in μm)
                    - initial_x - starting parameters vector (initial_x[end] - starting temperature,
                                other - emissivity approximation)
                    - polynomial_type - string of polynomial (this value governs the Vandermonde matrix form)
                      
"""
struct BandPyrometryPoint{Lx1,Px1,LxP,PxP,LxPm1} 
    # Stores data about the spectral band
    e_p::EmPoint
    # Additional data storages
    x::Px1 # Optimization variables vector
    # x[end] - temperature, x[1:end-1] - emissivity poynomial approximation
    Ic::Lx1 # sample spectral emittance
    Iₛᵤᵣ::Lx1 # surrounding radiation spectra
    ϵ::Lx1 # spectral emissivity in band
    r::Lx1 # residual vector
    jacobian::LxP # Jacobian matrix
    hessian_approx::PxP # approximate hessian matrix
    hessian::PxP # Hesse matrix
    vandermonde::LxPm1 # Vandermonde matrix type
    # internal usage
    x_em_vec::Px1 # vector of emissivity evaluation values
    x_jac_vec::Px1 # vector of jacobian claculation parameters
    x_hess_approx::Px1 # vector of the approximate hessian calculation
    x_hess_vec::Px1 # vector of hessian calculation parameters

    is_has_Iₛᵤᵣ::Bool # flag 
    """
    BandPyrometryPoint(measured_Intensity::AbstractVector,
                        λ::AbstractVector,
                        initial_x::AbstractVector;
                        polynomial_type::String="simple")

    Type for band pyrometry fitting
"""
BandPyrometryPoint(measured_Intensity::AbstractVector,
                        λ::AbstractVector,
                        initial_x::AbstractVector;
                        polynomial_type::String="simple",
                        I_sur::AbstractVector=[-1.0]) = begin 
       L = length(λ) #total number of spectral points
       P = length(initial_x) # full number of the optimization variables
       polynomial_degree =  P - 2 #degree of emissivity polynomial approximation
       # polynomial degree goes from 0,1... where 1 is linear approximation
       Lx1 = MVector{L,Float64} # independent data column
       Px1 = MVector{P,Float64} # optimization variables vector 
       LxP = MMatrix{L,P,Float64,L*P} # Jacobian type
       PxP = MMatrix{P,P,Float64,P*P} # Hessian type     
       LxPm1 = MMatrix{L,P-1,Float64,L*(P-1)} #Vandermonde matrix type
       is_has_Iₛᵤᵣ = length(I_sur)==length(λ)
       is_has_Iₛᵤᵣ ? Isr =Lx1(I_sur) : Isr =Lx1(undef)
       @assert length(measured_Intensity)==L
       new{Lx1,Px1,LxP,PxP,LxPm1}(
                EmPoint(measured_Intensity::AbstractVector{Float64},λ::AbstractVector{Float64}),# filling BB emission obj
                Px1(initial_x), #em_poly
                Lx1(undef), # emissivity
                Lx1(undef), # Ic corrected emission spectrum
                Isr, # Iₛᵤᵣ surrounding radiotion exclusion
                Lx1(undef), # r,residual vector function
                LxP(undef),# jacobian
                PxP(undef),# approximate hessian
                PxP(undef),# hessian
                Vandermonde(λ,polynomial_degree,LxPm1,polynomial_type), #vandermonde
                Px1(undef), # x_em_vec
                Px1(undef), # x_jac_vec
                Px1(undef), #x_hess_approx
                Px1(undef), # x_hess_vec
                is_has_Iₛᵤᵣ  # is_has_Iₛᵤᵣ
                )
    end
end

"""
    Vandermonde(λ,polynomial_degree,MatrixType,polynomial_type)

    Returns Vandermonde matrix

"""
function Vandermonde(λ,polynomial_degree,MatrixType,polynomial_type)
    V = MatrixType(repeat(λ,1,polynomial_degree+1))
    v_view = @views eachcol(V)
    if polynomial_type=="simple"
        for (i,col) in enumerate(v_view)
            col.^=(i-1)
        end
    end
    return V
end
    ## BAND PYROMETRY POINT METHODS
    """
        Evaluates lower box boundary
    """
    function lower_box_constraints(bp::BandPyrometryPoint) # methods calculates box constraints 
        # of the fesible region dumb version
        lower_bond_vector = copy(bp.x)
        a = @view lower_bond_vector[1:end-1]
        fill!(a,0.0)
        lower_bond_vector[end] =200
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
        vm = @view bp.vandermonde[max_ind,:]
        a .= 1.0 ./vm
        upper_bond_vector[end] =5000.0
        return upper_bond_vector
    end
    """
        Evaluates maximum emissivity in the whole wavelength range 
        This function is used in the constraints
    """
    function em_cons!(constraint_value::AbstractArray,x::AbstractVector, bp::BandPyrometryPoint)
        # evaluate the constraints on emissivity (it should not be greater than one in a whole spectra range)
        feval!(bp,x)  
        constraint_value.=extrema(bp.ϵ) # (minimum,maximum) values of the emissivity 
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
                bp.Ic .= (bp.e_p.Ib .- bp.e_p.Iₛᵤᵣ).*emissivity!(bp,x) # I=Ibb*ϵ
            else
                bp.Ic .= bp.e_p.Ib.*emissivity!(bp,x) # I=Ibb*ϵ
            end
            bp.x_em_vec.=x
        end
        return bp.Ic
    end
    """
        Fills emissivity  ,emission spectra and evaluates residual vector
    """
    function residual!(bp::BandPyrometryPoint,x::AbstractVector)
        feval!(bp,x)   
        bp.r .=bp.e_p.I_measured .- bp.Ic
        bp.e_p.r[] = 0.5*norm(bp.r)^2 # discrepancy value
        return bp.r # returns residual vector
    end
    """
        Calculates discrepancy function

    """    
    function  disc(x::AbstractVector,bp::BandPyrometryPoint)
        residual!(bp,x)
        return bp.e_p.r[]# returns current value of discrepancy
    end
    """
        Fills emissivity  ,emission spectra and evaluates residual vector
    """
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
    """
        Fills gradient
    """
    function grad!(g::AbstractVector,x::AbstractVector,bp::BandPyrometryPoint)
        residual!(bp,x)
        jacobian!(x,bp) # calculated Jₘ
        g .= -transpose(bp.jacobian)*bp.r # calculates gradient ∇f = -Jₘᵀ*r
        return nothing
    end
    """
        Calculates approximate hessian which is Hₐ = Jᵀ*J (J - Jacobian)

    """
    function hess_approx!(ha, x::AbstractVector,bp::BandPyrometryPoint)
        # calculates approximate hessian which is Hₐ = Jᵀ*J (J - Jacobian)
        if x!=bp.x_hess_approx
            jacobian!(x,bp)
            bp.hessian_approx .= transpose(bp.jacobian)*bp.jacobian
        end
        ha .= bp.hessian_approx
        return nothing
    end
    function hess!(h,x::AbstractVector,bp::BandPyrometryPoint)
        if x!=bp.x_hess_vec
            hess_approx!(h,x,bp)
            #error("NOT SUPPORTED")
        end
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
        optimizer_name="Default", 
        is_constraint=true)
        optimizer = optimizer_switch(optimizer_name,is_constraint)
        #e_p = EmPoint(data_to_fit,lam) # creating spectral struct
        #t0 = 
        fun=OptimizationFunction(disc,grad=grad!,hess=hess!) # by default all derivatives are supported
        if point isa EmPoint
            starting_vector = MVector{1}([235.0])
        else
            starting_vector = copy(point.x);
        end
        if is_constraint
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
        #display("NEWTON-TRUST-REGION")
        #display(@benchmark solve($prob_NEW_TR_AN,NewtonTrustRegion()))
        results = solve(probl,optimizer())
        #display(results.original)
        return (results.u, results)
    end
    function optimizer_switch(name::String,is_constraint::Bool)
        if is_constraint
            name = filter(x -> ==(name,x),support_constraint_optimizers)
        end
        return length(name)>=1 ? get(optim_dic,name[1],optim_dic["Default"]) :
                                 optim_dic["Default"]
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

