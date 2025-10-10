
module BandPyrometry
    using LinearAlgebra, #
    MKL, # using MKL turns default LinearAlgebra from library from openBLAS to mkl  
    Optimization,
    OptimizationOptimJL, 
    Interpolations,
    StaticArrays
    
    import PlanckFunctions as Planck
    include("BandPyrometryTypes.jl") # Brings types and functions for working with types
    #include("Pyrometers.jl") 
    

    export  BandPyrometryPoint,# type for least-square fitting 
            EmPoint, # type for BB temperature fitting
            fit_T!,# function to fit the BB and real surface temperature
            Pyrometers, # pyrometers module
            optimizer_switch, # supported optimizer selection function (returns optimizer constructor)
            optim_dic,
            support_constraint_optimizers,
            support_lagrange_constraints,
            em_cons!, # emissivity constraint function (returns maximum and minimum values of the emissivity in the whole spectral range)
            emissivity!, # fills the emissivity
            feval!,
            jacobian!,
            hess!,
            hess_approx!,
            grad!,
            residual!,
            disc
    """
    All supported optimizers
    """
    const optim_dic = Base.ImmutableDict("NelderMead"=>NelderMead,
                            "Newton"=>Newton,
                            "BFGS"=>BFGS,
                            "GradientDescent"=>GradientDescent,
                            "NewtonTrustRegion"=>NewtonTrustRegion,
                            "ParticleSwarm"=>ParticleSwarm,
                            "Default"=>NelderMead,
                            "LBFGS"=>LBFGS,
                            "IPNewton"=>IPNewton) # list of supported optimizers
    """
    Optimizers supporting box-constraint optimization
    """                        
    const support_constraint_optimizers = ["NelderMead", 
                                            "LBFGS",
                                            "IPNewton",
                                            "ParticleSwarm"]
    """
    Optimizer supporting lagrangian-constraint optimization
    """                                        
    const support_lagrange_constraints = ["IPNewton"]
    """
    optimizer_switch(name::String;is_constraint::Bool=false,
                is_lagrange_constraint::Bool=false)

Returns the appropriate optimizer constructor
Input:
    name - the name of th eoptimizer, 
    is_constraint - is box constraint problem formulation, 
    is_lagrange_constraint - use Lagrange constraints (supported only by IPNewton)

Returns the appropriate optimizer constructor
Input:
    name - the name of th eoptimizer, 
    is_constraint - is box constraint problem formulation, 
    is_lagrange_constraint - use Lagrange constraints (supported only by IPNewton)
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
    ## BAND PYROMETRY POINT METHODS
    """
    box_constraints(bp::BandPyrometryPoint)

Evaluates box-constraint of the problem
Evaluates box-constraint of the problem
"""
    function box_constraints(bp::BandPyrometryPoint) 
        # method calculates box constraints 
        # of the feasible region (dumb version)
            lb = copy(bp.x)
            ub = copy(bp.x)
            lb[1:end-1].=-1.0
            ub[1:end-1].=1.0
            lb[end] =20.0 # 20 Kelvins limit for the temperature
            ub[end] = 5000.0
        return (lb=lb,ub=ub)
    end
    
    """
    em_cons!(constraint_value::AbstractArray,
                            x::AbstractVector, 
                            bp::BandPyrometryPoint)


In-place filling of two-elemnt vector of [minimum,maximum] emissivity in the whole 
wavelength range  
This function is used in the constraints
Inputs:
    constraint_value - (modified)  two-element vector to be modified in-place
    x - optimization variables vector, x=[a1...an,T], where a1...an - emissivity approximations coefficients, T  - temperature 
    bp - (modified) 
"""
    function em_cons!(constraint_value::AbstractArray,
                            x::AbstractVector, 
                            bp::BandPyrometryPoint)
        # evaluate the constraints on emissivity (it should not be greater than one in a whole spectra range)
        feval!(bp,x)  
        constraint_value.=extrema(bp.ϵ) # (minimum,maximum) values of the emissivity 
        return constraint_value
        #   in a whole spectrum range
    end
    """
    emissivity!(bp::BandPyrometryPoint,x::AbstractVector)


Fills emissivity for the current BandPyrometry point
Input:
    bp - (modified) current spectral band pytometry point
    x - optimization variables vector, x=[a1...an,T],
    where a1...an - emissivity approximations coefficients, T  - temperature   
Fills emissivity for the current BandPyrometry point
Input:
    bp - (modified) current spectral band pytometry point
    x - optimization variables vector, x=[a1...an,T],
    where a1...an - emissivity approximations coefficients, T  - temperature   
"""
    function emissivity!(bp::BandPyrometryPoint,x::AbstractVector)
        a = @view x[1:end-1] #emissivity approximation variables
        bp.ϵ .= bp.vandermonde.v*a
        return bp.ϵ
    end
    """
    feval!(bp::BandPyrometryPoint,x::AbstractVector)

Fills both the emissivity and the thermal emission spectrum for the current BandPyrometry point

Input:
    bp - (modified) current spectral band pytometry point
    x - optimization variables vector, x=[a1...an,T],
    where a1...an - emissivity approximations coefficients, T  - temperature
"""
    function feval!(bp::BandPyrometryPoint,x::AbstractVector)
        # evaluates residual vector
        #a = @view x[1:end-1] #emissivity approximation variables
        feval!(bp.e_p,x[end]) # refreshes planck function values
        if x!=bp.x_em_vec # x_em_vec - emissivity calculation vector
            if bp.is_has_Iₛᵤᵣ # has surrounding radiation correction
                bp.Ic .= (bp.e_p.Ib .- bp.e_p.Iₛᵤᵣ).*emissivity!(bp,x) # I=(Ibb-Isur)*ϵ
            else
                bp.Ic .= bp.e_p.Ib.*emissivity!(bp,x) # I=Ibb*ϵ
            end
            bp.x_em_vec.=x
        end
        return bp.Ic
    end
    """
    residual!(bp::BandPyrometryPoint,x::AbstractVector)

Fills emissivity, thermal emission spectrum and evaluates the residuals vector 
for the current BandPyrometry point

Input:
    bp - (modified) current spectral band pytometry point
    x - optimization variables vector, x=[a1...an,T],
    where a1...an - emissivity approximations coefficients, T  - temperature    
"""
    function residual!(bp::BandPyrometryPoint,x::AbstractVector)
        feval!(bp,x)   # feval! calculates function value only if current x is not the same as 
        bp.r .=bp.e_p.I_measured .- bp.Ic # measured data - calculated 
        bp.e_p.r[] = 0.5*norm(bp.r)^2 # discrepancy value
        return bp.r # returns residual vector
    end
    
    """
    disc(x::AbstractVector,bp::BandPyrometryPoint)

Fills emissivity, thermal emission spectrum,evaluates the residuals vector
and calculates its norm for the current BandPyrometry point
Input:
    x - optimization variables vector, x=[a1...an,T],
    where a1...an - emissivity approximations coefficients, T  - temperature 
    bp - (modified) current spectral band pytometry point
Fills emissivity, thermal emission spectrum,evaluates the residuals vector
and calculates its norm for the current BandPyrometry point
Input:
    x - optimization variables vector, x=[a1...an,T],
    where a1...an - emissivity approximations coefficients, T  - temperature 
    bp - (modified) current spectral band pytometry point
"""
    function  disc(x::AbstractVector,bp::BandPyrometryPoint)
        residual!(bp,x)
        return bp.e_p.r[]# returns current value of discrepancy
    end

"""
    jacobian!(x::AbstractVector,bp::BandPyrometryPoint)

Fills the Jacobian matrix for current bandpyrometry point
Input:
    x - optimization variables vector, x=[a1...an,T],
    where a1...an - emissivity approximations coefficients, T  - temperature
    bp - (modified, stores Jacobian internally) current spectral band pytometry point 
Fills the Jacobian matrix for current bandpyrometry point
Input:
    x - optimization variables vector, x=[a1...an,T],
    where a1...an - emissivity approximations coefficients, T  - temperature
    bp - (modified, stores Jacobian internally) current spectral band pytometry point 
"""
function jacobian!(x::AbstractVector,bp::BandPyrometryPoint) # evaluates Planck function
        ∇!(x[end],bp.e_p) # refresh Planck function first derivative
        if x!=bp.x_jac_vec
            J1 = @view bp.jacobian[:,1:end-1] # Jacobian without temperature derivatives
            J2 = @view bp.jacobian[:,end] # Last column of the jacobian 
            #a  = @view (x,1,end-1)
            J1 .= bp.e_p.Ib.*bp.vandermonde.v # diag(ibb)*V
            J2 .= bp.e_p.∇I.*emissivity!(bp,x)# 
            bp.x_jac_vec .=x # refresh jacobian calculation vector
        end
    end   

    """
    grad!(g::AbstractVector,x::AbstractVector,bp::BandPyrometryPoint)

In-place filling of the gradient vector of BandPyrometryPoint at point x
Input:
    x - optimization variables vector, x=[a1...an,T],
    where a1...an - emissivity approximations coefficients, T  - temperature
    bp - (modified, recalculates residual vector and Jacobian if the 
    currently stored value was obtaibed for another optimization variables array)
    current spectral band pytometry point     
In-place filling of the gradient vector of BandPyrometryPoint at point x
Input:
    x - optimization variables vector, x=[a1...an,T],
    where a1...an - emissivity approximations coefficients, T  - temperature
    bp - (modified, recalculates residual vector and Jacobian if the 
    currently stored value was obtaibed for another optimization variables array)
    current spectral band pytometry point     
"""
    function grad!(g::AbstractVector,x::AbstractVector,bp::BandPyrometryPoint)
        residual!(bp,x)
        jacobian!(x,bp) # calculated Jₘ
        g .= -transpose(bp.jacobian)*bp.r # calculates gradient ∇f = -Jₘᵀ*r
        return nothing
    end

    """
    hess_approx!(ha, x::AbstractVector,bp::BandPyrometryPoint)

In-place filling of the approximate hessian (Hₐ = Jᵀ*J (J - Jacobian)) 
of BandPyrometryPoint at point vector x, approximate Hessian can be used 
in optimization methods to approximate the full Hessian (e.g. in Gauss-Newton
or Levenberg-Marquardt methods)

Input:
    ha - Hessian matrix to be filled in-place
    x - optimization variables vector, x=[a1...an,T],
    where a1...an - emissivity approximations coefficients, T  - temperature
    bp - (modified) current spectral band pytometry point     
"""
function hess_approx!(ha, x::AbstractVector,bp::BandPyrometryPoint)
        # calculates approximate hessian which is Hₐ = Jᵀ*J (J - Jacobian)
        if x!=bp.x_hess_approx
            jacobian!(x,bp)
            bp.hessian_approx .= transpose(bp.jacobian)*bp.jacobian 
            # this matrix is always symmetric positive definite
            bp.x_hess_approx .=x
        end
        ha .= bp.hessian_approx
        return nothing
    end
    """
    hess!(h,x::AbstractVector,bp::BandPyrometryPoint)
    
In-place filling of the whole hessian matrix for BandPyrometryPoint at point x
Input:
    ha - Hessian matrix to be filled in-place
    x - optimization variables vector, x=[a1...an,T],
    where a1...an - emissivity approximations coefficients, T  - temperature
    bp - (modified) current spectral band pytometry point      
In-place filling of the whole hessian matrix for BandPyrometryPoint at point x
Input:
    ha - Hessian matrix to be filled in-place
    x - optimization variables vector, x=[a1...an,T],
    where a1...an - emissivity approximations coefficients, T  - temperature
    bp - (modified) current spectral band pytometry point      
"""
function hess!(h,x::AbstractVector,bp::BandPyrometryPoint)
        if x!=bp.x_hess_vec
            hess_approx!(bp.hessian,x,bp) # refreech the approximate hessian 
            # and fill hessian with approximate hessian Jᵀ*J
            # refreshes second derivative of the Planck function
            ∇²!(x[end],bp.e_p) 
            # H = Ha - Hm, Ha is approximate Hessian
            # Hm_vec = Vᵀ*I'ᴰ*r - vector Hm,
            # V - Vandermonde matrix, I'ᴰ - first 
            # derivative diagonal matrix,
            # r - residual vector
            last_hess_col = @view bp.hessian[1:end-1,end] 
            # view of the last column of the hessian 
            # initial formula: Hm_vec = Vᵀ*I'ᴰ*r  => transpose(V)*diagm(I')*r 
            # A*diagm(b) <=> A.*transpose(b) <=> transpose(Aᵀ.*b) 
            # Hm_vec = (V.*I')ᵀ*r
            last_hess_col .-= transpose(bp.vandermonde.v.*bp.e_p.∇I)*bp.r
            bp.hessian[end,1:end-1] .= last_hess_col # the sample
            # only right-down corner of hessian contains the second derivative
            # hm = rᵀ*(∇²Ibb)ᴰ*V*a
            bp.hessian[end,end] =bp.hessian[end,end] - dot(bp.r.*bp.e_p.∇²I,bp.ϵ) # dot product
            bp.x_hess_vec .=x
        end
        h.=bp.hessian # filling external matrix with internally stored hessian
        return nothing
    end

    # EMISSION POINT METHODS
    box_constraints(::EmPoint) = (lb=[20.0],ub=[5000.0]) # limits on the BB temperature
    feval!(e::EmPoint,T::AbstractArray) = feval!(e,T[end])
    """
    feval!(e::EmPoint,t::Float64)

Evaluates bb intensity for temperature t
"""
function feval!(e::EmPoint,t::Float64) # fills planck spectrum
        if t!=e.Tib[] # if current temperature is the same as the last recorded, 
            #a₁₂₃!(e_obj.amat,e_obj.λ,t) # filling amat
            Planck.a₁₂₃!(e.amat,e.λ,t) #fills amatrix
            Planck.ibb!(e.Ib, e.λ, e.amat) #fills BB spectrum
            e.Tib[] = t # save the current temperature
        end
        return e.Ib
    end
    residual!(e::EmPoint,T::AbstractArray) = residual!(e::EmPoint,T[end])
"""
    residual!(e::EmPoint,t::Float64)

Evaluates the residual vector between calculated and measured bb thermal 
emission intensity spectrum
Evaluates the residual vector between calculated and measured bb thermal 
emission intensity spectrum
"""
function residual!(e::EmPoint,t::Float64)
        feval!(e,t)
        if t!=e.Tri[] # if current temperature is the same as the last recorded, 
            e.ri .= e.I_measured .- e.Ib# calculating discrepancy
            e.r[] =0.5* norm(e.ri)^2 # discrepancy value
            e.Tri[]=t# filling temperature of residual
        end
        return e.ri # returns residual vector
    end
    
    """
    disc(T,e::EmPoint)

Evaluates the least-square discrepancy between measured and calculates spectra
Input:
    x - optimization variables vector, x=[a1...an,T],
    where a1...an - emissivity approximations coefficients, T  - temperature 
    e - (modified) current bb thermal emission point 
Evaluates the least-square discrepancy between measured and calculates spectra
Input:
    x - optimization variables vector, x=[a1...an,T],
    where a1...an - emissivity approximations coefficients, T  - temperature 
    e - (modified) current bb thermal emission point 
"""
function  disc(T,e::EmPoint)
        residual!(e,T)# fills residuals
        return e.r[] # returns current value of discrepancy
    end

    ∇!(T::AbstractVector,e::EmPoint) = ∇!(T[end],e)

    """
    ∇!(t::Float64,e::EmPoint)

Fills the first derivative of Planck function 
Input:
    T  - temperature 
    e - (modified) current bb thermal emission point   
Fills the first derivative of Planck function 
Input:
    T  - temperature 
    e - (modified) current bb thermal emission point   
"""
function ∇!(t::Float64,e::EmPoint) # evaluates Planck function first derivative
        feval!(e,t)# refreshes amat and Ib
        if t!=e.T∇ib[] # current temperature is not equal to the temperature of gradient calculation
            Planck.∇ₜibb!(e.∇I,t, e.amat,e.Ib)# fills Planck first derivative
            e.T∇ib[] = t # refresh gradient calculation temperature
        end
        return e.∇I
    end
    grad!(g::AbstractVector,T::AbstractVector ,e::EmPoint)=grad!(g,T[end] ,e)
    """
    grad!(g::AbstractVector,t::Float64 ,e::EmPoint)

In-place filling of the gradient of least-square problem of bb thermal emission spectrum fitting 
Input:
    g - gradient vector to be filled
    t  - temperature 
    e - (modified) current bb thermal emission ppoint 
"""
function grad!(g::AbstractVector,t::Float64 ,e::EmPoint)
        ∇!(t,e)
        residual!(e,t)
        if t!=e.Tgrad[]
            g[end]= - dot(e.ri,e.∇I) # filling gradient vector
            e.Tgrad[] = t
        end
        return nothing
    end

    ∇²!(T::AbstractVector,e::EmPoint)=∇²!(T[end],e)

    """
    ∇²!(t::Float64,e::EmPoint)

Fills the second derivative of Planck function 
    Input:
        T  - temperature 
        e - (modified) current bb thermal emission point  
"""
function ∇²!(t::Float64,e::EmPoint)
        ∇!(t,e)# refreshes amat and Planck gradient
        if t!=e.T∇²ib[]
           Planck.∇²ₜibb!(e.∇²I,t,e.amat,e.∇I) 
           e.T∇²ib[] = t # ref value
        end
        return e.∇²I
    end
    hess!(h,T::AbstractVector,e::EmPoint) = hess!(h,T[end],e)
    """
    hess!(h,t::Float64,e::EmPoint)

In-place filling of least-square problem hessian matrix 
    Input:
        h - hessian 
        T - temperature 
        e - (modified) current bb thermal emission point  
In-place filling of least-square problem hessian matrix 
    Input:
        h - hessian 
        T - temperature 
        e - (modified) current bb thermal emission point  
"""
function hess!(h,t::Float64,e::EmPoint) # calculates hessian of a simple Planck function fitting
        ∇²!(t,e)
        if t!=e.Thess[]
            Base.@inbounds h[1]= dot(e.∇I,e.∇I) - dot(e.ri,e.∇²I) # fill hessian vector from current value
            e.Thess[] = t
        end
        return nothing
    end
    """
    fit_T!(point::Union{EmPoint,BandPyrometryPoint};
        optimizer_name="Default",is_constraint::Bool=false,
        is_lagrange_constraint::Bool=false)
    
Input:
    point - (modified) real suraface (BandPyrometryPoint) of blackbody thermal emission object
    (optional)
    optimizer_name - the name of optimizer must be the key of optim_dic
    is_constraint - is box-constraint flag
    is_lagrange_constraint - is Lagrange-constraint flag

Returns:
    if the point type is EmPoint, the output is:
        named tuple with (T - fitted temperature,
                            res - optimization output object,
                            optimizer - chosen optimizer)
    if the point is of the BandPyrometryPoint type, the output is:
        named tuple with (T - fitted temperature ,a - fitted emissivity approximation coefficients,
                            ϵ - emissivity spectrum in the whole wavelemgth range,
                            res - optimization output object,
                            optimizer - chosen optimizer)                
Input:
    point - (modified) real suraface (BandPyrometryPoint) of blackbody thermal emission object
    (optional)
    optimizer_name - the name of optimizer must be the key of optim_dic
    is_constraint - is box-constraint flag
    is_lagrange_constraint - is Lagrange-constraint flag
Returns:
    if the point type is EmPoint, the output is:
        named tuple with (T - fitted temperature,
                            res - optimization output object,
                            optimizer - chosen optimizer)
    if the point is of the BandPyrometryPoint type, the output is:
        named tuple with (T - fitted temperature ,a - fitted emissivity approximation coefficients,
                            ϵ - emissivity spectrum in the whole wavelemgth range,
                            res - optimization output object,
                            optimizer - chosen optimizer)                
"""
function fit_T!(point::Union{EmPoint,BandPyrometryPoint};
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
                            lcons = [0.0,0.0], # both min and max of emissivity should be not smaller than zero
                            ucons = [1.0,1.0]) # both min and max should be higher than one        
        elseif is_constraint
            bx = box_constraints(point)
            probl= OptimizationProblem(fun, 
                            starting_vector,
                            point,lb=bx.lb, ub=bx.ub)
        else
            probl= OptimizationProblem(fun, 
                                starting_vector,
                                point)           
        end
        results = solve(probl,optimizer())

        return  point isa BandPyrometryPoint ? 
                        (T=results.u[end],a=results.u[1:end-1],
                        ϵ=point.vandermonde*results.u[1:end-1],
                        res=results,
                        optimizer=optimizer) :
                        (T=results.u[end],res=results,optimizer=optimizer)
    end
end

