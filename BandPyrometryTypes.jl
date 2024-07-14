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