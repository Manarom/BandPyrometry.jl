
# BandPyrometryTypes should be included in the BandPyrometry module
include("PolynomialWrappers.jl")
#function 
"""
EmPoint type stores data on thermal emission spectrum and its 
first and second derivatives it also stores "Measurements " vector 
which further can be fitted? it also provides the constructor
EmPoint(I_measured,λ) -  I_measured is a measured spectrum
                      -  λ - wavelength vector (in μm)

"""
struct EmPoint{N,Nx3,T} 
    I_measured::MVector{N,T}# data to fit
    λ:: MVector{N,T}  # wavelength vector (it is constant during the optimization)
    Ib::MVector{N,T} # Planck function values vector ????
    ri::MVector{N,T} # residual vector
    r::Base.RefValue{T} # discrepancy value
    ∇I::MVector{N,T} # first derivative value
    ∇²I::MVector{N,T} # second derivative vector
    amat::MMatrix{N,3,T,Nx3} # intermediate private data used to speed up the planck function evaluation
    # temperatures of:
    Tib::Base.RefValue{T} # Planck intensity  
    Tri::Base.RefValue{T} # Residual vector  
    T∇ib::Base.RefValue{T} # Planck derivative  
    Tgrad::Base.RefValue{T} # Gradient of emission discrepancy function 
    T∇²ib::Base.RefValue{T} # Planck function second derivative 
    Thess::Base.RefValue{T} # Discrepancy function second derivative
    """
    EmPoint(I_measured::AbstractVector,λ::AbstractVector)

Constructor of the EmPoint object instance
Input: 
    I_measured - mesured blackbody spectral intensity
    λ - wavelength in μm
"""
function EmPoint(I_measured::StaticArray{Tuple{N},T,1},λ::StaticArray{Tuple{N},T,1}) where {N,T}

       #points_number = length(λ)
       SVectType = SVector{N,T}
       MVectType = MVector{N,T}
       MatType = MMatrix{N,3,T,3*N} 

       new{N,3*N,T}(
            SVectType(I_measured), # measured value
            MVectType(λ),# wavelength
            MVectType(undef),#Ib::AbstractVector# Planck function values vector
            MVectType(undef),#ri::AbstractVector  # discrepancy vector
            Ref(maxintfloat(Float64)),#r::Base.RefValue{Float64}# discrepancy value
            MVectType(undef),#∇I::AbstractVector # first derivative value
            MVectType(undef),#∇²I::AbstractVector # second derivative vector
            MatType(zeros(N,3)),#amat::Matrix{Float64} # intermediate private data 
            Ref(0.0),#  stores the temperature of Planck function evaluation
            Ref(0.0), # Tri
            Ref(0.0),# T∇ibb
            Ref(0.0),# Tgrad
            Ref(0.0),# T∇²ibb
            Ref(0.0),# Tsec          
       ) # calling the constructor

    end
end
VanderMatrix(em::EmPoint,vv::Val{CN};poly_type::Symbol = :stand) where CN = VanderMatrix(em.λ,vv,poly_type = poly_type)

"""
    BandPyrometryPoint type stores data of thermal emission spectrum of a real body with 
emissivity polynomial approximation, and  its first and second derivatives
it also stores "measurements" vector which further can be fitted, it also 
provides the constructor BandPyrometryPoint(I_measured,λ,initial_x,polynomial_type) 
where:
    -  I_measured is a measured spectrum
    -  λ - wavelength vector (in μm)
    - initial_x - starting parameters vector (initial_x[end] - starting temperature,
        x[1:end-1] - emissivity approximation coefficients)
    - polynomial_type - string of polynomial (this value governs the Vandermonde matrix form)

N - number ob wavelength points
Nx3 - 3*N used to store the intermediate data for planck function and derivarives evaluation
L - length of optimization variables vector
P -   number of variables to be optimized (can be less or equal to the number of )
NxP - N*P number of elements in jacobian 
PxP - P*P number of elements in hessian
Pm1 - P-1 number of parameters approximating 
NxPm1 - N*(P-1) number of vandermatrix elements
T - type of data
"""
struct BandPyrometryPoint{N , Nx3 ,P, NxP, PxP, Pm1 , NxPm1, Pm1xPm1 , T}#{N,P,T} # N - wavelength number, CN - parameters number + 1
    # N , Nx3 , P, NxP, PxP, Pm1 , NxCN, CNxCN , T
    # Stores data about the spectral band, BBemission spectrum and experimental measured spectrum
    e_p::EmPoint{N,Nx3,T} 
    # Additional data storages
    x::MVector{P,T}   #Px1 # Optimization variables vector
    # x[end] - temperature, x[1:end-1] - emissivity poynomial approximation
    Ic::MVector{N,T} #Lx1 # sample spectral emittance
    Iₛᵤᵣ::SVector{N,T} #::Lx1 # surrounding radiation spectra
    ϵ::MVector{N,T}  #Lx1 # spectral emissivity in band
    r::MVector{N,T} #Lx1 # residual vector
    jacobian::MMatrix{N,P,T,NxP}#LxP # Jacobian matrix
    hessian_approx::MMatrix{P,P,T,PxP} #PxP # approximate hessian matrix
    hessian::MMatrix{P,P,T,PxP} # Hesse matrix
    vandermonde::VanderMatrix{N, Pm1, T, NxPm1, Pm1xPm1} # Vandermonde matrix type VanderMatrix{N,CN,T,NxCN,CNxCN} - N - rows number, CN - columns of vander number ()
    # internal usage
    x_em_vec::MVector{P,T} # vector of emissivity evaluation values
    x_jac_vec::MVector{P,T} #Px1 # vector of jacobian claculation parameters
    x_hess_approx::MVector{P,T}#Px1 # vector of the approximate hessian calculation
    x_hess_vec::MVector{P,T} #Px1 # vector of hessian calculation parameters
    is_has_Iₛᵤᵣ::Bool # flag 

"""
    BandPyrometryPoint(measured_Intensity::AbstractVector,
                        λ::AbstractVector,
                        initial_x::AbstractVector;
                        polynomial_type::Symbol="stand",
                        I_sur::AbstractVector=[])

Constructor for band pyrometry fitting, 
    λ - wavelength vector, 
    initial_x - starting optimization vector 
    polynomial_type - type of polynomial for emissivity approximation
"""
function BandPyrometryPoint(measured_Intensity::StaticArray{Tuple{N},T,1},
                        λ::StaticArray{Tuple{N},T,1},
                        initial_x::StaticArray{Tuple{P},T,1};
                        polynomial_type::Symbol=:stand,
                        I_sur::Union{StaticArray{Tuple{N},T,1},Nothing}=nothing) where {N,P,T}

       PolyTypeAbs = haskey(SUPPORTED_POLYNOMIAL_TYPES,polynomial_type) ? SUPPORTED_POLYNOMIAL_TYPES[polynomial_type] : StandPolyWrapper

       # if entered polynomial type is not supported then it turns to "simple"
       #L = length(λ) #total number of spectral points
       #P = length(initial_x) # full number of the optimization variables
       #polynomial_degree =  P - 2 #degree of emissivity polynomial approximation
       # polynomial degree goes from 0,1... where 1 is linear approximation
       # {N , Nx3 , P, NxP, PxP, Pm1 , NxPm1, Pm1xPm1 , T}
       Nx3 = 3*N
       NxP = N*P
       PxP = P*P 
       Pm1 = P - 1
       NxPm1 = N*(P - 1)
       Pm1xPm1 = (P - 1)*(P - 1) 
       
       # @show isconcretetype(PolyTypeAbs{Pm1})

       Nx1_T = MVector{N,T} # independent data column
       Px1_T = MVector{P,T} # optimization variables vector 
       NxP_T = MMatrix{N,P,T,NxP} # Jacobian type
       PxP_T = MMatrix{P,P,T,PxP} # Hessian type     
       #LxPm1_T = MMatrix{N,Pm1,T,NxPm1} #Vandermonde matrix type
       is_has_Iₛᵤᵣ = !isnothing(I_sur) && length(I_sur)==length(λ)
       Isr =  is_has_Iₛᵤᵣ ? SVector{N}(I_sur) : SVector{N}(zeros(T,N))
       # {N , Nx3 , P, NxP, PxP, Pm1 , NxPm1, Pm1xPm1 , T}
       new{N , Nx3 , P, NxP, PxP, Pm1 , NxPm1, Pm1xPm1 , T}(
                EmPoint(measured_Intensity,λ),# filling BB emission obj
                Px1_T(initial_x), #em_poly
                Nx1_T(undef), # emissivity
                Nx1_T(undef), # Ic corrected emission spectrum
                Isr, # Iₛᵤᵣ surrounding radiation exclusion
                Nx1_T(undef), # r,residual vector function
                NxP_T(undef),# jacobian
                PxP_T(undef),# approximate hessian
                PxP_T(undef),# hessian
                VanderMatrix(SVector{N}(λ), # wavelength
                            PolyTypeAbs{Pm1}
                ),
                Px1_T(undef), # x_em_vec - emissivity evaluation vector
                Px1_T(undef), # x_jac_vec - Jacobian evaluation vector
                Px1_T(undef), # x_hess_approx - approximate Hessian evaluation vector
                Px1_T(undef), # x_hess_vec - rigorous hessian evaluation vector
                is_has_Iₛᵤᵣ  # is_has_Iₛᵤᵣ - flag is true if the point has surrounding radiation correction part
                )
    end
end

function temperature(emp::EmPoint)
    return emp.Tib[]
end
function temperature(bpp::BandPyrometryPoint)
    return bpp.e_p.Tib[]
end
pointsnumber(::Union{EmPoint{N},BandPyrometryPoint{N}}) where N = N
parnumber(::BandPyrometryPoint{N, Nx3, P}) where  {N, Nx3, P} = P
parnumber(::EmPoint) = 1
emissivity(p::BandPyrometryPoint) = copy(p.ϵ)
#emissivity(p::BandPyrometryPoint,λ::AbstractVector) = 