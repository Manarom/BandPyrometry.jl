using LinearAlgebra,Interpolations,Polynomials,LegendrePolynomials,StaticArrays
# BandPyrometryTypes should be included in the BandPyrometry module

abstract type AbstractPolyWrapper{P} end

struct TrigPolyWrapper{P} <: AbstractPolyWrapper{P}
    coeffs::MVector{P,Float64}
end
struct LegPolyWrapper{P}  <: AbstractPolyWrapper{P}
    coeffs::MVector{P,Float64}
end
struct StandPolyWrapper{P}  <: AbstractPolyWrapper{P} 
    coeffs::MVector{P,Float64}
end
struct ChebPolyWrapper{P}  <: AbstractPolyWrapper{P}
    coeffs::MVector{P,Float64}
end
(::Type{T})(x::Vector) where T<:AbstractPolyWrapper = T{length(x)}(MVector{length(x)}(x))
const SUPPORTED_POLYNOMIAL_TYPES = Base.ImmutableDict(
                :stand=> StandPolyWrapper,#Standard basis polynomials from Polynomials.Polynomial,
                :chebT=> ChebPolyWrapper,# Chebyshev polynomials from Polynomials.ChebyshevT,
                :trig => TrigPolyWrapper, # Self-made primitive type for trigonometric functions
                :leg =>  LegPolyWrapper # Legendre polynomials           
)
"""
    Function to evaluate polynomials
"""
eval_poly(poly::StandPolyWrapper,x) = Polynomials.Polynomial(poly.coeffs)(x)
eval_poly(poly::ChebPolyWrapper,x) = Polynomials.ChebyshevT(poly.coeffs)(x)
eval_poly(::LegPolyWrapper,degree,x) = LegendrePolynomials.Pl(x,degree)
function eval_poly(::TrigPolyWrapper,degree,x) 
     if degree==0
        return 1
     end  
     n = 1 + floor(degree/2) 
     if isodd(degree)
        return sin(n*pi*x)
     else
        return cos(n*pi*x)
     end
end
function (poly::Union{StandPolyWrapper,ChebPolyWrapper})(x::Number)
    return eval_poly(poly,x)
end

"""
(poly::Union{LegPolyWrapper,TrigPolyWrapper})(x::Number)

    Function returns the summation of polynomials values with coefficients 
    LegPolyWrapper wraps the LegendrePolynomials.jl library function to 
    make it consistent with Polynomials.jl 
    TrigPolyWrapper simple type for trigonometric function polynomils
"""
function (poly::Union{LegPolyWrapper,TrigPolyWrapper})(x::Number)
    #LegendrePolynomials.Pl(x,l) - computes Legendre polynomial of degree l at point x 
    res = 0.0
    for (i,coeff) ∈ enumerate(poly.coeffs)
        coeff != 0.0 || continue
        res += coeff*eval_poly(poly,i-1,x)
    end
    return res
end
"""
    VanderMatrix(λ,polynomial_degree,MatrixType,polynomial_type)
    
This type stores the Vandemonde matrix (fundamental matrix of basis functions),
supports various types of internal polynomials including provided externally
Structure VanderMatrix has the following fields:
    v - the matrix itself (each column of this matrix is the value of basis function)
    x_first -  first element of the initial vector 
    x_last  -  the last value of the initial vector
    xi - normalized vector 
    poly_type  - string of polynomial type name (nothing depends on this name)
    poly_constructor  - the constructor of the polynomial object, accepts the value of 
        polynomial coefficients and returns a callable object to evaluate the obtained 
        polynomial at a praticular point. V-matrix is filled by first creating the polynomial 
        obj with only one non-zero polynomial coefficient and then sending the values of xi
        to the created object 
This type stores the Vandemonde matrix (fundamental matrix of basis functions),
supports various types of internal polynomials including provided externally
Structure VanderMatrix has the following fields:
    v - the matrix itself (each column of this matrix is the value of basis function)
    x_first -  first element of the initial vector 
    x_last  -  the last value of the initial vector
    xi - normalized vector 
    poly_type  - string of polynomial type name (nothing depends on this name)
    poly_constructor  - the constructor of the polynomial object, accepts the value of 
        polynomial coefficients and returns a callable object to evaluate the obtained 
        polynomial at a praticular point. V-matrix is filled by first creating the polynomial 
        obj with only one non-zero polynomial coefficient and then sending the values of xi
        to the created object 
"""
struct VanderMatrix{M <: SMatrix , R <: SMatrix, V <: SVector}
    v::M # matrix of approximating functions 
    v_unnorm::M # unnormalized vandermatrix used to convert fitted parameters if necessarys
    # QR factorization matrices
    Q::M 
    R::R
    x_first::Float64 # first element of the initial array
    x_last::Float64 # normalizing coefficient 
    xi::V # normalized vector 
    poly_type::Symbol
end# struct spec

"""
    VanderMatrix(x::AbstractVector;
                    poly_degree::Int,
                    poly_type::Symbol = :stand)

Input: 
    x  - vector of independent variables (coordinates)
    poly_degree - degree of polynomial starting from zero, 
        (e.g. for standard basis constant value poly_degree=0)
        the number of V-matrix columns is poly_degree + 1
    poly_type - polynomial type name, must be member of SUPPORTED_POLYNOMIAL_TYPES

"""
function VanderMatrix(x::SVector{L},
                      ::Val{CN};
                      poly_type::Symbol = :stand,
                    ) where {L,CN}

            @assert issorted(x) "The x-vector must be sorted in ascending order"
            @assert allunique(x) "All x-values must be unique"
            # @assert P >= 0 "Degree of polynomial must be greater or equal zero"
            @assert haskey(SUPPORTED_POLYNOMIAL_TYPES,poly_type) "Polynomial type must be member of $(keys(SUPPORTED_POLYNOMIAL_TYPES))"# polynomial must be of supported type
            
            PolyWrapper = SUPPORTED_POLYNOMIAL_TYPES[poly_type]
            # CN = P + 1 # number of columns
            (_xi,x_first,x_last) = normalize_x(x)
           
            V = Matrix{Float64}(undef,L,CN) 
            Vunnorm = Matrix{Float64}(undef,L,CN)  
            poly_obj = PolyWrapper(zeros(CN))
            _fill_vander!(V, poly_obj,_xi)
            _fill_vander!(Vunnorm, poly_obj,x)

            MatrixType  = SMatrix{L, CN, Float64, L*CN}
            RMatrixType = SMatrix{CN, CN, Float64, CN * CN}
            VectorType = SVector{L,Float64}

            _V = MatrixType(V)
            (Q,R) = qr(_V)
            # {MatrixType,RMatrixType,VectorType,PolyWrapper}
            VanderMatrix{MatrixType,RMatrixType,VectorType}(_V,# Vandermonde matrix
                MatrixType(Vunnorm), #unnormalized vandermatrix
                MatrixType(Q),
                RMatrixType(R),
                x_first, # first element of the initial array
                x_last, # normalizing coefficient 
                VectorType(_xi),  
                poly_type
            )
end


function _fill_vander!(V, poly_obj::AbstractPolyWrapper,xi)
    VW = @views eachcol(V)
    for (i,col) ∈ enumerate(VW)
        poly_obj.coeffs[i] = 1.0                              
        @. col = poly_obj(xi)
        poly_obj.coeffs[i] = 0.0
    end  
    return poly_obj
end
"""
    is_the_same_x(v::VanderMatrix,x::AbstractVector)

Checks if input `x` is the same as the one used for VanderMatrix creation
"""
function is_the_same_x(vander::VanderMatrix,x::AbstractVector) 
    L = length(x)
    (L == length(vander.xi) && issorted(x) ) || return false
    x_f = first(x)
    vander.x_first == x_f || return false
    x_l = last(x)
    vander.x_last == x_l || return false
    for i in 1:L 
        x[i] == denormalize_x(vander.xi[i], x_f, x_l)  || return false
    end
    return true
end
"""
    *(V::VanderMatrix,a::AbstractVector)

VanderMatrix object can be directly multiplyed by a vector
"""
Base.:*(V::VanderMatrix,a::AbstractVector) =V.v*a 
"""
    polyfit(V::VanderMatrix,x::T,y::T) where T<:AbstractVector

Fits data x - coordinates, y - values using the VanderMatrix
basis function

Input:
    x - coordinates, [Nx0]
    y - values, [Nx0]
returns tuple with vector of polynomial coefficints, values of y_fitted at x points
and the norm of goodness of fit     
"""
function polyfit(V::VanderMatrix,x::T,y::T; normalized_x::Bool = true) where T<:AbstractVector
    yi =  !is_the_same_x(V,x) ? linear_interpolation(x,y)(denormalize_x(V)) : copy(y)
    a = V.R\(transpose(V.Q)*yi) # calculating pseudo-inverse
    y_fit = V*a
    goodness_fit = norm(yi .- y_fit)
    return  normalized_x ? (a, y_fit, goodness_fit) : (V.v_unnorm\(V*a), y_fit, goodness_fit)
end

"""
    normalize_x(x::AbstractVector)

Makes all elements of vector x to fit in range -1...1
returns normalized vector , xmin and xmax values
All elements of x must be unique
Makes all elements of vector x to fit in range -1...1
returns normalized vector , xmin and xmax values
All elements of x must be unique

"""
function normalize_x(x::AbstractVector)
    @assert allunique(x) "All elements of input vector must be unique"
    _x = !issorted(x) ? sort(x) : x

    x_min = _x[1]
    x_max = _x[end]
    a = 2.0/(x_max - x_min)
    xi_converted  = @. a*(_x - x_min) - 1.0
    return (xi_converted ,x_min,x_max)
end

"""
    denormalize_x(normalized_x::AbstractVector, x_min,x_max)

Creates normal vector from one created with [`normalize_x`](@ref)` function 
"""
function denormalize_x(normalized_x::Number, x_min, x_max)
    return 0.5*(normalized_x + 1.0)*(x_max - x_min) + x_min
end
denormalize_x(V::VanderMatrix) = Vector(denormalize_x.(V.xi, V.x_first,V.x_last))
#function 
"""
EmPoint type stores data on thermal emission spectrum and its 
first and second derivatives it also stores "Measurements " vector 
which further can be fitted? it also provides the constructor
EmPoint(I_measured,λ) -  I_measured is a measured spectrum
                      -  λ - wavelength vector (in μm)

"""
struct EmPoint{VectType,AmatType} 
    I_measured :: VectType#MVector{N,Float64}# data to fit
    λ:: VectType #MVector{N,Float64}  # wavelength vector (it is constant during the optimization)
    Ib::VectType#MVector{N,Float64} # Planck function values vector ????
    ri::VectType#MVector{N,Float64} # residual vector
    r::Base.RefValue{Float64} # discrepancy value
    ∇I::VectType#MVector{N,Float64} # first derivative value
    ∇²I::VectType#MVector{N,Float64} # second derivative vector
    amat::AmatType#MMatrix{N,3,Float64,L} # intermediate private data used to speed up the planck function evaluation
    # temperatures of:
    Tib::Base.RefValue{Float64} # Planck intensity  
    Tri::Base.RefValue{Float64} # Residual vector  
    T∇ib::Base.RefValue{Float64} # Planck derivative  
    Tgrad::Base.RefValue{Float64} # Gradient of emission discrepancy function 
    T∇²ib::Base.RefValue{Float64} # Planck function second derivative 
    Thess::Base.RefValue{Float64} # Discrepancy function second derivative
    points_number::Int64 # number of points in wavelength
    """
    EmPoint(I_measured::AbstractVector,λ::AbstractVector)

Constructor of the EmPoint object instance
Input: 
    I_measured - mesured blackbody spectral intensity
    λ - wavelength in μm
"""
function EmPoint(I_measured::AbstractVector,λ::AbstractVector)

       points_number = length(λ)
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
            Ref(0.0),#  stores the temperature of Planck function evaluation
            Ref(0.0), # Tri
            Ref(0.0),# T∇ibb
            Ref(0.0),# Tgrad
            Ref(0.0),# T∇²ibb
            Ref(0.0),# Tsec
            points_number#points_number::Int64 # number of points in wavelength            
       ) # calling the constructor

    end
end

function triplicate_columns(a::AbstractVector,T)
    return T(repeat(a,1,3))
end
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
                      
"""
struct BandPyrometryPoint{Lx1,Px1,LxP,PxP,LxPm1} 
    # Stores data about the spectral band, BBemission spectrum and experimental measured spectrum
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
    vandermonde::VanderMatrix # Vandermonde matrix type
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
                        polynomial_type::String="stand",
                        I_sur::AbstractVector=[])

Constructor for band pyrometry fitting, 
    λ - wavelength vector, 
    initial_x - starting optimization vector 
    polynomial_type - type of polynomial for emissivity approximation
"""
function BandPyrometryPoint(measured_Intensity::AbstractVector,
                        λ::AbstractVector,
                        initial_x::AbstractVector;
                        polynomial_type::Symbol=:stand,
                        I_sur::AbstractVector=[])

       polynomial_type = haskey(SUPPORTED_POLYNOMIAL_TYPES,polynomial_type) ? polynomial_type : :stand
       PolyWrapper = SUPPORTED_POLYNOMIAL_TYPES[polynomial_type]


       # if entered polynomial type is not supported then it turns to "simple"
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
       is_has_Iₛᵤᵣ ? Isr =Lx1(I_sur) : Isr = Lx1(undef)
       @assert length(measured_Intensity) == L
       new{Lx1,Px1,LxP,PxP,LxPm1}(
                EmPoint(measured_Intensity::AbstractVector{Float64},λ::AbstractVector{Float64}),# filling BB emission obj
                Px1(initial_x), #em_poly
                Lx1(undef), # emissivity
                Lx1(undef), # Ic corrected emission spectrum
                Isr, # Iₛᵤᵣ surrounding radiation exclusion
                Lx1(undef), # r,residual vector function
                LxP(undef),# jacobian
                PxP(undef),# approximate hessian
                PxP(undef),# hessian
                VanderMatrix(λ, # wavelength
                            polynomial_degree, # polynomial degree (constant polynomial is zero!)
                            MatrixType = LxPm1,# vandermonde matrix type
                            poly_constructor = polynomial_producing_function, # polynomial object Constructor
                            poly_type = polynomial_type # type of polynomial
                ),
                Px1(undef), # x_em_vec - emissivity evaluation vector
                Px1(undef), # x_jac_vec - Jacobian evaluation vector
                Px1(undef), # x_hess_approx - approximate Hessian evaluation vector
                Px1(undef), # x_hess_vec - rigorous hessian evaluation vector
                is_has_Iₛᵤᵣ  # is_has_Iₛᵤᵣ - flag is true if the point has surrounding radiation correction part
                )
    end
end


