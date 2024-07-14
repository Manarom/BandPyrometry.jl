# Planck function and its derivatives all tempeatures should be in Kelvins, all wavelengths in μm
module Planck
    export  Dₜibb!,
            ∫ibbₗ,
            Tₖ, 
            ibb,
            ibb!,
            ∇ₜibb,
            ∇ₜibb!, 
            ∇²ₜibb,
            ∇²ₜibb!,
            a₁₂₃!,
            power,
            band_power
             # first on is for matrix version Dibb! is simplified
    using StaticArrays # static arrays are used for a mutable c=onstant
    const ħ = 1.054_571_817E-34::Float64 # J*s
    const C₁   = 1.191043E8::Float64#(1.191043E8,"W*μm*/m²*sr"," ","Risch","2016"),
    const C₂ = 14387.752::Float64#(14387.752,"μm**K"," ","Risch","2016"),
    const C₃ = 2897.77::Float64#(2897.77,"μm*K"," ","Risch","2016"),
    const C₄ = 4.09567E-12::Float64#(4.09567E-12,"W/m^2*μm*sr*K^5"," ", "Risch","2016"),
    const σ  = 5.670400E-8::Float64 #(5.670400E-8,"W/(m²*K⁴)"," ", "Risch","2016"),
    const Tₖ = 273.15::Float64 #(273.15,"K"," ", "Risch","2016");
    const a = SizedVector{3,Float64}([0.0; 0.0; 0.0]) # intermediate vector
#const a = [0.0; 0.0; 0.0]
"""
Fill intermediate results to construct the Planck function and its derivatives 
    
"""
function  a₁₂₃(λ::Float64,T::Float64)
    a1=C₂/(λ*T)
    a2 = 1/expm1(a1)#1/eaxpm1(a)
    a3 = exp(a1)*a2#exp(a)/expm1(a)
    return (a1,a2,a3)
end
function _a₁₂₃!(λ::Float64,T::Float64) # filling constant vector 
    # instance version of constants filling
    a[1]=C₂/(λ*T)
    a[2] = 1.0/expm1(a[1])#1/expm1(a)
    a[3] = exp(a[1])*a[2]#exp(a)/expm1(a)
    return a
end
function a₁₂₃!(amat::AbstractMatrix,λ::AbstractVector,T::Float64)
    # TRY VIEW WITH broadcating
        a1,a2,a3   = @views eachcol(amat)
        a1 .= (C₂/T)./λ
        a3 .=exp.(a1)
        a2 .= 1.0 ./(-1.0 .+a3)
        #a2 .= 1 ./expm1.(a1)
        a3 .= a3.*a2
    return amat
end

# BLACKBODY INTENCITY
""" 
    Blackbody intencity , [W/m2-sr-mkm]
    Ibb =(λ⁻⁵)* C₁/(exp(a) - 1), where 
    a₁=C₂/(λ*T)
    a₂ = 1/(eᵃ¹-1)   
    a₃ = eᵃ¹/(eᵃ¹-1) 
    Ibb = (λ⁻⁵)* C₁*a₂
    Possible input args:
        λ - wavelength in μm
        T - tmperature in Kelvins
        amat - matrix of intermediate coefficients
"""
    function ibb(λ,T) # general version close to the symbolic
        return (C₁/expm1(C₂/(λ*T)))*λ^-5
    end
    function ibb(λ::AbstractVector,T) # this version is useful in Optimization with AutoDiff?
        return map((l)->(C₁/expm1(C₂/(l*T)))*l^-5,λ)
    end
    function ibb(λ::AbstractVector,T::Base.RefValue{Float64}) # temeperature is a reference 
        return map((l)->(C₁/expm1(C₂/(l*T[])))*l^-5,λ)
    end
    function ibb(λ::AbstractVector,amat::AbstractMatrix) # internal version with provided coefficients matrix
        a2 = view(amat,:,2)
        return C₁*a2.*((1 ./λ).^5)     
    end
    function ibb(λ::AbstractVector,T::AbstractVector)
        #out = Matrix{Float64}(undef,length(λ),length(T))
        return length(T)==1 ? ibb(λ,T[1]) : @. ($C₁/expm1($C₂/(λ*$transpose(T))))*λ^-5
    end
    
    function ibb!(i::AbstractVector,λ::AbstractVector,T::Float64)
        map!(l->(C₁/expm1(C₂/(l*T)))*l^-5,i,λ)
        return i
    end
    function ibb!(i::AbstractVector,λ::AbstractVector,amat::AbstractMatrix)::Nothing # this version is used in emissivity approximation 
        a2 = view(amat,:,2) # Ibb = (λ⁻⁵)* C₁*a₂
        i.=C₁*a2.*((1 ./λ).^5) 
        return  nothing   
    end   

"""
    BB intencity derivative with respect to temperature
    dIbb/dT = C₁*(eᵃ¹/(eᵃ¹-1)²)*(C₂/(λ⁶*T²))
    a₁=C₂/(λ*T)
    a₂ = 1/(eᵃ¹-1)   #  1/expm1(a1)
    a₃ = eᵃ¹/(eᵃ¹-1) #  exp(a)/expm1(a)
    dIbb/dT = C₁*a₃*a₂*a₁*(1/(λ⁵*T))
    as far as Ibb = C₁*a₂/λ⁵
    dIbb/dT = a₃*a₁*C₁*(a₂/λ⁵)*(1/T)=a₃*a₁*Ibb/T
    Possible input args:
        λ - wavelength in μm
        T - tmperature in Kelvins
        amat - matrix of intermediate coefficients
        i - BB intencity calculated elsewhere
"""
    function ∇ₜibb(λ,T)
        prod(_a₁₂₃!(λ,T))*C₁/(T*λ^5)
    end
    function ∇ₜibb(λ::AbstractVector,T)
        return C₁*_a₁₂₃!.(λ,T)./(T*λ.^5)
    end
    function ∇ₜibb(λ::AbstractVector,T,amat::AbstractMatrix)
        return C₁*prod(amat;dims=2)./(T*λ.^5)
    end
    function ∇ₜibb!(g::AbstractMatrix,λ::AbstractVector,T::AbstractVector)
        # instance version of Planck function first derivative with respect to T
        for (jjj,t) in enumerate(T)
            for (iii,l) in enumerate(λ) 
                g[iii,jjj] = ∇ₜibb(l,t)
            end
        end   
    end
    function ∇ₜibb!(g::AbstractVector,λ::AbstractVector,T)# for fixed value of temperature
        # instance version of Planck function first derivative with respect to T
        #a = zeros(3)
        for (iii,l) in enumerate(λ) 
            g[iii] = prod(_a₁₂₃!(l,T))*C₁/(T*l^5)
        end 
    end
    # this version uses Vector for T because in this way the handle to the optimization varibale is implemented
    # T should be one-element array!
    function ∇ₜibb!(g::AbstractVector,λ::AbstractVector,T,amat::AbstractMatrix)# for fixed value of temperature
        # instance version of Planck function first derivative with respect to T
        #a = zeros(3)
        #t = T[1]
        prod!(g,amat) # this puts all amat rows with column-wise product in g
        g./=(T*λ.^5) # a₃*a₁*a₂*C₁*(1/λ⁵)*(1/T)
        g.*=C₁
        return nothing
    end
    function ∇ₜibb!(g::AbstractVector,T, amat::AbstractMatrix,i::AbstractVector)::Nothing
        #t = T[1]
        a1 = view(amat,:,1)
        a3 = view(amat,:,3)
        g .= a1.*a3.*i/T  # dIbb/dT = a₃*a₁*Ibb/T
        return nothing
    end
"""
    BB intencity second order derivative with respect to temperature
    d²Ibb/dT² = C₁*(eᵃ¹/(eᵃ¹-1)²)*(C₂/(λ⁶*T³))*[(C₂/(λ*T))*(2*eᵃ¹/(eᵃ¹-1)-1)-2]
    d²Ibb/dT² = C₁*(eᵃ¹/(eᵃ¹-1)²)*(a₁/(λ⁵*T²))*[a₁*(2*eᵃ¹/(eᵃ¹-1) -1)-2]
    a₁=C₂/(λ*T)
    a₂ = 1/(eᵃ¹-1)   #  1/expm1(a1)
    a₃ = eᵃ¹/(eᵃ¹-1) #  exp(a)/expm1(a)
    d²Ibb/dT² = C₁*a₂*a₃*(a₁/(λ⁵*T²))*[a₁*(2*a₃ - 1))-2]
    as far as 
    Ibb = (λ⁻⁵)* C₁*a₂
    dIbb/dT = C₁*a₃*a₂*a₁*(1/(λ⁵*T)) = a₃*a₁*Ibb/T 
    d²Ibb/dT² = C₁*a₂*a₃*a₁*(1/(λ⁵*T²))*[a₁*(2*a₃ - 1))-2] 
              = [a₃*a₁*Ibb/T^2]*[a₁*(2*a₃ - 1))-2] 
              = [(dIbb/dT)/T]*[a₁*(2*a₃ - 1))-2] 
    Possible input args:
              λ - wavelength in μm
              T - tmperature in Kelvins
              amat - matrix of intermediate coefficients
              i - BB intencity calculated elsewhere              
"""
    function ∇²ₜibb(λ,T)
        _a₁₂₃!(λ,T)#ibb = C₁*a[2]*(l^-5)
        (a[1]*(2a[3]-1.0)-2.0)*a[1]*a[2]*a[3]*C₁/((T^3)*λ^5)
    end
    function ∇²ₜibb!(h::AbstractVector{Float64},λ::AbstractVector{Float64},T::Float64)# secpnd derivative for the fixed value of temperature
        # instance version of Planck function second derivative with respect to T
        for (iii,l) in enumerate(λ) 
            _a₁₂₃!(l,T) # ibb = C₁*a[2]*(l^-5)         
            h[iii] = (a[1]*(2.0*a[3]-1.0)-2.0)*a[1]*a[2]*a[3]*C₁/((T^2)*l^5)
            # a₂*a₃*a₁*[a₁*(2*a₃ - 1))-2]*(C₁/(λ⁵*T²))
        end 
        return h
    end
    function ∇²ₜibb!(h::AbstractMatrix{Float64},λ::AbstractVector{Float64},T::AbstractVector{Float64})
        # instance version of Planck function second derivative with respect to T
        for (jjj,t) in enumerate(T)
            for (iii,l) in enumerate(λ) 
                _a₁₂₃!(l,t) # ibb = C₁*a[2]*(l^-5)         
                h[iii,jjj] = (a[1]*(2a[3]-1.0)-2.0)*a[1]*a[2]*a[3]*C₁/((t^2)*l^5) 
                # d²Ibb/dT² = C₁*a₂*a₃*a₁*(1/(λ⁵*T²))*[a₁*(2*a₃ - 1))-2] 
            end
        end      
    end
    function ∇²ₜibb!(h::AbstractVector{Float64},λ::AbstractVector{Float64},T::Float64,amat::AbstractMatrix{Float64})::Nothing
        # instance version of Planck function second derivative with respect to T
        # with supplied coefficints matrix
            a1 = view(amat,:,1)
            a3 = view(amat,:,3)
            prod!(h,amat) # h = a₃*a₁*a₂
            h./=(T^2)*λ.^5 # h = a₃*a₁*a₂*(1/λ⁵)*(1/T²)
            h.*=C₁ # h = C₁*a₃*a₁*a₂*(1/λ⁵)*(1/T²)
            h .*= a1.*(2.0*a3.-1) .-2.0  # C₁*a₂*a₃*a₁*(1/(λ⁵*T²))*[a₁*(2*a₃ - 1))-2] 
        return nothing        
    end
    function ∇²ₜibb!(h::AbstractVector{Float64},T::Float64,amat::AbstractMatrix{Float64},∇i::AbstractVector{Float64})::Nothing
        # instance version of Planck function second derivative with respect to T
        # with supplied coefficients matrix
            a1 = view(amat,:,1)
            a3 = view(amat,:,3)
            h.=∇i/T # h = (∇Ibb/T)
            h .*= a1.*(2.0*a3.-1.0) .-2.0  # h = (∇Ibb/T)*[a₁*(2*a₃ - 1))-2] 
        return nothing       
    end
    """
    The BB  with temperature T  spectral maximum 

    """
    function λₘ(T)
        # maximum wavelength of BB intencity in μm at temperature T (in Kelvins)
        C₃/T
    end
    """
     The temperature of BB with maximum wavelength λ
    """
    function tₘ(λ)
        # temperature (in Kelvins) of BB with intencity maximum at λ μm  
        C₃/λ
    end
"""
    BB intencity derivative with respect to wavelength

"""
    function ∇ₗibb(λ,T)
        # first derivative of Planck function with respect to wavelength
        #double a = C2/(lam*T);
        _a₁₂₃!(λ,T)
        (a[1]/λ)*(a[3]-5/λ)*(C₁*a[2]*(λ^-5)) #(C₁*a₁₂₃(λ,T)[2])*λ^-5
    end
    function  ∇²ₗibb(λ,T)
        # second derivative of Planck function with respect to wavelength
        #local a,e2,e3
        _a₁₂₃!(λ,T)
        C₁*a[2]*(a[1]*a[3]*(2a[1]*a[3]-a[1]-12)+30.0)/(λ^7)
    end
    function Dₗibb(λ,T)
        # methods returns PLanck function and its derivatives with respect to the wavelength
        # output is a tuple with (Planckfunction, Its first derivatiev with respect to the wavelength, Its second derivative with respect to the wavelength)
        _a₁₂₃!(λ,T)
        return (
            (C₁*a[2])*λ^-5,(a[1]/λ)*(a[3]-5/λ)*(C₁*a[2])*λ^-5, # first derivative
            (C₁/(λ^7) )*a[2]*(a[1]*a[3]*(2a[1]*a[3]-a[1]-12)+30.0) # second derivative
        )
    end
    function Dₗibb(λ::AbstractVector,T::AbstractVector)
        # returns spectral intencity and its first and second derivatives with respect to the wavelength
        i = fill(0.0,length(λ), length(T))
        d1i = fill(0.0,length(λ), length(T))
        d2i = fill(0.0,length(λ), length(T))
        for (jjj,t) in enumerate(T)
            for (iii,l) in enumerate(λ) 
                _a₁₂₃!(l,t)
                i[iii,jjj] = C₁*a[2]*(l^-5)
                d1i[iii,jjj] = (a[1]/l)*(a[3]-5/l)*i[iii,jjj]
                d2i[iii,jjj] = (i[iii,jjj]/(l^2))*(a[1]*a[3]*(2a[1]*a[3]-a[1]-12.0)+30.0)
                # (C₁/(λ^7))*a[2]*(a[1]*a[3]*(2a[1]*a[3]-a[1]-12)+30)
            end
        end
        return (i,d1i,d2i)
    end
    function power(T)
        # integral intencity of BB at temperature T
        return σ*(T^4)/π
    end
    """ 
        functionDₜibb! fills input tuple contains three vectors to be filled:
        (Planck_function_vector,Planck_function_first derivative,Planck_function_second_derivative)

    """ 
    # ::Tuple{Union{Array,Nothing,SizedArray},Union{Array,Nothing,SizedArray},Union{Array,Nothing,SizedArray}}
    # FILLS THE VALUE, FIRST AND SECOND DERIVATIVES
    function Dₜibb!(input_tuple::Tuple{AbstractVector,AbstractVector,AbstractVector}, λ::AbstractVector,T)
        # returns spectral intencity and its first and second derivatives with 
        # respect to the temperature
        # check if any of the output elements are ignored
        #for (iii,l) in enumerate(λ) 
        for (iii,l) in enumerate(λ) 
            _a₁₂₃!(l,T) #this function mutates global variable
            input_tuple[1][iii] = C₁*a[2]*(l^-5)
            input_tuple[2][iii] = a[1]*a[3]*input_tuple[1][iii]/T #a[1]*a[2]*a[3]*C₁/(T*l^5)
            input_tuple[3][iii] = (a[1]*(2a[3]-1.0) -2.0)*input_tuple[2][iii]/T# = [(dIbb/dT)/T]*[a₁*(2*a₃ - 1))-2] 
        end
        return input_tuple
    end
    # FILLS BOTH FIRST AND SECOND DERIVATIVES
    function Dₜibb!(input_tuple:: Tuple{nothing,AbstractVector,AbstractVector}, λ::AbstractVector,T)
        # returns spectral intencity and its first and second derivatives with 
        # respect to the temperature
        # check if any of the output elements are ignored
        #for (iii,l) in enumerate(λ) 
        for (iii,l) in enumerate(λ) 
            _a₁₂₃!(l,T) #this function mutates global variable
            input_tuple[2][iii] = a[1]*a[3]*C₁*a[2]*(l^-5)/T #a[1]*a[2]*a[3]*C₁/(T*l^5)
            input_tuple[3][iii] = (a[1]*(2a[3]-1.0) -2.0)*input_tuple[2][iii]/T# = [(dIbb/dT)/T]*[a₁*(2*a₃ - 1))-2] 
        end
        return input_tuple
    end
    # FILLS THE VALUE AND THE SE
    # FILLS ONLY SECOND DERIVATIVE
    function Dₜibb!(input_tuple:: Tuple{Nothing,Nothing,AbstractVector}, λ::AbstractVector,T::Float64)
        # returns spectral intencity and its first and second derivatives with 
        # respect to the temperature
        # check if any of the output elements are ignored
        #for (iii,l) in enumerate(λ) 
        ∇²ₜibb!(input_tuple[3],λ,T)
        return input_tuple
    end
    # FILLS ONLY VALUE
    function Dₜibb!(input_tuple:: Tuple{AbstractVector,Nothing,Nothing}, λ::AbstractVector,T::Float64)
        # returns spectral intencity and its first and second derivatives with 
        # respect to the temperature
        # check if any of the output elements are ignored
        #for (iii,l) in enumerate(λ) 
        ibb!(input_tuple[2],λ,T)
        return input_tuple
    end
    # FILLS ONLY FIRST DERIVATIVE
    function Dₜibb!(input_tuple::Tuple{Nothing,AbstractVector,Nothing},λ::AbstractVector,T::Float64)
        ∇ₜibb!(input_tuple[2],λ,T)
        return input_tuple
    end
    # HERE WE HAVE VERSIONS WHEN T is a Vector 
    function Dₜibb!(input_tuple::Tuple{Matrix{Float64},Matrix{Float64},Matrix{Float64}}, λ::AbstractVector,T::AbstractVector)
        # returns spectral intencity and its first and second derivatives with 
        # respect to the temperature
        # check if any of the output elements are ignored
        for (jjj,t) in enumerate(T)
            for (iii,l) in enumerate(λ) 
                _a₁₂₃!(l,t) 
                input_tuple[1][iii,jjj] = C₁*a[2]*(l^-5)
                input_tuple[2][iii,jjj] = a[1]*a[3]*input_tuple[1][iii,jjj]/t #a[1]*a[2]*a[3]*C₁/(T*l^5)
                input_tuple[3][iii,jjj] = (a[1]*(2a[3]-1.0) -2.0)*input_tuple[2][iii,jjj]/t^2#(a[1]*(2a[3]-1)-2)*a[1]*a[2]*a[3]*C₁/((T^3)*l^5)
            end
        end
        return input_tuple
    end
    function Dₜibb(λ::AbstractVector,T::AbstractVector)
        # returns spectral intencity and its first and second derivatives with respect to the temperature
        i = fill(0.0,length(λ), length(T))
        d1i = fill(0.0,length(λ), length(T))
        d2i = fill(0.0,length(λ), length(T))
        for (jjj,t) in enumerate(T)
            for (iii,l) in enumerate(λ) 
                _a₁₂₃!(l,t);
                i[iii,jjj] = C₁*a[2]*(l^-5)
                d1i[iii,jjj] = a[1]*a[3]*i[iii,jjj]/t #a[1]*a[2]*a[3]*C₁/(T*l^5)
                d2i[iii,jjj] = (a[1]*(2a[3]-1.0) -2.0)*d1i[iii,jjj]/(t^2)#(a[1]*(2a[3]-1)-2)*a[1]*a[2]*a[3]*C₁/((T^3)*l^5)
            end
        end
        return (i,d1i,d2i)
    end
    function band_power(T;λₗ=0.0,λᵣ=Inf,tol=1e-6)
        return power(T)*∫ibbₗ(T;λₗ=λₗ,λᵣ=λᵣ,tol=tol)
    end
    function ∫ibbₗ(T;λₗ=0.0,λᵣ=Inf,tol=1e-6)
        # calculates the integral of spectral intencity over the wavelength
        if λₗ>=λᵣ
            throw("The left wavelength boundary should be smaller than the right on");
        end
        if ~isfinite(λᵣ)# the right boundary is infinite
            if λₗ==0.0
                return 1
            else #integration from fixed wavelength to infinity
                return ∫ibbₗ(T)- ∫ibbₗ(T,λᵣ=λₗ) 
            end
        else# righ wavelength boundary is finite
            if λₗ==0.0# integration from zero to fixed wavelength
                n=1
                ϵ=tol*100
                sum=0
                a = _a₁₂₃!(λᵣ,T)[1];
                while  (ϵ>tol)&&(n<1e4) 
                    etan = a*n
                    sum+=(exp(-etan)/n)*(etan*(etan*(etan + 3) + 6) + 6)/(n^4)
                    n+=1;
                end
                return 15*sum/(pi^4)
            else# both wavelength resions are limited 
                return ∫ibbₗ(T,λᵣ=λᵣ) - ∫ibbₗ(T,λᵣ=λₗ)
            end
        end
    end
end