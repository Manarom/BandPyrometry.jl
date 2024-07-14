


include("JDXreader.jl")
include("PLanck.jl")
include("BBSource.jl")
include("BandPyrometry.jl")
module Epoint
        using Interpolations
        import ..JDXreader,
        ..Planck, 
        ..BBSource  
        include("EmissivityPointTypes.jl")
        const ConfigurationNames = Dict(
            "AP_sample" => SPD((:raw_Sₛ,# FTIR signal in AP_sample - sample raw spectrum
                                :dark_Wₛ, # PM signal in AP_sample - sample power dark signal
                                "aperture, sample radiation directed to the FTIR"
            )),
            "AP_ethalon" =>SPD((:raw_Sₑ, # FTIR signal in AP_ethalon - ethalon raw spectrum
                                :dark_Wₑ, # PM signal in AP_ethalon - ethalon power dark signal
                                "aperture, ethalon radiation directed to the FTIR"
            )),
            "AB_sample" => SPD((:dark_Sₛ,# FTIR signal in AB_sample - sample dark spectrum
                                :raw_Wₛ,# PM signal in AB_sample - sample power signal
                                "aperture, sample radiation directed to the  PM"
            )),
            "AB_ethalon" => SPD((:dark_Sₑ, # FTIR signal in AB_ethalon - ethalon dark signal
                                 :raw_Wₑ,  # PM signal in AB_ethalon - ethalon power
                                "aperture, ethalon radiation directed to the PM"
            )),
            "SP_sample" => SPD((:cold_raw_Sₛ,# FTIR signal in SP_sample - cold sample raw spectrum (reflected flux)
                                :cold_dark_Wₛ,# PM signal in SP_sample - cold sample power dark signal 
                                "aperture, sample radiation directed to FTIR (cold sample)")),
            "SP_ethalon" => SPD((:cold_raw_Sₑ, # FTIR signal in SP_ethalon - ethalon raw spectrum (cold sample)
                                 :cold_dark_Wₑ,# PM signal in SP_ethalon - ethalon power dark signal (cold sample)            
                                "aperture, ethalon radiation directed to FTIR (cold sample)"
            )),
            "BL_sample" => SPD((:cold_dark_Sₛ,# FTIR signal in BL_sample - cold sample dark spectrum
                                :cold_raw_Wₛ,# PM signal in BL_sample - cold sample power signal 
                                "aperture, sample radiation directed to PM (cold sample)"
            )),
            "BL_ethalon" => SPD((:cold_dark_Sₑ, # FTIR signal in AB_ethalon - ethalon dark signal
                                 :cold_raw_Wₑ,  # PM signal in AB_ethalon - ethalon power =#
                                 "aperture, ethalon radiation directed to PM (cold sample)"
            ))
        )
        const DarkDictSignal = Dict(
                            :Sₑ => (:raw_Sₑ,:dark_Sₑ),                            
                            :Sₛ => (:raw_Sₛ,:dark_Sₛ),                            
                            :cold_Sₑ => (:cold_raw_Sₑ,:cold_dark_Sₑ),
                            :cold_Sₛ => (:cold_raw_Sₛ,:cold_dark_Sₛ)                          
        ) 
        const DarkDictPower = Dict(
                            :Wₑ => (:raw_Wₑ,:dark_Wₑ),
                            :Wₛ =>  (:raw_Wₛ, :dark_Wₛ),
                            :cold_Wₛ => (:cold_raw_Wₛ, :cold_dark_Wₛ), 
                            :cold_Wₑ => (:cold_raw_Wₑ,:cold_dark_Wₑ) 
        )                                     
        """
            Stores all raw and calculated data for particular point

        """
        mutable struct EmissionPoint
            # project is a folder name
            # strores data of emissivity spectrum
            # RAW data right from the file
            # device configurations:
            #                       AP_sample - aperture, sample radiation directed to the FTIR
            #                       AP_ethalon - aperture, ethalon radiation directed to the FTIR
            #                       AB_sample  - aperture, sample radiation directed to the  PM (power meter)
            #                       AB_ethalon - aperture, ethalon radiation directed to the PM
            #                       SP_sample - aperture, sample radiation directed to FTIR (cold sample)
            #                       SP_ethalon - aperture, ethalon radiation directed to FTIR (cold sample)
            #                       BL_sample - aperture, sample radiation directed to PM (cold sample)
            #                       BL_ethalon - aperture, ethalon radiation directed to PM (cold sample)
            # MEASURED DATA - raw data without any interpolation 
            # AP_sample configuration:
                raw_Sₛ::XYdata # FTIR signal in AP_sample - sample raw spectrum
                dark_Wₛ::Float64 # PM signal in AP_sample - sample power dark signal
            # AP_ethalon configuration:
                raw_Sₑ::XYdata # FTIR signal in AP_ethalon - ethalon raw spectrum
                dark_Wₑ::Float64 # PM signal in AP_ethalon - ethalon power dark signal
            # AB_sample configuration:
                dark_Sₛ::XYdata # FTIR signal in AB_sample - sample dark spectrum
                raw_Wₛ::Float64 # PM signal in AB_sample - sample power signal
            # AB_ethalon configuration
                dark_Sₑ::XYdata # FTIR signal in AB_ethalon - ethalon dark signal
                raw_Wₑ::Float64 # PM signal in AB_ethalon - ethalon power
            # SP_sample sonfiguration 
                cold_raw_Sₛ::XYdata # FTIR signal in SP_sample - cold sample raw spectrum (reflected flux)
                cold_dark_Wₛ::Float64# PM signal in SP_sample - cold sample power dark signal     
            # SP_ethalon configuration:
                cold_raw_Sₑ::XYdata # FTIR signal in SP_ethalon - ethalon raw spectrum (cold sample)
                cold_dark_Wₑ::Float64# PM signal in SP_ethalon - ethalon power dark signal (cold sample)           
            # BL_sample configuration:
                cold_dark_Sₛ::XYdata# FTIR signal in BL_sample - cold sample dark spectrum
                cold_raw_Wₛ::Float64# PM signal in BL_sample - cold sample power signal                
            # BL_ethalon configuration
                cold_dark_Sₑ::XYdata # FTIR signal in AB_ethalon - ethalon dark signal
                cold_raw_Wₑ::Float64 # PM signal in AB_ethalon - ethalon power
            # Spectra and power interpolated after dark signal extraction
                λ_interp::Vector{Float64} # wavelengths should be the same for all data and spectral points
                # Fields to be filled 
                Sₛ::Vector{Float64} # sample spectrum after extraction of the dark signal and interpolation
                Wₛ::Float64 # sample power after extraction of the dark signal
                Sₑ::Vector{Float64} # ethalon spectrum after extraction of the dark signal and interpolation
                Wₑ::Float64 # ethalon power after extraction of the dark signal
                cold_Sₛ::Vector{Float64} # cold sample spectrum after extraction of the dark signal and interpolation
                cold_Wₛ::Float64 # cold sample power after extraction of the dark signal
                cold_Sₑ::Vector{Float64} # ethalon spectrum after extraction of the dark signal and interpolation (cold sample)
                cold_Wₑ::Float64 # ethalon power after extraction of the dark signal  (cold sample)
            # Flags
                has_cold_signal::Bool # true if we use the cold signal (or just as is)
                has_ethalon_signal::Bool # true if we use the ethalon spectrum , otherwise all ethalon raw data 
                # remains unfefined and Cᵢⱼ should be provided externally
                is_ready_to_fit_flag::Bool
                is_all_ready_flag::Bool


            #Correction spectra
            Cᵢⱼ::Vector{Float64} # instrumental function of the unit
            cold_Cᵢⱼ::Vector{Float64} # instrumental function of the unit with cold sepctrum correction
            RT_emissivity::XYdata
            ϵ₀::Vector{Float64} # room-temperature spectral emissivity
            ϵₑ::Vector{Float64} # ethalon spectral emissivity
            Iₛᵤᵣ::Vector{Float64} # external radiation vector

            Yᵢⱼ::Vector{Float64} # measured signal to be fitted using spectral pyrometry

            ϵ_int::Float64 # integral emissivity    
            T_header::String # data header (measurement temperature)
            T_sample::Float64# sample temperature measured by the pyrometer            
            T_ethalon::Float64

            band_pyrometry_range::Vector{Float64} # band pyrometry range
            T_real::Float64 # temperature fitted from the emission spectrum
            # CONSTRUCTORS
            EmissionPoint()=new()
        end
        
        """
    check_files(dir_name::String,title_temp::String;
            has_cold_signal::Bool=false,file_extention::String=".jdx")::file_checher_return_type

            Return NamedTuple with the following fields:
                        :is_ready::Bool  - flag is point ready to be  loaded (all necessary files are in folder)
                        :full_file_names::Vector{String} - full file names of files to be loaded

"""
    function check_files(dir_name::String,
            title_temp::String;
            has_cold_signal::Bool=false,
            file_extention::String=".jdx",
            has_ethalon_signal::Bool=false)::file_checher_return_type
            D = ConfigurationNames
            if !has_cold_signal
                D = defilter_dict(D,"SP","BL")
            end
            if !has_ethalon_signal
                D = defilter_dict(D,"ethalon")
            end
            dict_keys_Set = keys(D)
            if !isdir(dir_name)
                return file_checher_return_type((false,
                    Vector{String}(),
                    Vector{String}(),
                    dict_keys_Set,
                    ))
            end
            file_name_tail = "_"*title_temp*file_extention
            # filters all files in dir containing the title_temp
            file_names = filter( fl->contains(fl,title_temp),readdir(dir_name)) 
            # creating set from file names, we should check if all files we need to fill are 
            file_names_Set = Set(replace.(file_names,file_name_tail => ""))
            # keys which we need to fully create the 
            keys_with_no_files_set = filter(ky->!in(ky,file_names_Set),dict_keys_Set)
            keys_with_files_set = intersect(dict_keys_Set,file_names_Set)
            keys_with_files_vec = collect(keys_with_files_set)
            return  file_checher_return_type((
                                isempty(keys_with_no_files_set),
                                joinpath.(dir_name,keys_with_files_vec .* file_name_tail),
                                keys_with_files_vec,
                                keys_with_no_files_set
                                ))
        end

    """

    EmissionPoint(dir_name::String,title_temp::String;
                  has_cold_signal::Bool=false,
                  pyrometer_index::Int=1,
                  has_ethalon_signal::Bool=true)

    External Constructor for the emission point,
        dir_name is the project directory name
        title_temp is the temperature of the current point which should be included in 
        file names, 
        (optional):
        has_cold_signal - if true cold signal was measured for this point, 
                program searches the directory for the following files:
                    ["SP_ethalon_", "SP_sample","BL_ethalon","BL_sample"].*title_temp
        pyrometer_index - int from 1 to 3 which of three pyrometer signal is considered as the 
                    temperature of a sample
        has_ethalon_signal - if true seraches for files containing "_ethalon_", if false it is 
                    supposed that calubration spectrum is provided externaly 

"""
        function EmissionPoint(dir_name::String,title_temp::String;
            has_cold_signal::Bool=false,pyrometer_index::Int=1,
            has_ethalon_signal::Bool=false)
            em_point = EmissionPoint()
            em_point.T_header = title_temp
            em_point.is_all_ready_flag = false
            em_point.is_ready_to_fit_flag = false
            check_files_res = check_files(dir_name,title_temp,
                                            has_cold_signal=has_cold_signal,
                                            has_ethalon_signal=has_ethalon_signal)
 #=         file_checher_return_type = NamedTuple{
                    (:is_ready,:full_file_names,:checked_tags,:needed_tags),
                    Tuple{Bool,Vector{String},Vector{String},AbstractSet}
            } =#
            em_point.is_all_ready_flag = check_files_res.is_ready
            em_point.has_cold_signal = has_cold_signal
            em_point.has_ethalon_signal = has_ethalon_signal
            if !check_files_res.is_ready
                @warn "Please provide the file for: "*string(collect(check_files_res.needed_file_names))
                return (check_files_res.is_ready,em_point)
            end
            for (file,tag) in zip(check_files_res.full_file_names,check_files_res.checked_tags) 
                _fill_raw_data(em_point,file,tag)
            end
            selected_pyrometer = [:pyr_T1,:pyr_T2,:pyr_T3][pyrometer_index]
            em_point.T_sample = getfield(em_point.raw_Sₛ.state,selected_pyrometer) + Planck.Tₖ
            if has_ethalon_signal
                em_point.T_ethalon = em_point.raw_Sₑ.state.m330_T+ Planck.Tₖ # all temperatures should be in Kelvins!!!
            else
                em_point.T_ethalon = -1.0
            end

            return (check_files_res.is_ready,em_point)
        end

        """
        _fill_raw_data(e::EmissionPoint,full_file::String,tag::String)

        Internal function fills raw (initial data by ConfigurationNames key)
        no checks version 
"""
        function _fill_raw_data(e::EmissionPoint,full_file::String,tag::String)
            # (x=jdx.x_data,
            # y=jdx.y_data, 
            # headers=jdx.data_headers) x_units are in 1/cm-1 by default
            conf_names = ConfigurationNames[tag] # returns SPD type NamedTuple 
                                                    # SPD = NamedTuple{
                                                    #(:signal,:power,:description),
                                                    # Tuple{Symbol,Symbol,String}}
            setfield!(e,
                conf_names.signal,
                makeXYdata_from_file(full_file))
            pow = getfield(e,conf_names.signal).state.power
            setfield!(e,
                conf_names.power,
                pow
            ) # set the value of power to the power field of e-struct
            return nothing
        end

        """
        add_RT_emissivity(e::EmissionPoint,em_source::Union{String,AbstractMatrix})

        Adds RT emissivity to the em::EmissionPoint, em_source can be full file name or
        matrix of float
"""
        function add_RT_emissivity(em::EmissionPoint,em_source::Union{String,Matrix{Float64}})
            if em_source isa String
                em.RT_emissivity = makeXYdata_from_file(em_source)
                
            else
                em.RT_emissivity = XYdata((em_source[:,1],
                em_source[:,1],
                Dict("COMMENT"=>"from matrix"),
                InstrumentState("")))
            end
        end
    """
        evaluate_correction(e::EmissionPoint,
                λ_interp::Vector{Float64},
                m330::BBSource.M330;
                Cij::Vector{Float64}=[-1],Iₛ::Vector{Float64}=[-1])

        Internal function, evaluates correction vector of the emission point e, 
                λ_interp - wavelength interpolation vector (for all spectra) 
                m330 - BlackBody reference
            (Optional):
                Cij - optional correction signal (this can be used in case if we have only one correction function 
                for all temperature points)
                cold_Cij - optional cold correction signal
                Iₛ
"""
        function evaluate_correction(e::EmissionPoint,
                λ_interp::Vector{Float64},
                m330::BBSource.M330;
                Cij::Vector{Float64}=[-1.0],
                Iₛ::Vector{Float64}=[-1.0])
            e.λ_interp = λ_interp 
            has_RT_emissivity = isdefined(e,:RT_emissivity)
            if has_RT_emissivity
                e.ϵ₀ = e.RT_emissivity(λ_interp,out_units="") # wavelength of RT emissivity 
                #   should be in microns
            end
            _extract_dark_signal(e) # dark signal extraction
            if isnothing(m330) # if ethalon object is not provided, than we suppose the emissivity is 1
                e.ϵₑ=ones(length(e.λ_interp))
            else
                if e.has_ethalon_signal
                    e.ϵₑ = BBSource.interp(m330,e.T_ethalon,e.λ_interp) # evaluates interpolated 
                    # BBSource spectral emissivity
                else # if we dont have ethalon spectra, emissivity is calculated from the sample temperature
                    e.ϵₑ = BBSource.interp(m330,e.T_sample,e.λ_interp) # evaluates interpolated 
                    # BBSource spectral emissivity
                end
            end
            if !e.has_cold_signal &&  length(Iₛ)!=length(λ_interp) 
                # if does  not has cold signal 
                # and surrounding radiation 
                # correction is provided externally
                e.Iₛᵤᵣ=zeros(length(λ_interp))  
            elseif !e.has_cold_signal
                e.Iₛᵤᵣ = Iₛ
            else 
                # has cold signal data
                Tbb_cold = e.cold_raw_Sₑ.state.m330_T + Planck.Tₖ 
                # temperature of the cold ethalon spectrum
                e.cold_Cᵢⱼ=e.cold_Sₑ./(e.ϵₑ.*Planck.ibb(λ_interp,Tbb_cold))
                if has_RT_emissivity
                    e.Iₛᵤᵣ = (e.cold_Sₛ./e.cold_Cᵢⱼ - e.ϵ₀.*Planck.ibb(λ_interp,300.0))./(1-e.ϵ₀)
                else
                    e.Iₛᵤᵣ = Planck.ibb(λ_interp,300.0) # if we dont know the RT emissivity we just 
                end

            end 
            if !e.has_ethalon_signal && length(λ_interp)!=length(Cij)
                error("Current point marked as 
                        not having self ethalon spectrum, 
                        thus Cij sould be provided externally")
            elseif !e.has_ethalon_signal
                e.Cᵢⱼ = Cij
            else # emissivity point has ethalon
                 # signal and intrument function should be calculated for this particular point
                e.Cᵢⱼ = e.Sₑ./(e.ϵₑ .* Planck.ibb(λ_interp,e.T_ethalon))
            end
            e.Yᵢⱼ = e.Sₛ./e.Cᵢⱼ - e.Iₛᵤᵣ
            return e
        end
        """
        _extract_dark_signal(e::EmissionPoint)

            Internal function. Extracts dark signals from the measured signal, interpolates the results 
            and fills the corresponding fields of EmissionPoint struct

        
"""        
        function _extract_dark_signal(e::EmissionPoint)
            DS = DarkDictSignal
            DP = DarkDictPower
            if !e.has_cold_signal # if EmissionPoint does not has cold Measurements
                # we need to filter fields stored in fields dictionaries
                DS = defilter_dict(DS,"cold")
                DP = defilter_dict(DP,"cold")
            end
            if !e.has_ethalon_signal # remove all ethalon spectra from the dictionary
                DS = defilter_dict(DS,"Sₑ")
                DP = defilter_dict(DP,"Wₑ")                   
            end
            for d in DS # for signals we need to interpolate and extract
                setfield!(e,d[1], .-(
                                    interp(e,getfield(e,d[2][1])), 
                                    interp(e,getfield(e,d[2][2])))
                )                            
            end
            for d in DP # for power we need just to extract on from another
                # dictionary iterator returns (key-value) pair
                setfield!(e,d[1],-(
                                    getfield(e,d[2][1]),
                                    getfield(e,d[2][2])
                                    )
                )
            end
            return e
        end
        """

        Interpolates data in the form of Matrix at the point EmissionPoint.λ_interp

        """
        function interp(e::EmissionPoint,data::XYdata)
            return data(e.λ_interp) # flat extrapolation gives const for extrapolation
        end
end