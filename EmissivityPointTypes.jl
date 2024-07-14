# Script stores type declaration and utilities to work with this types
"""
Single spectra return type (internal)

"""
const SPD = NamedTuple{(:signal,:power,:description),
Tuple{Symbol,Symbol,String}}
"""

Type return by the file_checker function 

"""
const file_checher_return_type = NamedTuple{
    (:is_ready, # is file ready to be loaded
    :full_file_names, # Vector of full file names
    :checked_tags, # Vector of checked file type tags
    :needed_tags # Tags which we need to add
    ),
    Tuple{Bool,Vector{String},Vector{String},AbstractSet}
}
"""
Stores instrument position and measured temperatures

"""
struct InstrumentState
    pyr_T1::Float64
    pyr_T2::Float64
    pyr_T3::Float64
    m330_T::Float64
    power::Float64
    mir1_rotation::Float64
    mir2_position::Float64
    sample_coordinate::Float64
    InstrumentState(state_string::String) = begin
        data_intermediate = map((X)->Base.tryparse(Float64,X), split(state_string))[2:end]
        if length(data_intermediate)!=length(fieldnames(InstrumentState))
            return new(-1,-1,-1,-1,-1,-1,-1,-1)
        end
        if any(isnothing.(data_intermediate))
            return new(foreach(x-> isnothing(x) ? -1 : x
                            ,data_intermediate)...)
        else
            return new(data_intermediate...)
        end    
#=              
        code from Qt project
        oldData.insert(index, tr("\t%1\t%2\t%3\t%4\t%5\t%6\t%7\t%8")
                        .arg((s_t1 / iter), 0, 'f', 3)
                        .arg(s_t2 / iter) 
                        .arg(s_t3 / iter) pyrometer_temperature
                        .arg(tMode)  - m330 mode temperature
                        .arg(s_pw / iter) power
                        .arg(d_func()->motion->getCurrentPositions().value(1), 0, 'f', 1)
                        .arg(d_func()->motion->getCurrentPositions().value(2), 0, 'f', 1)
                        .arg(d_func()->motion->getCurrentPositions().value(3), 0, 'f', 1));
=#
    end
end
""" 
    main type to store two column data
"""
const XYdata = NamedTuple{(:x,:y,:headers,:state),Tuple{
                                            Vector{Float64},
                                            Vector{Float64},
                                            Dict{String,Union{String,Number}},
                                            InstrumentState
                                            }}

"""
    (xydata::XYdata)(x_interp::Vector{Float64})
    Calling xydata(x_interp) returns interpolated xydata.y

"""
function (xydata::XYdata)(x_interp::Vector{Float64};out_units::String="MKM")
    if out_units!="MKM"
        linear_interpolation(xydata.x,xydata.y,extrapolation_bc = Flat())(x_interp)
    end
    xunits = get(xydata.headers,"XUNITS","NONE")
    if (xunits=="1/CM")
        if issorted(xydata.x) # ascending order
            first_ind = findfirst(i->(xydata.x[i]>0),eachindex(xydata.x)) # first non-zero ind
            x_data = 10000.0 ./reverse(xydata.x[first_ind:end])
            y_data = reverse(xydata.y[first_ind:end])
            return linear_interpolation(x_data,y_data,extrapolation_bc = Flat())(x_interp)
        else
            return linear_interpolation(xydata.x,xydata.y,extrapolation_bc = Flat())(x_interp)
        end
    else
        return linear_interpolation(xydata.x,xydata.y,extrapolation_bc = Flat())(x_interp)
    end
    
end

"""
    (xydata::XYdata)(key_name::String)
    Returns headers element attached to the xydata
"""
function (xydata::XYdata)(key_name::String)
    return get(xydata.headers,key_name,"")
end
"""
    makeXYdata_from_file(xydata::XYdata,full_file::String)

    Fills XYdata obj from file 

"""
function makeXYdata_from_file(full_file::String)
    #reader_out = JDXreader.read_jdx_file(full_file) # returns named tuple 
    # (x=jdx.x_data,
    # y=jdx.y_data, 
    # headers=jdx.data_headers) x_units are in 1/cm-1 by default
    try
    reader_out= JDXreader.read_jdx_file(full_file) 
        return XYdata(
            (reader_out...,InstrumentState(
                            get(reader_out.headers
                                    ,"TITLE",""
                                )#endof "TITLE" header getter
                                )#endof InstrumentsState constructor
                                )#endof tuple
        )# endof XYdata condtructor
    catch
        @warn "Unable to open file as JCAMP"
        return XYdata(([-1.0], [-1.0],Dict{String,Float64}(),InstrumentState("")))
    end
end



"""
    (xydata::XYdata)(state_field::Symbol)

    when calling xydata(Symbol) returns instruments state position value
"""
function (xydata::XYdata)(state_field::Symbol)
    return  ∈(state_field,fieldnames(InstrumentState)) ? 
                    getfield(xydata.state,state_field) :
                    NaN64
end
"""
    defilter_dict(D::AbstractDict{T,V} where {V,T<:Union{Symbol, String}}, 
    filter_by::Union{Symbol, String}... )

    Simple function which return subdictionary of D where  elements containing  filter_by
    are removed 
"""
function defilter_dict(D::AbstractDict{T,V} where {V,T<:Union{Symbol, String}}, 
    filter_by::Union{Symbol, String}... ) 
    keys_ = keys(D) # returns keys iterator
    for f in filter_by
        keys_ = filter(k-> !contains(string(f))(string(k)),keys_)
    end
    return Dict(zip(keys_,[D[k] for k in keys_])) 
end     