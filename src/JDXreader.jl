
# module to read JCAMP-DX=4.24 file formats

module JDXreader
using Dates,Interpolations
"""
JCAMP file content example: 

            ##TITLE=1 
            ##JCAMP-DX=4.24
            ##DATATYPE=INFRARED SPECTRUM
            ##DATE=2021/10/17
            ##TIME=11:30
            ##XUNITS=1/CM
            ##YUNITS=TRANSMITTANCE
            ##XFACTOR=1
            ##YFACTOR=0.00699183
            ##FIRSTX=0
            ##LASTX=15801.4671743
            ##FIRSTY=11746.3412893072
            ##MAXX=15801.4671743
            ##MINX=0
            ##MAXY=3.75371e+006
            ##MINY=4040.25
            ##NPOINTS=16384
            ##XYDATA=(X++(Y..Y))
            0 1680010 821286 2148133 1505245 1537124 1367661 1147725 1134981
            7.71603 1166853 1213186 1029828 1067595 1135904 1195128 1157134 1150556
            15.4321 1266743 1164401 1014224 1022338 999780 1138781 1208950 1161258
            .
            .
            ##END

Header lines start with ##
Data start after headers lines

"""
    JDXreader
    const YMAX_INT = round(Int64,typemax(Int32)/4)-1
    const default_headers = Dict{String,Union{Float64,String}}(
            "TITLE"=>"NO TITLE",
            "JCAMP-DX"=>4.24,
            "DATATYPE"=>"INFRARED SPECTRUM",
            "DATE"=>"2021/10/17",
            "TIME"=>"11:30",
            "XUNITS"=>"1/CM",
            "YUNITS"=>"TRANSMITTANCE",
            "XFACTOR"=>1.0,
            "YFACTOR"=>0.00699183,
            "FIRSTX"=>0.0,
            "LASTX"=>15801.4671743,
            "FIRSTY"=>11746.3412893072,
            "MAXX"=>15801.4671743,
            "MINX"=>0.0,
            "MAXY"=>3.75371e+006,
            "MINY"=>4040.25,
            "NPOINTS"=>16384.0,
            "XYDATA"=>"(X++(Y..Y))"
    )
    const sym_default_dict = Dict(Symbol(k[1])=>k[2] for k in default_headers)
    const default_keys = collect(default_headers)
    export JDXfile,read!,read_jdx_file,write_jdx_file

    function fill_headers(x::Vector{Float64},y::Vector{Float64}; kwargs...)
        @assert length(x)==length(y)

        headers = copy(default_headers)
        for (k,v) in kwargs
            k_str = string(k)
            if k_str == "XUNITS"
                 @assert v=="1/CM" "x-values should be in 1/CM, other units are unsupported"  
            end
            headers[k_str] = v 
        end
        x_factor = headers["XFACTOR"] 
        n_points = length(x)
        if is_linspaced(x)# checks if all coordinates are equally spaced        
            x_copy = copy(x)
            x./=x_factor
            y_copy = copy(y)
            if !issorted(x_copy)
                y_int = sortperm(x)
                @. x=x[y_int]
                @. y_copy=y_copy[y_int]
            else
                y_int = Vector{Int}(undef,n_points) 
            end
        else # if x is not equally spaced we perform linear interpolation
            x_copy = collect(range(minimum(x),maximum(x),n_points)) 
            if !issorted(x)
                y_int = sortperm(x)
                y_copy = linear_interpolation(x[y_int],y[y_int])(x_copy)
            else
                y_int = Vector{Int}(undef,n_points)
                y_copy = linear_interpolation(x,y)(x_copy)
            end
            x_copy./=x_factor
        end
        y_columns_number = 1
        cur_date_time = string(now())
        ind = findfirst("T",cur_date_time)[1]
        headers["NPOINTS"] = Float64(n_points)
        headers["DATE"] = cur_date_time[1:ind-1]
        headers["TIME"] = cur_date_time[ind+1:end]
        headers["FIRSTX"] = x_copy[begin]
        headers["FIRSTY"] = y_copy[begin]
        (headers["MINY"],headers["MAXY"]) = extrema(y_copy)
        (headers["MINX"],headers["MAXX"]) = extrema(x_copy)
        headers["FIRSTX"] = headers["MINX"]
        headers["LASTX"] = headers["MAXX"]
        y_factor = (headers["MAXY"]/YMAX_INT)
        @. y_int = round(Int,y_copy/y_factor)
        headers["YFACTOR"] = y_factor
        return (x_copy,y_copy,y_int,headers,y_columns_number)
    end
    function is_linspaced(x::Vector{T}) where T<:Number
        if length(x)<=2
            return true
        end
        dx1 = x[2] - x[1]
        for i in eachindex(x)[2:end-1]
           dx = x[i+1] - x[i]
           isapprox(dx,dx1,rtol=1e-8) ? continue : return false
        end
        return true
    end
    mutable struct JDXfile
        # Main struct, prepares loads file name, parses data and data headers
        file_name::String
        data_headers::Dict{String,Union{String,Float64}}
        x_data::Vector{Float64}
        y_data::Vector{Float64}
        is_violated_flag::Vector{Bool}
        """
        JDXfile()
        JDXreader obj constructor, creates empty object with no data
    """
        JDXfile()=begin
            new("",
                Dict{String,Union{String,Float64}}(),
                Vector{Float64}(),
                Vector{Float64}(),
                Vector{Bool}())
        end
    end
    """
    read_jdx_file(file_name::String)

    Read JCAMP format file file_name - full file name,
    Input: 
        file_name - full file name
    returns named tuple with fields :
           x - coordinate (wavelength, wavenumber or other)
           y - data
           headers - dictionary in "String => value" format with 
           headers values  
"""
    function read_jdx_file(file_name::String)
        return file_name |> JDXfile |> read!
    end
    """
    JDXfile(file_name::String)

    Creates JDXreader object from full file name 
"""
    function JDXfile(file_name::String) # external constructor
        if !isfile(file_name)
            error("Input filename must be the fullfile name to the existing file")
        end
        jdx = JDXfile()
        jdx.file_name = file_name
        return jdx
    end
    function parseJDXheaders(jdx::JDXfile,headers::Vector{String})
        for head in headers
            name = strip(head,'#')
            name = strip(name,'$')
            splitted_string = split(name,"=")
            if length(splitted_string)!=2
                jdx.data_headers[name] = name
                continue
            end
            jdx.data_headers[string(splitted_string[1])] = 
                isnothing(tryparse(Float64,splitted_string[2])) ? string(splitted_string[2]) : parse(Float64,splitted_string[2])
        end
        #if haskey(jdx.data_headers,"NPOINTS")
        #    jdx.data_headers["NPOINTS"] = round(Int64, jdx.data_headers["NPOINTS"])
        #end
    end
    function addYline!(jdx::JDXfile, current_line::String,number_of_y_point_per_chunk,chunk_index)
        data_intermediate = map((X)->Base.parse(Float64,X), split(current_line))
        starting_index =1+ (chunk_index-1)*number_of_y_point_per_chunk
        jdx.y_data[starting_index:starting_index+number_of_y_point_per_chunk-1] = data_intermediate[2:end]
        return data_intermediate[1]
    end
    function addYline!(jdx::JDXfile, current_line::String)
        data_intermediate = map((X)->Base.parse(Float64,X), split(current_line))
        resize!(jdx.y_data,0)
        append!(jdx.y_data,data_intermediate[2:end])
        return data_intermediate[1]
    end
    function generateXvector!(jdx::JDXfile)
            point_number = round(Int64,jdx.data_headers["NPOINTS"])
            starting_X = haskey(jdx.data_headers,"FIRSTX") ? jdx.data_headers["FIRSTX"] : 0.0
            if haskey(jdx.data_headers,"DELTAX")
                step_value =  jdx.data_headers["DELTAX"]  
            else
                if haskey(jdx.data_headers,"XFACTOR")
                    step_value =  ((jdx.data_headers["LASTX"] -  starting_X)/jdx.data_headers["XFACTOR"])/(point_number-1)
                else
                    step_value =  (jdx.data_headers["LASTX"] -  starting_X)/(point_number-1)
                end
            end
            resize!(jdx.x_data,point_number)
            jdx.x_data .= [starting_X + i*step_value for i in 0:point_number-1]
    end
    """
    read!(jdx::JDXfile)

fills precreated JDXfile object

"""
function read!(jdx::JDXfile)
        if !isfile(jdx.file_name)
            return nothing
        end
        header_lines = Vector{String}()
        total_number_Of_lines = countlines(jdx.file_name)
        x_point =0.0;
        jdx.y_data = Vector{Float64}()
        open(jdx.file_name) do io_file
            x_point=0.0
            for ln in eachline(io_file)
                if ln[1]!='#'
                    x_point = addYline!(jdx, ln)
                    break
                else
                    push!(header_lines,ln)
                end
            end
            number_of_y_point_per_chunk = length(jdx.y_data)
            parseJDXheaders(jdx,header_lines)
            data_lines_number = total_number_Of_lines - length(header_lines)-1
            if haskey(jdx.data_headers,"NPOINTS") # correct JDX file
                total_point_number = round(Int64,jdx.data_headers["NPOINTS"])
                generateXvector!(jdx)
                resize!(jdx.y_data,total_point_number)
                for i in 2:data_lines_number
                    x_point = addYline!(jdx,readline(io_file),number_of_y_point_per_chunk,i)
                end
            else number_of_y_point_per_chunk==1 # file with two columns like CSV
                total_point_number = data_lines_number
                resize!(jdx.y_data,total_point_number)
                jdx.x_data = similar(jdx.y_data)
                jdx.x_data[1] = x_point
                for i in 2:data_lines_number
                    ln = readline(io_file)
                    jdx.x_data[i] = addYline!(jdx,ln,1,i)
                end
            end

        end 
        if haskey(jdx.data_headers,"YFACTOR") 
            y_factor = jdx.data_headers["YFACTOR"] 
            jdx.y_data .*= y_factor
        end        
        return (x=jdx.x_data,
                y=jdx.y_data, 
                headers=jdx.data_headers)
    end
end