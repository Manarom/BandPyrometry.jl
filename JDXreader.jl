using DelimitedFiles #,CSV

module JDXreader
# reader for JDXformat files
    export JDXfile,read!
    using ..DelimitedFiles
    mutable struct JDXfile
        # Main struct, prepares loads file name, parses data and data headers
        file_name::String
        data_headers::Dict{String,Union{String,Number}}
        x_data::Vector{Float64}
        y_data::Vector{Float64}
        is_violated_flag::Vector{Bool}
        JDXfile()=begin
            new("",
                Dict{String,Union{String,Number}}(),
                Vector{Float64}(),
                Vector{Float64}(),
                Vector{Bool}())
        end
    end
    function read_jdx_file(file_name::String)
        return file_name |> JDXfile |> read!
    end
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
        if haskey(jdx.data_headers,"NPOINTS")
            jdx.data_headers["NPOINTS"] = round(Int64, jdx.data_headers["NPOINTS"])
        end
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
            #xydata.xValues = new double[xydata.dataspec.NPoints];// инициализируем массивы, который будут выдаваться в качестве результата чтения файла
            #xydata.yValues = new double[xydata.dataspec.NPoints];
            #double step = ((xydata.dataspec.LastX - xydata.dataspec.FirstX)/xydata.dataspec.XFactor)/(xydata.dataspec.NPoints-1);
            ##DATE=2021/10/06
            ##TIME=13:42
            ##XUNITS=1/CM
            ##YUNITS=TRANSMITTANCE
            ##XFACTOR=1
            ##YFACTOR=0.00702876
            ##FIRSTX=0
            ##LASTX=15794.7156470938
            ##FIRSTY=8921.3875152695
            ##MAXX=15794.7156470938
            ##MINX=0
            ##MAXY=3.77354e+006
            ##MINY=1894.21
            ##NPOINTS=2048
            ##XYDATA=(X++(Y..Y))
            point_number = jdx.data_headers["NPOINTS"]
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
                total_point_number = jdx.data_headers["NPOINTS"]
                generateXvector!(jdx)
                resize!(jdx.y_data,total_point_number)
                for i in 2:data_lines_number
                    ln = readline(io_file)
                    x_point = addYline!(jdx,ln,number_of_y_point_per_chunk,i)
                end
            else number_of_y_point_per_chunk==1
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