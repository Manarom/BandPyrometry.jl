push!(LOAD_PATH,"../src/")

include("../src/Planck.jl")
include("../src/BandPyrometry.jl")
include("../src/JDXreader.jl")
include("../src/Pyrometers.jl")

using Documenter,.BandPyrometry,.Planck, .JDXreader, .Pyrometers
#using Base
println(pwd())
#mathengine = Documenter.MathJax3()
makedocs(
         #checkdocs =:none,
        sitename = "BandPyrometry.jl",
        highlightsig = false,
        checkdocs = :none,
        remotes=nothing,
       pages=[
                "Home" => "index.md"
                "Modules" => [
                    "Planck" =>"planck.md"
                    "BandPyrometry" => "bandpyrometry.md"
                    "Pyrometers" => "pyrometers.md"
                    "JDXreader" => "jcamp-reader.md"
                 ] 
               ]#
			   )