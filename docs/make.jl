push!(LOAD_PATH,"../src/")

include("../src/Planck.jl")
include("../src/BandPyrometry.jl")
include("../src/JDXreader.jl")
include("../src/Pyrometers.jl")

using Documenter,.BandPyrometry,.Planck, .JDXreader, .Pyrometers
makedocs(
        sitename = "BandPyrometry.jl",
        highlightsig = false,
        checkdocs = :none,
        format=Documenter.HTML(),
        repo = "https://github.com/Manarom/BandPyrometry/blob/{commit}{path}#{line}",
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
deploydocs(;
         repo="github.com/Manarom/BandPyrometry"
)