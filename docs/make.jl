push!(LOAD_PATH,"../src/")
include("../src/BandPyrometry.jl")
using Documenter,.BandPyrometry
mathengine = Documenter.MathJax3()
makedocs(
        sitename = "BandPyrometry.jl",
        highlightsig = false,
        checkdocs = :none,
        format=Documenter.HTML(size_threshold = 2000 * 2^10),
        pages=[
                "Home" => "index.md"
                "Examples"=>["BandPyrometry"=>"pluto_tests_git.md"
                              "Planck" =>"pluto_tests_git.md"
                ]
                "Modules" => [
                    "Planck" =>"planck.md"
                    "BandPyrometry" => "bandpyrometry.md"
                    "Pyrometers" => "pyrometers.md"
                    "JDXreader" => "jcamp-reader.md"
                 ] 
               ]#
			   )
#deploydocs(;
#         repo="github.com/Manarom/BandPyrometry"
#)