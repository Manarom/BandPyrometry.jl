push!(LOAD_PATH,"../src/")
using Documenter
using BandPyrometry
mathengine = Documenter.MathJax3()
makedocs(
        sitename = "BandPyrometry.jl",
        highlightsig = false,
        checkdocs = :none,
        format=Documenter.HTML(size_threshold = 2000 * 2^10),
        pages=[
                "Home" => "index.md"
                #"Examples"=>["BandPyrometry"=>"pluto_tests_git.md"
                #]
                "API" => [
                    "BandPyrometry" => "bandpyrometry.md"
                 ] 
               ]#
			   )