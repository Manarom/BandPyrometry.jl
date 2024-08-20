```@raw html
<p>The static version (with no interactive elements) of Pluto notebook with <b>BandPyrometry.jl</b> usage and main ideas description is available at: <a href="../bandpyrometry_test_git.html">bandpyrometry_test_git.html</a>,
pdf version is available at:<a href="../assets/bandpyrometry_test_git.pdf">bandpyrometry_test_git.pdf</a>
.</p



<p>The static version of <b>Planck.jl</b> testing Pluto notebook is available at :<a href="../planck_test_git.html">planck_test_git.html</a>,
pdf version is available at:<a href="../assets/planck_test_git.pdf">planck_test_git.pdf</a>
.</p
```


Notebook file with `BandPyrometry.jl` usage examples are available at [Pluto notebooks](https://github.com/Manarom/BandPyrometry.jl/blob/main/tests).

To run these notebooks, you need:
1) Install `julia` language itself from its official [download page](https://julialang.org/downloads) 
2) Install [Pluto](https://plutojl.org/) notebook from `julia` REPL by entering the following commands line-by-line:
```julia
import Pkg
Pkg.add("Pluto")
using Pluto
Pluto.run()
```
The last line will launch the Pluto starting page in your default browser 
3) Copy/Clone the entire GitHub [repository](https://github.com/Manarom/BandPyrometry.jl.git) to your local folder. As far as `BandPyrometry.jl` is not a registered package, all files needed to run the notebooks must be in the `../src` folder with respect to  `note_book_name.jl` file location.
4) Open notebook's `.jl`- file in `Pluto` by providing the full path to the *"Open a notebook"* text field on `Pluto`'s starting page. As far as `Pluto` has its own package manager, it will automatically install all necessary dependancies, which are marked in `using` cell of this file. 
