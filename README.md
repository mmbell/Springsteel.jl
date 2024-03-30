# Springsteel.jl

Springsteel is a semi-spectral grid engine that uses a mixture of cubic B-spline, Fourier, and Chebyshev basis functions to represent physical variables and their spatial derivatives. Springsteel currently supports 1-D radius (R), 2-D polar radius-azimuth (RL, where L denotes azimuthal angle lamba), 2-D radius-height (RZ), and 3-D cylindrical radius-azimuth-height (RLZ) grids. The name comes from an amalgamation of "spectral grid engine" and the idea that this package provides the "steel" for the Scythe.jl numerical model and the Daisho.jl data analysis and assimilation software.

Springsteel has some limited distributed grid capabilities. The grid can be decomposed in the radial direction into a user-specified number of "tiles" with a halo at the tile interfaces.  

### Installation

After cloning this repository, start Julia using Springsteel.jl as the project directory. This can be done on the command line using `julia --project` or set the JULIA_PROJECT environmental variable:

`export JULIA_PROJECT=/path/to/Springsteel.jl`

To install Springsteel, in the REPL, go into Package mode by pressing `]`. You will see the REPL change color and indicate `pkg` mode. 

If you are actively developing or modifying Springsteel then you can install the module using `dev /path/to/Springsteel.jl` in `pkg` mode. This will update the module as changes are made to the code. You should see the dependencies being installed, and then the Springsteel package will be precompiled. Exit Package mode with ctrl-C.

If you wish to just install a static version of the latest code, run `activate` to activate the package environment. Then, run `instantiate` to install the necessary dependencies. Exit Package mode with ctrl-C.

Test to make sure the precompilation was successful by running `using Springsteel` in the REPL. If everything is successful then you should get no errors and it will just move to a new line.

### Springsteel API

To create a grid, a GridParameters structure needs to be declared first. A sample declaration for a 1-D radius grid with a single variable is:
```
grid_params = GridParameters(
        geometry = "R",
        xmin = -50.0,
        xmax = 50.0,
        num_cells = 100,
        BCL = Dict(
            "u" => CubicBSpline.PERIODIC),
        BCR = Dict(
            "u" => CubicBSpline.PERIODIC),
        vars = Dict(
            "u" => 1))
```

After defining the parameters, you can call `grid = createGrid(grid_params)` to return a grid object. It is often useful to call `getGridpoints(grid)` to get a multi-dimensional array of the underlying physical gridpoints. 

The supertype of all grids is `AbstractGrid` which can be used to define functional interfaces that don't depend on the underlying grid geometry or basis functions. Multiple dispatch is used to define the various grid types, their transforms, and I/O.

The primary transforms are `spectralTransform!(grid)` which converts from physical to spectral space, and `gridTransform!(grid)` which converts back to physical space. Additional partial transforms are defined in the code and will documented better in future code releases.

### Future plans
Support for CF-compliant NetCDF input and output will be added in the near future. Future versions will also support different geometries that can be constructed from the underlying basis functions (e.g. Cartesian and spherical). Support for grid nesting using the cubic B-splines will be added in future versions. Interested users are welcome to contribute to improve the grid engine. Stay tuned for more functionality!
