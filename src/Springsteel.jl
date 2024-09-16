module Springsteel

# Functions to define a spectral grid
abstract type AbstractGrid end

using CSV
using DataFrames
using SharedArrays
using SparseArrays

#Define some convenient aliases
const real = Float64
const int = Int64
const uint = UInt64

# These are declared as submodules to avoid namespace clashes with each other and other packages
include("CubicBSpline.jl")
include("Fourier.jl")
include("Chebyshev.jl")
using .CubicBSpline, .Fourier, .Chebyshev

export AbstractGrid, GridParameters
export CubicBSpline, Chebyshev
export SplineParameters, Spline1D
export createGrid, getGridpoints, calcTileSizes
export read_physical_grid, write_grid
export calcPatchMap, calcHaloMap, allocateSplineBuffer, num_columns
export spectralTransform!, splineTransform!, tileTransform!, gridTransform!

Base.@kwdef struct GridParameters
    geometry::String = "R"
    xmin::real = 0.0
    xmax::real = 0.0
    num_cells::int = 0
    rDim::int = num_cells * CubicBSpline.mubar
    b_rDim::int = num_cells + 3
    l_q::real = 2.0
    BCL::Dict = CubicBSpline.R0
    BCR::Dict = CubicBSpline.R0
    ymin::real = 0.0
    ymax::real = 2 * π
    lDim::int = 0
    b_lDim::int = 0
    BCU::Dict = Fourier.PERIODIC
    BCB::Dict = Fourier.PERIODIC
    zmin::real = 0.0
    zmax::real = 0.0
    zDim::int = 0
    b_zDim::int = min(zDim, floor(((2 * zDim) - 1) / 3) + 1)
    BCB::Dict = Chebyshev.R0
    BCT::Dict = Chebyshev.R0
    vars::Dict = Dict("u" => 1)
    # Patch indices
    spectralIndexL::int = 1
    spectralIndexR::int = spectralIndexL + b_rDim - 1
    patchOffsetL::int = (spectralIndexL - 1) * 3
    patchOffsetR::int = patchOffsetL + rDim
    tile_num::int = 0
    r_incr_out::real = (xmax - xmin) / num_cells
    l_incr_out::real = (ymax - ymin) / (rDim*2+1)
    z_incr_out::real = (zmax - zmin) / zDim
end

# Include functions for implemented grids
include("r_grid.jl")
include("rz_grid.jl")
include("rl_grid.jl")
include("rlz_grid.jl")

# Not yet implemented
struct Z_Grid <: AbstractGrid
    params::GridParameters
    columns::Array{Chebyshev1D}
    spectral::Array{Float64}
    physical::Array{Float64}
end

# I/O routines
include("io.jl")

function createGrid(gp::GridParameters)

    # Call the respective grid factory
    if gp.geometry == "R"
        # R grid
        grid = create_R_Grid(gp)
        return grid

    elseif gp.geometry == "RZ"
        # RZ grid
        grid = create_RZ_Grid(gp)
        return grid

    elseif gp.geometry == "RL"
        # RL grid
        grid = create_RL_Grid(gp)
        return grid

    elseif gp.geometry == "RLZ"
        # RLZ grid
        grid = create_RLZ_Grid(gp)
        return grid
        
    elseif gp.geometry == "Z"
        # Z grid
        throw(DomainError(0, "Z column model not implemented yet"))
    else
        # Unknown grid
        throw(DomainError(0, "Unknown geometry"))
    end
    
end

# Module end
end
