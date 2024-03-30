module CubicBSpline

#= This module holds the functions for manipulating cubic B-splines.
The math and terminology closely follow Ooyama, K. V., 2002: The cubic-spline transform method: Basic definitions and tests in a 1d single domain. Mon. Wea. Rev., 130, 2392–2415. =#

using LinearAlgebra
using SparseArrays
using SuiteSparse

export SplineParameters, Spline1D
#export R0, R1T0, R1T1, R1T2, R2T10, R2T20, R3, PERIODIC
export SBtransform, SBtransform!, SAtransform!, SItransform!
export SAtransform, SBxtransform, SItransform, SIxtransform, SIxxtransform
export setMishValues

#Define some convenient aliases
const real = Float64
const int = Int64
const uint = UInt64
const ONESIXTH = 1.0/6.0
const FOURSIXTH = 4.0/6.0
const sqrt35 = sqrt(3.0/5.0)

# Define homogeneous boundary conditions
# Inhomgoneous conditions to be implemented later
const R0 = Dict("R0" => 0)
const R1T0 = Dict("α1" => -4.0, "β1" => -1.0)
const R1T1 = Dict("α1" =>  0.0, "β1" =>  1.0)
const R1T2 = Dict("α1" =>  2.0, "β1" => -1.0)
const R2T10 = Dict("α2" => 1.0, "β2" => -0.5)
const R2T20 = Dict("α2" => -1.0, "β2" => 0.0)
const R3 = Dict("R3" => 0)
const PERIODIC = Dict("PERIODIC" => 0)

# Fix the mish to 3 points
const mubar = 3
const gaussweight = [8.0/18.0, 5.0/18.0, 8.0/18.0]

# Define the spline parameters
Base.@kwdef struct SplineParameters
    xmin::real = 0.0
    xmax::real = 0.0
    num_cells::int = 1
    l_q::real = 2.0
    BCL::Dict = R0
    BCR::Dict = R0
    DX::real = (xmax - xmin) / num_cells
    DXrecip::real = 1.0/DX
end

# This structure holds a 1D Cubic B-spline type
struct Spline1D
    params::SplineParameters
    gammaBC::Matrix{real}
    pq
    pqFactor::SuiteSparse.CHOLMOD.Factor{Float64}
    mishDim::int
    bDim::int
    mishPoints::Vector{real}
    uMish::Vector{real}
    b::Vector{real}
    a::Vector{real}
end

function basis(sp::SplineParameters, m::int, x::real, derivative::int)

    # Calculate the cubic b-spline basis function for a given set of spline parameters
    # m = node index
    # x = physical location
    # derivative = 0, 1, or 2 for none, first, or second derivative
    
    b = 0.0
    if (x < sp.xmin) || (x > sp.xmax)
        throw(DomainError(x, "x outside spline domain"))
    end
    xm = sp.xmin + (m * sp.DX)
    delta = (x - xm) * sp.DXrecip
    z = abs(delta)
    if (z < 2.0)
        if (derivative == 0)
            z = 2.0 - z
            b = (z*z*z) * ONESIXTH
            z -= 1.0
            if (z > 0)
                b -= (z*z*z) * FOURSIXTH
           end
        elseif (derivative == 1)
            z = 2.0 - z
            b = (z*z) * ONESIXTH
            z -= 1.0
            if (z > 0)
                b -= (z*z) * FOURSIXTH
            end
            b *= ((delta > 0) ? -1.0 : 1.0) * 3.0 * sp.DXrecip
        elseif (derivative == 2)
            z = 2.0 - z
            b = z
            z -= 1.0
            if (z > 0)
                b -= z * 4
            end
            b *= sp.DXrecip * sp.DXrecip
        elseif (derivative == 3)
            if (z > 1.0)
                b = 1.0
            elseif (z < 1.0)
                b = -3.0
            end
            b *= ((delta > 0) ? -1.0 : 1.0) * sp.DXrecip * sp.DXrecip * sp.DXrecip
        end
    end
    return b
end

function calcGammaBC(sp::SplineParameters)

    # Calculate the boundary condition matrix
    if haskey(sp.BCL,"α1")
        rankL = 1
    elseif haskey(sp.BCL,"α2")
        rankL = 2
    elseif sp.BCL == R0
        rankL = 0
    elseif sp.BCL == R3
        rankL = 3
    elseif sp.BCL == PERIODIC
        rankL = 1
    end

    if haskey(sp.BCR,"α1")
        rankR = 1
    elseif haskey(sp.BCR,"α2")
        rankR = 2
    elseif sp.BCR == R0
        rankR = 0
    elseif sp.BCR == R3
        rankR = 3
    elseif sp.BCR == PERIODIC
        rankR = 2
    end

    Mdim = sp.num_cells + 3
    Minterior_dim = Mdim - rankL - rankR

    # Create the BC matrix
    gammaBC = zeros(real, Minterior_dim, Mdim)
    if haskey(sp.BCL,"α1")
        gammaBC[1,1] = sp.BCL["α1"]
        gammaBC[2,1] = sp.BCL["β1"]
        gammaBC[:,2:(Mdim-rankR)] = Matrix(1.0I, Minterior_dim, Minterior_dim)
    elseif haskey(sp.BCL,"α2")
        gammaBC[1,1] = sp.BCL["α2"]
        gammaBC[1,2] = sp.BCL["β2"]
        gammaBC[:,3:(Mdim-rankR)] = Matrix(1.0I, Minterior_dim, Minterior_dim)
    elseif sp.BCL == PERIODIC
        gammaBC[Minterior_dim,1] = 1.0
        gammaBC[:,2:(Mdim-rankR)] = Matrix(1.0I, Minterior_dim, Minterior_dim)
    elseif sp.BCL == R3
        gammaBC[:,4:(Mdim-rankR)] = Matrix(1.0I, Minterior_dim, Minterior_dim)
    else
        gammaBC[:,1:(Mdim-rankR)] = Matrix(1.0I, Minterior_dim, Minterior_dim)
    end

    if haskey(sp.BCR,"α1")
        gammaBC[Minterior_dim,Mdim] = sp.BCR["α1"]
        gammaBC[Minterior_dim-1,Mdim] = sp.BCR["β1"]
    elseif haskey(sp.BCR,"α2")
        gammaBC[Minterior_dim,Mdim] = sp.BCR["α2"]
        gammaBC[Minterior_dim,Mdim-1] = sp.BCR["β2"]
    elseif sp.BCR == PERIODIC
        gammaBC[1,Mdim-1] = 1.0
        gammaBC[2,Mdim] = 1.0
    end
    return gammaBC
end

function calcPQfactor(sp::SplineParameters, gammaBC::Matrix{real})

    # Calculate the P + Q matrix, convert to sparse representation, 
    # and factor with Cholesky decomposition
    eps_q = ((sp.l_q*sp.DX)/(2*π))^6
    Mdim = sp.num_cells + 3

    # Create the P and Q matrices
    P = zeros(real, Mdim, Mdim)
    Q = zeros(real, Mdim, Mdim)

    for mi1 = 1:Mdim
        for mi2 = 1:Mdim
            if abs(mi1 - mi2) > 3
                continue
            end
            m1 = mi1 - 2
            m2 = mi2 - 2
            for mc = 0:(sp.num_cells-1)
                if (mc < (m1 - 2)) || (mc > (m1 + 1))
                    continue
                end
                if (mc < (m2 - 2)) || (mc > (m2 + 1))
                    continue
                end
                for mu = 1:mubar
                    i = mu + (mubar * mc)
                    x = sp.xmin + (mc * sp.DX) + sp.DX * ((mu/2.0 - 1.0) * sqrt35) + sp.DX * 0.5
                    pm1 = basis(sp, m1, x, 0)
                    qm1 = basis(sp, m1, x, 3)
                    pm2 = basis(sp, m2, x, 0)
                    qm2 = basis(sp, m2, x, 3)
                    P[mi1,mi2] += sp.DX * gaussweight[mu] * pm1 * pm2
                    Q[mi1,mi2] += sp.DX * gaussweight[mu] * eps_q * qm1 * qm2
                end
            end
        end
    end

    # Fold in the BCs to get open form
    PQ = Symmetric(P + Q)
    PQopen = Symmetric((gammaBC * P * gammaBC') + (gammaBC * Q * gammaBC'))
    PQsparse = sparse(PQopen)
    PQfactor = (cholesky(PQsparse))
    return PQ, PQfactor
end

function calcMishPoints(sp::SplineParameters)

    # Calculate the Gaussian quadrature "mish" points in between the nodes
    x = zeros(real,sp.num_cells*mubar)
    for mc = 0:(sp.num_cells-1)
        for mu = 1:mubar
            i = mu + (mubar * mc)
            x[i] = sp.xmin + (mc * sp.DX) + sp.DX * ((mu/2.0 - 1.0) * sqrt35) + sp.DX * 0.5
        end
    end
    return x
end

function Spline1D(sp::SplineParameters)

    # Constructor for the 1D Spline structure
    gammaBC = calcGammaBC(sp)
    pq, pqFactor = calcPQfactor(sp, gammaBC)

    mishDim = sp.num_cells*mubar
    mishPoints = calcMishPoints(sp)
    uMish = zeros(real,sp.num_cells*mubar)

    bDim = sp.num_cells + 3
    b = zeros(real,bDim)
    a = zeros(real,bDim)

    spline = Spline1D(sp,gammaBC,pq,pqFactor,mishDim,bDim,mishPoints,uMish,b,a)
    return spline
end

function setMishValues(spline::Spline1D, uMish::Vector{real})

    # Convenience function to assign a vector to the spline mish values
    spline.uMish .= uMish
end

function SBtransform(sp::SplineParameters, uMish::Vector{real})

    # SB transform from physical space (u) to B vector
    # The various versions of this function take different type inputs, which is useful for different calling scenarios
    Mdim = sp.num_cells + 3
    b = zeros(real,Mdim)

    for mi = 1:Mdim
        m = mi - 2
        for mc = 0:(sp.num_cells-1)
            if (mc < (m - 2)) || (mc > (m + 1))
                continue
            end
            for mu = 1:mubar
                i = mu + (mubar * mc)
                x = sp.xmin + (mc * sp.DX) + sp.DX * ((mu/2.0 - 1.0) * sqrt35) + sp.DX * 0.5
                bm = basis(sp, m, x, 0)
                b[mi] += sp.DX * gaussweight[mu] * bm * uMish[i]
            end
        end
    end

    # Don't border fold it here, only in SA
    return b
end

function SBtransform(spline::Spline1D, uMish::Vector{real})

    b = SBtransform(spline.params,uMish)
    return b
end

function SBtransform!(spline::Spline1D)

    b = SBtransform(spline.params,spline.uMish)
    spline.b .= b
end

#= function SBxtransform(sp::SplineParameters, uMish::Vector{real}, BCL, BCR)

    # This function does the SBx transform but is not working
    # It is essentially integration by parts, but still needs some bug fixes evidently
    # Do not use at this time!
    Mdim = sp.num_cells + 3
    b = zeros(real,Mdim)

    for mi = 1:Mdim
        m = mi - 2
        for mc = 0:(sp.num_cells-1)
            if (mc < (m - 2)) || (mc > (m + 1))
                continue
            end
            for mu = 1:mubar
                i = mu + (mubar * mc)
                x = sp.xmin + (mc * sp.DX) + sp.DX * ((mu/2.0 - 1.0) * sqrt35) + sp.DX * 0.5
                bm = basis(sp, m, x, 1)
                b[mi] += sp.DX * gaussweight[mu] * bm * uMish[i]
            end
        end
        bl = basis(sp, m, sp.xmin, 0) * BCL
        br = basis(sp, m, sp.xmax, 0) * BCR
        b[mi] = br - bl - b[mi] 
    end
    
    return b
end 

function SBxtransform(spline::Spline1D, uMish::Vector{real}, BCL, BCR)

    bx = SBxtransform(spline.params,uMish,BCL,BCR)
    return bx
end =#

function SAtransform(sp::SplineParameters, gammaBC::Matrix{Float64}, pqFactor, b::Vector{real})

    # The SA transform converts a B matrix into A coefficients, applying the boundary conditions
    # The various versions of this function take different type inputs, which is useful for different calling scenarios
    a = gammaBC' * (pqFactor \ (gammaBC * b))
    return a
end

function SAtransform(spline::Spline1D, b::AbstractVector)

    a = spline.gammaBC' * (spline.pqFactor \ (spline.gammaBC * b))
    return a
end

function SAtransform!(spline::Spline1D)

    # In-place version of the SA transform
    spline.a .= spline.gammaBC' * (spline.pqFactor \ (spline.gammaBC * spline.b))
end

function SAtransform(spline::Spline1D, b::Vector{real}, ahat::Vector{real})

    btilde = spline.gammaBC * (b - (spline.pq * ahat))
    a = (spline.gammaBC' * (spline.pqFactor \ btilde)) + ahat
    return a
end

function SItransform(sp::SplineParameters, a::Vector{real}, x::real, derivative::int = 0)

    # The SI transform converts A coefficients back to physical space (u)
    # The various versions of this function take different type inputs, which is useful for different calling scenarios
    u = 0.0
    xm = ceil(int,(x - sp.xmin - (2.0 * sp.DX)) * sp.DXrecip)
    for m = xm:(xm + 3)
        if (m >= -1) && (m <= (sp.num_cells+1))
            mi = m + 2
            u += basis(sp, m, x, derivative) * a[mi]
        end
    end
    return u
end

function SItransform(sp::SplineParameters, a::Vector{real}, derivative::int = 0)

    u = zeros(real,sp.num_cells*mubar)
    for mc = 0:(sp.num_cells-1)
        for mu = 1:mubar
            i = mu + (mubar * mc)
            x = sp.xmin + (mc * sp.DX) + sp.DX * ((mu/2.0 - 1.0) * sqrt35) + sp.DX * 0.5
            for m = (mc-1):(mc+2)
                if (m >= -1) && (m <= (sp.num_cells+1))
                    mi = m + 2
                    u[i] += basis(sp, m, x, derivative) * a[mi]
                end
            end
        end
    end
    return u
end

function SItransform(sp::SplineParameters, a::Vector{real}, points::Vector{real}, derivative::int = 0)

    u = zeros(real,length(points))
    for i in eachindex(points)
        xm = ceil(int,(points[i] - sp.xmin - (2.0 * sp.DX)) * sp.DXrecip)
        for m = xm:(xm + 3)
            if (m >= -1) && (m <= (sp.num_cells+1))
                mi = m + 2
                u[i] += basis(sp, m, points[i], derivative) * a[mi]
            end
        end
    end
    return u
end

function SItransform(sp::SplineParameters, a::Vector{real}, points::Vector{real}, u::AbstractVector, derivative::int = 0)

    for i in eachindex(points)
        u[i] = 0.0
        xm = ceil(int,(points[i] - sp.xmin - (2.0 * sp.DX)) * sp.DXrecip)
        for m = xm:(xm + 3)
            if (m >= -1) && (m <= (sp.num_cells+1))
                mi = m + 2
                u[i] += basis(sp, m, points[i], derivative) * a[mi]
            end
        end
    end
    return u
end

function SItransform!(spline::Spline1D)

    # In-place SI transform
    u = SItransform(spline.params,spline.a,spline.mishPoints,spline.uMish)
    return u
end

function SItransform(spline::Spline1D, u::AbstractVector)

    u = SItransform(spline.params,spline.a,spline.mishPoints,u)
    return u
end

function SItransform(spline::Spline1D, points::Vector{real}, u::AbstractVector)

    u = SItransform(spline.params,spline.a,points,u)
    return u
end

function SItransform_matrix(spline::Spline1D, points::Vector{Float64}, derivative::Int64 = 0)

    # This function creates a matrix for the SI transform, which is mostly useful for debugging
    # but could have value in linear equation solutions later on
    sp = spline.params
    u = zeros(Float64,sp.num_cells*mubar,spline.bDim)
    for i in eachindex(points)
        xm = ceil(Int64,(points[i] - sp.xmin - (2.0 * sp.DX)) * sp.DXrecip)
        for m = xm:(xm + 3)
            if (m >= -1) && (m <= (sp.num_cells+1))
                mi = m + 2
                u[i,mi] = basis(sp, m, points[i], derivative)
            end
        end
    end
    return u
end

function SIxtransform(spline::Spline1D)

    # SIx transform gives back the 1st derivative of the function in physical space
    uprime = SItransform(spline.params,spline.a,spline.mishPoints,1)
    return uprime
end

function SIxtransform(spline::Spline1D, uprime::AbstractVector)

    uprime = SItransform(spline.params,spline.a,spline.mishPoints,uprime,1)
    return uprime
end

function SIxtransform(spline::Spline1D, points::Vector{real}, uprime::AbstractVector)

    uprime = SItransform(spline.params,spline.a,points,uprime,1)
    return uprime
end

function SIxtransform(sp::SplineParameters, a::Vector{real}, points::AbstractVector)

    uprime = SItransform(sp,a,points,1)
    return uprime
end

function SIxxtransform(spline::Spline1D)

    # SIxx transform gives back the 2nd derivative of the function in physical space
    uprime2 = SItransform(spline.params,spline.a,spline.mishPoints,2)
    return uprime2
end

function SIxxtransform(spline::Spline1D, uprime2::AbstractVector)

    uprime2 = SItransform(spline.params,spline.a,spline.mishPoints,uprime2,2)
    return uprime2
end

function SIxxtransform(spline::Spline1D, points::Vector{real}, uprime2::AbstractVector)

    uprime2 = SItransform(spline.params,spline.a,points,uprime2,2)
    return uprime2
end

function SIxxtransform(sp::SplineParameters, a::Vector{real}, points::Vector{real})

    uprime2 = SItransform(sp,a,points,2)
    return uprime2
end

#Module end
end