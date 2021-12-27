# Calculation of the combined WENO interpolant evaluated at x.
function weno6_interpolant(δ, f1, f2, f3, f4, f5, f6)

    # Interpolants
    p = weno6_candidates(δ, f1, f2, f3, f4, f5, f6)

    # weight for each interpolant
    w = weno6_weights(δ, f1, f2, f3, f4, f5, f6)

    return sum(@. p * w)
end

function weno6_candidates(δ, f1, f2, f3, f4, f5, f6)

    # Candidate interpolants over stencils
    # S1 = {x1, ..., x4},
    # S2 = {x2, ..., x5}, and
    # S3 = {x3, ..., x6},
    # respectively

    p = @fastmath (
        f1 + (f2 - f1) * (δ + 3)
        + 0.5 * (f3 - 2*f2 + f1) * (δ + 3) * (δ + 2)
        + (f4 - 3*f3 + 3*f2 - f1) * (δ + 3) * (δ + 2) * (δ + 1) / 6,

        f2 + (f3 - f2) * (δ + 2)
        + 0.5 * (f4 - 2*f3 + f2) * (δ + 2) * (δ + 1)
        + (f5 - 3*f4 + 3*f3 - f2) * (δ + 2) * (δ + 1) * (  δ  ) / 6,

        f3 + (f4 - f3) * (δ + 1)
        + 0.5 * (f5 - 2*f4 + f3) * (δ + 1) * (  δ  )
        + (f6 - 3*f5 + 3*f4 - f3) * (δ + 1) * (  δ  ) * (δ - 1) / 6,
    )

    return p
end

function weno6_weights(δ, f1, f2, f3, f4, f5, f6)

    # Ideal weights
    C = @fastmath (
        0.05 * (1 - δ) * (2 - δ),
        0.1  * (2 - δ) * (3 + δ),
        0.05 * (3 + δ) * (2 + δ),
    )

    # Smoothness indicators. Their general form is shown on Eq. (7)
    S = @fastmath (
        (814*f4^2 + 4326*f3^2 + 2976*f2^2 + 244*f1^2
         - 3579*f3*f4 - 6927*f3*f2 + 1854*f3*f1
         + 2634*f4*f2 -  683*f4*f1 - 1659*f2*f1 ) / 180,

        (1986*f4^2 + 1986*f3^2 + 244*f2^2 + 244*f5^2
         + 1074*f3*f5 - 3777*f3*f4 - 1269*f3*f2
         + 1074*f4*f2 - 1269*f5*f4 -  293*f5*f2 ) / 180,

        (814*f3^2 + 4326*f4^2 + 2976*f5^2 + 244*f6^2
         -  683*f3*f6 + 2634*f3*f5 - 3579*f3*f4
         - 6927*f4*f5 + 1854*f4*f6 - 1659*f5*f6 ) / 180,
    )

    # Smoothness-corrected, nonlinear weights
    # eps() to avoid divission by zero
    α = @fastmath @. C / (eps() + S)^2

    # Normalization
    α = @fastmath α ./ sum(α)

    return α
end

"""
    advect_weno6(f::AbstractVector, dx::Real, shift::Real)

Perform a 6th-order WENO reconstruction to advect the values of
the 1-dimensional array `f` according to a shifting `shift` on space.

Note that `dx` is the (constant) step used to discretize the space points
where `f` is sampled.

"""
function advect_weno6(f::AbstractVector, dx::Real, shift::Real)
    N = length(f)

    xs = shift/dx + 1     # Shifted, normalized position
    xn = floor(Int64, xs) # Number of the first node to the left of xs
    δ = xn - xs           # Position of the reference node relative to xs

    # Prepare a shifted copy of f such that fs[i] = f[mod1(i-xn, N)]
    fs = circshift(f, xn)
    
    advected = similar(f)
    @tturbo for i in 3:(N-3)
        advected[i] = weno6_interpolant(δ, fs[i-2], fs[i-1], fs[i], fs[i+1], fs[i+2], fs[i+3])
    end
    
    ############ Periodic BC's
    # i = N-2
    advected[N-2] = weno6_interpolant(δ, fs[N-4], fs[N-3], fs[N-2], fs[N-1], fs[ N ], fs[ 1 ])
    # i = N-1
    advected[N-1] = weno6_interpolant(δ, fs[N-3], fs[N-2], fs[N-1], fs[ N ], fs[ 1 ], fs[ 2 ])
    # i = N
    advected[ N ] = weno6_interpolant(δ, fs[N-2], fs[N-1], fs[ N ], fs[ 1 ], fs[ 2 ], fs[ 3 ])
    # i = 1
    advected[ 1 ] = weno6_interpolant(δ, fs[N-1], fs[ N ], fs[ 1 ], fs[ 2 ], fs[ 3 ], fs[ 4 ])
    # i = 2
    advected[ 2 ] = weno6_interpolant(δ, fs[ N ], fs[ 1 ], fs[ 2 ], fs[ 3 ], fs[ 4 ], fs[ 5 ])

    return advected
end

function advect_weno6_2d(f::AbstractMatrix{Float64}, dx::Real, shift::AbstractVector;
                         dim = 1)
    # For dim ∈ {1, 2}
    other_dim = 3-dim

    # Checks
    @assert dim<3 "The value of 'dim' must be 0 or 1. You entered $dim."
    @assert length(shift)==size(f, other_dim) "The shift should have the sane number of elements as the non-advected dimension."

    # Preallocate arrays and views
    adv = similar(f)
    A = PermutedDimsArray(adv, (dim, other_dim)) # Transposed view
    F = PermutedDimsArray(f,   (dim, other_dim)) # Transposed view

    # Perform 1d advections
    @batch for i in axes(f, other_dim)
        A[:, i] = advect_weno6(view(F, :, i), dx, shift[i])
    end

    return adv
end
