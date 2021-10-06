
# Calculation of the combined WENO interpolant eval. at x.
# i is the position of the node immediately at the left of x
function weno_interpolant(f1, f2, f3, f4, f5, f6,
                          i, x)

    # 3 interpolants
    p = interpolants(f1, f2, f3, f4, f5, f6, i, x)

    # weight for each interpolant
    w = weights(f1, f2, f3, f4, f5, f6, i, x)

    return sum(@. p * w)
end

function interpolants(f1, f2, f3, f4, f5, f6,
                      i, x)

    # Candidate interpolants over stencils
    # S1 = {x[i-2], ..., x[i+1]},
    # S2 = {x[i-1], ..., x[i+2]}, and
    # S3 = {x[i], ..., x[i+3]},
    # respectively

    # The evaluation is done for x ∈ [x3, x4]

    # NOTE WENO order may be generalized by pre-calculating
    # interpolant shapes and weights
    p = SVector(
        f1 + (f2 - f1) * (x - i + 3)
        + (f3 - 2*f2 + f1) * (x - i + 3) * (x - i + 2) / 2
        + (f4 - 3*f3 + 3*f2 - f1) * (x - i + 3) * (x - i + 2) * (x - i + 1) / 6,

        f2 + (f3 - f2) * (x - i + 2)
        + (f4 - 2*f3 + f2) * (x - i + 2) * (x - i + 1) / 2
        + (f5 - 3*f4 + 3*f3 - f2) * (x - i + 2) * (x - i + 1) * (x - i) / 6,

        f3 + (f4 - f3) * (x - i + 1)
        + (f5 - 2*f4 + f3) * (x - i + 1) * (x - i) / 2
        + (f6 - 3*f5 + 3*f4 - f3) * (x - i + 1) * (x - i) * (x - i - 1) / 6,
    )

    return p
end

function weights(f1, f2, f3, f4, f5, f6,
                 i, x)

    # Ideal weights
    C = SVector(
        (i + 1 - x) * (i + 2 - x) / 20,
        (i + 2 - x) * (x - i + 3) / 10,
        (x - i + 3) * (x - i + 2) / 20,
    )

    # Smoothness indicators. Their general form is shown on Eq. (7)
    S = SVector(
        (814*f4^2 + 4326*f3^2 + 2976*f2^2 + 244*f1^2
         - 3579*f3*f4 - 6927*f3*f2 + 1854*f3*f1
         + 2634*f4*f2 -  683*f4*f1 - 1659*f2*f1 ) / 180,

        (1986*f4^2 + 1986*f3^2 + 244*f2^2 + 244*f5^2
         + 1074*f3*f5 - 3777*f3*f4 - 1269*f3*f2
         + 1074*f4*f2 - 1269*f5*f4 -  293*f5*f2 )/180,

        (814*f3^2 + 4326*f4^2 + 2976*f5^2 + 244*f6^2
         -  683*f3*f6 + 2634*f3*f5 - 3579*f3*f4
         - 6927*f4*f5 + 1854*f4*f6 - 1659*f5*f6 )/180,
    )

    # Smoothness-corrected, nonlinear weights
    # eps() to avoid divission by zero
    α = @. C / (eps() + S)^2

    # Normalization
    α = α / sum(α)

    return α
end

"""
    advect_weno6(f, dx, shift)

Perform a 6th-order WENO reconstruction to advect the values of
the 1-dimensional array `f` according to a shifting `shift` on space.

Note that `dx` is the (constant) step used to sample the space points
where `f` is sampled.

"""
function advect_weno6(f, dx, shift)
    N = length(f)

    # Normalize shift and index-shift
    normshift = shift/dx + 1
    fnormshift = floor(Int64, normshift)

    ff = circshift(f, fnormshift)
    advected = similar(f)
    for i in 3:(N-3)
        xs = i - normshift      # Normalized and shifted position
        xl = i - fnormshift     # Position of the left node from xs
        advected[i] = weno_interpolant(ff[i-2], ff[i-1], ff[i],
                                       ff[i+1], ff[i+2], ff[i+3],
                                       xl, xs)
    end

    ############ Periodic BC's
    # i = 1
    advected[1] = weno_interpolant(ff[N-1], ff[N], ff[1], ff[2], ff[3], ff[4],
                                   1-fnormshift, 1-normshift)
    # i = 2
    advected[2] = weno_interpolant(ff[N], ff[1], ff[2], ff[4], ff[4], ff[5],
                                   2-fnormshift, 2-normshift)
    # i = N-2
    advected[N-2] = weno_interpolant(ff[N-4], ff[N-3], ff[N-2], ff[N-1], ff[N], ff[1],
                                     N-2-fnormshift, N-2-normshift)
    # i = N-1
    advected[N-1] = weno_interpolant(ff[N-3], ff[N-2], ff[N-1], ff[N], ff[1], ff[2],
                                     N-1-fnormshift, N-1-normshift)
    # i = N
    advected[N] = weno_interpolant(ff[N-2], ff[N-1], ff[N], ff[1], ff[2], ff[3],
                                   N-fnormshift, N-normshift)

    return advected
end

# # This version looks simpler but is slower because of the mod1
# function advect_weno6(f, dx, shift)
#     N = length(f)

#     # Normalize shift and node positions
#     normshift = shift/dx + 1
#     fnormshift = floor(Int64, normshift)

#     advected = similar(f)
#     for i in 1:N
#         xs = i - normshift  # x normalized and shifted
#         j  = i - fnormshift # Position of the left node from xs

#         indices = SVector(j-2, j-1, j, j+1, j+2, j+3)

#         # Periodic BC's
#         indices = mod1.(indices, N)

#         f1, f2, f3, f4, f5, f6 = @views f[indices]
#         advected[i] = weno_interpolant(f1, f2, f3, f4, f5, f6, j, xs)
#     end

#     return advected
# end
