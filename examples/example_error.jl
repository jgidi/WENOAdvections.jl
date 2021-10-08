using WENOAdvection

# Test function
f(x) = @. exp(-0.5*(x/0.1)^2) + sin(x)^2 #+ float(-pi/2 < x < pi/2)

# L2 norm of the errors
L2(fref, fnum, dx) = dx * sqrt( sum(abs2, fref .- fnum) )


# Shift and space domain
shift = rand()
xmin, xmax = -pi, pi

# Number of nodes to consider
Nnodes = @. 2^(3:16)

errors = Float64[]
for Nx in Nnodes
    # Make space
    dx = (xmax-xmin)/Nx
    x = range(xmin, step=dx, length=Nx)

    # Compute shifted functions
    # fref: analytycal
    # fnum: WENO-shifted
    fref = f(x .- shift)
    fnum = advect_weno6(f(x), dx, shift)

    # Compute and save error
    append!(errors, L2(fref, fnum, dx))
end

# Make a plot with both initial and shifted functions
using Plots

p = plot(legend = :false,
         box = :false,
         scale = :log10,
         xlabel = "dx",
         ylabel = "L2(Err)",
         xticks = 10.0 .^ (-6:0),
         yticks = 10.0 .^ (-16:2:0),
         title = "shift = $shift",
         )

scatter!(p, 2pi ./ Nnodes, errors)

display(p)
