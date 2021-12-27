using WENOAdvections

# Construct uniform space
Nx, Nv = 64, 128
vmin, vmax = -2.0, 2.0
xmin, xmax = -2.0, 2.0

# Make space
dx = (xmax-xmin)/Nx
x  = range(xmin, step=dx, length=Nx)
#
dv = (vmax-vmin)/Nv
v  = range(vmin, step=dv, length=Nv)

# Test function

# Discontinuous step function
f = float(@. abs(x) < 0.5 ) * ones(Nv)'
# f = ones(Nx) * float(@. abs(v) < 0.5 )'

# Shift
shiftx = 2v
# shiftv = @. x^2

# Make advection
fadv = WENOAdvections.advect_weno6_2d(f, dx, shiftx, dim=1);
# fadv = WENOAdvections.advect_weno6_2d(f, dv, shiftv, dim=2);


# Make plot of te original and WENO-shifted function
using Plots

p1 = heatmap(x, v, f', title = "\$ f(x, v) \$")
p2 = heatmap(x, v, fadv', title = "\$ f(x-2v, v) \$")
# p2 = heatmap(x, v, fadv', title = "\$ f(x, v-x^2) \$")

p = plot(p1, p2,
         xlabel = "\$ x \$",
         ylabel = "\$ v \$",
         dpi = 300,
         )

# Save plot
savefig(p, "example2d.png")

# Show plot
display(p)
