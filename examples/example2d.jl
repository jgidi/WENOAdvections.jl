using WENOAdvection

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
fadv = WENOAdvection.advect_weno6_2d(f, dx, shiftx, dim=1);
# fadv = WENOAdvection.advect_weno6_2d(f, dv, shiftv, dim=2);

using Plots

p = plot(xlabel = "\$ x \$",
         ylabel = "\$ v \$",
         )

p1 = heatmap(p, x, v, f', title = "\$ f(x, v) \$")
p2 = heatmap(p, x, v, fadv', title = "\$ f(x-2v, v) \$")

display(plot(p1, p2))
