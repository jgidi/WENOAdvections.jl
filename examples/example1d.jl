using WENOAdvections

# Construct uniform space
Nv = 128
vmin, vmax = -4.0, 4.0

# Make space
dv = (vmax-vmin)/Nv
v = range(vmin, vmax, length=Nv)

# Test function

# # Cold Gaussian
# vth = 0.1
# f = exp.(-0.5*v.^2/vth^2)

# Discontinuous step function
f = float(@. abs(v) < 0.5 )

shift = 3.65

# Make advection
fadv = advect_weno6(f, dv, shift);

# Make a plot with both initial and shifted functions
using Plots

p = plot(
    framestyle = :box,
    legend = :false,
    box = :false,
    xlabel = "\$ v \$"
)

plot!(p, v, f,    label = "\$ f(v)   \$", l = (:blue, 2))
plot!(p, v, fadv, label = "\$ f(v-s) \$", l = (:red, 1.5))

# Save plot
savefig(p, "example1d.png")

# Show plot
display(p)
