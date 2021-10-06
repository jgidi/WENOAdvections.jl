using WENOAdvection

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

shift = -3.7

fadv = advect_weno6(f, dv, shift);

# Make a plot with both initial and shifted functions
using PyPlot

plot(v, f,    color = "red",  linewidth=1,   label="\$ f(v) \$")
plot(v, fadv, color = "blue", linewidth=0.5, label="\$ f(v - s) \$")
legend(loc="best")
title("s = $shift")

show()
