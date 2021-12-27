# Based upon Sec. 3 from:
# 'Level set equations on surfaces via de Closest Point Method' (2008),
# by Colin Macdonald & Steven Ruuth

module WENOAdvections

# External dependencies
using LoopVectorization         # Fast loops
using Polyester                 # Cheap threads

# Exported functions
export advect_weno6, advect_weno6_2d

include("include/advect_weno6.jl")

end # module
