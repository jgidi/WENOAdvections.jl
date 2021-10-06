# Based upon Sec. 3 from:
# 'Level set equations on surfaces via de Closest Point Method' (2008),
# by Colin Macdonald & Steven Ruuth

module WENOAdvection

# External dependencies
using StaticArrays, LinearAlgebra

# Exported functions
export advect_weno6

include("include/advect_weno6.jl")

end # module
