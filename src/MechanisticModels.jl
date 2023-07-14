include("KineticsLibrary.jl")
include("BioreactorDynamics.jl")

dX_monod(in, out, transfer) = change(in, out, transfer, monod(12, 1.2, 0.04))

# Add other dynamics here
