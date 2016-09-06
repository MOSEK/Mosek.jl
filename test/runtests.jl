using FactCheck

include("apitest.jl")
include("mathprogtest.jl")
include("testexamples.jl")
include("mathprogtestextra.jl")
include("testconvex.jl")
#include("test_consts.jl") # Something in Julia 0.5 make this test go nuts

FactCheck.exitstatus()
