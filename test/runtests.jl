using FactCheck

include("apitest.jl")
include("mathprogtest.jl")
include("testexamples.jl")
include("mathprogtestextra.jl")
include("testconvex.jl")
include("test_consts.jl")

FactCheck.exitstatus()
