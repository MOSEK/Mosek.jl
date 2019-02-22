using Test

include("apitest.jl")
#include("mathprogtest.jl")
include("testexamples.jl")
#include("mathprogtestextra.jl")

#import Pkg
#if (try Pkg.installed("Convex") != nothing catch; false; end)
#    include("testconvex.jl")
#end

include("test_mathprogbase.jl")
include("test_optimizer.jl")
