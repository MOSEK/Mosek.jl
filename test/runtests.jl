using Base.Test

include("apitest.jl")
#include("mathprogtest.jl")
include("testexamples.jl")
#include("mathprogtestextra.jl")


if (try Pkg.installed("Convex") != nothing catch false end)
    include("testconvex.jl")
end

if (try Pkg.installed("MathOptInterface") != nothing catch false end)
    include("test_mathoptinterface.jl")
end

if (try Pkg.installed("MathProgBase") != nothing catch false end)
    include("test_mathprogbase.jl")
end
