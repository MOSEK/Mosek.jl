using Base.Test

include("apitest.jl")
#include("mathprogtest.jl")
include("testexamples.jl")
#include("mathprogtestextra.jl")


if false
  # Not yet compatible with the switch to MOI 
  include("testconvex.jl")
end

include("test_mathoptinterface.jl")
