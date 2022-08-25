# To ensure that tests can run on travis we have to do a little
# hackadoodle here. The tests require a license file. We include
# a license file that is only valid for one day (the day when
# change is submitted).
# If there is no valid license file, we default to that file.

using Test

include("apitest.jl")
include("testexamples.jl")
#include("testshow.jl")
#include("test_mathprogbase.jl")
include("issues.jl")
