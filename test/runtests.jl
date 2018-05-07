using Base.Test

# To ensure that tests can run on travis we have to do a little
# hackadoodle here. The tests require a license file. We include
# a license file that is only valid for one day (the day when
# change is submitted).
# If there is no valid license file, we default to that file.

if haskey(ENV,"MOSEKLM_LICENSE_FILE")
    # that's nice
elseif haskey(ENV,"HOME")
    if isfile(joinpath(ENV["HOME"],"mosek","mosek.lic"))
        # our lucky day!
    else
        import Mosek
        Mosek.putlicensepath(Mosek.msk_global_env,joinpath(@__DIR__,"..","test",".dontuse-probablyexpired.lic"))
    end
end

include("apitest.jl")
#include("mathprogtest.jl")
include("testexamples.jl")
#include("mathprogtestextra.jl")


if (try Pkg.installed("Convex") != nothing catch false end)
    include("testconvex.jl")
end

if (try Pkg.installed("MathProgBase") != nothing catch false end)
    include("test_mathprogbase.jl")
end
