using Mosek

try
    Pkg.installed("Convex")
    using Convex
    using Base.Test

    set_default_solver(MosekSolver(QUIET=true));
    include(joinpath(Pkg.dir("Convex"),"test","test_sdp.jl"))
catch
    # If convex not installed, skip this test
end
