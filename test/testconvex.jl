using Mosek

try
    Pkg.installed("Convex")
    Pkg.installed("FactCheck")
    using Convex
    using FactCheck

    set_default_solver(MosekSolver());
    include(joinpath(Pkg.dir("Convex"),"test","test_sdp.jl"))
    FactCheck.exitstatus()
catch
    # If convex+factcheck not installed, skip this test
end
