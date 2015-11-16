using Mosek

try
    Pkg.installed("Convex")
    using Convex
    using FactCheck

    set_default_solver(MosekSolver());
    include(joinpath(Pkg.dir("Convex"),"test","test_sdp.jl"))
    FactCheck.exitstatus()
catch
    # If convex not installed, skip this test
end
