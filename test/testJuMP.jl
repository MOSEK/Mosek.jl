import Mosek
try
    Pkg.installed("JuMP")
    using JuMP
    using Base.Test
    using Base.Test

    solver = Mosek.MosekSolver(QUIET=false)
    lp_solvers         = Any[solver]
    ip_solvers         = Any[solver]
    quad_solvers       = Any[solver]
    quad_mip_solvers   = Any[solver]
    convex_nlp_solvers = Any[solver]
    sdp_solvers        = Any[solver]
    semi_solvers       = Any[]
    sos_solvers        = Any[]
    lazy_solvers, cut_solvers, heur_solvers, info_solvers = Any[], Any[], Any[], Any[]
    minlp_solvers      = Any[]

    # Solver-dependent tests
    include(joinpath(Pkg.dir("JuMP"),"test","model.jl"))
    include(joinpath(Pkg.dir("JuMP"),"test","probmod.jl"))
    include(joinpath(Pkg.dir("JuMP"),"test","qcqpmodel.jl"))
    include(joinpath(Pkg.dir("JuMP"),"test","nonlinear.jl"))
    include(joinpath(Pkg.dir("JuMP"),"test","sdp.jl"))
catch
    # If convex not installed, skip this test
end
