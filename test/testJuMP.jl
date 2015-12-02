import Mosek
try
    Pkg.installed("JuMP")
    using JuMP
    using FactCheck
    using Base.Test

    solver = Mosek.MosekSolver(QUIET=false)
    lp_solvers         = Any[solver]
    ip_solvers         = Any[solver]
    quad_solvers       = Any[solver]
    convex_nlp_solvers = Any[solver]
    sdp_solvers        = Any[solver]
    semi_solvers       = Any[]
    sos_solvers        = Any[]
    lazy_solvers, cut_solvers, heur_solvers, info_solvers = Any[], Any[], Any[], Any[]
    minlp_solvers      = Any[]

    include(joinpath(Pkg.dir("JuMP"),"test","print.jl"))
    include(joinpath(Pkg.dir("JuMP"),"test","variable.jl"))
    include(joinpath(Pkg.dir("JuMP"),"test","expr.jl"))
    include(joinpath(Pkg.dir("JuMP"),"test","operator.jl"))
    include(joinpath(Pkg.dir("JuMP"),"test","macros.jl"))

    # Fuzzer of macros to build expressions
    include(joinpath(Pkg.dir("JuMP"),"test","fuzzer.jl"))

    # Load solvers
    include(joinpath(Pkg.dir("JuMP"),"test","solvers.jl"))

    # Solver-dependent tests
    include(joinpath(Pkg.dir("JuMP"),"test","model.jl"))
    include(joinpath(Pkg.dir("JuMP"),"test","probmod.jl"))
    #include(joinpath(Pkg.dir("JuMP"),"test","callback.jl"))
    include(joinpath(Pkg.dir("JuMP"),"test","qcqpmodel.jl"))
    include(joinpath(Pkg.dir("JuMP"),"test","nonlinear.jl"))
    include(joinpath(Pkg.dir("JuMP"),"test","sdp.jl"))

    FactCheck.exitstatus()
catch
    # If convex not installed, skip this test
end
