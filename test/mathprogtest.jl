using Mosek

include(joinpath(Pkg.dir("MathProgBase"),"test","linprog.jl"))
include(joinpath(Pkg.dir("MathProgBase"),"test","linproginterface.jl"))
include(joinpath(Pkg.dir("MathProgBase"),"test","quadprog.jl"))
include(joinpath(Pkg.dir("MathProgBase"),"test","mixintprog.jl"))
linprogtest(MosekSolver())
linprogsolvertest(MosekSolver())
quadprogtest(MosekSolver())
socptest(MosekSolver())
mixintprogtest(MosekSolver())
