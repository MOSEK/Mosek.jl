using Mosek

include(joinpath(Pkg.dir("MathProgBase"),"test","linprog.jl"))
linprogtest(MosekSolver())
