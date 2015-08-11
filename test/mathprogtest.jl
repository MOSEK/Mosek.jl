using Mosek
using MathProgBase
using MathProgBase.SolverInterface

include(joinpath(Pkg.dir("MathProgBase"),"test","linprog.jl"))
include(joinpath(Pkg.dir("MathProgBase"),"test","quadprog.jl"))
include(joinpath(Pkg.dir("MathProgBase"),"test","mixintprog.jl"))
include(joinpath(Pkg.dir("MathProgBase"),"test","nlp.jl"))
include(joinpath(Pkg.dir("MathProgBase"),"test","conicinterface.jl"))



linprogtest(MosekSolver())
quadprogtest(MosekSolver())
socptest(MosekSolver())
mixintprogtest(MosekSolver())
convexnlptest(MosekSolver())
coniclineartest(MosekSolver(),duals=true)
conicSOCtest(MosekSolver(),duals=true)


