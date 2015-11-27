module MosekMathProgTest

import Mosek
using MathProgBase
using MathProgBase.SolverInterface

include(joinpath(Pkg.dir("MathProgBase"),"test","linprog.jl"))
include(joinpath(Pkg.dir("MathProgBase"),"test","quadprog.jl"))
include(joinpath(Pkg.dir("MathProgBase"),"test","mixintprog.jl"))
include(joinpath(Pkg.dir("MathProgBase"),"test","conicinterface.jl"))
include(joinpath(Pkg.dir("MathProgBase"),"test","nlp.jl"))



linprogtest(Mosek.MosekSolver())
quadprogtest(Mosek.MosekSolver())
#socptest(MosekSolver())
mixintprogtest(Mosek.MosekSolver())
convexnlptest(Mosek.MosekSolver())
coniclineartest(Mosek.MosekSolver(),duals=true)
conicSOCtest(Mosek.MosekSolver(),duals=true)

end
