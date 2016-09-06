module MosekMathProgTest

import Mosek
using MathProgBase
using MathProgBase.SolverInterface

include(joinpath(Pkg.dir("MathProgBase"),"test","linprog.jl"))
include(joinpath(Pkg.dir("MathProgBase"),"test","quadprog.jl"))
include(joinpath(Pkg.dir("MathProgBase"),"test","mixintprog.jl"))
include(joinpath(Pkg.dir("MathProgBase"),"test","conicinterface.jl"))
include(joinpath(Pkg.dir("MathProgBase"),"test","linproginterface.jl"))
include(joinpath(Pkg.dir("MathProgBase"),"test","nlp.jl"))


solver = Mosek.MosekSolver(QUIET=true)

linprogtest(solver)
mixintprogtest(solver)
quadprogtest(solver)
convexnlptest(solver)

#socptest(MosekSolver())
coniclineartest(solver,duals=true)
conicSOCtest(solver,duals=true)
conicSOCRotatedtest(solver,duals=true)
conicSDPtest(solver,duals=true)

    
linprogsolvertest(solver)

# linprogsolvertestextra(solver)
end
