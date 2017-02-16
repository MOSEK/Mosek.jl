import Mosek

include(joinpath(Pkg.dir("MathProgBase"),"test","linprog.jl"))
linprogtest(Mosek.MosekSolver())

include(joinpath(Pkg.dir("MathProgBase"),"test","mixintprog.jl"))
mixintprogtest(Mosek.MosekSolver())

include(joinpath(Pkg.dir("MathProgBase"),"test","quadprog.jl"))
quadprogtest(Mosek.MosekSolver())
#qpdualtest(Mosek.MosekSolver())
#socptest(Mosek.MosekSolver())

include(joinpath(Pkg.dir("MathProgBase"),"test","conicinterface.jl"))
coniclineartest(Mosek.MosekSolver(), duals=true)
# Test conic fallback for LPs
coniclineartest(Mosek.MosekSolver(), duals=true)

include(joinpath(Pkg.dir("MathProgBase"),"test","linproginterface.jl"))
linprogsolvertest(Mosek.MosekSolver())
linprogsolvertestextra(Mosek.MosekSolver())
# Test LP fallback for conics
linprogsolvertest(Mosek.MosekSolver())
