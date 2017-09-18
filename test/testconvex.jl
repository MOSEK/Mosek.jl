using Convex
using FactCheck

solvers = Any[]

if isdir(Pkg.dir("Mosek"))
    using Mosek
    push!(solvers, MosekSolver(LOG=0))
end

for solver in solvers
    println("Running tests with $(solver):")
    set_default_solver(solver)
    println(get_default_solver())
    include(joinpath(Pkg.dir("Convex"),"test","runtests_single_solver.jl"))
end


FactCheck.exitstatus()
