
module TestMiscConvex

import Mosek
using Convex, FactCheck
if false
    facts("[misc tests]") do
        context("convex psd") do
            # Borderline case:
            x = Variable(Positive())
            y = Variable((3, 3))
            p = minimize(x + y[1, 1],
                         isposdef(y),
                         x       >= 1,
                         y[2, 1] == 1)
            try
                solve!(p, Mosek.MosekSolver())
            catch
            end
            @fact p.status --> :Optimal
        end
    end
end
end # module


module TestMiscJuMP

import Mosek
using JuMP,FactCheck

TOL = 1e-8

facts("[model] Relaxation keyword argument to solve") do
    m = Model()
    @defVar(m, 1.5 <= y <= 2, Int)
    @defVar(m, z, Bin)
    @defVar(m, 0.5 <= w <= 1.5, Int)
    @defVar(m, 1 <= v <= 2)

    @setObjective(m, Min, y + z + w + v)
    
    # Force LP solver since not all MIP solvers
    # return duals (i.e. Cbc)
    setSolver(m, Mosek.MosekSolver())
    @fact solve(m, relaxation=true) --> :Optimal
    @fact getValue(y) --> 1.5
    @fact getValue(z) --> 0
    @fact getValue(w) --> 0.5
    @fact getValue(v) --> 1
    @fact getDual(y) --> 1
    @fact getDual(z) --> 1
    @fact getDual(w) --> 1
    @fact getDual(v) --> 1
    @fact getObjectiveValue(m) --> 1.5 + 0 + 0.5 + 1

    # Let JuMP choose solver again
    setSolver(m, JuMP.UnsetSolver())
    @fact solve(m) --> :Optimal
    @fact getValue(y) --> 2
    @fact getValue(z) --> 0
    @fact getValue(w) --> 1
    @fact getValue(v) --> 1
    @fact getObjectiveValue(m) --> 2 + 0 + 1 + 1

    @defVar(m, 1 <= x <= 2, SemiCont)
    @defVar(m, -2 <= t <= -1, SemiInt)

    addSOS1(m, [x, 2y, 3z, 4w, 5v, 6t])
    @setObjective(m, Min, x + y + z + w + v - t)

    @fact solve(m, relaxation=true) --> :Optimal

    @fact getValue(x) --> 0
    @fact getValue(y) --> 1.5
    @fact getValue(z) --> 0
    @fact getValue(w) --> 0.5
    @fact getValue(v) --> 1
    @fact getValue(t) --> 0
    @fact getObjectiveValue(m) --> 0 + 1.5 + 0 + 0.5 + 1 + 0
end

end # module
