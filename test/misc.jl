
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

facts("[misc tests]") do
    if  false
        context("Test fixed variables don't leak through MPB") do
            mod = Model(solver=Mosek.MosekSolver())
            @defVar(mod, x >= 3, SemiInt)
            @defVar(mod, y >= 2, SemiInt)
            @addConstraint(mod, x + y >= 2.5)
            @setObjective(mod, Min, x+1.1y)
            solve(mod)

            @fact getValue(x) --> roughly(3.0, TOL)
            @fact getValue(y) --> 0.0
        end

    end

    facts("[probmod] Test buildInternalModel") do
        solver = Mosek.MosekSolver()
        context("With solver $(typeof(solver))") do
            m = Model(solver=solver)
            @defVar(m, x >= 0)
            @defVar(m, y >= 0)
            @addConstraint(m, x + y == 1)
            @setObjective(m, Max, y)
            buildInternalModel(m)
            @fact getInternalModel(m) --> not(nothing)
            @fact m.internalModelLoaded --> true
            stat = solve(m)
            Mosek.writedata(getInternalModel(m).task,"problem.opf")
            @fact stat --> :Optimal
            @fact getValue(x) --> roughly( 0.0, TOL)
            @fact getValue(y) --> roughly( 1.0, TOL)
            @fact getObjectiveValue(m) --> roughly(1.0, TOL)
            @fact getDual(x)  --> roughly(-1.0, TOL)
            @fact getDual(y)  --> roughly( 0.0, TOL)
        end

        facts("[qcqpmodel] Test simple normed problem") do
            context("With solver $(typeof(solver))") do
                m = Model(solver=solver);
                @defVar(m, x[1:3]);
                @addConstraint(m, 2norm2{x[i]-1, i=1:3} <= 2)
                @setObjective(m, Max, x[1]+x[2])

                @fact solve(m) --> :Optimal
                @fact getObjectiveValue(m) --> roughly(2+sqrt(2), 1e-5)
                @fact norm(getValue(x)-[1+sqrt(1/2),1+sqrt(1/2),1]) --> roughly(0, 1e-6)
            end
        end
    end
end

end # module
