
module TestMisc

import Mosek
using Convex, FactCheck

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
