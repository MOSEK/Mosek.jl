#   Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
#
#   File :      solutionquality.jl
#
#   Purpose :   To demonstrate how to examine the quality of a solution.

using Mosek

if length(ARGS) < 1
    println("Usage: solutionquality FILENAME")
else
    filename = ARGS[1]

    maketask() do task
        #putstreamfunc(task,MSK_STREAM_LOG,msg -> print(msg))

        # We assume that a problem file was given as the first command
        # line argument (received in `args')
        readdata(task,filename)
        # Solve the problem
        optimize(task)
        # System.Out.Println (a summary of the solution
        solutionsummary(task,MSK_STREAM_LOG)

        solsta = getsolsta(task,MSK_SOL_BAS)

        (pobj,
         pviolcon,
         pviolvar,
         pviolbarvar,
         pviolcones,
         pviolitg,
         dobj,
         dviolcon,
         dviolvar,
         dviolbarvar,
         dviolcones) = getsolutioninfo(task,MSK_SOL_BAS)

        if solsta == MSK_SOL_STA_OPTIMAL
            abs_obj_gap     = abs(dobj - pobj)
            rel_obj_gap     = abs_obj_gap / (1.0 + min(abs(pobj), abs(dobj)))
            max_primal_viol = max(pviolcon, pviolvar)
            max_primal_viol = max(max_primal_viol, pviolbarvar)
            max_primal_viol = max(max_primal_viol, pviolcones)

            max_dual_viol = max(dviolcon,      dviolvar)
            max_dual_viol = max(max_dual_viol, dviolbarvar)
            max_dual_viol = max(max_dual_viol, dviolcones)

            # Assume the application needs the solution to be within
            #    1e-6 ofoptimality in an absolute sense. Another approach
            #   would be looking at the relative objective gap

            println("Customized solution information.")
            println("  Absolute objective gap: $abs_obj_gap")
            println("  Relative objective gap: $rel_obj_gap")
            println("  Max primal violation  : $max_primal_viol")
            println("  Max dual violation    : $max_dual_viol")

            accepted = true

            if rel_obj_gap > 1e-6
                println("Warning: The relative objective gap is LARGE.")
                accepted = false
            end

            # We will accept a primal infeasibility of 1e-8 and
            # dual infeasibility of 1e-6. These number should chosen problem
            # dependent.
            if max_primal_viol > 1e-8
                println("Warning: Primal violation is too LARGE")
                accepted = false
            end

            if max_dual_viol > 1e-6
                println("Warning: Dual violation is too LARGE.")
                accepted = false
            end

            if accepted
                numvar = getnumvar(task)
                println("Optimal primal solution")

                xx = getxxslice(task,MSK_SOL_BAS,1,numvar+1)
                println("  xx = $xx")
            else
                # print etailed information about the solution
                analyzesolution(task,MSK_STREAM_LOG, MSK_SOL_BAS)
            end
        elseif solsta == MSK_SOL_STA_DUAL_INFEAS_CER || solsta == MSK_SOL_STA_PRIM_INFEAS_CER
            println("Primal or dual infeasibility certificate found.")
        elseif solsta == MSK_SOL_STA_UNKNOWN
            println("The status of the solution is unknown.")
        else
            println("Other solution status")
        end
    end
end
