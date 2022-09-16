##
#  Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
#
#  File :      ceo1.jl
#
#  Purpose :   Demonstrates how to solve small conic exponential
#              optimization problem using the MOSEK Python API.
##

using Mosek
using Printf, SparseArrays



# Create a task
maketask() do task
    # Attach a printer to the task
    putstreamfunc(task,MSK_STREAM_LOG,msg -> print(msg))

    c = [1.0, 1.0, 0.0]
    a = [1.0, 1.0, 1.0]
    numvar = 3
    numcon = 1
    
    # Append 'numcon' empty constraints.
    # The constraints will initially have no bounds.
    appendcons(task,numcon)

    # Append 'numvar' variables.
    # The variables will initially be fixed at zero (x=0).
    appendvars(task,numvar)

    # Set up the linear part of the problem
    putcslice(task,1, numvar+1, c)
    putarow(task,1, [1, 2, 3], a)
    putvarboundsliceconst(task,1, numvar+1, MSK_BK_FR, -Inf, Inf)
    putconbound(task,1, MSK_BK_FX, 1.0, 1.0)
    # Add a conic constraint
    # Create a 3x3 identity matrix F
    appendafes(task,3)
    putafefentrylist(task,
                     [1, 2, 3],         # Rows
                     [1, 2, 3],         # Columns
                     ones(3))

    # Exponential cone (x(0),x(1),x(2)) \in EXP
    expdomain = appendprimalexpconedomain(task)
    appendacc(task,
              expdomain,               # Domain
              [1, 2, 3],               # Rows from F
              zeros(3))                # Unused

    # Input the objective sense (minimize/maximize)
    putobjsense(task,MSK_OBJECTIVE_SENSE_MINIMIZE)

    # Optimize the task
    optimize(task)
    # Print a summary containing information
    # about the solution for debugging purposes
    solutionsummary(task,MSK_STREAM_MSG)
    prosta = getprosta(task,MSK_SOL_ITR)
    solsta = getsolsta(task,MSK_SOL_ITR)

    # Output a solution
    xx = getxx(task,MSK_SOL_ITR)

    if solsta == MSK_SOL_STA_OPTIMAL
        println("Optimal solution: $xx")
    elseif solsta == MSK_SOL_STA_DUAL_INFEAS_CER
        println("Primal or dual infeasibility.")
    elseif solsta == MSK_SOL_STA_PRIM_INFEAS_CER
        println("Primal or dual infeasibility.")
    elseif solsta == MSK_SOL_STA_UNKNOWN
        println("Unknown solution status")
    else
        println("Other solution status")
    end

    @assert maximum(abs.(xx-[0.6117882543880403, 0.17040004803746528, 0.21781169885758184])) < 1e-7

end
