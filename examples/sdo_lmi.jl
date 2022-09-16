#
# Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
#
# File :      sdo_lmi.jl
#
# Purpose :   To solve a problem with an LMI and an affine conic constrained problem with a PSD term
#
#                 minimize    Tr [1, 0; 0, 1]*X + x(1) + x(2) + 1
#
#                 subject to  Tr [0, 1; 1, 0]*X - x(1) - x(2) >= 0
#                             x(1) [0, 1; 1, 3] + x(2) [3, 1; 1, 0] - [1, 0; 0, 1] >> 0
#                             X >> 0
#

using Mosek

let numafe      = 4,  # Number of affine expressions.
    numvar      = 2,  # Number of scalar variables
    dimbarvar = [2],         # Dimension of semidefinite cone
    lenbarvar = [2 * (2 + 1) / 2], # Number of scalar SD variables
    barc_j  = [1, 1],
    barc_k  = [1, 2],
    barc_l  = [1, 2],
    barc_v  = [1.0, 1.0],

    afeidx  = [1, 1, 2, 3, 3, 4],
    varidx  = [1, 2, 2, 1, 2, 1],
    f_val   = Float64[-1, -1, 3, sqrt(2.0), sqrt(2.0), 3],
    g       = Float64[0, -1, 0, -1],

    barf_i = [1, 1],
    barf_j = [1, 1],
    barf_k = [1, 2],
    barf_l = [1, 1],
    barf_v = [0.0, 1.0]

    maketask() do task
        # Append 'NUMAFE' empty affine expressions.
        appendafes(task,numafe)

        # Append 'NUMVAR' variables.
        # The variables will initially be fixed at zero (x=0).
        appendvars(task,numvar)

        # Append 'NUMBARVAR' semidefinite variables.
        appendbarvars(task,dimbarvar)

        # Optionally add a constant term to the objective.
        putcfix(task,1.0)

        # Set the linear term c_j in the objective.
        putcj(task,1, 1.0)
        putcj(task,2, 1.0)

        for j in 1:numvar
            putvarbound(task,j, MSK_BK_FR, -0.0, 0.0)
        end
        # Set the linear term barc_j in the objective.
        putbarcblocktriplet(task, barc_j, barc_k, barc_l, barc_v)

        # Set up the affine conic constraints

        # Construct the affine expressions
        # F matrix
        putafefentrylist(task,afeidx, varidx, f_val)
        # g vector
        putafegslice(task,1, 5, g)

        # barF block triplets
        putafebarfblocktriplet(task,barf_i, barf_j, barf_k, barf_l, barf_v)

        # Append R+ domain and the corresponding ACC
        appendacc(task,appendrplusdomain(task,1), [1],[0.0])
        # Append SVEC_PSD domain and the corresponding ACC
        appendacc(task,appendsvecpsdconedomain(task,3), [2,3,4], [0.0,0.0,0.0])

        # Run optimizer
        optimize(task)

        # Print a summary containing information
        # about the solution for debugging purposes
        solutionsummary(task,MSK_STREAM_MSG)

        solsta = getsolsta(task,MSK_SOL_ITR)

        if solsta == MSK_SOL_STA_OPTIMAL
            xx = getxx(task,MSK_SOL_ITR)
            barx = getbarxj(task,MSK_SOL_ITR,1);    # Request the interior solution.
            println("Optimal primal solution, x = $xx, barx = $barx")
            @assert maximum(abs.(xx-[1.0, 1.0])) < 1e-6
            @assert maximum(abs.(barx- [1.0, 1.0, 1.0])) < 1e-6
        elseif solsta == MSK_SOL_STA_PRIM_INFEAS_CER || solsta == MSK_SOL_STA_DUAL_INFEAS_CER
            println("Primal or dual infeasibility certificate found.")
        elseif solsta == MSK_SOL_STA_UNKNOWN
            println("The status of the solution could not be determined.")
        else
            println("Other solution status.")
        end
    end
end
