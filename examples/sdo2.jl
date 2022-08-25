#  Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
#
#  File :      sdo2.jl
#
#  Purpose :   Solves the semidefinite problem with two symmetric variables:
#
#                 min   <C1,X1> + <C2,X2>
#                 st.   <A1,X1> + <A2,X2> = b
#                             (X2)_{1,2} <= k
#
#                 where X1, X2 are symmetric positive semidefinite,
#
#                 C1, C2, A1, A2 are assumed to be constant symmetric matrices,
#                 and b, k are constants.
#

using Mosek

# Input data
let numcon    = 2,              # Number of constraints.
    numbarvar = 2,
    dimbarvar = Int32[3, 4],         # Dimension of semidefinite variables

    # Objective coefficients concatenated
    Cj = Int32[ 1, 1, 2, 2, 2, 2 ],   # Which symmetric variable (j)
    Ck = Int32[ 1, 3, 1, 2, 2, 3 ],   # Which entry (k,l)->v
    Cl = Int32[ 1, 3, 1, 1, 2, 3 ],
    Cv = [ 1.0, 6.0, 1.0, -3.0, 2.0, 1.0 ],

    # Equality constraints coefficients concatenated
    Ai = Int32[ 1, 1, 1, 1, 1, 1 ],   # Which constraint (i = 0)
    Aj = Int32[ 1, 1, 1, 2, 2, 2 ],   # Which symmetric variable (j)
    Ak = Int32[ 1, 3, 3, 2, 2, 4 ],   # Which entry (k,l)->v
    Al = Int32[ 1, 1, 3, 1, 2, 4 ],
    Av = [ 1.0, 1.0, 2.0, 1.0, -1.0, -3.0 ],

    # The second constraint - one-term inequality
    A2i = Int32[ 2 ],                        # Which constraint (i = 1)
    A2j = Int32[ 2 ],                        # Which symmetric variable (j = 1)
    A2k = Int32[ 2 ],                        # Which entry A(1,0) = A(0,1) = 0.5
    A2l = Int32[ 1 ],
    A2v = [ 0.5 ],

    bkc = [ MSK_BK_FX,
            MSK_BK_UP ],
    blc = [ 23.0, 0.0 ],
    buc = [ 23.0, -3.0 ]

    maketask() do task
        putstreamfunc(task,MSK_STREAM_LOG,msg -> print(msg))

        # Append numcon empty constraints.
        # The constraints will initially have no bounds.
        appendcons(task,numcon)

        # Append numbarvar semidefinite variables.
        appendbarvars(task,dimbarvar)

        # Set objective (6 nonzeros).
        putbarcblocktriplet(task,Cj, Ck, Cl, Cv)

        # Set the equality constraint (6 nonzeros).
        putbarablocktriplet(task,Ai, Aj, Ak, Al, Av)

        # Set the inequality constraint (1 nonzero).
        putbarablocktriplet(task,A2i, A2j, A2k, A2l, A2v)

        # Set constraint bounds
        putconboundslice(task,1, 3, bkc, blc, buc)

        # Run optimizer
        optimize(task)
        solutionsummary(task,MSK_STREAM_MSG)

        solsta = getsolsta(task,MSK_SOL_ITR)

        if solsta == MSK_SOL_STA_OPTIMAL
            # Retrieve the soution for all symmetric variables
            println("Solution (lower triangular part vectorized):")
            for i in 1:numbarvar
                barx = getbarxj(task,MSK_SOL_ITR, i)
                println("X$i: $barx")
                Xexpect = [[21.04706098136004, 0.0, 4.077117604213827, 5.5337653353154215, 0.0, 0.7897961639459908],
                           [5.053657059322597, -2.9999999957221952, 0.0, 0.0, 1.7808885486033874, 0.0, 0.0, 1.1353820371132691e-08, 0.0, 4.3074656422047645e-09]]
                @assert maximum(abs.(barx- Xexpect[i])) < 1e-3
            end
        elseif solsta == MSK_SOL_STA_DUAL_INFEAS_CER || if solsta == MSK_SOL_STA_PRIM_INFEAS_CER
            println("Primal or dual infeasibility certificate found.")
        elseif solsta == MSK_SOL_STA_UNKNOWN
            println("The status of the solution could not be determined.")
        else
            println("Other solution status.")
        end
        end
    end
end
