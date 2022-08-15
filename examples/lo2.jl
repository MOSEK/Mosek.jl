# Copyright : Copyright (c) 2022 MOSEK ApS
#
# File :      lo2.jl
#
# Purpose :   Demonstrates how to solve a small linear
#             optimization problem using the MOSEK Java API.

##TAG:begin-code
using Mosek


let numcon = 3,
    numvar = 4,
    NUMANZ = 9,

    c = [3.0, 1.0, 5.0, 1.0],
##TAG:begin-A
    aptrb = [1,4,8],
    aptre = [4,8,10],
    asub = [ 1, 2, 3,
             1, 2, 3, 4,
             2, 4 ],
    aval = [ 3.0, 1.0, 2.0 ,
             2.0, 1.0, 3.0, 1.0,
             2.0, 3.0 ],
##TAG:end-A
    bkc  = [ MSK_BK_FX,
             MSK_BK_LO,
             MSK_BK_UP ],
    blc  = [30.0,
            15.0,
            -Inf
            ],
    buc  = [30.0,
            +Inf,
            25.0
            ],
##TAG:begin-bx
    bkx  = [ MSK_BK_LO,
             MSK_BK_RA,
             MSK_BK_LO,
             MSK_BK_LO ],
    blx  = [ 0.0,
             0.0,
             0.0,
             0.0 ],
    bux  = [ Inf,
             10.0,
             Inf,
             Inf ]
##TAG:end-bx

    maketask() do task
        putstreamfunc(task,MSK_STREAM_LOG,msg -> print(msg))

##TAG:begin-append
        # Append 'numcon' empty constraints.
        # The constraints will initially have no bounds.
        appendcons(task,numcon)

        # Append 'numvar' variables.
        # The variables will initially be fixed at zero (x=0).
        appendvars(task,numvar)
##TAG:end-append

        for j in 1:numvar
            # Set the linear term c_j in the objective.
            putcj(task,j,c[j])
            # Set the bounds on variable j.
            # blx[j] <= x_j <= bux[j]
            putvarbound(task,j, bkx[j], blx[j], bux[j])
        end

        # Set the bounds on constraints.
        # for i=1, ...,numcon : blc[i] <= constraint i <= buc[i]
        for i in 1:numcon
            putconbound(task,i, bkc[i], blc[i], buc[i])
        end

        #Input rows of A #
        putarowslice(task,
                     1,numcon+1,
                     aptrb,aptre,
                     asub,
                     aval)

        # A maximization problem
        putobjsense(task,MSK_OBJECTIVE_SENSE_MAXIMIZE)

        # Solve the problem
        trm = optimize(task)

        # Print a summary containing information
        #   about the solution for debugging purposes
        solutionsummary(task,MSK_STREAM_MSG)

##TAG:begin-getsolutionstatus
        # Get status information about the solution 
        solsta = getsolsta(task,MSK_SOL_BAS)
##TAG:end-getsolutionstatus

        x = getxx(task,MSK_SOL_BAS) # Basic solution.
        ##TAG:ASSERT:begin-check-solution
        @assert maximum(abs.(x-[0.0, 0.0, 15.0, 8.333333333333334])) < 1e-7
        ##TAG:ASSERT:end-check-solution

        if solsta == MSK_SOL_STA_OPTIMAL
            println("Optimal primal solution")
            println("  x = $x")
        elseif ( solsta == MSK_SOL_STA_DUAL_INFEAS_CER ||
                 solsta == MSK_SOL_STA_PRIM_INFEAS_CER)
            println("Primal or dual infeasibility.")
        elseif solsta == MSK_SOL_STA_UNKNOWN
            println("Unknown solution status.")
        else
            println("Other solution status")
        end
    end
end
