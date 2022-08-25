##
#  Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
#
#  File :      lo1.jl
#
#  Purpose :   Demonstrates how to solve small linear
#              optimization problem using the MOSEK Python API.
##
using Mosek
using Printf, SparseArrays

############################
## Define problem data

bkc = [MSK_BK_FX
       MSK_BK_LO
       MSK_BK_UP]

# Bound values for constraints
blc = [30.0, 15.0, -Inf]
buc = [30.0, +Inf, 25.0]

# Bound keys for variables
bkx = [ MSK_BK_LO
        MSK_BK_RA
        MSK_BK_LO
        MSK_BK_LO ]

# Bound values for variables
blx = [   0.0,  0.0,    0.0,    0.0]
bux = [+Inf, 10.0, +Inf, +Inf]

numvar = length(bkx)
numcon = length(bkc)

# Objective coefficients
c = [ 3.0, 1.0, 5.0, 1.0 ] 

# Below is the sparse representation of the A
# matrix stored by column. 
A = sparse([1, 2, 1, 2, 3, 1, 2, 2, 3], 
           [1, 1, 2, 2, 2, 3, 3, 4, 4], 
           [3.0, 2.0, 1.0, 1.0, 2.0, 2.0, 3.0, 1.0, 3.0 ],
           numcon,numvar)

############################

maketask() do task
    putstreamfunc(task,MSK_STREAM_LOG,msg -> print(msg))

    putobjname(task,"lo1")

    # Append 'numcon' empty constraints.
    # The constraints will initially have no bounds. 
    appendcons(task,numcon)
    for i=1:numcon
        putconname(task,i,@sprintf("c%02d",i))
    end

    # Append 'numvar' variables.
    # The variables will initially be fixed at zero (x=0). 
    appendvars(task,numvar)
    for j=1:numvar
        putvarname(task,j,@sprintf("x%02d",j))
    end

    putclist(task,[1,2,3,4], c)
    putacolslice(task,1,numvar+1,A)

    putvarboundslice(task, 1, numvar+1, bkx,blx,bux)

    # Set the bounds on constraints.
    # blc[i] <= constraint_i <= buc[i]
    putconboundslice(task,1,numcon+1,bkc,blc,buc)

    # Input the objective sense (minimize/maximize)
    putobjsense(task,MSK_OBJECTIVE_SENSE_MAXIMIZE)

    # Solve the problem
    optimize(task,"mosek://solve.mosek.com:30080")

    # Print a summary containing information
    # about the solution for debugging purposes
    solutionsummary(task,MSK_STREAM_MSG)

    # Get status information about the solution
    solsta = getsolsta(task,MSK_SOL_BAS)

    if solsta == MSK_SOL_STA_OPTIMAL
        xx = getxx(task,MSK_SOL_BAS)
        print("Optimal solution:")
        println(xx)

        @assert maximum(abs.(xx-[0, 0, 15, 8.333333333333334])) < 1e-7
    elseif solsta in [ MSK_SOL_STA_DUAL_INFEAS_CER,
                       MSK_SOL_STA_PRIM_INFEAS_CER ]
        println("Primal or dual infeasibility certificate found.\n")
    elseif solsta == MSK_SOL_STA_UNKNOWN
        println("Unknown solution status")
    else
        @printf("Other solution status (%d)\n",solsta)
    end
end
