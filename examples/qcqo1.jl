#
#  File:    qcqo1.jl
#
#  Purpose: Demonstrates how to solve small quadratic and quadratically 
#           constrained optimization problem using the MOSEK Python API.
##

using Mosek
using Printf
# Since the actual value of Infinity is ignores, we define it solely
# for symbolic purposes:


# Set up and input bounds and linear coefficients
bkc   = [ MSK_BK_LO ]
blc   = [ 1.0 ]
buc   = [ Inf ]
  
bkx   = [ MSK_BK_LO
          MSK_BK_LO
          MSK_BK_LO ]
blx   = [ 0.0,  0.0, 0.0 ]
bux   = [ Inf,  Inf, Inf ]

c     = [ 0.0, -1.0, 0.0 ]

asub  = [ 1 ,2, 3 ]
aval  = [ 1.0, 1.0, 1.0 ]

numvar = length(bkx)
numcon = length(bkc)


# Create a task
maketask() do task
    # Append 'numcon' empty constraints.
    # The constraints will initially have no bounds. 
    appendcons(task,numcon)
    
    #Append 'numvar' variables.
    # The variables will initially be fixed at zero (x=0). 
    appendvars(task,numvar)

    #Optionally add a constant term to the objective. 
    putcfix(task,0.0)
    # Set the linear term c_j in the objective.
    putclist(task,[1:numvar;],c)

    # Set the bounds on variable j
    # blx[j] <= x_j <= bux[j] 
    putvarboundslice(task,1,numvar+1,bkx,blx,bux)
    # Input column j of A 
    putarow(task,1,asub,aval)

    putconbound(task,1,bkc[1],blc[1],buc[1])
    
    # Set up and input quadratic objective

    qsubi = [ 1,   2,    3,   3   ]
    qsubj = [ 1,   2,    1,   3   ]
    qval  = [ 2.0, 0.2, -1.0, 2.0 ]

    putqobj(task,qsubi,qsubj,qval)

    # The lower triangular part of the Q^0
    # matrix in the first constraint is specified.
    # This corresponds to adding the term
    # - x0^2 - x1^2 - 0.1 x2^2 + 0.2 x0 x2

    qsubi = [  1,    2,    3,   3   ]
    qsubj = [  1,    2,    3,   1   ]
    qval  = [ -2.0, -2.0, -0.2, 0.2 ]

    # put Q^0 in constraint with index 0. 

    putqconk(task,1, qsubi,qsubj, qval) 

    # Input the objective sense (minimize/maximize)
    putobjsense(task,MSK_OBJECTIVE_SENSE_MINIMIZE)
    
    # Optimize the task
    optimize(task,"mosek://solve.mosek.com:30080")
    # Print a summary containing information
    # about the solution for debugging purposes
    solutionsummary(task,MSK_STREAM_MSG)
    prosta = getprosta(task,MSK_SOL_ITR)
    solsta = getsolsta(task,MSK_SOL_ITR)

    if solsta == MSK_SOL_STA_OPTIMAL
        # Output a solution
        xx = getxx(task,MSK_SOL_ITR)
        @printf("Optimal solution: %s\n", xx')
    elseif solsta == MSK_SOL_STA_DUAL_INFEAS_CER
        println("Primal or dual infeasibility.\n")
    elseif solsta == MSK_SOL_STA_PRIM_INFEAS_CER
        println("Primal or dual infeasibility.\n")
    elseif solsta == MSK_SOL_STA_UNKNOWN
        println("Unknown solution status")
    else
        println("Other solution status")
    end
end
