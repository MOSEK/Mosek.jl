##
#   Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
#
#   File :      qo1.jl
#
#   Purpose: Demonstrate how to solve a quadratic
#            optimization problem using the MOSEK Python API.
##

using Mosek
using Printf, SparseArrays

# Define a stream printer to grab output from MOSEK

bkc   = [ MSK_BK_LO ]
blc   = [ 1.0 ]
buc   = [ Inf ]
  
bkx   = [ MSK_BK_LO, MSK_BK_LO, MSK_BK_LO ]
blx   = [ 0.0,  0.0, 0.0 ]
bux   = [ Inf,  Inf, Inf ]
numvar = length(bkx)
numcon = length(bkc)

c     = [ 0.0, -1.0, 0.0 ]
A     = sparse( [ 1, 1, 1 ], 
                [ 1, 2, 3 ], 
                [ 1.0, 1.0, 1.0 ],
                numcon, numvar )

maketask() do task
    putstreamfunc(task,MSK_STREAM_LOG,msg -> print(msg))

    # Append 'numcon' empty constraints.
    # The constraints will initially have no bounds.  
    appendcons(task,numcon)

    # Append 'numvar' variables.
    # The variables will initially be fixed at zero (x=0). 
    appendvars(task,numvar)

    # Set the linear term c_j in the objective.
    putclist(task,[1:numvar;],c)

    # Set the bounds on variable j
    # blx[j] <= x_j <= bux[j] 
    putvarboundslice(task,1,numvar+1,bkx,blx,bux)

    putacolslice(task,1,numvar+1,
                 A.colptr[1:numvar],A.colptr[2:numvar+1],
                 A.rowval,A.nzval)

    # Set up and input quadratic objective
    qsubi = [ 1,   2,    3,   3   ]
    qsubj = [ 1,   2,    1,   3   ]
    qval  = [ 2.0, 0.2, -1.0, 2.0 ]

    putqobj(task,qsubi,qsubj,qval)

    putobjsense(task,MSK_OBJECTIVE_SENSE_MINIMIZE)

    # Optimize
    r = optimize(task,"mosek://solve.mosek.com:30080")
    # Print a summary containing information
    # about the solution for debugging purposes
    solutionsummary(task,MSK_STREAM_MSG)

    prosta = getprosta(task,MSK_SOL_ITR)
    solsta = getsolsta(task,MSK_SOL_ITR)

    if solsta == MSK_SOL_STA_OPTIMAL
        xx = getxx(task,MSK_SOL_ITR)
        println("Optimal solution:")
        println(xx)
        @assert maximum(abs.(xx-[0.00015777655846067916, 4.999999955843769, 0.00015779898859216525])) < 1e-5
    elseif solsta in [ MSK_SOL_STA_DUAL_INFEAS_CER,
                       MSK_SOL_STA_PRIM_INFEAS_CER ]
        println("Primal or dual infeasibility certificate found.\n")
    elseif solsta == MSK_SOL_STA_UNKNOWN
        println("Unknown solution status")
    else
        @printf("Other solution status (%d)\n",solsta)
    end

end
