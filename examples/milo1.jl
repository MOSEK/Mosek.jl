##
#    File:    milo1.jl
#
#    Purpose:  Demonstrates how to solve a small mixed
#              integer linear optimization problem using the MOSEK Python API.
##

using Mosek
using Printf, SparseArrays

# Define a stream printer to grab output from MOSEK
printstream(msg::String) = print(msg)


bkc = [ MSK_BK_UP, MSK_BK_LO  ]
blc = [      -Inf,      -4.0  ]
buc = [     250.0,       Inf  ]

bkx = [ MSK_BK_LO, MSK_BK_LO  ]
blx = [       0.0,       0.0  ]
bux = [       Inf,       Inf  ]

c   = [       1.0,      0.64 ]

A    = sparse( [ 1, 1, 2, 2],
               [ 1, 2, 1, 2],
               [ 50.0, 31.0,
                 3.0, -2.0] )

numvar = length(bkx)
numcon = length(bkc)


# Create a task
maketask() do task
    # Attach a printer to the task
    putstreamfunc(task,MSK_STREAM_LOG,printstream)

    # Append 'numcon' empty constraints.
    # The constraints will initially have no bounds.
    appendcons(task,numcon)

    #Append 'numvar' variables.
    # The variables will initially be fixed at zero (x=0).
    appendvars(task,numvar)

    # Set the linear term c_j in the objective.
    putclist(task,[1:numvar;],c)

    # Set the bounds on variables
    # blx[j] <= x_j <= bux[j]
    putvarboundslice(task,1,numvar+1,bkx,blx,bux)

    # Input columns of A
    putacolslice(task,1,numvar+1, A.colptr[1:numvar],A.colptr[2:numvar+1],A.rowval,A.nzval)

    putconboundslice(task,1,numcon+1,bkc,blc,buc)

    # Input the objective sense (minimize/maximize)
    putobjsense(task,MSK_OBJECTIVE_SENSE_MAXIMIZE)

    # Define variables to be integers
    putvartypelist(task,[ 1, 2 ],
                   [ MSK_VAR_TYPE_INT, MSK_VAR_TYPE_INT ])

    # Set max solution time
    putdouparam(task,MSK_DPAR_MIO_MAX_TIME, 60.0)

    # Optimize the task
    optimize(task,"mosek://solve.mosek.com:30080")

    writedata(task,"milo1.ptf")

    # Print a summary containing information
    # about the solution for debugging purposes
    solutionsummary(task,MSK_STREAM_MSG)

    prosta = getprosta(task,MSK_SOL_ITG)
    solsta = getsolsta(task,MSK_SOL_ITG)

    if solsta == MSK_SOL_STA_INTEGER_OPTIMAL
        # Output a solution
        xx = getxx(task,MSK_SOL_ITG)
        @printf("Optimal solution: %s\n", xx')
        @assert maximum(abs.(xx-[5.0, 0.0])) < 1e-7
    elseif solsta == MSK_SOL_STA_UNKNOWN
        println("Unknown solution status")
    else
        println("Other solution status")
    end

end
