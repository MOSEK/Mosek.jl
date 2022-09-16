#   Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
#
#   File :      pow1.jl
#
#   Purpose: Demonstrates how to solve the problem
#
#     maximize x^0.2*y^0.8 + z^0.4 - x
#           st x + y + 0.5z = 2
#              x,y,z >= 0
#

using Mosek
printstream(msg::AbstractString) = print(msg)


csub = [ 4, 5, 1 ]
cval = [ 1.0, 1.0, -1.0]
asub = [ 1, 2, 3]
aval = [ 1.0, 1.0, 0.5]
numvar = 5
numcon = 1

# Create a task
maketask() do task
    putstreamfunc(task,MSK_STREAM_LOG,msg -> print(msg))

    appendcons(task,numcon)
    appendvars(task,numvar)

    # Set up the linear part of the problem
    putclist(task,csub, cval)
    putarow(task,1, asub, aval)
    putconbound(task,1, MSK_BK_FX, 2.0, 2.0)

    putvarboundsliceconst(task,1, numvar+1,MSK_BK_FR,-Inf,Inf)

    # Input the cones
    pc1 = appendprimalpowerconedomain(task,3,[0.2, 0.8])
    pc2 = appendprimalpowerconedomain(task,3,[4.0, 6.0])

    appendafes(task,6)
    putafefentrylist(task,
                     [1, 2, 3, 4, 6], # Rows
                     [1, 2, 4, 3, 5],  #Columns
                     [1.0, 1.0, 1.0, 1.0, 1.0])
    putafeg(task,5,1.0)

    # Append the two conic constraints
    appendacc(task,
              pc1,           # Domain
              [1, 2, 3],     # Rows from F
              [0.0,0.0,0.0]) # rhs offset
    appendacc(task,
              pc2,           # Domain
              [4, 5, 6],     # Rows from F
              [0.0,0.0,0.0]) # rhs offset

    # Input the objective sense (minimize/maximize)
    putobjsense(task,MSK_OBJECTIVE_SENSE_MAXIMIZE)

    # Optimize the task
    optimize(task,"mosek://solve.mosek.com:30080")
    writedata(task,"pow1.ptf")
    # Print a summary containing information
    # about the solution for debugging purposes
    solutionsummary(task,MSK_STREAM_MSG)
    prosta = getprosta(task,MSK_SOL_ITR)
    solsta = getsolsta(task,MSK_SOL_ITR)

    if solsta == MSK_SOL_STA_OPTIMAL
        # Output a solution
        xx = getxx(task,MSK_SOL_ITR)
        println("Optimal solution: $(xx[1:3])")
        @assert maximum(abs.(xx[1:3]-[0.063938, 0.78328, 2.305562])) < 1e-4
    elseif solsta == MSK_SOL_STA_DUAL_INFEAS_CER
        println("Primal or dual infeasibility.")
    elseif solsta == MSK_SOL_STA_PRIM_INFEAS_CER
        println("Primal or dual infeasibility.")
    elseif solsta == MSK_SOL_STA_UNKNOWN
        println("Unknown solution status")
    else
        println("Other solution status")
    end
end
