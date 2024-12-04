#
#  Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
#
#  File :      acc2.jl
#
#  Purpose :   Tutorial example for affine conic constraints.
#              Models the problem:
#
#              maximize c^T x
#              subject to  sum(x) = 1
#                          gamma >= |Gx+h|_2
#
#              This version inputs the linear constraint as an affine conic constraint.
##
using Mosek


# Define problem data
n = 3
k = 2

# Create a task
maketask() do task
    # Use remote server: putoptserverhost(task,"http://solve.mosek.com:30080")
    # Attach a printer to the task
    putstreamfunc(task,MSK_STREAM_LOG,msg -> print(msg))

    # Create n free variables
    appendvars(task,n)
    putvarboundsliceconst(task,1, n+1, MSK_BK_FR, -Inf, Inf)

    # Set up the objective
    c = Float64[2, 3, -1]
    putobjsense(task,MSK_OBJECTIVE_SENSE_MAXIMIZE)
    putclist(task,[1:n...], c)

    # Set AFE rows representing the linear constraint
    appendafes(task,1)
    putafefrow(task,1, [1:n...], ones(n))
    putafeg(task,1, -1.0)

    # Set AFE rows representing the quadratic constraint
    appendafes(task,k + 1)
    putafefrow(task,3,          # afeidx, row number
                    [1, 2],     # varidx, column numbers
                    [1.5, 0.1]) # values
    putafefrow(task,4,          # afeidx, row number
                    [1, 3],     # varidx, column numbers
                    [0.3, 2.1]) # values

    h     = [0, 0.1]
    gamma = 0.03
    putafeg(task,2, gamma)
    putafegslice(task,3, k+2+1, h)

    # Define domains
    zeroDom = appendrzerodomain(task,1)
    quadDom = appendquadraticconedomain(task,k + 1)

    # Append affine conic constraints
    appendacc(task,zeroDom,    # Domain index
              [1],             # Indices of AFE rows
              nothing)         # Ignored
    appendacc(task,quadDom,    # Domain index
              [2,3,4],         # Indices of AFE rows
              nothing)         # Ignored

    # Solve and retrieve solution
    optimize(task)
    writedata(task,"acc2.ptf")
    @assert getsolsta(task,MSK_SOL_ITR) == MSK_SOL_STA_OPTIMAL
    xx = getxx(task,MSK_SOL_ITR)
    if getsolsta(task,MSK_SOL_ITR) == MSK_SOL_STA_OPTIMAL
        println("Solution: $xx")
    end

    # Demonstrate retrieving activity of ACC
    activity = evaluateacc(task,MSK_SOL_ITR,2)
    println("Activity of quadratic ACC:: $activity")

    # Demonstrate retrieving the dual of ACC
    doty = getaccdoty(task,MSK_SOL_ITR,2)
    println("Dual of quadratic ACC:: $doty")

end
