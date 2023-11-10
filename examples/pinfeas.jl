#  File : pinfeas.jl
#
#  Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
#
#  Purpose: Demonstrates how to fetch a primal infeasibility certificate
#           for a linear problem
#
using Mosek

# Set up a simple linear problem from the manual for test purposes
function testProblem(func :: Function)
    maketask() do task
        appendvars(task,7)
        appendcons(task,7);
        putclist(task,
                 Int32[1,2,3,4,5,6,7],
                 Float64[1,2,5,2,1,2,1])
        putaijlist(task,
                   Int32[1,1,2,2,3,3,3,4,4,5,6,6,7,7],
                   Int32[1,2,3,4,5,6,7,1,5,2,3,6,4,7],
                   Float64[1,1,1,1,1,1,1,1,1,1,1,1,1,1])
        putconboundslice(task,
                         1, 8,
                         [MSK_BK_UP,MSK_BK_UP,MSK_BK_UP,MSK_BK_FX,MSK_BK_FX,MSK_BK_FX,MSK_BK_FX],
                         Float64[-Inf, -Inf, -Inf, 1100, 200, 500, 500],
                         Float64[200, 1000, 1000, 1100, 200, 500, 500])
        putvarboundsliceconst(task,1, 8, MSK_BK_LO, 0.0, +Inf)

        func(task)
    end
end

"""
Analyzes and prints infeasibility contributing elements
sl - dual values for lower bounds
su - dual values for upper bounds
eps - tolerance for when a nunzero dual value is significant
"""
function analyzeCertificate(sl :: Vector{Float64}, su :: Vector{Float64}, eps :: Float64)
    for i in 1:length(sl)
        if abs(sl[i]) > eps
            println("#$i, lower,  dual = $(sl[i])")
        end
        if abs(su[i]) > eps
            println("#$i, upper,  dual = $(su[i])")
        end
    end
end

# In this example we set up a simple problem
# One could use any task or a task read from a file
testProblem() do task
    # Useful for debugging
    writedata(task,"pinfeas.ptf");                          # Write file in human-readable format
    # Attach a log stream printer to the task
    putstreamfunc(task,MSK_STREAM_LOG,msg -> print(msg))

    # Perform the optimization.
    optimize(task)
    solutionsummary(task,MSK_STREAM_LOG)

    # Check problem status, we use the interior point solution
    if getprosta(task,MSK_SOL_ITR) == MSK_PRO_STA_PRIM_INFEAS
        # Set the tolerance at which we consider a dual value as essential
        eps = 1e-7
        println("Variable bounds important for infeasibility: ");
        analyzeCertificate(getslx(task,MSK_SOL_ITR), getsux(task,MSK_SOL_ITR), eps)

        println("Constraint bounds important for infeasibility: ")
        analyzeCertificate(getslc(task,MSK_SOL_ITR), getsuc(task,MSK_SOL_ITR), eps)
    else
        println("The problem is not primal infeasible, no certificate to show")
        @assert false
    end
end
