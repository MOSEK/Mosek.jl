#
# Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
#
# File :      mico1.jl
#
# Purpose :   Demonstrates how to solve a small mixed
#             integer conic optimization problem.

#             minimize    x^2 + y^2
#             subject to  x >= e^y + 3.8
#                         x, y - integer

using Mosek

maketask() do task
    # Directs the log task stream to the user specified
    # method task_msg_obj.stream
    appendvars(task,3);   # x, y, t
    x=1
    y=2
    t=3
    putvarboundsliceconst(task,1, 4, MSK_BK_FR, -Inf, Inf)

    # Integrality constraints for x, y
    putvartypelist(task,
                   [x,y],
                   [MSK_VAR_TYPE_INT,MSK_VAR_TYPE_INT])

    # Set up the affine expressions
    # x, x-3.8, y, t, 1.0
    appendafes(task,5)
    putafefentrylist(task,
                     [1,2,3,4],
                     [x,x,y,t],
                     [1,1,1,1])
    putafegslice(task,1, 6, Float64[0, -3.8, 0, 0, 1.0])

    # Add constraint (x-3.8, 1, y) \in \EXP
    appendacc(task,appendprimalexpconedomain(task), Int64[2, 5, 3], zeros(3))

    # Add constraint (t, x, y) \in \QUAD
    appendacc(task,appendquadraticconedomain(task,3), Int64[4, 1, 3], zeros(3))

    # Objective
    putobjsense(task,MSK_OBJECTIVE_SENSE_MINIMIZE)
    putcj(task,t, 1.0)

    # Optimize the task
    optimize(task)
    solutionsummary(task,MSK_STREAM_MSG)

    xx = getxxslice(task,MSK_SOL_ITG, 1, 3)
    println("x = $(xx[1]) y = $(xx[2])")

    @assert maximum(abs.(xx-[4.0, -2.0])) < 1e-7
    
end
