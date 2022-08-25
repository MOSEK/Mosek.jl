##
#  Copyright: Copyright (c) MOSEK ApS, Denmark. All rights reserved.
#
#  File:    helloworld.jl
#
#  The most basic example of how to get started with MOSEK.
##


using Mosek

maketask() do task
    appendvars(task,1)                             # 1 variable x
    putcj(task,1, 1.0)                             # c_0 = 1.0
    putvarbound(task,1,MSK_BK_RA, 2.0, 3.0)        # 2.0 <= x <= 3.0
    putobjsense(task,MSK_OBJECTIVE_SENSE_MINIMIZE) # minimize

    optimize(task)                                 # Optimize

    x = getxx(task,MSK_SOL_ITR)                    # Get solution
    println("Solution x = $(x[1])")                # Print solution
end
