##
#  Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
#
#  File :      response.jl
#
#  Purpose :   This example demonstrates proper response handling
#              for problems solved with the interior-point optimizers.
#
##

using Mosek


cqo1_ptf = "
Task ''
Objective obj
    Minimize + x4 + x5 + x6
Constraints
    c1 [1] + x1 + x2 + 2 x3
    k1 [QUAD(3)]
        @ac1: + x4
        @ac2: + x1
        @ac3: + x2
    k2 [RQUAD(3)]
        @ac4: + x5
        @ac5: + x6
        @ac6: + x3
Variables
    x4
    x1 [0;+inf]
    x2 [0;+inf]
    x5
    x6
    x3 [0;+inf]
"

maketask() do task
    putstreamfunc(task,MSK_STREAM_LOG,msg -> print(msg))

    if length(ARGS) < 1
        readptfstring(task,cqo1_ptf)
    else
        readdata(task,ARGS[1])
    end

    # (Optional) uncomment to see what happens when solution status is unknown
    # putintparam(task,MSK_IPAR_INTPNT_MAX_ITERATIONS, 1)


    # Optimize
    trmcode = optimize(task)
    solutionsummary(task,MSK_STREAM_LOG)

    # We expect solution status OPTIMAL
    solsta = getsolsta(task,MSK_SOL_ITR)

    if solsta == MSK_SOL_STA_OPTIMAL
        # Optimal solution. Fetch and print it.
        println("An optimal interior-point solution is located.")
        numvar = getnumvar(task)
        xx = getxx(task,MSK_SOL_ITR)
        println("x = $xx")
    elseif solsta == MSK_SOL_STA_DUAL_INFEAS_CER
        println("Dual infeasibility certificate found.")
    elseif solsta == MSK_SOL_STA_PRIM_INFEAS_CER
        println("Primal infeasibility certificate found.")
    elseif solsta == MSK_SOL_STA_UNKNOWN
        # The solutions status is unknown. The termination code
        # indicates why the optimizer terminated prematurely.
        println("The solution status is unknown.")
        (symname, desc) = getcodedesc(trmcode)
        println("   Termination code: $symname $desc")
    else
        println("An unexpected solution status $solsta is obtained.")
    end
end
