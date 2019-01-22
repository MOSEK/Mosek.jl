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
numvar = 6
numcon = 1

# Create a task
maketask() do task
    putstreamfunc(task,MSK_STREAM_LOG,msg -> print(msg))

    appendcons(task,numcon)
    appendvars(task,numvar)

    # Set up the linear part of the problem
    putclist(task,csub, cval)
    putarow(task,1, asub, aval)
    putvarboundslice(task,1, 6,
                     [ MSK_BK_FR, MSK_BK_FR, MSK_BK_FR, MSK_BK_FR, MSK_BK_FR, MSK_BK_FX ],
                     [ -Inf,      -Inf,      -Inf,      -Inf,      -Inf,      1.0],
                     [  Inf,       Inf,       Inf,       Inf,       Inf,      1.0])
    putconbound(task,1, MSK_BK_FX, 2.0, 2.0)

    # Input the cones
    appendcone(task,MSK_CT_PPOW, 0.2, [1, 2, 4])
    appendcone(task,MSK_CT_PPOW, 0.4, [3, 6, 5])

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
        x = getxx(task,MSK_SOL_ITR)
        println("Optimal solution: $x")
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
