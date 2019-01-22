#
# Solves the exponential-cone problem
#
#    minimize    x1+x2
#    subject to  x1 + x2 + x3 = 1
#                (x1,x2,x3) in Kexp
#
#    where Kexp = { (x1,x2,x3) | x1 >= x2*exp(x3/x2) }.
#

using Mosek
printstream(msg::AbstractString) = print(msg)

bkc = [ MSK_BK_FX ]
blc = [ 1.0 ]
buc = [ 1.0 ]

c   = [ 1.0, 1.0, 0.0 ]
bkx = [ MSK_BK_FR,MSK_BK_FR,MSK_BK_FR ]
blx = [ 0.0, 0.0, 0.0 ]
bux = [ 0.0, 0.0, 0.0 ]

asub  = [ 1, 2, 3 ]
aval  = [ 1.0, 1.0, 1.0 ]

numvar = length(bkx)
numcon = length(bkc)


# Create a task
maketask() do task
    putstreamfunc(task,MSK_STREAM_LOG,msg -> print(msg))

    appendcons(task,numcon)
    appendvars(task,numvar)

    # Set the linear term c_j in the objective.
    putclist(task,[1:3;],c)

    # Set the bounds on variable j
    putvarboundslice(task,1,numvar+1,bkx,blx,bux)

    putarow(task,1,asub,aval)
    putconbound(task,1,bkc[1],blc[1],buc[1])

    # Input the cones
    appendcone(task,MSK_CT_PEXP, 0.0, [1, 2, 3 ])

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
        x1,x2,x3 = getxx(task,MSK_SOL_ITR)
        println("Optimal solution: (x1,x2,x3)=($([x1,x2,x3]))")
        println("Solution is on boundary: x1-x2*exp(x3/x2)=$(x1-x2*exp(x3/x2))")
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
