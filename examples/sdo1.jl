##
#  Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
#
#  File :      sdo1.jl
#
#   Purpose:   Demonstrates how to solve a small mixed semidefinite and conic quadratic
#              optimization problem using the MOSEK Julia API.
##

using Mosek
using Printf, SparseArrays

printstream(msg::String) = print(msg)

# Bound keys for constraints
bkc = [ MSK_BK_FX
        MSK_BK_FX]

# Bound values for constraints
blc = [1.0, 0.5]
buc = [1.0, 0.5]


A = sparse( [1,2,2],[1,2,3],[1.0, 1.0, 1.0])
conesub = [1, 2, 3]

barci = [1, 2, 2, 3, 3]
barcj = [1, 1, 2, 2, 3]
barcval = [2.0, 1.0, 2.0, 1.0, 2.0]

barai   = Any[ [1, 2, 3],
               [1, 2, 3, 2, 3, 3] ]
baraj   = Any[ [1, 2, 3],
               [1, 1, 1, 2, 2, 3] ]
baraval = Any[ [1.0, 1.0, 1.0],
               [1.0, 1.0, 1.0, 1.0, 1.0, 1.0] ]

numvar = 3
numcon = length(bkc)
barvardim = [3]

# Create a task object and attach log stream printer
maketask() do task
    putstreamfunc(task,MSK_STREAM_LOG,printstream)

    # Append 'numvar' variables.
    # The variables will initially be fixed at zero (x=0).
    appendvars(task,numvar)

    # Append 'numcon' empty constraints.
    # The constraints will initially have no bounds.
    appendcons(task,numcon)

    # Append matrix variables of sizes in 'BARVARDIM'.
    # The variables will initially be fixed at zero.
    appendbarvars(task,barvardim)

    # Set the linear term c_0 in the objective.
    putcj(task, 1, 1.0)

    # Set the bounds on variable j
    # blx[j] <= x_j <= bux[j]
    putvarboundsliceconst(task,1,numvar+1,
                          MSK_BK_FR,
                          -Inf,
                          +Inf)

    # Set the bounds on constraints.
    # blc[i] <= constraint_i <= buc[i]
    putconboundslice(task,1,numcon+1, bkc,blc,buc)


    
    # Append the conic quadratic cone
    let afei = getnumafe(task)+1,
        dom = appendquadraticconedomain(task,3)
        
        appendafes(task,3)
        putafefentrylist(task,[1,2,3],
                              [1,2,3],
                              [1.0,1.0,1.0])
        appendaccseq(task,dom,afei,[0.0,0.0,0.0])
    end

    
    # Input row i of A
    putacolslice(task,1,numvar+1,
                 A.colptr[1:numvar], A.colptr[2:numvar+1],
                 A.rowval,A.nzval)


    symc  = appendsparsesymmat(task,barvardim[1],
                               barci,
                               barcj,
                               barcval)

    syma0 = appendsparsesymmat(task,barvardim[1],
                               barai[1],
                               baraj[1],
                               baraval[1])

    syma1 = appendsparsesymmat(task,barvardim[1],
                               barai[2],
                               baraj[2],
                               baraval[2])

    putbarcj(task,1, [symc], [1.0])

    putbaraij(task,1, 1, [syma0], [1.0])
    putbaraij(task,2, 1, [syma1], [1.0])

    # Input the objective sense (minimize/maximize)
    putobjsense(task,MSK_OBJECTIVE_SENSE_MINIMIZE)

    # Solve the problem and print summary
    optimize(task,"mosek://solve.mosek.com:30080")
    writedata(task,"sdo1.ptf")
    solutionsummary(task,MSK_STREAM_MSG)

    # Get status information about the solution
    prosta = getprosta(task,MSK_SOL_ITR)
    solsta = getsolsta(task,MSK_SOL_ITR)


    if solsta == MSK_SOL_STA_OPTIMAL
        # Output a solution
        xx = getxx(task,MSK_SOL_ITR)
        barx = getbarxj(task,MSK_SOL_ITR, 1)
        @assert maximum(abs.(xx-[0.254404851118613, 0.17989139511438135, 0.17989139511438135])) < 1e-6
        @assert maximum(abs.(barx-[0.21725335998058032, -0.25997116466128245, 0.21725335966506842, 0.3110884301498656, -0.2599711646612823, 0.21725335998058004])) < 1e-6

        @printf("Optimal solution: \n  xx = %s\n  barx = %s\n", xx',barx')
    elseif solsta == MSK_SOL_STA_DUAL_INFEAS_CER
        println("Primal or dual infeasibility.\n")
    elseif solsta == MSK_SOL_STA_PRIM_INFEAS_CER
        println("Primal or dual infeasibility.\n")
    elseif  solsta == MSK_SOL_STA_UNKNOWN
        println("Unknown solution status")
    else
        println("Other solution status")
    end
end
