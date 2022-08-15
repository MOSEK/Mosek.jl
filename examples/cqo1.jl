#
#  Copyright : Copyright (c) 2022 MOSEK ApS
#
#  File :      cqo1.jl
#
#  Purpose: Demonstrates how to solve small conic
#           optimization problem using the MOSEK Python API.
##

##TAG:begin-code
using Mosek
using Printf

printstream(msg::AbstractString) = print(msg)
callback(where,dinf,iinf,liinf) = 0 

# Since the actual value of Infinity is ignores, we define it solely
# for symbolic purposes:

bkc = [ MSK_BK_FX ]
blc = [ 1.0 ]
buc = [ 1.0 ]

c   = [               0.0,              0.0,              0.0,
                      1.0,              1.0,              1.0 ]
bkx = [ MSK_BK_LO,MSK_BK_LO,MSK_BK_LO,
        MSK_BK_FR,MSK_BK_FR,MSK_BK_FR ]
blx = [               0.0,              0.0,              0.0,
                     -Inf,             -Inf,             -Inf ]
bux = [               Inf,              Inf,              Inf,
                      Inf,              Inf,              Inf ]

asub  = [ 1 ,2, 3 ]
aval  = [ 1.0, 1.0, 2.0 ]

numvar = length(bkx)
numcon = length(bkc)


# Create a task
maketask() do task
    putstreamfunc(task,MSK_STREAM_LOG,printstream)
    putcallbackfunc(task,callback)

##TAG:begin-append
    # Append 'numcon' empty constraints.
    # The constraints will initially have no bounds.
    appendcons(task,numcon)

    #Append 'numvar' variables.
    # The variables will initially be fixed at zero (x=0).
    appendvars(task,numvar)
##TAG:end-append

##TAG:end-putcj
##TAG:begin-putbound-var
    # Set the linear term c_j in the objective.
    putclist(task,[1:6;],c)
##TAG:end-putcj

##TAG:begin-putbound-var
    # Set the bounds on variable j
    # blx[j] <= x_j <= bux[j]
    putvarboundslice(task,1,numvar+1,bkx,blx,bux)
##TAG:end-putbound-var

##TAG:begin-putavec
    putarow(task,1,asub,aval)
    putconbound(task,1,bkc[1],blc[1],buc[1])
##TAG:end-putavec

##TAG:begin-appendcone
    # Input the cones
    appendafes(task,6)
    putafefentrylist(task,
                     [1, 2, 3, 4, 5, 6],         # Rows
                     [4, 1, 2, 5, 6, 3],         # Columns */
                     [1.0, 1.0, 1.0, 1.0, 1.0, 1.0])

    # Quadratic cone (x(3),x(0),x(1)) \in QUAD_3
    quadcone  = appendquadraticconedomain(task,3)
    appendacc(task,
              quadcone,  # Domain
              [1, 2, 3], # Rows from F
              [0.0,0.0,0.0])

    # Rotated quadratic cone (x(4),x(5),x(2)) \in RQUAD_3
    rquadcone = appendrquadraticconedomain(task,3)
    appendacc(task,
              rquadcone, # Domain
              [4, 5, 6], # Rows from F
              [0.0,0.0,0.0]);
##TAG:end-appendcone

    # Input the objective sense (minimize/maximize)
    putobjsense(task,MSK_OBJECTIVE_SENSE_MINIMIZE)

##TAG:begin-optimize
    # Optimize the task
    #optimize(task,"mosek://solve.mosek.com:30080")
    optimize(task)
##TAG:end-optimize
    writedata(task,"cqo1.ptf")
    # Print a summary containing information
    # about the solution for debugging purposes
    solutionsummary(task,MSK_STREAM_MSG)
##TAG:begin-getsolutionstatus
    prosta = getprosta(task,MSK_SOL_ITR)
    solsta = getsolsta(task,MSK_SOL_ITR)
##TAG:end-getsolutionstatus

    if solsta == MSK_SOL_STA_OPTIMAL
        # Output a solution
        xx = getxx(task,MSK_SOL_ITR)
        @printf("Optimal solution: %s\n", xx')
    elseif solsta == MSK_SOL_STA_DUAL_INFEAS_CER
        println("Primal or dual infeasibility.\n")
    elseif solsta == MSK_SOL_STA_PRIM_INFEAS_CER
        println("Primal or dual infeasibility.\n")
    elseif solsta == MSK_SOL_STA_UNKNOWN
        println("Unknown solution status")
    else
        println("Other solution status")
    end

    ##TAG:ASSERT:begin-check-solution
    @assert maximum(abs.(xx - [0.2609204081408032, 0.2609204081408032, 0.23907959185918956, 0.36899717989264824, 0.1690548006469457, 0.1690548006469457])) < 1e-7
    ##TAG:ASSERT:end-check-solution
end
