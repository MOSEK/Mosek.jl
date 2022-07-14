# Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
#
# File :      reoptimization.jl
#
# Purpose :   Demonstrates how to solve a  linear
#             optimization problem using the MOSEK API
#             and modify and re-optimize the problem.

using Mosek

let numcon = 3,
    numvar = 3,
    c      = [1.5, 2.5, 3.0 ],
    bkc    = [ MSK_BK_UP,
               MSK_BK_UP,
               MSK_BK_UP
               ],
    blc    = [ -Inf,
               -Inf,
               -Inf ],
    buc    = [ 100000.0,
               50000.0,
               60000.0
            ],
    bkx    = [ MSK_BK_LO,
               MSK_BK_LO,
               MSK_BK_LO ],
    blx    = [ 0.0, 0.0, 0.0 ],
    bux    = [ +Inf,
               +Inf,
               +Inf ],
    aptrb  = Int64[ 1,4,7 ],
    aptre  = Int64[ 4,7,10 ],
    asub   = Int32[ 1, 2, 3,
                    1, 2, 3,
                    1, 2, 3 ],
    aval   = [ 2.0, 3.0, 2.0,
               4.0, 2.0, 3.0,
               3.0, 3.0, 2.0 ]

    maketask() do task          # 126
        putstreamfunc(task,MSK_STREAM_LOG,msg -> print(msg))
        # Append the constraints.
        appendcons(task,numcon)

        # Append the variables.
        appendvars(task,numvar)

        # Put C.
        for j in 1:numvar
            putcj(task,j, c[j])
        end

        # Put constraint bounds
        for i in 1:numcon
            putconbound(task,i, bkc[i], blc[i], buc[i])
        end

        # Put variable bounds.
        for j in 1:numvar
            putvarbound(task,j, bkx[j], blx[j], bux[j])
        end

        # Put A.
        putacolslice(task,1,numvar+1,aptrb,aptre,asub,aval)

        # A maximization problem
        putobjsense(task,MSK_OBJECTIVE_SENSE_MAXIMIZE)
        # Solve the problem
        optimize(task)

        xx = getxx(task,MSK_SOL_BAS) # Request the basic solution.

        println("x = $xx")

        #***************** Make a change to the A matrix *****************
        putaij(task,1,1, 3.0)

        optimize(task);
        xx = getxx(task,MSK_SOL_BAS)

        println("x = $xx")

        #**************** Add a new variable *****************************
        # Get index of new variable.
        varidx = getnumvar(task)+1

        # Append a new variable x_3 to the problem
        appendvars(task,1)
        numvar += 1

        # Set bounds on new varaible
        putvarbound(task,
                    varidx,
                    MSK_BK_LO,
                    0.0,
                    Inf);

        # Change objective
        putcj(task,varidx, 1.0)

        # Put new values in the A matrix
        let acolsub = Int32[1,3], 
            acolval = Float64[4.0, 1.0]

            putacol(task,varidx, # column index
                    acolsub,
                    acolval)
        end

        # Change optimizer to simplex free and reoptimize
        putintparam(task,MSK_IPAR_OPTIMIZER, MSK_OPTIMIZER_FREE_SIMPLEX)
        optimize(task)

        xx = getxx(task,MSK_SOL_BAS) # Request the basic solution.
        println("x = $xx")

        #********************* Add a new constraint **************************
        # Get index of new constraint.
        conidx = getnumcon(task)+1

        # Append a new constraint
        appendcons(task,1)
        numcon += 1;

        # Set bounds on new constraint
        putconbound(task,
                    conidx,
                    MSK_BK_UP,
                    -Inf,
                    30000.0)

        # Put new values in the A matrix
        let arowsub = [1,   2,   3,   4  ],
            arowval = [1.0, 2.0, 1.0, 1.0]

            putarow(task,conidx, # row index
                    arowsub,
                    arowval)
        end

        optimize(task)

        xx = getxx(task,MSK_SOL_BAS)

        println("x = $xx")

        #********************* Change constraint bounds *******************
        let newbkc  = [MSK_BK_UP,
                       MSK_BK_UP,
                       MSK_BK_UP,
                       MSK_BK_UP],
            newblc = [ -Inf,
                       -Inf,
                       -Inf,
                       -Inf],
            newbuc = [ 80000.0, 40000.0, 50000.0, 22000.0 ]

            putconboundslice(task,1, numcon+1, newbkc, newblc, newbuc)
        end
        optimize(task)

        xx = getxx(task,MSK_SOL_BAS) # Request the basic solution.

        println("x = $xx")
    end
end
