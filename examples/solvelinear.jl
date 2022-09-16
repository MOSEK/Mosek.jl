#
#   Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
#
#   File :      solvelinear.jl
##
#   Purpose :   To demonstrate the usage of MSK_solvewithbasis
#               when solving the linear system:
#
#               1.0  x1             = b1
#               -1.0  x0  +  1.0  x1 = b2
#
#               with two different right hand sides
#
#               b = (1.0, -2.0)
#
#               and
#
#               b = (7.0, 0.0)

using Mosek

function setup(task   :: Mosek.Task,
               aval   :: Vector{Float64},
               asub   :: Vector{Int32},
               ptrb   :: Vector{Int64},
               ptre   :: Vector{Int64},
               numvar :: Int32)

    appendvars(task,numvar)
    appendcons(task,numvar)

    putacolslice(task,1,numvar+1,ptrb,ptre,asub,aval)

    putconboundsliceconst(task,1,numvar+1,MSK_BK_FX,0.0,0.0)
    putvarboundsliceconst(task,1,numvar+1,MSK_BK_FR,-Inf,Inf)

    # Define a basic solution by specifying status keys for variables
    # & constraints.
    deletesolution(task,MSK_SOL_BAS)

    putskcslice(task,MSK_SOL_BAS, 1, numvar+1, fill(MSK_SK_FIX,numvar))
    putskxslice(task,MSK_SOL_BAS, 1, numvar+1, fill(MSK_SK_BAS,numvar))

    return initbasissolve(task)
end

let numcon = Int32(2),
    numvar = Int32(2),

    aval = [ -1.0 ,
             1.0, 1.0 ],
    asub = Int32[ 2,
                  1, 2 ],
    ptrb  = Int64[1, 2],
    ptre  = Int64[2, 4]

    # bsub  = new int[numvar];
    # b     = new double[numvar];
    # basis = new int[numvar];

    maketask() do task
        putstreamfunc(task,MSK_STREAM_LOG,msg -> print(msg))

        # Put A matrix and factor A. Call this function only once for a
        # given task.

        basis = setup(task,
                      aval,
                      asub,
                      ptrb,
                      ptre,
                      numvar)

        # now solve rhs
        let b    = Float64[1,-2],
            bsub = Int32[1,2],
            nz = solvewithbasis(task,false, 2, bsub, b)
            println("Solution to Bx = b:")

            # Print solution and show correspondents to original variables in the problem
            for i in 1:nz
                if basis[bsub[i]] <= numcon
                    println("This should never happen")
                else
                    println("x $(basis[bsub[i]] - numcon) = $(b[bsub[i]])")
                end
            end
        end

        let b    = Float64[7,0],
            bsub = Int32[1,0],
            nz = solvewithbasis(task,false,1, bsub, b)

            println("Solution to Bx = b:")
            # Print solution and show correspondents to original variables in the problem
            for i in 1:nz
                if (basis[bsub[i]] <= numcon)
                    println("This should never happen")
                else
                    println("x $(basis[bsub[i]] - numcon) = $(b[bsub[i]])")
                end
            end
        end
    end
end
