#
#  Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
#
#  File :      mioinitsol.jl
#
# Purpose :   Demonstrates how to solve a MIP with a start guess.


using Mosek

let numvar = 4,
    numcon = 1,
    c      = [ 7.0, 10.0, 1.0, 5.0 ],
    bkc    = [MSK_BK_UP],
    blc    = [ -Inf],
    buc    = [2.5],
    bkx    = [ MSK_BK_LO,
               MSK_BK_LO,
               MSK_BK_LO,
               MSK_BK_LO],
    blx    = [ 0.0,
               0.0,
               0.0,
               0.0 ],
    bux    = [ Inf,
               Inf,
               Inf,
               Inf ],
    ptrb   = Int64[1, 2, 3, 4],
    ptre   = Int64[2, 3, 4, 5],
    aval   = Float64[1.0, 1.0, 1.0, 1.0],
    asub   = Int32[1,   1,   1,   1  ],
    intsub = Int32[1, 2, 3],
    inttype = [ MSK_VAR_TYPE_INT,
                MSK_VAR_TYPE_INT,
                MSK_VAR_TYPE_INT ]

    maketask() do task
        putstreamfunc(task,MSK_STREAM_LOG,msg -> print(msg))

        inputdata(task,numcon, numvar,
                  c, 0.0,
                  ptrb, ptre,
                  asub, aval,
                  bkc, blc, buc,
                  bkx, blx, bux)

        putvartypelist(task,intsub, inttype)

        # A maximization problem
        putobjsense(task,MSK_OBJECTIVE_SENSE_MAXIMIZE)

        # Assign values to integer variables
        # We only set that slice of xx
        putxxslice(task,MSK_SOL_ITG, 1, 4, [1.0, 1.0, 0.0])

        # Request constructing the solution from integer variable values
        putintparam(task,MSK_IPAR_MIO_CONSTRUCT_SOL, MSK_ON)

        # solve
        optimize(task)
        writedata(task,"mioinitsol.ptf")
        solutionsummary(task,MSK_STREAM_LOG)

        # Read and print solution
        if solutiondef(task,MSK_SOL_ITG)
            # Output a solution
            xx = getxx(task,MSK_SOL_ITG)
            println("Integer optimal solution:")
            for i in 1:numvar
                println("  x[$i] = $(xx[i])")
            end

            # Was the initial guess used?
            constr = getintinf(task,MSK_IINF_MIO_CONSTRUCT_SOLUTION)
            constrVal = getdouinf(task,MSK_DINF_MIO_CONSTRUCT_SOLUTION_OBJ)
            println("Construct solution utilization: $constr")
            println("Construct solution objective: $constrVal")

            @assert maximum(abs.(xx-[0.0, 2.0, 0.0, 0.5])) < 1e-7
            @assert abs(constrVal-19.5) < 1e-7
            @assert constr == 1
        else
            println("No integer solution is available.")
        end
        # Was the initial solution used?
        constr = getintinf(task,MSK_IINF_MIO_CONSTRUCT_SOLUTION)
        constrVal = getdouinf(task,MSK_DINF_MIO_CONSTRUCT_SOLUTION_OBJ)
        println("Construct solution utilization: $constr")
        println("Construct solution objective: $constrVal")
    end
end
