#
#  Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
#
#  File :      djc1.jl
#
#  Purpose: Demonstrates how to solve the problem with two disjunctions:
#
#      minimize    2x0 + x1 + 3x2 + x3
#      subject to   x0 + x1 + x2 + x3 >= -10
#                  (x0-2x1<=-1 and x2=x3=0) or (x2-3x3<=-2 and x1=x2=0)
#                  x0=2.5 or x1=2.5 or x2=2.5 or x3=2.5
#

using Mosek


maketask() do task
    # Append free variables
    numvar = 4
    appendvars(task,numvar)
    putvarboundsliceconst(task,1, numvar+1, MSK_BK_FR, -Inf, Inf)

    # The linear part: the linear constraint
    appendcons(task,1)
    putarow(task,1,[1:numvar...], ones(numvar))
    putconbound(task,1,MSK_BK_LO, -10.0, -10.0)

    # The linear part: objective
    putobjsense(task,MSK_OBJECTIVE_SENSE_MINIMIZE)
    putclist(task,[1:numvar...], Float64[2, 1, 3, 1])

    # Fill in the affine expression storage F, g
    numafe = 10
    appendafes(task,numafe)

    fafeidx = Int64[1, 1, 2, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    fvaridx = Int32[1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4]
    fval    = Float64[1.0, -2.0, 1.0, -3.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
    g       = Float64[1.0, 2.0, 0.0, 0.0, 0.0, 0.0, -2.5, -2.5, -2.5, -2.5]

    putafefentrylist(task,fafeidx, fvaridx, fval)
    putafegslice(task,1, numafe+1, g)

    # Create domains
    zero1   = appendrzerodomain(task,1)
    zero2   = appendrzerodomain(task,2)
    rminus1 = appendrminusdomain(task,1)

    # Append disjunctive constraints
    numdjc = 2
    appenddjcs(task,numdjc)

    # First disjunctive constraint
    putdjc(task,1,                                        # DJC index
                [rminus1, zero2, rminus1, zero2],         # Domains     (domidxlist)
                [1, 5, 6, 2, 3, 4],                       # AFE indices (afeidxlist)
                zeros(6),                                 # Unused
                [2, 2] )                                  # Term sizes  (termsizelist)

    # Second disjunctive constraint
    putdjc(task,2,                                        # DJC index
                [zero1, zero1, zero1, zero1],             # Domains     (domidxlist)
                [7, 8, 9, 10],                            # AFE indices (afeidxlist)
                zeros(4),                                 # Unused
                [1, 1, 1, 1] )                            # Term sizes  (termidxlist)

    # Useful for debugging
    writedata(task,"djc1.ptf")                         # Write file in human-readable format

    # Solve the problem
    optimize(task)

    # Print a summary containing information
    # about the solution for debugging purposes
    solutionsummary(task,MSK_STREAM_MSG)

    # Get status information about the solution
    sta = getsolsta(task,MSK_SOL_ITG)

    println(sta)
    @assert sta == MSK_SOL_STA_INTEGER_OPTIMAL

    xx = getxx(task,MSK_SOL_ITG)

    println("Optimal solution:")
    for i in 1:numvar
        println("x[$i]=$(xx[i])")
    end
    @assert maximum(abs.(xx - [0.0, 0.0, -12.5, 2.5])) < 1e-7

    @assert maximum(abs.(xx-[0.0, 0.0, -12.5, 2.5])) < 1e-7
end
