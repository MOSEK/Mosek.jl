
# Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
#
# File :      portfolio_1_basic.jl
#
# Description :  Implements a basic portfolio optimization model.

using Mosek

function portfolio( mu :: Vector{Float64},
                    x0 :: Vector{Float64},
                    w  :: Float64,
                    gamma :: Float64,
                    GT :: Array{Float64,2})
    (k,n) = size(GT)
    maketask() do task
        # Directs the log task stream
        putstreamfunc(task,MSK_STREAM_LOG,msg -> print(msg))

        totalBudget = sum(x0)+w

        #Offset of variables into the API variable.
        x_ofs = 0

        # Constraints offsets
        budget_ofs = 0

        # Holding variable x of length n
        # No other auxiliary variables are needed in this formulation
        appendvars(task,n)

        # Setting up variable x
        for j in 1:n
            # Optionally we can give the variables names
            putvarname(task, x_ofs+j, "x[$(j)]");
            # No short-selling - x^l = 0, x^u = inf
            putvarbound(task,x_ofs+j, MSK_BK_LO, 0.0, Inf);
        end

        # One linear constraint: total budget
        appendcons(task,1);
        putconname(task,1,"budget");
        for j in 1:n
            # Coefficients in the first row of A
            putaij(task,budget_ofs+1, x_ofs+j, 1.0)
        end

        putconbound(task, budget_ofs+1, MSK_BK_FX, totalBudget, totalBudget)

        # Input (gamma, GTx) in the AFE (affine expression) storage
        # We need k+1 rows
        appendafes(task,k + 1)
        # The first affine expression = gamma
        putafeg(task,1, gamma)
        # The remaining k expressions comprise GT*x, we add them row by row
        # In more realisic scenarios it would be better to extract nonzeros and input in sparse form


        subj = [1:n...]
        for i in 1:k
            putafefrow(task,i + 1, subj, GT[i,:])
        end

        # Input the affine conic constraint (gamma, GT*x) \in QCone
        # Add the quadratic domain of dimension k+1
        qdom = appendquadraticconedomain(task,k + 1)
        # Add the constraint
        appendaccseq(task,qdom,1,zeros(k+1))
        putaccname(task,1, "risk")



        # Objective: maximize expected return mu^T x
        putclist(task,[x_ofs+1:x_ofs+n...],mu)
        putobjsense(task,MSK_OBJECTIVE_SENSE_MAXIMIZE)

        optimize(task)

        # Display solution summary for quick inspection of results
        solutionsummary(task,MSK_STREAM_LOG)

        writedata(task,"portfolio_1_basic.ptf");

        # Read the results
        xx = getxxslice(task,MSK_SOL_ITR, x_ofs+1,x_ofs+n+1)
        expret = mu' * xx

        (xx,expret)
    end
end # portfolio()

let w      = 59.0,
    mu     = [0.07197349, 0.15518171, 0.17535435, 0.0898094 , 0.42895777, 0.39291844, 0.32170722, 0.18378628],
    x0     = [8.0, 5.0, 3.0, 5.0, 2.0, 9.0, 3.0, 6.0],
    gamma  = 36.0,
    GT     = [ 0.30758 0.12146 0.11341 0.11327 0.17625 0.11973 0.10435 0.10638
               0.      0.25042 0.09946 0.09164 0.06692 0.08706 0.09173 0.08506
               0.      0.      0.19914 0.05867 0.06453 0.07367 0.06468 0.01914
               0.      0.      0.      0.20876 0.04933 0.03651 0.09381 0.07742
               0.      0.      0.      0.      0.36096 0.12574 0.10157 0.0571
               0.      0.      0.      0.      0.      0.21552 0.05663 0.06187
               0.      0.      0.      0.      0.      0.      0.22514 0.03327
               0.      0.      0.      0.      0.      0.      0.      0.2202 ],
    (k,n) = size(GT)

    let (xx,expret) = portfolio(mu,x0,w,gamma,GT)
        println("Expected return $(expret) for gamma $(gamma)")
        println("Solution vector = $(xx)")

            @assert abs(expret - 4.1922467685e+01) < 1e-7
            @assert sum(abs.(xx - [7.7127e-08, 1.0252e-07, 9.9557e-08, 7.1640e-08, 7.2993e+01, 2.7007e+01, 1.8063e-07, 9.9400e-08])) < 1e-3
    end
end
