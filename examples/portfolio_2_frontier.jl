
# Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
#
# File :      portfolio_2_frontier.jl
#
# Description :  Implements a basic portfolio optimization model.

using Mosek
using Printf

function portfolio( mu :: Vector{Float64},
                    x0 :: Vector{Float64},
                    w  :: Float64,
                    alphas :: Vector{Float64},
                    GT :: Array{Float64,2})
    (k,n) = size(GT)
    maketask() do task
        # Directs the log task stream
        #putstreamfunc(task,MSK_STREAM_LOG,msg -> print(msg))

        totalBudget = sum(x0)+w

        #Offset of variables into the API variable.
        x_ofs = 0
        s_ofs = n

        # Constraints offsets
        budget_ofs = 0

        # Holding variable x of length n
        # No other auxiliary variables are needed in this formulation
        appendvars(task,n+1)

        # Setting up variable x
        for j in 1:n
            # Optionally we can give the variables names
            putvarname(task, x_ofs+j, "x[$(j)]")
            # No short-selling - x^l = 0, x^u = inf
            putvarbound(task,x_ofs+j, MSK_BK_LO, 0.0, Inf)
        end
        putvarname(task, s_ofs+1, "s")
        putvarbound(task,s_ofs+1, MSK_BK_FR, -Inf, Inf)

        # One linear constraint: total budget
        appendcons(task,1)
        putconname(task,1,"budget")
        for j in 1:n
            # Coefficients in the first row of A
            putaij(task,budget_ofs+1, x_ofs+j, 1.0)
        end

        putconbound(task, budget_ofs+1, MSK_BK_FX, totalBudget, totalBudget)

        # Input (gamma, GTx) in the AFE (affine expression) storage
        # We build the following F and g for variables [x, s]:
        #     [0, 1]      [0  ]
        # F = [0, 0], g = [0.5]
        #     [GT,0]      [0  ]
        # We need k+2 rows
        appendafes(task,k + 2)
        # The first affine expression = alpha
        putafefentry(task,1,s_ofs+1,1.0)
        putafeg(task,2,0.5)
        # The remaining k expressions comprise GT*x, we add them row by row
        # In more realisic scenarios it would be better to extract nonzeros and input in sparse form
        subj = [1:n...]
        for i in 1:k
            putafefrow(task,i + 2, subj, GT[i,:])
        end

        # Input the affine conic constraint (alpha, GT*x) \in QCone
        # Add the quadratic domain of dimension k+1
        qdom = appendrquadraticconedomain(task,k + 2)
        # Add the constraint
        appendaccseq(task,qdom,1,zeros(k+2))
        putaccname(task,1, "risk")


        # Objective: maximize expected return mu^T x
        putclist(task,[x_ofs+1:x_ofs+n...],mu)
        putobjsense(task,MSK_OBJECTIVE_SENSE_MAXIMIZE)

        expret = Float64[]
        stddev = Float64[]
        for alpha in alphas
            putcj(task,s_ofs+1,-alpha)
            optimize(task)
            writedata(task,"portfolio_2_frontier-$alpha.ptf")

            # Display solution summary for quick inspection of results
            solutionsummary(task,MSK_STREAM_LOG)

            # Read the results
            r = mu' * getxxslice(task,MSK_SOL_ITR, x_ofs+1,x_ofs+n+1)
            slvl =    getxxslice(task,MSK_SOL_ITR, s_ofs+1,s_ofs+2)

            push!(expret,r)
            push!(stddev,sqrt(slvl[1]))
        end
        (expret,stddev)
    end
end # portfolio()

let w    = 1.0,
    mu   = [0.07197, 0.15518, 0.17535, 0.08981, 0.42896, 0.39292, 0.32171, 0.18379],
    x0   = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    GT   = [0.30758 0.12146 0.11341 0.11327 0.17625 0.11973 0.10435 0.10638
            0.      0.25042 0.09946 0.09164 0.06692 0.08706 0.09173 0.08506
            0.      0.      0.19914 0.05867 0.06453 0.07367 0.06468 0.01914
            0.      0.      0.      0.20876 0.04933 0.03651 0.09381 0.07742
            0.      0.      0.      0.      0.36096 0.12574 0.10157 0.0571
            0.      0.      0.      0.      0.      0.21552 0.05663 0.06187
            0.      0.      0.      0.      0.      0.      0.22514 0.03327
            0.      0.      0.      0.      0.      0.      0.      0.2202 ],
    alphas = [0.0, 0.01, 0.1, 0.25, 0.30, 0.35, 0.4, 0.45, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 10.0],
    (k,n) = size(GT),

    N = length(alphas)
    let (expret,stddev) = portfolio(mu,x0,w,alphas,GT)
        for i in 1:N
            @printf("alpha = %2.2e, exp. ret. = %2.3e, std. dev. = %2.3e\n",alphas[i],expret[i],stddev[i])
        end
    end
end
