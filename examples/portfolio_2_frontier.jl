
# File : portfolio_1_basic
#
# Copyright : Mosek ApS
#
# Description :  Implements a basic portfolio optimization model.

include("portfolio_data.jl")

using Mosek
using Printf
#TAG:begin-code
#TAG:begin-basic-markowitz

function portfolio( mu :: Vector{Float64},
                    x0 :: Vector{Float64},
                    w  :: Float64,
                    gamma :: Vector{Float64},
                    GT :: Array{Float64,2})
    (k,n) = size(GT)
    maketask() do task
        # Directs the log task stream
        #putstreamfunc(task,MSK_STREAM_LOG,msg -> print(msg))

        totalBudget = sum(x0)+w

        #TAG:begin-offsets
        #Offset of variables into the API variable.
        x_ofs = 0
        s_ofs = n

        # Constraints offsets
        #TAG:end-offsets
        budget_ofs = 0

        # Holding variable x of length n
        # No other auxiliary variables are needed in this formulation
        appendvars(task,n+1)

        # Setting up variable x
        for j in 1:n
            # Optionally we can give the variables names
            putvarname(task, x_ofs+j, "x[$(j)]");
            # No short-selling - x^l = 0, x^u = inf
            putvarbound(task,x_ofs+j, MSK_BK_LO, 0.0, Inf);
        end
        putvarname(task, s_ofs+1, "s");
        putvarbound(task,s_ofs+1, MSK_BK_FR, -Inf, Inf);

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
        appendafes(task,k + 2)
        # The first affine expression = gamma
        putafefentry(task,1,s_ofs+1,1.0)
        putafeg(task,2,0.5)
        # The remaining k expressions comprise GT*x, we add them row by row
        # In more realisic scenarios it would be better to extract nonzeros and input in sparse form
        subj = [1:n...]
        for i in 1:k
            putafefrow(task,i + 2, subj, GT[i,:])
        end

        #TAG:begin-basic-markowitz-appendaccseq
        # Input the affine conic constraint (gamma, GT*x) \in QCone
        # Add the quadratic domain of dimension k+1
        qdom = appendrquadraticconedomain(task,k + 2)
        # Add the constraint
        appendaccseq(task,qdom,1,zeros(k+2))
        #TAG:end-basic-markowitz-appendaccseq
        putaccname(task,1, "risk")


        # Objective: maximize expected return mu^T x
        putclist(task,[x_ofs+1:x_ofs+n...],mu)
        putobjsense(task,MSK_OBJECTIVE_SENSE_MAXIMIZE)

        expret = Float64[]
        stddev = Float64[]
        for gamma in gammas
            putcj(task,s_ofs+1,-gamma)
            optimize(task)
            writedata(task,"pf-frontier-$gamma.ptf")

            #TAG:begin-solutionsummary
            # Display solution summary for quick inspection of results
            solutionsummary(task,MSK_STREAM_LOG)
            #TAG:end-solutionsummary

            # Read the results
            r = mu' * getxxslice(task,MSK_SOL_ITR, x_ofs+1,x_ofs+n+1)
            slvl =    getxxslice(task,MSK_SOL_ITR, s_ofs+1,s_ofs+2)

            push!(expret,r)
            push!(stddev,sqrt(slvl[1]))
        end
        (expret,stddev)
    end
end # portfolio()
#TAG:end-code
#TAG:end-basic-markowitz

gammas = [0.0, 0.01, 0.1, 0.25, 0.30, 0.35, 0.4, 0.45, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 10.0]
N = length(gammas)
let (expret,stddev) = portfolio(mu,x0,w,gammas,GT)
    for i in 1:N
        @printf("gamma = %2.2e, exp. ret. = %2.3e, std. dev. = %2.3e\n",gammas[i],expret[i],stddev[i])
    end
end
