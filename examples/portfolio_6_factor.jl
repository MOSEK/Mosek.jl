#
#  Copyright : Copyright (c) 2022 MOSEK ApS
#
#  File :      portfolio_6_factor.jl
#
#  Description :  Implements a portfolio optimization model using factor model

using Mosek
using LinearAlgebra
#TAG:begin-code
##TAG:begin-factor-markowitz
function portfolio( mu :: Vector{Float64},
                    x0 :: Vector{Float64},
                    w  :: Float64,
                    gammas :: Vector{Float64},
                    GT :: Array{Float64,2})
    (k,n) = size(GT)
    maketask() do task
        # Directs the log task stream
        putstreamfunc(task,MSK_STREAM_LOG,msg -> print(msg))

        totalBudget = sum(x0)+w

        #TAG:begin-offsets
        #Offset of variables into the API variable.
        x_ofs = 0

        # Constraints offsets
        #TAG:end-offsets
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
        # The remaining k expressions comprise GT*x, we add them row by row
        # In more realisic scenarios it would be better to extract nonzeros and input in sparse form


        subj = [1:n...]
        for i in 1:k
            putafefrow(task,i + 1, subj, GT[i,:])
        end

        #TAG:begin-basic-markowitz-appendaccseq
        # Input the affine conic constraint (gamma, GT*x) \in QCone
        # Add the quadratic domain of dimension k+1
        qdom = appendquadraticconedomain(task,k + 1)
        # Add the constraint
        appendaccseq(task,qdom,1,zeros(k+1))
        #TAG:end-basic-markowitz-appendaccseq
        putaccname(task,1, "risk")

        # Objective: maximize expected return mu^T x
        putclist(task,[x_ofs+1:x_ofs+n...],mu)
        putobjsense(task,MSK_OBJECTIVE_SENSE_MAXIMIZE)


        res = Tuple{Float64,Float64}[]
        for gamma in gammas
            putafeg(task,1, gamma)

            optimize(task)

            #TAG:begin-writedata
            writedata(task,"portfolio_6_factor-$(gamma).ptf");
            #TAG:end-writedatahol

            # Read the results
            xx = getxxslice(task,MSK_SOL_ITR, x_ofs+1,x_ofs+n+1)
            expret = mu' * xx

            push!(res,(gamma,expret))
        end
        
        res
    end
end # portfolio()
#TAG:end-code
##TAG:begin-factor-markowitz

let w  = 1.0,
    mu = [0.07197, 0.15518, 0.17535, 0.08981, 0.42896, 0.39292, 0.32171, 0.18379],
    x0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    gammas = [0.24, 0.28, 0.32, 0.36, 0.4, 0.44, 0.48],
    #TAG:begin-factor-model-vars
    # Factor exposure matrix, n x 2
    B = [ 0.4256  0.1869
          0.2413  0.3877
          0.2235  0.3697
          0.1503  0.4612
          1.5325 -0.2633
          1.2741 -0.2613
          0.6939  0.2372
          0.5425  0.2116 ],
    # Factor covariance matrix, 2x2
    S_F = [ 0.0620 0.0577
            0.0577 0.0908 ],
    # Specific risk components
    theta = [0.0720, 0.0508, 0.0377, 0.0394, 0.0663, 0.0224, 0.0417, 0.0459],
    S_sqrt_theta = Matrix(Diagonal(sqrt.(theta))),
    #TAG:end-factor-model-vars
    #TAG:begin-factor-model-G
    
    P = Matrix(cholesky(S_F).L),
    BP = B * P,
    G = [ BP S_sqrt_theta ],
    GT = Matrix(G') # 10 x 8

    
    
    for (gamma,expret) in portfolio(mu,x0,w,gammas,GT)
        println("Expected return $expret for gamma $gamma");
    end

    print(P)
end
