##
# Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
#
# File :      portfolio_3_impact.jl
#
# Description :  Implements a basic portfolio optimization model.
##

using Mosek
function portfolio( mu :: Vector{Float64},
                    x0 :: Vector{Float64},
                    w  :: Float64,
                    gamma :: Float64,
                    GT :: Array{Float64,2},
                    m  :: Vector{Float64})
    (k,n) = size(GT)
    maketask() do task
        # Directs the log task stream
        putstreamfunc(task,MSK_STREAM_LOG,msg -> print(msg))

        totalBudget = sum(x0)+w

        #Offset of variables into the API variable.
        x_ofs = 0
        c_ofs = n
        z_ofs = 2*n
        # Constraints offsets
        budget_ofs = 0

        # Holding variable x of length n
        # No other auxiliary variables are needed in this formulation
        appendvars(task,3*n)

        # Setting up variable x
        for j in 1:n
            # Optionally we can give the variables names
            putvarname(task, x_ofs+j, "x[$(j)]");
            putvarname(task, c_ofs+j, "c[$(j)]");
            putvarname(task, z_ofs+j, "z[$(j)]");
        end
        # No short-selling - x^l = 0, x^u = inf
        putvarboundsliceconst(task,x_ofs+1,x_ofs+n+1, MSK_BK_LO,  0.0, Inf);
        putvarboundsliceconst(task,c_ofs+1,c_ofs+n+1, MSK_BK_FR, -Inf, Inf);
        putvarboundsliceconst(task,z_ofs+1,z_ofs+n+1, MSK_BK_FR, -Inf, Inf);

        # Linear constraint: total budget
        let budget_ofs = getnumcon(task)
            appendcons(task,1);
            putconname(task,1,"budget");
            for j in 1:n
                # Coefficients in the first row of A
                putaij(task,budget_ofs+1, x_ofs+j, 1.0)
                putaij(task,budget_ofs+1, c_ofs+j, m[j])
            end
        end

        putconbound(task, budget_ofs+1, MSK_BK_FX, totalBudget, totalBudget)

        # - Absolute value
        let zabs1_ofs = getnumcon(task),
            zabs2_ofs = zabs1_ofs+n

            appendcons(task,2*n)
            for i in 1:n
                putconname(task,zabs1_ofs + i, "zabs1[$i]")
                putaij(task,zabs1_ofs + i, x_ofs + i, -1.0)
                putaij(task,zabs1_ofs + i, z_ofs + i, 1.0)
                putconbound(task,zabs1_ofs + i, MSK_BK_LO, -x0[i], Inf)
                putconname(task,zabs2_ofs + i, "zabs2[$i]")
                putaij(task,zabs2_ofs + i, x_ofs + i, 1.0)
                putaij(task,zabs2_ofs + i, z_ofs + i, 1.0)
                putconbound(task,zabs2_ofs + i, MSK_BK_LO, x0[i], Inf)
            end
        end

        let qdom = appendquadraticconedomain(task,k + 1)
            afei = getnumafe(task)
            acci = getnumacc(task)
            # Input (gamma, GTx) in the AFE (affine expression) storage
            # We need k+1 rows
            appendafes(task,k + 1)
            # The first affine expression = gamma
            putafeg(task,afei+1, gamma)
            # The remaining k expressions comprise GT*x, we add them row by row
            # In more realisic scenarios it would be better to extract nonzeros and input in sparse form

            subj = [1:n...]
            for i in 1:k
                putafefrow(task,afei+i+1, subj, GT[i,:])
            end

            # Input the affine conic constraint (gamma, GT*x) \in QCone
            # Add the quadratic domain of dimension k+1

            # Add the constraint
            appendaccseq(task,qdom,1,zeros(k+1))
            putaccname(task,acci+1, "risk")
        end

        let dom = appendprimalpowerconedomain(task,3,[2.0,1.0])
            afei = getnumafe(task)
            acci = getnumacc(task)
            afe0 = afei

            appendafes(task,2*n+1)
            putafeg(task,afe0+1,1.0)

            afei += 1

            for i in 1:n
                putafefentry(task,afei + i,     c_ofs + i, 1.0);
                putafefentry(task,afei + n + i, z_ofs + i, 1.0);
            end

            accafes = Int64[ k for i in 1:n for k in [ afei + i, afe0+1, afei + n + i ] ]
            accb    = zeros(n*3)
            accdom  = Int64[ dom for i in 1:n ]

            appendaccs(task,accdom,accafes,accb)

            for i in 1:n
                putaccname(task,acci+i,"market_impact[$i]")
            end
        end

        # Objective: maximize expected return mu^T x
        putclist(task,[x_ofs+1:x_ofs+n...],mu)
        putobjsense(task,MSK_OBJECTIVE_SENSE_MAXIMIZE)

        optimize(task)

        # Display solution summary for quick inspection of results
        solutionsummary(task,MSK_STREAM_LOG)

        writedata(task,"portfolio_3_impact.ptf");

        # Read the results
        xx = getxxslice(task,MSK_SOL_ITR, x_ofs+1,x_ofs+n+1)
        expret = mu' * xx

        (xx,expret)
    end
end # portfolio()

let w    = 1.0,
    mu = [0.07197, 0.15518, 0.17535, 0.08981, 0.42896, 0.39292, 0.32171, 0.18379],
    x0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    GT = [0.30758 0.12146 0.11341 0.11327 0.17625 0.11973 0.10435 0.10638
          0.      0.25042 0.09946 0.09164 0.06692 0.08706 0.09173 0.08506
          0.      0.      0.19914 0.05867 0.06453 0.07367 0.06468 0.01914
          0.      0.      0.      0.20876 0.04933 0.03651 0.09381 0.07742
          0.      0.      0.      0.      0.36096 0.12574 0.10157 0.05710
          0.      0.      0.      0.      0.      0.21552 0.05663 0.06187
          0.      0.      0.      0.      0.      0.      0.22514 0.03327
          0.      0.      0.      0.      0.      0.      0.      0.22020],
    gamma = 0.36,
    (k,n) = size(GT),
    market_impact = 0.01 * ones(n)

    let (xx,expret) = portfolio(mu,x0,w,gamma,GT,market_impact)
        println("Solution x = $(xx)")
        println("Expected return $(expret) for gamma $(gamma)")
    end
end
