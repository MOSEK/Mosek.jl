#
#  Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
#
#  File :      portfolio_4_transcost.jl
#
#  Description :  Implements a basic portfolio optimization model
#                 with fixed setup costs and transaction costs
#                 as a mixed-integer problem.

using Mosek

function portfolio( mu :: Vector{Float64},
                    x0 :: Vector{Float64},
                    w  :: Float64,
                    gamma :: Float64,
                    GT :: Array{Float64,2},
                    f  :: Vector{Float64},
                    g  :: Vector{Float64})
    (k,n) = size(GT)
    maketask() do task
        # Directs the log task stream
        putstreamfunc(task,MSK_STREAM_LOG,msg -> print(msg))

        totalBudget = sum(x0)+w

        #Offset of variables into the API variable.
        x_ofs = 0
        y_ofs = n
        z_ofs = 2*n

        # Constraints offsets

        # Holding variable x of length n
        # No other auxiliary variables are needed in this formulation
        appendvars(task,3*n)

        # Setting up variable x
        for j in 1:n
            # Optionally we can give the variables names
            putvarname(task, x_ofs+j, "x[$(j)]");
            putvarname(task, z_ofs+j, "z[$(j)]");
            putvarname(task, y_ofs+j, "y[$(j)]");
            # No short-selling - x^l = 0, x^u = inf
            putvartype(task, y_ofs+j, MSK_VAR_TYPE_INT);
        end
        putvarboundsliceconst(task,x_ofs+1,x_ofs+n+1, MSK_BK_LO, 0.0, Inf);
        putvarboundsliceconst(task,y_ofs+1,y_ofs+n+1, MSK_BK_RA, 0.0, 1.0);
        putvarboundsliceconst(task,z_ofs+1,z_ofs+n+1, MSK_BK_FR, -Inf, Inf);

        # One linear constraint: total budget
        let coni = getnumcon(task)
            appendcons(task,1);
            putconname(task,coni+1,"budget")

            consub = Int32[ coni+1 for i in 1:n ]
            putaijlist(task,consub,[x_ofs+1:x_ofs+n...], ones(n))
            putaijlist(task,consub,[z_ofs+1:z_ofs+n...],g)
            putaijlist(task,consub,[y_ofs+1:y_ofs+n...],f)

            putconbound(task, coni+1, MSK_BK_FX, totalBudget, totalBudget)
        end

        let coni = getnumcon(task)
            appendcons(task,2*n)
            for i in 1:n
                putconname(task,coni+i,"zabs1[$i]")
                putconname(task,coni+n+i,"zabs2[$i]")

                putaij(task,coni + i, x_ofs + i, -1.0)
                putaij(task,coni + i, z_ofs + i, 1.0)
                putconbound(task,coni + i, MSK_BK_LO, -x0[i], Inf)

                putaij(task,coni + n + i, x_ofs + i, 1.0)
                putaij(task,coni + n + i, z_ofs + i, 1.0)
                putconbound(task,coni + n + i, MSK_BK_LO, x0[i], Inf)
            end
        end

        # - Switch
        let coni = getnumcon(task)
            appendcons(task,n)
            for i in 1:n
                putconname(task,coni+i, "switch[$i]")
                putaij(task,coni + i, z_ofs + i, 1.0)
                putaij(task,coni + i, y_ofs + i, -totalBudget)
                putconbound(task,coni + i, MSK_BK_UP, -Inf, 0.0)
            end
        end


        let afei = getnumafe(task),
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
            qdom = appendquadraticconedomain(task,k + 1)
            # Add the constraint
            appendaccseq(task,qdom,1,zeros(k+1))
            putaccname(task,acci+1, "risk")
        end


        # Objective: maximize expected return mu^T x
        putclist(task,[x_ofs+1:x_ofs+n...],mu)
        putobjsense(task,MSK_OBJECTIVE_SENSE_MAXIMIZE)

        optimize(task)

        writedata(task,"portfolio_4_transcost.ptf")

        # Display solution summary for quick inspection of results
        solutionsummary(task,MSK_STREAM_LOG)

        writedata(task,"portfolio_4_transcost.ptf");

        # Read the results
        xx = getxxslice(task,MSK_SOL_ITG, x_ofs+1,x_ofs+n+1)
        expret = mu' * xx

        (xx,expret)
    end
end # portfolio()

let w = 1.0,
    mu = [0.07197, 0.15518, 0.17535, 0.08981, 0.42896, 0.39292, 0.32171, 0.18379],
    x0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    GT = [ 0.30758 0.12146 0.11341 0.11327 0.17625 0.11973 0.10435 0.10638
           0.      0.25042 0.09946 0.09164 0.06692 0.08706 0.09173 0.08506
           0.      0.      0.19914 0.05867 0.06453 0.07367 0.06468 0.01914
           0.      0.      0.      0.20876 0.04933 0.03651 0.09381 0.07742
           0.      0.      0.      0.      0.36096 0.12574 0.10157 0.0571
           0.      0.      0.      0.      0.      0.21552 0.05663 0.06187
           0.      0.      0.      0.      0.      0.      0.22514 0.03327
           0.      0.      0.      0.      0.      0.      0.      0.2202 ],
    (k,n) = size(GT),
    f = 0.01 * ones(n),
    g = 0.001 * ones(n),
    gamma = 0.36

    let (xx,expret) = portfolio(mu,x0,w,gamma,GT,f,g)
        println("Expected return $(expret) for gamma $(gamma)")
        println("Solution vector = $(xx)")
    end
end
