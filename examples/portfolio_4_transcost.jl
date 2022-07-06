
# File : portfolio_4_transcost
#
# Copyright : Mosek ApS
#
# Description :  Implements a basic portfolio optimization model.

include("portfolio_data.jl")

using Mosek
#TAG:begin-code
#TAG:begin-basic-markowitz

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

        #TAG:begin-offsets
        #Offset of variables into the API variable.
        x_ofs = 0
        z_ofs = n
        y_ofs = 2*n

        # Constraints offsets
        #TAG:end-offsets

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
            putvarbound(task,x_ofs+j, MSK_BK_LO, 0.0, Inf);
            putvarbound(task,z_ofs+j, MSK_BK_FR, -Inf, Inf);
            putvarbound(task,y_ofs+j, MSK_BK_RA, 0.0, 1.0);
            putvartype(task, y_ofs+j, MSK_VAR_TYPE_INT);
        end

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

            #TAG:begin-basic-markowitz-appendaccseq
            # Input the affine conic constraint (gamma, GT*x) \in QCone
            # Add the quadratic domain of dimension k+1
            qdom = appendquadraticconedomain(task,k + 1)
            # Add the constraint
            appendaccseq(task,qdom,1,zeros(k+1))
            #TAG:end-basic-markowitz-appendaccseq
            putaccname(task,acci+1, "risk")
        end


        # Objective: maximize expected return mu^T x
        putclist(task,[x_ofs+1:x_ofs+n...],mu)
        putobjsense(task,MSK_OBJECTIVE_SENSE_MAXIMIZE)

        optimize(task)

        writedata(task,"portfolio_4_transcost.ptf")
        
        #TAG:begin-solutionsummary
        # Display solution summary for quick inspection of results
        solutionsummary(task,MSK_STREAM_LOG)
        #TAG:end-solutionsummary

        #TAG:begin-writedata
        writedata(task,"portfolio_1_basic.ptf");
        #TAG:end-writedata

        # Read the results
        xx = getxxslice(task,MSK_SOL_ITG, x_ofs+1,x_ofs+n+1)
        expret = mu' * xx

        (xx,expret)
    end
end # portfolio()
#TAG:end-code
#TAG:end-basic-markowitz

gamma = 0.36
let (xx,expret) = portfolio(mu,x0,w,gamma,GT,fixed_tcost,proportional_tcost)
    println("Expected return $(expret) for gamma $(gamma)")
    println("Solution vector = $(xx)")
end
