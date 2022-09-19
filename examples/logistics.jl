#
# Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
#
# File :      logistics.jl
#
# Purpose: Implements logistic regression with regulatization.
#
#          Demonstrates using the exponential cone and log-sum-exp in Optimizer API.

using Mosek

"""
Adds ACCs for t_i >= log ( 1 + exp((1-2*y[i]) * theta' * X[i]) )
Adds auxiliary variables, AFE rows and constraints
"""
function softplus(task :: Mosek.Task, d::Int, n::Int, theta::Int, t::Int, X::Matrix{Float64}, y::Vector{Bool})
    nvar = getnumvar(task)
    ncon = getnumcon(task)
    nafe = getnumafe(task)
    z1       = nvar
    z2       = z1+n
    zcon     = ncon
    thetaafe = nafe
    tafe     = thetaafe+n
    z1afe    = tafe+n
    z2afe    = z1afe+n

    appendvars(task,2*n)   # z1, z2
    appendcons(task,n)     # z1 + z2 = 1
    appendafes(task,4*n)   #theta * X[i] - t[i], -t[i], z1[i], z2[i]

    # Linear constraints
    let subi = Vector{Int32}(undef,2*n),
        subj = Vector{Int32}(undef,2*n),
        aval = Vector{Float64}(undef,2*n),
        k = 1

        for i in 1:n
            # z1 + z2 = 1
            subi[k] = zcon+i;  subj[k] = z1+i;  aval[k] = 1;  k += 1
            subi[k] = zcon+i;  subj[k] = z2+i;  aval[k] = 1;  k += 1
        end

        putaijlist(task,subi, subj, aval)
    end
    putconboundsliceconst(task,zcon+1, zcon+n+1, MSK_BK_FX, 1, 1)
    putvarboundsliceconst(task,nvar+1, nvar+2*n+1, MSK_BK_FR, -Inf, Inf)

    # Affine conic expressions
    let afeidx = Vector{Int64}(undef,d*n+4*n),
        varidx = Vector{Int32}(undef,d*n+4*n),
        fval   = Vector{Float64}(undef,d*n+4*n),
        k = 1
        # Thetas
        for i in 1:n
            for j in 1:d
                afeidx[k] = thetaafe + i; varidx[k] = theta + j
                fval[k]   = if y[i] X[i,j]*(-1.0) else X[i,j] end
                k += 1
            end
        end

        # -t[i]
        for i in 1:n
            afeidx[k] = thetaafe + i; varidx[k] = t + i; fval[k] = -1; k += 1
            afeidx[k] = tafe + i;     varidx[k] = t + i; fval[k] = -1; k += 1
        end

        # z1, z2
        for i in 1:n
            afeidx[k] = z1afe + i; varidx[k] = z1 + i; fval[k] = 1; k += 1
            afeidx[k] = z2afe + i; varidx[k] = z2 + i; fval[k] = 1; k += 1
        end

        # Add the expressions
        putafefentrylist(task,afeidx, varidx, fval)
    end

    # Add a single row with the constant expression "1.0"

    let oneafe = getnumafe(task)+1,
        # Add an exponential cone domain
        expdomain = appendprimalexpconedomain(task)

        appendafes(task,1)
        putafeg(task,oneafe,1.0)

        # Conic constraints
        for i in 1:n
            appendacc(task,expdomain, [z1afe+i, oneafe, thetaafe+i], zeros(3))
            appendacc(task,expdomain, [z2afe+i, oneafe, tafe+i],     zeros(3))
        end
    end
end # softplus


"""
Model logistic regression (regularized with full 2-norm of theta)
X - n x d matrix of data points
y - length n vector classifying training points
lamb - regularization parameter
"""
function logisticRegression(X    :: Matrix{Float64},
                            y    :: Vector{Bool},
                            lamb :: Float64)
    (n,d) = size(X)

    maketask() do task
        # Variables [r; theta; t]
        nvar = 1+d+n
        appendvars(task,nvar)
        putvarboundsliceconst(task,1, nvar+1, MSK_BK_FR, -Inf, Inf)

        r     = 0
        theta = r+1
        t     = theta+d

        # Objective lambda*r + sum(t)
        putcj(task,r+1,lamb)
        for i in 1:n
            putcj(task,t+i, 1.0)
        end

        # Softplus function constraints
        softplus(task, d, n, theta, t, X, y)

        # Regularization
        # Append a sequence of linear expressions (r, theta) to F
        numafe = getnumafe(task)
        appendafes(task,1+d)
        putafefentry(task,numafe+1, r+1, 1.0)

        for i in 1:d
            putafefentry(task,numafe+i+1, theta+i, 1.0)
        end

        # Add the constraint
        appendaccseq(task,appendquadraticconedomain(task,1+d), numafe+1, zeros(d+1))

        # Solution
        optimize(task)

        getxxslice(task,MSK_SOL_ITR, theta+1, theta+d+1)
    end
end # logisticRegression




# Test: detect and approximate a circle using degree 2 polynomials
let n = 30
    X = Matrix{Float64}(undef,(n*n,6))
    Y = Vector{Bool}(undef,n*n)

    k = 1
    for i in 1:n
        for j in 1:n
            x = -1 + 2.0*i/(n-1)
            y = -1 + 2.0*j/(n-1)

            X[k,:] = [1.0, x, y, x*y, x*x, y*y]

            # X[k,0] = 1.0; X[k,1] = x; X[k,2] = y; X[k,3] = x*y
            # X[k,4] = x*x; X[k,5] = y*y
            Y[k]   = x*x+y*y>=0.69

            k += 1
        end
    end

    theta = logisticRegression(X, Y, 0.1)
    println("theta = $theta")

end
