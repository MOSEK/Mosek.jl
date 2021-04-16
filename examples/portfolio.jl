using Mosek

"""
    basicPortfolio(n::Int, x0 :: Vector{Float64})

Implements the simple Portfolio Selection example maximizing return subject to a limited risk.

    maximize
        mu' * x
    such that
        sum(x) = budget
        GT * x - t = 0
        s >= || t ||
        s = gamma

Parameters:

- `w`     Initial uninvested wealth
- `x0`    Initial invested wealth for each asset
- `mu`    Expected returns for each asset
- `gamma` Risk limit as maximum standard deviation
- `GT`    Factorized covariance matrix
"""
function basicPortfolio(w :: Float64, x0 :: Vector{Float64}, mu :: Vector{Float64}, GT :: Array{Float64,2}, gamma :: Float64)
    n = length(mu)
    maketask() do t
        putstreamfunc(t,MSK_STREAM_LOG,msg -> print(msg))

        budget = w + sum(x0)

        # variables
        appendvars(t,2*n + 1)
        x_start  = 1
        x_end    = n
        s_idx  = n+1

        putobjsense(t,MSK_OBJECTIVE_SENSE_MAXIMIZE)
        putclist(t,Int32[Int32(x_start):Int32(x_end) ...], mu)

        putvarboundslice(t,Int32(x_start), Int32(x_end+1),
                         fill(MSK_BK_LO,n),
                         zeros(Float64,n),
                         fill(Inf,n))

        for j in 1:n
            putvarname(t,Int32(x_start+j-1), "x$j")
        end
        putvarname(t,Int32(s_idx),"s")

        appendcons(t,1)

        # sum(x) = w+sum(x0)
        putconbound(t,Int32(1),MSK_BK_FX,w + sum(x0),w + sum(x0))
        putconname(t,Int32(1),"budget")
        putarow(t,Int32(1),Int32[Int32(x_start):Int32(x_end)...],ones(Float64,n))

        # [ s ; GT * x ] in C_r
        domidx = appendquadraticconedomain(t,n+1)
        appendafes(t,1+n)

        appendacc(t,domidx,Int64[1:n+1...],zeros(Float64,4))
        putaccname(t,1,"GT")
        putafefrow(t,1,Int32[s_idx],Float64[-1])
        for i in 1:3
            putafefrow(t,i+1,
                       Int32[ x_start:x_end... ],
                       GT[i,:])
        end

        # s = gamma
        putvarbound(t,Int32(s_idx),MSK_BK_UP,gamma,gamma)

        # Turn all log output off.
        putintparam(t,MSK_IPAR_LOG, Int32(0))

        optimize(t,"mosek://solve.mosek.com:30080")

        # Display the solution summary for quick inspection of results.
        solutionsummary(t, MSK_STREAM_MSG)


        xxvector = getxx(t,MSK_SOL_ITR)

        #expectedReturn = sum(xxvector[x_start:x_end] .* mu)
        #standardDeviation = xxvector[s_start]

        xxvector[x_start:x_end],xxvector[s_idx]
    end
end

# default data
gamma = 0.05
mu = Float64[0.1073, 0.0737, 0.0627]
GT = Float64[0.1667 0.0232  0.0013
             0.0000 0.1033 -0.0022
             0.0000 0.0000  0.0338]
x0 = Float64[0.0, 0.0, 0.0]
w = 1.0

(x,s) = basicPortfolio(w,x0,mu,GT,gamma)

expectedReturn = sum(x .* mu)
println("------------------------------------------------------------")
println("Expected return $expectedReturn for gamma $s")
println("Optimal investment: $x")
println("------------------------------------------------------------")
