using Mosek

"""
    basicPortfolio(n::Int, x0 :: Vector{Float64})

Implements the simple Portfolio Selection example maximizing return subject to a limited risk.

    maximize
        mu' * x        
    such that
        sum(x) = budget
        GT * x - s = 0
        s >= || t ||
        s <= gamma

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
        s_start  = n+1
        s_end    = n+1
        t_start  = n+2
        t_end    = t_start+n-1
        
        putobjsense(t,MSK_OBJECTIVE_SENSE_MAXIMIZE)
        putclist(t,Int32[Int32(x_start):Int32(x_end) ...], mu)
        
        putvarboundslice(t,Int32(x_start), Int32(x_end+1),
                         fill(MSK_BK_LO,n),
                         zeros(Float64,n),
                         fill(Inf,n))

        putvarboundslice(t,Int32(t_start), Int32(t_end+1),
                         fill(MSK_BK_LO,n),
                         zeros(Float64,n),
                         fill(Inf,n))
        
        for j in 1:n
            putvarname(t,Int32(x_start+j-1), "x$j")
            putvarname(t,Int32(t_start+j-1), "t$j")
        end
        putvarname(t,Int32(s_start),"s")
                        
        # constraints
        appendcons(t,Int32(1+n))

        # sum(x) = w+sum(x0)
        putconbound(t,Int32(1),MSK_BK_FX,w + sum(x0),w + sum(x0))
        putconname(t,Int32(1),"budget")
        putarow(t,Int32(1),Int32[Int32(x_start):Int32(x_end)...],ones(Float64,n))

        # GT * x - s = 0
        putconboundslice(t,2,n+2,
                         fill(MSK_BK_FX, n),
                         zeros(Float64,n),
                         zeros(Float64,n))

        for i in 1:n
            putconname(t,Int32(i+1),"GT[$i]")            
            putarow(t,Int32(i+1),Int32[Int32(x_start):Int32(x_end)...],GT[i,:])
            putaij(t,Int32(i+1),s_start+i-1, -1.0)
        end
        
        # s > || t ||
        appendcone(t,MSK_CT_QUAD, 0.0, [ s_start, t_start:t_end... ])
        putconename(t,Int32(1),"stddev")
        
        # s = gamma
        putvarbound(t,Int32(s_start),MSK_BK_FX,gamma,gamma)


        # Turn all log output off.
        putintparam(t,MSK_IPAR_LOG, Int32(0))

        optimize(t,"mosek://solve.mosek.com:30080")

        # Display the solution summary for quick inspection of results.
        solutionsummary(t, MSK_STREAM_MSG) 


        xxvector = getxx(t,MSK_SOL_ITR)
        
        #expectedReturn = sum(xxvector[x_start:x_end] .* mu)
        #standardDeviation = xxvector[s_start]

        xxvector[x_start:x_end],xxvector[s_start]
    end
end


# default data
gamma = 0.05
mu = [0.1073, 0.0737, 0.0627]
GT = [0.1667 0.0232  0.0013
      0.0000 0.1033 -0.0022
      0.0000 0.0000  0.0338]
x0 = [0.0, 0.0, 0.0]
w = 1.0

x,s = basicPortfolio(w,x0,mu,GT,gamma)

expectedReturn = sum(x .* mu)
println("------------------------------------------------------------")
println("Expected return $expectedReturn for gamma $s")
println("Optimal investment: $x")
println("------------------------------------------------------------")
