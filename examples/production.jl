using Mosek
using Printf

# Demonstrates how to solve a linear optimization problem using the
# MOSEK API and modify and re-optimize the problem.


function production(productNames     :: Vector{String},
                    processNames     :: Vector{String},
                    timeRequirements :: Array{Float64,2},
                    timeResources    :: Vector{Float64},
                    profit           :: Vector{Float64},
                    # Modify A entry
                    mod_i               :: Int,
                    mod_j               :: Int,
                    mod_cof             :: Float64,
                    
                    # Add product Name
                    addProductName             :: String,
                    addProductTimeRequirements :: Vector{Float64},
                    addProfit                  :: Float64,
                    # Add product Name
                    addProcessName      :: String,
                    addTimeRequirements :: Vector{Float64},
                    addTimeResources    :: Float64 )
    n = length(productNames)
    m = length(processNames)
    maketask() do t
        appendvars(t,n)
        appendcons(t,m)

        
        for j in 1:n
            putvarname(t,Int32(j),productNames[j])
            putvarbound(t,Int32(j),MSK_BK_LO,0.0,Inf)
        end

        for i in 1:m
            putconname(t,Int32(i),processNames[i])
            putconbound(t,i, MSK_BK_UP, -Inf,timeResources[i])
            putarow(t,i,Int32[1:n...],timeRequirements[:,i])
        end

        putobjsense(t,MSK_OBJECTIVE_SENSE_MAXIMIZE)
        putclist(t,Int32[1:n...],profit)

        r = optimize(t,"mosek://solve.mosek.com:30080")
        xx = getxx(t,MSK_SOL_BAS)
        println("Initial solution:")
        for i in 1:n
            @printf("    %-15s: %10.1f\n",productNames[i],xx[i])
        end

        ############### Make a change to the A matrix ###############
        putaij(t,Int32(mod_i),Int32(mod_j),mod_cof)

        optimize(t,"mosek://solve.mosek.com:30080")
        xx = getxx(t,MSK_SOL_BAS)
        println("After modified coefficient:")
        for i in 1:n
            @printf("    %-15s: %10.1f\n",productNames[i],xx[i])
        end

        ############## Add a new variable / column ##################
        appendvars(t,Int32(1))
        putvarname(t,Int32(n+1),addProductName)
        putcj(t,Int32(n+1),addProfit)
        putacol(t,Int32(n+1),Int32[1:n+1...],addProductTimeRequirements)
        putvarbound(t,Int32(n+1),MSK_BK_LO,0.0,Inf)


        # reoptimize using simplex
        putintparam(t,MSK_IPAR_OPTIMIZER, MSK_OPTIMIZER_FREE_SIMPLEX)

        optimize(t,"mosek://solve.mosek.com:30080")
        xx = getxx(t,MSK_SOL_BAS)
        println("After added product (column):")
        for i in 1:n
            @printf("    %-15s: %10.1f\n",productNames[i],xx[i])
        end
        @printf("    %-15s: %10.1f\n",addProductName,xx[n+1])

        ################ Add a new constraint / row #################
        appendcons(t,Int32(1))
        putconname(t,Int32(m+1),addProcessName)
        putarow(t,Int32(m+1),Int32[1:n...],addTimeRequirements)
        putconbound(t,Int32(n+1),MSK_BK_UP,-Inf,addTimeResources)

        optimize(t,"mosek://solve.mosek.com:30080")
        xx = getxx(t,MSK_SOL_BAS)
        println("After added process (row):")
        for i in 1:n
            @printf("    %-15s: %10.1f\n",productNames[i],xx[i])
        end
        @printf("    %-15s: %10.1f\n",addProductName,xx[n+1])
    end
end


timeRequirements = Float64[ 2.0  3.0  2.0 
                            4.0  2.0  3.0 
                            3.0  3.0  2.0 ]
timeResources    = Float64[100000.0, 50000.0, 60000.0]
profit           = Float64[1.5, 2.5, 3.0]
processNames     = String["Assembly", "Polish","Packing" ]
productNames     = String["Chairs", "Tables", "Beds"]


production(productNames,processNames,timeRequirements,timeResources,profit,
           1,1,3.0, # modify a coefficient
           "Cutting_board", Float64[4.0,0.0,1.0], 1.0, # additional product
           "Varnish", Float64[1.0,2.0,1.0,1.0], 30000.0) # additional process
