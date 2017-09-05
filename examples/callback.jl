using Mosek


function intProduction(productNames     :: Vector{String},
                       processNames     :: Vector{String},
                       timeRequirements :: Array{Float64,2},
                       timeResources    :: Vector{Float64},
                       profit           :: Vector{Float64})
    n = length(productNames)
    m = length(processNames)
    maketask() do t
        appendvars(t,n)
        appendcons(t,m)
        
        for j in 1:n
            putvarname(t,Int32(j),productNames[j])
            putvarbound(t,Int32(j),MSK_BK_LO,0.0,Inf)
            putvartype(t,Int32(j),MSK_VAR_TYPE_INT)
        end

        for i in 1:m
            putconname(t,Int32(i),processNames[i])
            putconbound(t,i, MSK_BK_UP, -Inf,timeResources[i])
            putarow(t,i,Int32[1:n...],timeRequirements[:,i])
        end

        putobjsense(t,MSK_OBJECTIVE_SENSE_MAXIMIZE)
        putclist(t,Int32[1:n...],profit)


        function callbackfunc(code::Callbackcode, dinf, iinf, liinf)
            if code == MSK_CALLBACK_NEW_INT_MIO && solutiondef(t,MSK_SOL_ITG)
                xx = getxx(t,MSK_SOL_ITG)
                println("New solution: $xx")
            end
            0
        end        

        putcallbackfunc(t,callbackfunc)

        writedata(t,"test.opf")
        optimize(t)
        
        xx = getxx(t,MSK_SOL_ITG)
        println("Final solution: $xx")
    end
end


timeRequirements = Float64[ 2.0  3.0  2.0  4.0
                            4.0  2.0  3.0  0.0
                            1.0  2.0  1.0  1.0 
                            3.0  3.0  2.0  1.0 ]
timeResources    = Float64[100.0, 50.0, 30.0, 60.0]
profit           = Float64[1.7, 2.3, 1.0, 1.6]
processNames     = String["Assembly", "Polish","Varnish", "Packing" ]
productNames     = String["Chairs", "Tables", "Beds", "CuttingBoards"]


intProduction(productNames,processNames,timeRequirements,timeResources,profit)
