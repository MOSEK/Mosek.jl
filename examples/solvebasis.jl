##
#
#  File: solvebasis.jl
#
#  Purpose :   To demonstrate the usage of
#              MSK_solvewithbasis on the problem:
#
#              maximize  x0 + x1
#              st.
#                      x0 + 2.0 x1 <= 2
#                      x0  +    x1 <= 6
#                      x0 >= 0, x1>= 0
#
#               The problem has the slack variables
#               xc0, xc1 on the constraints
#               and the variabels x0 and x1.
#
#               maximize  x0 + x1
#               st.
#                  x0 + 2.0 x1 -xc1       = 2
#                  x0  +    x1       -xc2 = 6
#                  x0 >= 0, x1>= 0,
#                  xc1 <=  0 , xc2 <= 0
##

using Mosek

maketask() do task


    putobjname(task,"solvebasis")
    

    putstreamfunc(task,MSK_STREAM_LOG,msg -> print(msg))
    
    numcon = 2
    numvar = 2

    c = Float64[1.0, 1.0]
    ptrb = Int64[1,3]
    ptre = Int64[3,4]
    asub = Int32[1,2,
                 1,2]
    aval = Float64[1.0, 1.0,
                   2.0, 1.0]
    bkc = Boundkey[MSK_BK_UP
                   MSK_BK_UP]

    blc = Float64[-Inf
                  -Inf]
    buc = Float64[2.0
                  6.0]

    bkx = Boundkey[MSK_BK_LO
                   MSK_BK_LO]
    blx = Float64[0.0
                  0.0]

    bux = Float64[Inf
                  Inf]
    w1 = Float64[2.0, 6.0]
    w2 = Float64[1.0, 0.0]

    inputdata(task,
              Int32(numcon), Int32(numvar),
              c,
              0.0,
              ptrb,
              ptre,
              asub,
              aval,
              bkc,
              blc,
              buc,
              bkx,
              blx,
              bux)
    
    putobjsense(task,MSK_OBJECTIVE_SENSE_MAXIMIZE)
    
    
    r = optimize(task,"mosek://solve.mosek.com:30080")
    if r != MSK_RES_OK
        println("Mosek warning: $r")
    end
    
    basis = initbasissolve(task)
    
    #List basis variables corresponding to columns of B
    varsub = Int32[1, 2]
    
    for i in 1:numcon
        if basis[varsub[i]] < numcon
            println("Basis variable no $i is xc$(basis[i])")
        else
            println("Basis variable no $i is x$(basis[i]-numcon)")
            
            # solve Bx = w1
            # varsub contains index of non-zeros in b.
            #  On return b contains the solution x and
            # varsub the index of the non-zeros in x.
            nz = 2

            nz = solvewithbasis(task,Int32(0), nz, varsub, w1)
            println("nz = $nz")
            println("Solution to Bx = $w1")

            for i in 1:nz
                if basis[varsub[i]] < numcon
                    println("xc $(basis[varsub[i]]) = $(w1[varsub[i]])")
                else
                    println("x$(basis[varsub[i]] - numcon) = $(w1[varsub[i]])")
                end
            end

            # Solve B^Tx = w2
            nz = 1
            varsub[1] = 1

            
            nz = solvewithbasis(task,1,nz,varsub,w2)
            println(nz)

            println("Solution to B^Tx = $(w2)")

            for i in 1:nz
                if basis[varsub[i]] < numcon
                    print("xc$(basis[varsub[i]]) = $(w2[varsub[i]])")
                else
                    print("x$(basis[varsub[i]] - numcon) = $(w2[varsub[i]])")
                end
            end
        end
    end

end
