#
#  Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
#
#  File :      acc1.jl
#
#  Purpose :   Tutorial example for affine conic constraints.
#              Models the problem:
#
#              maximize c^T x
#              subject to  sum(x) = 1
#                          gamma >= |Gx+h|_2
##

using Mosek

# Define problem data
n, k = 3, 2


let Gsubi = Int64[1, 1, 2, 2],
    Gsubj = Int32[1, 2, 1, 3],
    Gval  = Float64[1.5, 0.1, 0.3, 2.1]

    # Make a MOSEK environment
    maketask() do task
        # Attach a printer to the task
        putstreamfunc(task,MSK_STREAM_LOG,msg -> print(msg))

        # Create n free variables
        appendvars(task,n)
        putvarboundsliceconst(task,1, n+1, MSK_BK_FR, -Inf, Inf)

        # Set up the objective
        putobjsense(task,MSK_OBJECTIVE_SENSE_MAXIMIZE)
        putclist(task,[1,2,3],[2.0,3.0,-1.0])

        # One linear constraint - sum(x) = 1
        appendcons(task,1)
        putarow(task,1, [1,2,3], [1.0,1.0,1.0])
        putconbound(task,1, MSK_BK_FX, 1.0, 1.0)

        # Append empty AFE rows for affine expression storage
        appendafes(task,k + 1)

        # G matrix in sparse form
        # Other data
        h     = Float64[0, 0.1]
        gamma = 0.03

        # Construct F matrix in sparse form
        Fsubi = Int32[i + 1 for i in Gsubi]   # G will be placed from row number 1 in F
        Fsubj = Gsubj
        Fval  = Gval

        # Fill in F storage
        putafefentrylist(task,Fsubi, Fsubj, Fval)

       # Fill in g storage
        putafeg(task,1, gamma)
        putafegslice(task,2, k+2, h)

        # Define a conic quadratic domain
        quadDom = appendquadraticconedomain(task,k + 1)

        # Create the ACC
        appendaccseq(task,
                     quadDom,    # Domain index
                     1,          # Indices of AFE rows [0,...,k]
                     zeros(k+1)) # Ignored

        # Solve and retrieve solution
        optimize(task)
        writedata(task,"acc1.ptf")

        xx = getxx(task,MSK_SOL_ITR)
        println("Solution: $xx")
        # Demonstrate retrieving activity of ACC
        activity = evaluateacc(task,MSK_SOL_ITR,
                               1)     # ACC index
        println("Activity of ACC:: $activity")

        # Demonstrate retrieving the dual of ACC
        doty = getaccdoty(task,MSK_SOL_ITR,
                          1)          # ACC index
        println("Dual of ACC: $doty")

        compl = sum(activity' * doty)

        @assert (abs(compl) < 1e-7) "Complementarity is invalid"
        @assert (maximum(abs.(xx      -[-0.07838011145615721, 1.1289128998004547, -0.0505327883442975])) < 1e-7) "Variable solution is incorrect"
        @assert (maximum(abs.(doty    -[-1.9429680870375095, -0.30303030303030304, -1.9191919191919191])) < 1e-7) "Constraint dual solution is incorrect"
        @assert (maximum(abs.(activity-[0.03, -0.004678877204190343, -0.029632888959872067])) < 1e-7) "Constraint level solution is incorrect"
        println("Complementarity $compl")

    end
end
