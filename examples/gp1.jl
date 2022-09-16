#
#   Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
#
#   File :      gp1.jl
#
#   Purpose:   Demonstrates how to solve a simple Geometric Program (GP)
#              cast into conic form with exponential cones and log-sum-exp.
#
#              Example from
#                https:#gpkit.readthedocs.io/en/latest/examples.html#maximizing-the-volume-of-a-box
#

using Mosek

function max_volume_box(Aw    :: Float64,
                        Af    :: Float64,
                        alpha :: Float64,
                        beta  :: Float64,
                        gamma :: Float64,
                        delta :: Float64)
    numvar = 3 # Variables in original problem
        # Create the optimization task. 

    maketask() do task
        putstreamfunc(task,MSK_STREAM_LOG,msg -> print(msg))

        # Add variables and constraints
        appendvars(task,numvar)

        x = 1
        y = 2
        z = 3

        # Objective is the sum of three first variables
        putobjsense(task,MSK_OBJECTIVE_SENSE_MAXIMIZE)
        putcslice(task,1, numvar+1, [1.0,1.0,1.0])

        putvarboundsliceconst(task,1, numvar+1, MSK_BK_FR, -Inf, Inf)

        appendcons(task,3)
        # s0+s1 < 1 <=> log(s0+s1) < 0
        putaijlist(task,
                   [1,1,2,2,3,3],
                   [y, z, x, y, z, y],
                   [1.0, 1.0, 1.0, -1.0, 1.0, -1.0])

        putconbound(task,1,MSK_BK_UP,-Inf,log(Af))
        putconbound(task,2,MSK_BK_RA,log(alpha),log(beta))
        putconbound(task,3,MSK_BK_RA,log(gamma),log(delta))


        let afei = getnumafe(task)+1,
            u1 = getnumvar(task)+1,
            u2 = u1+1,
            afeidx = [1, 2, 3, 3, 4, 4, 6, 6],
            varidx = [u1, u2, x, y, x, z, u1, u2],
            fval   = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            gfull  = [0.0, 0.0, log(2.0/Aw), log(2.0/Aw), 1.0, -1.0]

            appendvars(task,2)
            appendafes(task,6)

            putvarboundsliceconst(task,u1, u1+2, MSK_BK_FR, -Inf, Inf)

            # Affine expressions appearing in affine conic constraints
            # in this order:
            # u1, u2, x+y+log(2/Awall), x+z+log(2/Awall), 1.0, u1+u2-1.0
            putafefentrylist(task,afeidx, varidx, fval)
            putafegslice(task,afei, afei+6, gfull)

            let dom = appendprimalexpconedomain(task)

                # (u1, 1, x+y+log(2/Awall)) \in EXP
                appendacc(task,dom, [1, 5, 3], [0.0,0.0,0.0])

                # (u2, 1, x+z+log(2/Awall)) \in EXP
                appendacc(task,dom, [2, 5, 4], [0.0,0.0,0.0])
            end
            let dom = appendrzerodomain(task,1)
                # The constraint u1+u2-1 \in \ZERO is added also as an ACC
                appendacc(task,dom, [6], [0.0])
            end
        end

        optimize(task)
        writedata(task,"gp1.ptf")

        exp.(getxxslice(task,MSK_SOL_ITR, 1, numvar+1))
    end # maketask
end # max_volume_box

# maximize     h*w*d
# subjecto to  2*(h*w + h*d) <= Awall
#              w*d <= Afloor
#              alpha <= h/w <= beta
#              gamma <= d/w <= delta
#
# Variable substitutions:  h = exp(x), w = exp(y), d = exp(z).
#
# maximize     x+y+z
# subject      log( exp(x+y+log(2/Awall)) + exp(x+z+log(2/Awall)) ) <= 0
#                              y+z <= log(Afloor)
#              log( alpha ) <= x-y <= log( beta )
#              log( gamma ) <= z-y <= log( delta )
#
# Finally, the model we will implement:
#
# maximize   x+y+z
# subject to s0 > exp(x+y+log(2/Awall); (s0,1,x+y+log(2/Awall)) in PEXP
#            s1 > exp(x+z+log(2/Awall); (s1,1,x+z+log(2/Awall)) in PEXP
#            s0+s1 < 1
#
#            y+z < log Afloor
#
#            x-y in [log alpha; log beta]
#            z-y in [log gamma; log delta]
#
#            (x,y,z) in pexp : x0 > x1 * exp(x2/x1)

hwd = let Aw    = 200.0,
    Af    = 50.0,
    alpha = 2.0,
    beta  = 10.0,
    gamma = 2.0,
    delta = 10.0

    max_volume_box(Aw, Af, alpha, beta, gamma, delta)
end
println("h=$(hwd[1]) w=$(hwd[2]) d=$(hwd[3])\n");

@assert maximum(abs.(hwd-[8.164, 4.082, 8.167])) < 1e-3

