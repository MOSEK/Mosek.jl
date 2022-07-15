#   Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
#
#   File :      sensitivity.jl
#
#   Purpose :   To demonstrate how to perform sensitivity
#               analysis from the API on a small problem:
#
#               minimize
#
#               obj: +1 x11 + 2 x12 + 5 x23 + 2 x24 + 1 x31 + 2 x33 + 1 x34
#               st
#               c1:   +  x11 +   x12                                           <= 400
#               c2:                  +   x23 +   x24                           <= 1200
#               c3:                                  +   x31 +   x33 +   x34   <= 1000
#               c4:   +  x11                         +   x31                   = 800
#               c5:          +   x12                                           = 100
#               c6:                  +   x23                 +   x33           = 500
#               c7:                          +   x24                 +   x34   = 500
#
#               The example uses basis type sensitivity analysis.

using Mosek

prob_ptf = "Task ''
Objective
  Minimize  + 1 x11 + 2 x12 + 5 x23 + 2 x24 + 1 x31 + 2 x33 + 1 x34
Constraints
  c1 [-inf;400 ]  +  x11 +   x12
  c2 [-inf;1200]                +   x23 +   x24
  c3 [-inf;1000]                                 +   x31 +   x33 +   x34
  c4 [800]        +  x11                         +   x31
  c5 [100]               +   x12
  c6 [500]                       +   x23                 +   x33
  c7 [500]                               +   x24                 +   x34
Variables
  x11 [0;+inf]
  x12 [0;+inf]
  x23 [0;+inf]
  x24 [0;+inf]
  x31 [0;+inf]
  x33 [0;+inf]
  x34 [0;+inf]
"

maketask() do task
    putstreamfunc(task,MSK_STREAM_LOG,msg -> print(msg))

    readptfstring(task,prob_ptf)

    # A maximization problem

    optimize(task)
    checkmem(task,@__FILE__,@__LINE__)
    # Analyze upper bound on c1 and the equality constraint on c4
    let subi = [1, 4],
        marki = [MSK_MARK_UP, MSK_MARK_UP],
        # Analyze lower bound on the variables x12 and x31
        subj  = [2, 5],
        markj = [MSK_MARK_LO,MSK_MARK_LO],
        (leftpricei,
         rightpricei,
         leftrangei,
         rightrangei,
         leftpricej,
         rightpricej,
         leftrangej,
         rightrangej) = primalsensitivity(task,
                                          subi,
                                          marki,
                                          subj,
                                          markj)
        checkmem(task,@__FILE__,@__LINE__)
        println("Results from sensitivity analysis on bounds:  ")

        println("For constraints:")
        for i in 1:2
            println("leftprice = $(leftpricei[i]) | rightprice = $(rightpricei[i]) | leftrange = $(leftrangei[i]) | rightrange = $(rightrangei[i])");
        end
        println("For variables:")
        for i in 1:2
            println("leftprice = $(leftpricej[i]) rightprice = $(rightpricej[i]) leftrange = $(leftrangej[i]) rightrange = $(rightrangej[i])")
        end
    end

    let subc  = Int32[2, 5],
        (leftprice,
         rightprice,
         leftrange,
         rightrange) = dualsensitivity(task,subc)

        println("Results from sensitivity analysis on objective coefficients:")

        for i in 1:2
            println("leftprice = $(leftprice[i]) rightprice = $(rightprice[i]) leftrange = $(leftrange[i]) rightrange = $( rightrange[i])")
        end
    end
end
