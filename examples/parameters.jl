##
#   Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
#
#   File :      parameters.jl
#
#   Purpose :   Demonstrates a very simple example about how to get/set
#               parameters with MOSEK Julia API
#

using Mosek

maketask() do task
    println("Test MOSEK parameter get/set functions");

    # Set log level (integer parameter)
    putintparam(task,MSK_IPAR_LOG, 1)
    # Select interior-point optimizer... (integer parameter)
    putintparam(task,MSK_IPAR_OPTIMIZER, MSK_OPTIMIZER_INTPNT)
    # ... without basis identification (integer parameter)
    putintparam(task,MSK_IPAR_INTPNT_BASIS,MSK_BI_NEVER)
    # Set relative gap tolerance (double parameter)
    putdouparam(task,MSK_DPAR_INTPNT_CO_TOL_REL_GAP, 1.0e-7)

    # The same using explicit string names
    putparam(task,"MSK_DPAR_INTPNT_CO_TOL_REL_GAP", "1.0e-7");
    putnadouparam(task,"MSK_DPAR_INTPNT_CO_TOL_REL_GAP",  1.0e-7 )

    # Incorrect value
    try
        putdouparam(task,MSK_DPAR_INTPNT_CO_TOL_REL_GAP, -1.0)
    catch
        println("Wrong parameter value")
    end


    param = getdouparam(task,MSK_DPAR_INTPNT_CO_TOL_REL_GAP)
    println("Current value for parameter intpnt_co_tol_rel_gap = $param")

    # Define and solve an optimization problem here
    # optimize(task,)
    # After optimization:

    println("Get MOSEK information items")

    tm = getdouinf(task,MSK_DINF_OPTIMIZER_TIME)
    iter = getintinf(task,MSK_IINF_INTPNT_ITER)

    println("Time: $tm");
    println("Iterations: $iter");

end
