#     
#   File:      feasrepairex1.jl
#
#   Purpose:    To demonstrate how to use the MSK_relaxprimal function to
#               locate the cause of an infeasibility.
#
#   Syntax: On command line
#           julia-basic feasrepairex1.jl feasrepair.lp

using Mosek

printstream(msg::String) = print(msg)

maketask() do task
    putstreamfunc(task,MSK_STREAM_LOG,printstream)

    inputfile = ARGS[1]

    # Read data 
    readdata(task,inputfile)

    nvars = getnumvar(task)
    ncons = getnumcon(task)

    putintparam(task,MSK_IPAR_LOG_FEAS_REPAIR,3)

    try 
        primalrepair(task,ones(ncons),ones(ncons),ones(nvars),ones(nvars))

        sum_viol = getdouinf(task,MSK_DINF_PRIMAL_REPAIR_PENALTY_OBJ)
        println("Minimized sum of violations = $sum_viol\n")
        
        optimize(task,"mosek://solve.mosek.com:30080")

        solutionsummary(task,MSK_STREAM_MSG)
    catch e
        println("Failed to solve: $e")
    end
        
end
