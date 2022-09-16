
#  Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
#
#  File :      feasrepairex1.jl
#
#  Purpose :   To demonstrate how to use the MSK_relaxprimal function to
#              locate the cause of an infeasibility.
#
#  Syntax :     On command line
#
#                  feasrepairex1.jl [ feasrepair.lp | - ]
#
#               feasrepair.lp is located in mosek\<version>\tools\examples.

feasrepair_lp = """
minimize
 obj: - 10 x1 - 9 x2
st
 c1: + 7e-01 x1 + x2 <= 630
 c2: + 5e-01 x1 + 8.333333333e-01 x2 <= 600
 c3: + x1 + 6.6666667e-01 x2 <= 708
 c4: + 1e-01 x1 + 2.5e-01 x2 <= 135
bounds
x2 >= 650
end
"""

using Mosek

if length(ARGS) < 1
    println("Syntax: feasrepairex1 [ FILENAME | - ]")
else
    filename = ARGS[1]
    maketask() do task
        putstreamfunc(task,MSK_STREAM_LOG,msg -> print(msg))

        if filename != "-"
            readdata(task,filename)
        else
            readlpstring(task,feasrepair_lp)
        end

        putintparam(task,MSK_IPAR_LOG_FEAS_REPAIR, 3)

        numvar = getnumvar(task)
        numcon = getnumcon(task)
        println(numvar,",",numcon)
        primalrepair(task,zeros(numcon),zeros(numcon),zeros(numvar),zeros(numvar))

        sum_viol = getdouinf(task,MSK_DINF_PRIMAL_REPAIR_PENALTY_OBJ);

        println("Minimized sum of violations = $(sum_viol)");

        optimize(task)
        solutionsummary(task,MSK_STREAM_LOG)
    end
end
