using Mosek


##
#  Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
#
#  File :      callback.jl
#
#  Purpose :   To demonstrate how to use the progress
#              callback.
#
#              Use this script as follows:
#
#               callback.py psim  25fv47.mps
#               callback.py dsim  25fv47.mps
#               callback.py intpnt 25fv47.mps
#
#              The first argument tells which optimizer to use
#              i.e. psim is primal simplex, dsim is dual simplex
#              and intpnt is interior-point.
##
using Mosek
using Printf

function callback(task    :: Mosek.Task,
                  maxtime :: Float64,
                  caller  :: Callbackcode,
                  douinf  :: Vector{Float64},
                  intinf  :: Vector{Int32},
                  lintinf :: Vector{Int64})
        opttime = 0.0

        if caller == MSK_CALLBACK_BEGIN_INTPNT
            println("Starting interior-point optimizer")
        elseif caller == MSK_CALLBACK_INTPNT
            itrn    = intinf[MSK_IINF_INTPNT_ITER]
            pobj    = douinf[MSK_DINF_INTPNT_PRIMAL_OBJ]
            dobj    = douinf[MSK_DINF_INTPNT_DUAL_OBJ]
            stime   = douinf[MSK_DINF_INTPNT_TIME]
            opttime = douinf[MSK_DINF_OPTIMIZER_TIME]

            println("Iterations: $itrn")
            @printf("  Elapsed time: %6.2f(%.2f) ",opttime, stime)
            @printf("  Primal obj.: %-18.6e  Dual obj.: %-18.6e",pobj, dobj)
        elseif caller == MSK_CALLBACK_END_INTPNT
            println("Interior-point optimizer finished.")
        elseif caller == MSK_CALLBACK_BEGIN_PRIMAL_SIMPLEX
            println("Primal simplex optimizer started.")
        elseif caller == MSK_CALLBACK_UPDATE_PRIMAL_SIMPLEX
            itrn = intinf[MSK_IINF_SIM_PRIMAL_ITER]
            pobj = douinf[MSK_DINF_SIM_OBJ]
            stime = douinf[MSK_DINF_SIM_TIME]
            opttime = douinf[MSK_DINF_OPTIMIZER_TIME]

            println("Iterations: %-3d", itrn)
            println("  Elapsed time: %6.2f(%.2f)",opttime, stime)
            println("  Obj.: %-18.6e", pobj)
        elseif caller == MSK_CALLBACK_END_PRIMAL_SIMPLEX
            println("Primal simplex optimizer finished.")
        elseif caller == MSK_CALLBACK_BEGIN_DUAL_SIMPLEX
            println("Dual simplex optimizer started.")
        elseif caller == MSK_CALLBACK_UPDATE_DUAL_SIMPLEX
            itrn = intinf[MSK_IINF_SIM_DUAL_iter]
            pobj = douinf[MSK_DINF_SIM_OBJ]
            stime = douinf[MSK_DINF_SIM_TIME]
            opttime = douinf[MSK_DINF_OPTIMIZER_TIME]
            println("Iterations: %-3d",itrn)
            println("  Elapsed time: %6.2f(%.2f)",opttime, stime)
            println("  Obj.: %-18.6e",pobj)
        elseif caller == MSK_CALLBACK_END_DUAL_SIMPLEX
            println("Dual simplex optimizer finished.")
        elseif caller == MSK_CALLBACK_NEW_INT_MIO
            println("New integer solution has been located.")
            xx = getxx(task,MSK_SOL_ITG)
            println("  x = $xx")
            println("Obj.: %f", douinf[MSK_DINF_MIO_OBJ_INT])
        end

        if opttime >= maxtime
            # mosek is spending too much time. Terminate it.
            println("Terminating.")
            1
        else
            0
        end
end


# To run a continuous problem example

if length(ARGS) < 2
    println("Usage: callback ( psim | dsim | intpnt ) filename")
else
    slvr = ARGS[1]
    filename = ARGS[2]

    maketask() do task
        readdata(task,filename)

        if     slvr == "psim"
            putintparam(task,MSK_IPAR_OPTIMIZER, MSK_OPTIMIZER_PRIMAL_SIMPLEX)
        elseif slvr == "dsim"
            putintparam(task,MSK_IPAR_OPTIMIZER, MSK_OPTIMIZER_DUAL_SIMPLEX)
        elseif slvr == "intpnt"
            putintparam(task,MSK_IPAR_OPTIMIZER, MSK_OPTIMIZER_INTPNT)
        end
        putcallbackfunc(task,(caller,dinf,iinf,linf) -> callback(task,0.05,caller,dinf,iinf,linf))

        optimize(task)
    end
end

