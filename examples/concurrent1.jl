#
#  File:    concurrent1.jl
#
#  Purpose: Demonstrates a simple implementation of a concurrent optimizer.
#
#           The concurrent optimizer starts a few parallel optimizations
#           of the same problem using different algorithms, and reports
#           a solution when the first optimizer is ready.
#
#           This example also demonstrates how to define a simple callback handler
#           that stops the optimizer when requested.


using Mosek

# Defines a Mosek callback function whose only function
# is to indicate if the optimizer should be stopped.
stop = false
firstStop = 0
function callback(caller  :: Callbackcode,
                  douinf  :: Vector{Float64},
                  intinf  :: Vector{Int32},
                  lintinf :: Vector{Int64})
    if stop
        1
    else
        0
    end
end

# firstOK, res, trm = optimize(tasks)
#
# Takes a list of tasks and optimizes then in parallel. The
# response code and termination code from each optimization is
# returned in ``res`` and ``trm``.
#
# When one task completes with rescode == ok, others are terminated.
#
# firstOK is the index of the first task that returned
# with rescode == ok. Whether or not this task contains the
# most valuable answer, is for the caller to decide. If none
# completed without error returns -1.
function runTask(num, task)
    global stop
    global firstStop

    ttrm =
        try
            optimize(task)
        catch e
            if isa(e,MosekError)
                return r.rcode,MSK_RES_ERR_UNKNOWN
            else
                rethrow()
            end
        end

    # If this finished with success, inform other tasks to interrupt
    # Note that data races around stop/firstStop are irrelevant
    if ! stop
        stop = true
        firstStop = num
    end

    return MSK_RES_OK,ttrm
end

function optimizeconcurrent(tasks::Vector{Mosek.Task})
    res = [ MSK_RES_ERR_UNKNOWN for t in tasks ]
    trm = [ MSK_RES_ERR_UNKNOWN for t in tasks ]

    # Set a callback function
    for t in tasks
        putcallbackfunc(t, callback)
    end

    # Start parallel optimizations, one per task

    
    Threads.@threads for i in 1:length(tasks)
        (tres,ttrm) = runTask(i,tasks[i])
        res[i] = tres
        trm[i] = ttrm
    end

    # For debugging, print res and trm codes for all optimizers
    for (i,(tres,ttrm)) in enumerate(zip(res,trm))
        println("Optimizer  $i   res $tres   trm $ttrm")
    end

  return firstStop, res, trm
end

#
# idx, winTask, winTrm, winRes = optimizeconcurrent(task, optimizers)
#
# Given a continuous task, set up jobs to optimize it
# with a list of different solvers.
#
# Returns an index, corresponding to the optimization
# task that is returned as winTask. This is the task
# with the best possible status of those that finished.
# If none task is considered successful returns -1.
function optimizeconcurrent(task, optimizers)
    # Choose various optimizers for cloned tasks
    tasks = Mosek.Task[ let t = maketask()
                            putintparam(t,MSK_IPAR_OPTIMIZER, opt)
                            t
                        end for opt in optimizers ]


    # Solve tasks in parallel
    firstOK, res, trm = optimizeconcurrent(tasks)

    if firstOK > 0
        return firstOK, tasks[firstOK], trm[firstOK], res[firstOK]
    else
        return 0, Nothing, Nothing, Nothing
    end
end

#
# idx, winTask, winTrm, winRes = optimizeconcurrent(task, optimizers)
#
# Given a mixed-integer task, set up jobs to optimize it
# with different values of seed. That will lead to
# different execution paths of the optimizer.
#
# Returns an index, corresponding to the optimization
# task that is returned as winTask. This is the task
# with the best value of the objective function.
# If none task is considered successful returns -1.
#
# Typically, the input task would contain a time limit. The two
# major scenarios are:
# 1. Some clone ends before time limit - then it has optimum.
# 2. All clones reach time limit - pick the one with best objective.
function optimizeconcurrentMIO(task, seeds)
    # Choose various seeds for cloned tasks
    tasks = Mosek.Task[ let t = maketask(task)
                            putintparam(MSK_IPAR_MIO_SEED, seed)
                            t
                        end for seed in seeds ]

    # Solve tasks in parallel
    (firstOK, res, trm) = optimizeconcurrent(tasks)

    sense = getobjsense(task)
    bestObj = if sense == MSK_OBJECTIVE_SENSE_MINIMIZE 1.0e+10 else -1.0e+10 end
    bestPos = -1

    if firstOK >= 0
        # Pick the task that ended with res = ok
        # and contains an integer solution with best objective value

        for (i,t) in enumerate(tasks)
            pobj = getprimalobj(t,MSK_SOL_ITG)
            print("$i   $pobj")
        end

        for (i,(tres,ttrm,t)) in enumerate(zip(res,trm,tasks))
            solsta = getsolsta(t,MSK_SOL_ITG)
            if tres == MSK_RES_OK &&
                ( solsta == MSK_SOL_STA_PRIM_FEAS ||
                  solsta == MSK_SOL_STA_INTEGER_OPTIMAL)
                pobj = getprimalobj(t,MSK_SOL_ITG)
                if ( ( sense == MSK_OBJECTIVE_SENSE_MINIMIZE &&
                       getprimalobj(t,MSK_SOL_ITG) < bestObj ) ||
                     ( sense == MSK_OBJECTIVE_SENSE_MAXIMIZE &&
                       getprimalobj(t,MSK_SOL_ITG) > bestObj ) )
                    bestObj = pobj
                    bestPos = i
                end
            end
        end
    end

    if bestPos > 0
        return bestPos, tasks[bestPos], trm[bestPos], res[bestPos]
    else
        return 0, Nothing, Nothing, Nothing
    end
end

# This is an example of how one can use the methods
#       optimizeconcurrent
#       optimizeconcurrentMIO
#
#   argv[0] : name of file with input problem
#   argv[1]: (optional) time limit
function main(fname::String,tlimit)
    maketask() do task
        putstreamfunc(task,MSK_STREAM_LOG,msg -> print(msg))
        if fname != "-"
            readdata(task,fname)
        else
            readptfstring(task,lo1_ptf)
        end

        putobjname(task,"concurrent1")
        # Optional time limit
        if tlimit !== Nothing
            putdouparam(task,MSK_DPAR_OPTIMIZER_MAX_TIME, tlimit)
        end

        let (idx,t,trm,res) = (
            if getnumintvar(task) == 0
                # If the problem is continuous
                # optimize it with three continuous optimizers.
                # (Simplex will fail for non-linear problems)
                optimizers = [ MSK_OPTIMIZER_CONIC,
                               MSK_OPTIMIZER_DUAL_SIMPLEX,
                               MSK_OPTIMIZER_PRIMAL_SIMPLEX ]
                optimizeconcurrent(task, optimizers)
            else
                # Mixed-integer problem.
                # Try various seeds.
                seeds = [ 42, 13, 71749373 ]

                optimizeconcurrentMIO(task, seeds)
            end )

            # Check results and print the best answer
            if idx > 0
                println("Result from optimizer with index $idx:  res $res  trm $trm")
                putstreamfunc(t,MSK_STREAM_LOG,msg -> print(msg))
                optimizersummary(t,MSK_STREAM_LOG)
                solutionsummary(t,MSK_STREAM_LOG)
            else
                println("All optimizers failed.")
            end
        end
    end
end

const lo1_ptf = """Task
Objective
    Maximize + 3 @x0 + @x1 + 5 @x2 + @x3
Constraints
    @c0 [30] + 3 @x0 + @x1 + 2 @x2
    @c1 [15;+inf] + 2 @x0 + @x1 + 3 @x2 + @x3
    @c2 [-inf;25] + 2 @x1 + 3 @x3
Variables
    @x0 [0;+inf]
    @x1 [0;10]
    @x2 [0;+inf]
    @x3 [0;+inf]
"""

let fname = if length(ARGS) < 1 "-" else ARGS[1] end,
    tlimit = if length(ARGS) < 2 Nothing else parse(Float64,ARGS[2]) end

    if false
        main(fname,tlimit)
    else
        println!("Disabled: concurrent1.jl. Example is broken.")
    end
end
