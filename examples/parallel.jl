
#   Copyright: Copyright (c) MOSEK ApS, Denmark. All rights reserved.
#
#   File:      parallel.jl
#
#   Purpose: Demonstrates parallel optimization using optimizebatch()


using Mosek

# Example of how to use env.optimizebatch().
# Optimizes tasks whose names were read from command line.
if length(ARGS) < 2
    println("Usage: parallel FILENAME FILENAME [ FILENAME ... ]")
else
    n = length(ARGS)
    makeenv() do env
        maketask() do task
            tasks = [ maketask(filename=f) for f in ARGS ]

            # Size of thread pool available for all tasks
            threadpoolsize = 6

            for t in tasks
                putintparam(t,MSK_IPAR_NUM_THREADS, 2)
            end

            # Optimize all the given tasks in parallel
            (trm,res) = optimizebatch(env,
                                      false,          # No race
                                      -1.0,           # No time limit
                                      threadpoolsize,
                                      tasks)          # Array of tasks to optimize

            for (i,t) in enumerate(tasks)
                println("Task  $i  res $(res[i])   trm $(trm[i])   obj_val  $(getdouinf(t,MSK_DINF_INTPNT_PRIMAL_OBJ))  time $(getdouinf(t,MSK_DINF_OPTIMIZER_TIME))")
            end
        end
    end
end
