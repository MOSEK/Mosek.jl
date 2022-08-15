#   Copyright : Copyright (c) 2022 MOSEK ApS
#
#   File :      simple.jl
#
#   Purpose :   Demonstrates a very simple example using MOSEK by
#               reading a problem file, solving the problem and
#               writing the problem+solution to a file.

##TAG:begin-code
using Mosek

if length(ARGS) < 1
    println("Syntax: simple FILENAME [ OUTFILE ]")
else
    let filename = ARGS[1],
        outfile = if length(ARGS) > 1 ARGS[2] else Nothing end

        maketask() do task
            putstreamfunc(task,MSK_STREAM_LOG,msg -> print(msg))

##TAG:begin-readdata
            # We assume that a problem file was given as the first command
            # line argument (received in `args')
            readdata(task,filename)
##TAG:end-readdata

            # Solve the problem
            optimize(task)

            # Print a summary of the solution
            solutionsummary(task,MSK_STREAM_LOG)

            # If an output file was specified, save problem to file
            if outfile != Nothing
                # If using OPF format, these parameters will specify what to include in output
                put_int_param(task,MSK_IPAR_OPF_WRITE_SOLUTIONS,  MSK_ON)
                put_int_param(task,MSK_IPAR_OPF_WRITE_PROBLEM,    MSK_ON)
                put_int_param(task,MSK_IPAR_OPF_WRITE_HINTS,      MSK_OFF)
                put_int_param(task,MSK_IPAR_OPF_WRITE_PARAMETERS, MSK_OFF)
                put_int_param(task,MSK_IPAR_PTF_WRITE_SOLUTIONS,  MSK_ON)

##TAG:begin-writedata
                write_data(task,outfile)
##TAG:end-writedata
            end
        end
    end
end
##TAG:end-code
