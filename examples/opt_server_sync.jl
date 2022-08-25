##
#  Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
#
#  File :      opt_server_sync.jl
#
#  Purpose :   Demonstrates how to use MOSEK OptServer
#              to solve optimization problem synchronously
##

using Mosek

if length(ARGS) < 2
    println("Missing argument, syntax is:")
    println("  opt_server_sync inputfile addr [certpath]")
else
    maketask() do task
        putstreamfunc(task,MSK_STREAM_LOG,msg -> print(msg))

        inputfile  = ARGS[1]
        serveraddr = ARGS[2]
        tlscert    = if length(ARGS) < 3 Nothing else ARGS[3] end

        # We assume that a problem file was given as the first command
        # line argument (received in `argv')
        readdata(task,inputfile)

        # Set OptServer URL
        putoptserverhost(task,serveraddr)

        # Path to certificate, if any
        if tlscert !== Nothing
            putstrparam(task,MSK_SPAR_REMOTE_TLS_CERT_PATH, tlscert)

            # Solve the problem remotely, no access token
            trm = optimize(task)

            # Print a summary of the solution
            solutionsummary(task,MSK_STREAM_LOG)
        end
    end
end
