##
#  Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
#
#  File :      opt_server_async.jl
#
#  Purpose :   Demonstrates how to use MOSEK OptServer
#              to solve optimization problem asynchronously
##
using Mosek

if length(ARGS) < 2
    println("Missing argument, syntax is:")
    println("  opt_server_async inputfile host:port numpolls [cert]")
else
    filename   = ARGS[1]
    serveraddr = ARGS[2]
    numpolls   = parse(Int,ARGS[3])
    cert       = (if length(ARGS) > 3
                      ARGS[4]
                  else
                      Nothing
                  end)

    token = maketask() do task
        putstreamfunc(task,MSK_STREAM_LOG,msg -> print(msg))

        token = Nothing

        print("reading task from file")
        readdata(task,filename)

        if cert !== Nothing
            putstrparam(task,MSK_SPAR_REMOTE_TLS_CERT_PATH,cert)
        end

        println("Solve the problem remotely (async)")
        asyncoptimize(task,serveraddr,"")
    end

    println("Task token: '$token'")

    maketask() do task
        putstreamfunc(task,MSK_STREAM_LOG,msg -> print(msg))
        readdata(task,filename)

        if cert !== Nothing
            putstrparam(task,MSK_SPAR_REMOTE_TLS_CERT_PATH,cert)
        end


        i = 0
        while i < numpolls
            sleep(0.1)

            println("poll $i...")
            respavailable, res, trm = asyncpoll(task, serveraddr, "", token)

            println("done!")

            if respavailable
                println("solution available!")
                if res!=MSK_RES_OK
                    println("Wrong response code from remote server: expected OK, got $res")
                    @assert false
                end

                respavailable, res, trm = asyncgetresult(task, serveraddr, "", token)

                solutionsummary(task,MSK_STREAM_LOG)
                break
            end
            i = i + 1

            if i == numpolls
                println("max number of polls reached, stopping host.")
                asyncstop(task,serveraddr,"", token)
            end
        end
    end
end
