##
#  Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
#
#  File :      $${file}
#
#  Purpose :   Demonstrates how to use MOSEK OptServer
#              to solve optimization problem asynchronously
##
using Mosek

if length(ARGS) < 2
    println("Missing argument, syntax is:")
    println("  opt_server_async inputfile host:port numpolls [cert]")
else
    token = maketask() do task
        putstreamfunc(task,MSK_STREAM_LOG,msg -> print(msg))

        inputfile  = ARGS[1]
        serveraddr = ARGS[2]
        numpolls   = Int(ARGS[3])

        tlscert    = if length(ARGS) < 4 Nothing else ARGS[4] end

        token = Nothing

        print("reading task from file")
        readdata(task,filename)

        if cert !== Nothing:
            putstrparam(task,MSK_SPAR_REMOTE_TLS_CERT_PATH,cert)
        end

        println("Solve the problem remotely (async)")
        asyncoptimize(task,addr,"")
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
            respavailable, res, trm = asyncpoll(task, addr, "", token)

            println("done!")

            if respavailable:
                println("solution available!")

                respavailable, res, trm = asyncgetresult(task, addr, "", token)

                solutionsummary(task,MSK_STREAM_LOG)
                break
            end
            i = i + 1

            if i == numpolls:
                println("max number of polls reached, stopping host.")
                asyncstop(task,addr,"", token)
            end
        end
    end
end
