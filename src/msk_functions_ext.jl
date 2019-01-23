optimize(task:: MSKtask,server:: AbstractString,port:: AbstractString) = optimizermt(task,server,port)
optimize(task:: MSKtask,server:: AbstractString,port:: Int) = optimizermt(task,server,"$port")

function optimize(task:: MSKtask, fallback :: String)
    licfile = (if     haskey(ENV,"MOSEKLM_LICENSE_FILE") && isfile(ENV["MOSEKLM_LICENSE_FILE"]) ENV["MOSEKLM_LICENSE_FILE"]
               elseif haskey(ENV,"HOME") && isfile(joinpath(ENV["HOME"],"mosek","mosek.lic"))       joinpath(ENV["HOME"],"mosek","mosek.lic")
               elseif haskey(ENV,"PROFILE") && isfile(joinpath(ENV["PROFILE"],"mosek","mosek.lic")) joinpath(ENV["PROFILE"],"mosek","mosek.lic")
               else
                 nothing
               end)

    if licfile == nothing
        m = match(r"mosek://([a-zA-Z-]+(\.[a-zA-Z-]+)*)(:([0-9]+))?$",fallback)
        if m == nothing
            error("Invalid server specification $fallback")
        else
            server = m.captures[1]
            port   = if m.captures[4] == nothing "30080" else m.captures[4] end
            optimizermt(task,server,port)
        end
    else
        return optimize(task)
    end
end
