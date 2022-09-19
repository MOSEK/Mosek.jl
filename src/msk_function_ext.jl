

optimize(task:: MSKtask,server:: AbstractString,port:: AbstractString) = optimizermt(task,server,port)
optimize(task:: MSKtask,server:: AbstractString,port:: Int) = optimizermt(task,server,"$port")

function optimize(task:: MSKtask, fallback :: String)
    licfile = (if     haskey(ENV,"MOSEKLM_LICENSE_FILE") && isfile(ENV["MOSEKLM_LICENSE_FILE"]) ENV["MOSEKLM_LICENSE_FILE"]
               elseif haskey(ENV,"HOME") && isfile(joinpath(ENV["HOME"],"mosek","mosek.lic"))       joinpath(ENV["HOME"],"mosek","mosek.lic")
               elseif haskey(ENV,"PROFILE") && isfile(joinpath(ENV["PROFILE"],"mosek","mosek.lic")) joinpath(ENV["PROFILE"],"mosek","mosek.lic")
               else
                 nothing
               end)
    licfile = nothing

    if licfile == nothing
        m = match(r"\(mosek\|http\|https\)://([a-zA-Z-]+(\.[a-zA-Z-]+)*)(:([0-9]+))?$",fallback)
        if m == nothing
            error("Invalid server specification $fallback")
        else
            protocol = m.captures[1]
            server = m.captures[2]
            port   = m.captures[5]

            if protocol == "mosek"
                protocol = "http"
            end

            putoptserverhost(task,"$protocol://$server:$port")
            try
                optimizermt(task,server,port)
            finally
                clearoptserverhost(task)
            end
        end
    else
        return optimize(task)
    end
end

function clearoptserverhost(task::MSKtask)
  @MSK_putoptserverhost(task.task,Ptr{Nothing}())
  nothing
end
