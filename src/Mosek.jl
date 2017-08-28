__precompile__()
module Mosek

if isfile(joinpath(dirname(@__FILE__),"..","deps","deps.jl"))
    include("../deps/deps.jl")
else
    error("Mosek not properly installed. Please run Pkg.build(\"Mosek\")")
end

export
    makeenv, maketask, maketask_ptr,
    MosekError

# A macro to make calling C API a little cleaner
macro msk_ccall(func, args...)
    f = Base.Meta.quot(Symbol("MSK_$(func)"))
    args = [esc(a) for a in args]
    quote
        ccall(($f,libmosek), $(args...))
    end
end

# -----
# Types
# -----
type MosekError <: Exception
    rcode :: Int32
    msg   :: String
end



# Environment: typedef void * Env_t;
type Env
    env::Ptr{Void}
    streamcallbackfunc::Any
end


# Task: typedef void * Task_t;
type Task
    env::Env
    task::Ptr{Void}
    borrowed::Bool
    # need to keep a reference to callback funcs for GC
    streamcallbackfunc:: Any
    userstreamcallbackfunc:: Any
    callbackfunc:: Any
    usercallbackfunc:: Any
    nlinfo:: Any

    function Task(env::Env)
        temp = Array{Ptr{Void}}(1)
        res = @msk_ccall(maketask, Int32, (Ptr{Void}, Int32, Int32, Ptr{Void}), env.env, 0, 0, temp)

        if res != MSK_RES_OK.value
            throw(MosekError(res,""))
        end

        task = new(env,temp[1],false,nothing,nothing,nothing,nothing,nothing)
        task.callbackfunc = cfunction(msk_info_callback_wrapper, Cint, (Ptr{Void}, Ptr{Void}, Int32, Ptr{Float64}, Ptr{Int32}, Ptr{Int64}))
        task.usercallbackfunc = nothing
        finalizer(task,deletetask)

        r = @msk_ccall(putcallbackfunc, Cint, (Ptr{Void}, Ptr{Void}, Any), task.task, task.callbackfunc, task)
        if r != MSK_RES_OK.value
            throw(MosekError(r,getlasterror(t)))
        end

        task
    end

    function Task(t::Task)
        temp = Array{Ptr{Void}}(1)
        res = @msk_ccall(clonetask, Int32, (Ptr{Void}, Ptr{Void}), t.task, temp)

        if res != MSK_RES_OK.value
            throw(MosekError(res,""))
        end

        task = new(env,temp[1],false,nothing,nothing,nothing,nothing,nothing)
        task.callbackfunc = cfunction(msk_info_callback_wrapper, Cint, (Ptr{Void}, Ptr{Void}, Int32, Ptr{Float64}, Ptr{Int32}, Ptr{Int64}))
        task.usercallbackfunc = nothing

        finalizer(task,deletetask)

        r = @msk_ccall(putcallbackfunc, Cint, (Ptr{Void}, Ptr{Void}, Any), task.task, task.callbackfunc, task)
        if r != MSK_RES_OK.value
            throw(MosekError(r,getlasterror(t)))
        end

        task
    end

    function Task(t::Ptr{Void},borrowed::Bool)
        task = new(msk_global_env,t,borrowed,nothing,nothing,nothing,nothing,nothing)
        task.callbackfunc = cfunction(msk_info_callback_wrapper, Cint, (Ptr{Void}, Ptr{Void}, Int32, Ptr{Float64}, Ptr{Int32}, Ptr{Int64}))
        task.usercallbackfunc = nothing

        finalizer(task,deletetask)

        r = @msk_ccall(putcallbackfunc, Cint, (Ptr{Void}, Ptr{Void}, Any), task.task, task.callbackfunc, task)
        if r != MSK_RES_OK.value
            throw(MosekError(r,getlasterror(t)))
        end

        task
    end
end
const MSKtask = Task
const MSKenv  = Env

# ------------
# API wrappers
# ------------
# TODO: Support other argument
"""
    makeenv()

Create a MOSEK environment.
"""
function makeenv()
    temp = Array{Ptr{Void}}(1)
    res = @msk_ccall(makeenv, Int32, (Ptr{Ptr{Void}}, Ptr{UInt8}), temp, C_NULL)
    if res != 0
        # TODO: Actually use result code
        error("MOSEK: Error creating environment")
    end
    Env(temp[1],nothing)
end

"""
    makeenv(func::Function)

Create a MOSEK environment for use with `do`-syntax.
"""
function makeenv(func::Function)
    temp = Array{Ptr{Void}}(1)
    res = @msk_ccall(makeenv, Int32, (Ptr{Ptr{Void}}, Ptr{UInt8}), temp, C_NULL)
    if res != 0
        # TODO: Actually use result code
        error("MOSEK: Error creating environment")
    end
    env = Env(temp[1],nothing)

    try
        func(env)
    finally
        deleteenv(env)
    end
end

# Note on initialization of msk_global_env:
#
#  When loading Mosek from source this works fine, but when loading
#  precompiled module, makeenv() is not called (and some garbage
#  value is put in msk_global_env). It appears that static
#  initializers must be called from __init__(). That is called a bad solution here:
#    https://github.com/JuliaLang/julia/issues/12010
msk_global_env = makeenv() :: Env
__init__() = (global msk_global_env = makeenv())


"""
    maketask(env:Env)
    maketask()

Create a task that belongs to either the given `env` or to the global `env`.

    maketask(task::Task)

Create a clone of `task`.

    maketask(func::Function, env::Env)
    maketask(func::Function)

Create a task that either longs to the given `env` or to the global `env`. The `func` parameter is a function

```julia
func(task :: Task) :: Any
```

This can be used with the `do`-syntax:

```julia
maketask(env) do task
    readdata(task,"MyFile.task")
end
```
"""
function maketask end
maketask(env::Env) = Task(env)
maketask(task::Task) = Task(task)
maketask() = Task(msk_global_env)
function maketask(func::Function, env::Env)
    t = Task(env)
    try
        func(t)
    finally
        deletetask(t)
    end
end
function maketask(func::Function)
    t = Task(msk_global_env)
    try
        func(t)
    finally
        deletetask(t)
    end
end


function maketask_ptr(t::Ptr{Void},borrowed::Bool)
    Task(t,borrowed)
end


"""
    deletetask(t::Task)

Destroy the task object.
"""
function deletetask(t::Task)
    if t.task != C_NULL
        if ! t.borrowed
            temp = Array{Ptr{Void}}(1)
            temp[1] = t.task
            @msk_ccall(deletetask,Int32,(Ptr{Ptr{Void}},), temp)
        end
        t.task = C_NULL
    end
end

"""
    deleteenv(t::Env)

Destroy the task object.
"""
function deleteenv(e::Env)
    if e.env != C_NULL
        temp = Array{Ptr{Void}}(1)
        temp[1] = t.env
        @msk_ccall(deleteenv,Int32,(Ptr{Ptr{Void}},), temp)
        e.env = C_NULL
    end
end

function getlasterror(t::Task)
    lasterrcode = Array{Cint}(1)
    lastmsglen = Array{Cint}(1)

    @msk_ccall(getlasterror,Cint,(Ptr{Void},Ptr{Cint},Cint,Ptr{Cint},Ptr{UInt8}),
               t.task, lasterrcode, 0, lastmsglen, C_NULL)
    lastmsg = Array{UInt8}(lastmsglen[1])
    @msk_ccall(getlasterror,Cint,(Ptr{Void},Ptr{Cint},Cint,Ptr{Cint},Ptr{UInt8}),
               t.task, lasterrcode, lastmsglen[1], lastmsglen, lastmsg)
    convert(String,lastmsg[1:lastmsglen[1]-1])
end

include("msk_enums.jl")
include("msk_functions.jl")
include("msk_callback.jl")
include("msk_geco.jl")


end
