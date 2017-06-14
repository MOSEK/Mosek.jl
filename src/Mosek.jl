__precompile__()
module Mosek
if isfile(joinpath(dirname(@__FILE__),"..","deps","deps.jl"))
    include("../deps/deps.jl")
else
    error("Mosek not properly installed. Please run Pkg.build(\"Mosek\")")
end

import MathProgBase


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



# Environment: typedef void * MSKenv_t;
type MSKenv
    env::Ptr{Void}
    streamcallbackfunc::Any
end


# Task: typedef void * MSKtask_t;
type MSKtask
    env::MSKenv
    task::Ptr{Void}
    borrowed::Bool
    # need to keep a reference to callback funcs for GC
    streamcallbackfunc:: Any
    userstreamcallbackfunc:: Any
    callbackfunc:: Any
    usercallbackfunc:: Any
    nlinfo:: Any

    function MSKtask(env::MSKenv)
        temp = Array{Ptr{Void}}(1)
        res = @msk_ccall(maketask, Int32, (Ptr{Void}, Int32, Int32, Ptr{Void}), env.env, 0, 0, temp)

        if res != MSK_RES_OK
            throw(MosekError(res,""))
        end

        task = new(env,temp[1],false,nothing,nothing,nothing,nothing,nothing)
        task.callbackfunc = cfunction(msk_info_callback_wrapper, Cint, (Ptr{Void}, Ptr{Void}, Int32, Ptr{Float64}, Ptr{Int32}, Ptr{Int64}))
        task.usercallbackfunc = nothing
        finalizer(task,deletetask)

        r = @msk_ccall(putcallbackfunc, Cint, (Ptr{Void}, Ptr{Void}, Any), task.task, task.callbackfunc, task)
        if r != MSK_RES_OK
            throw(MosekError(r,getlasterror(t)))
        end

        task
    end

    function MSKtask(t::MSKtask)
        temp = Array{Ptr{Void}}(1)
        res = @msk_ccall(clonetask, Int32, (Ptr{Void}, Ptr{Void}), t.task, temp)

        if res != MSK_RES_OK
            throw(MosekError(res,""))
        end

        task = new(env,temp[1],false,nothing,nothing,nothing,nothing,nothing)
        task.callbackfunc = cfunction(msk_info_callback_wrapper, Cint, (Ptr{Void}, Ptr{Void}, Int32, Ptr{Float64}, Ptr{Int32}, Ptr{Int64}))
        task.usercallbackfunc = nothing

        finalizer(task,deletetask)

        r = @msk_ccall(putcallbackfunc, Cint, (Ptr{Void}, Ptr{Void}, Any), task.task, task.callbackfunc, task)
        if r != MSK_RES_OK
            throw(MosekError(r,getlasterror(t)))
        end

        task
    end

    function MSKtask(t::Ptr{Void},borrowed::Bool)
        task = new(msk_global_env,t,borrowed,nothing,nothing,nothing,nothing,nothing)
        task.callbackfunc = cfunction(msk_info_callback_wrapper, Cint, (Ptr{Void}, Ptr{Void}, Int32, Ptr{Float64}, Ptr{Int32}, Ptr{Int64}))
        task.usercallbackfunc = nothing

        finalizer(task,deletetask)

        r = @msk_ccall(putcallbackfunc, Cint, (Ptr{Void}, Ptr{Void}, Any), task.task, task.callbackfunc, task)
        if r != MSK_RES_OK
            throw(MosekError(r,getlasterror(t)))
        end

        task
    end
end


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
    MSKenv(temp[1],nothing)
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
    env = MSKenv(temp[1],nothing)

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
msk_global_env = makeenv() :: MSKenv
__init__() = (global msk_global_env = makeenv())

"""
    mosektask(env:MSKenv)

Create a MOSEK task belonging to `env`.
"""
function maketask(env::MSKenv)

    MSKtask(env)
end

"""
    mosektask(task::MSKtask)

Create a clone of a MOSEK task.
"""
function maketask(task::MSKtask)
    MSKtask(task)
end

"""
    mosektask()

Create a new MOSEK task belonging to the global environment.
"""
function maketask()
    MSKtask(msk_global_env)
end

function maketask_ptr(t::Ptr{Void},borrowed::Bool)
    MSKtask(t,borrowed)
end

"""
    maketask(fnuc::Function, env::MSKenv)

Create a new MOSEK task for use with `do`-notation.
"""
function maketask(func::Function, env::MSKenv)
    t = MSKtask(env)
    try
        func(t)
    finally
        deletetask(t)
    end
end

"""
    maketask(func::Function, env::MSKenv)

Create a new MOSEK task belonging to the global environment for use
with `do`-notation.
"""
function maketask(func::Function)
    t = MSKtask(msk_global_env)
    try
        func(t)
    finally
        deletetask(t)
    end
end

"""
    deletetask(t::MSKtask)

Destroy the task object.
"""
function deletetask(t::MSKtask)
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
    deleteenv(t::MSKenv)

Destroy the task object.
"""
function deleteenv(e::MSKenv)
    if e.env != C_NULL
        temp = Array{Ptr{Void}}(1)
        temp[1] = t.env
        @msk_ccall(deleteenv,Int32,(Ptr{Ptr{Void}},), temp)
        e.env = C_NULL
    end
end

function getlasterror(t::MSKtask)
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

include("MosekSolverInterface.jl")
using Mosek.MosekMathProgSolverInterface
export MosekSolver


end
