__precompile__()
module Mosek
using SparseArrays

if isfile(joinpath(dirname(@__FILE__),"..","deps","deps.jl"))
    include("../deps/deps.jl")
else
    error("Mosek not properly installed. Please run Pkg.build(\"Mosek\")")
end

# The MathOptInterface wrapper is in MosekTools.jl but it can be created with
# `using MosekTools; Mosek.Optimizer()`.
function Optimizer(args...; kwargs...)
    error("To use Mosek with JuMP (or MathOptInterface), you need to use ",
          "the package `MosekTools` (via `using MosekTools`). You may need ",
          "to first install it via `import Pkg; Pkg.add(\"MosekTools\")`.")
end

export
    makeenv, maketask, maketask_ptr,
    MosekError,
    writedatastream

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
struct MosekError <: Exception
    rcode :: Int32
    msg   :: String
end



# Environment: typedef void * Env_t;
mutable struct Env
    env::Ptr{Nothing}
    streamcallbackfunc::Any
    userstreamcallbackfunc:: Any
end

function msk_info_callback_wrapper end

# Task: typedef void * Task_t;
mutable struct Task
    env::Env
    task::Ptr{Nothing}
    borrowed::Bool
    # need to keep a reference to callback funcs for GC
    streamcallbackfunc:: Any
    userstreamcallbackfunc:: Any
    callbackfunc:: Any
    usercallbackfunc:: Any
    nlinfo:: Any

    function Task(env::Env)
        temp = Array{Ptr{Nothing}}(undef,1)
        res = @msk_ccall(maketask, Int32, (Ptr{Nothing}, Int32, Int32, Ptr{Nothing}), env.env, 0, 0, temp)

        if res != MSK_RES_OK.value
            throw(MosekError(res,""))
        end

        task = new(env,temp[1],false,nothing,nothing,nothing,nothing,nothing)
        task.callbackfunc = @cfunction(msk_info_callback_wrapper, Cint, (Ptr{Nothing}, Ptr{Nothing}, Int32, Ptr{Float64}, Ptr{Int32}, Ptr{Int64}))
        task.usercallbackfunc = nothing
        finalizer(deletetask,task)

        r = @msk_ccall(putcallbackfunc, Cint, (Ptr{Nothing}, Ptr{Nothing}, Any), task.task, task.callbackfunc, task)
        if r != MSK_RES_OK.value
            throw(MosekError(r,getlasterrormsg(t)))
        end

        task
    end

    function Task(t::Task)
        temp = Array{Ptr{Nothing}}(undef,1)
        res = @msk_ccall(clonetask, Int32, (Ptr{Nothing}, Ptr{Nothing}), t.task, temp)

        if res != MSK_RES_OK.value
            throw(MosekError(res,""))
        end

        task = new(t.env,temp[1],false,nothing,nothing,nothing,nothing,nothing)
        task.callbackfunc = @cfunction(msk_info_callback_wrapper, Cint, (Ptr{Nothing}, Ptr{Nothing}, Int32, Ptr{Float64}, Ptr{Int32}, Ptr{Int64}))
        task.usercallbackfunc = nothing

        finalizer(deletetask,task)

        r = @msk_ccall(putcallbackfunc, Cint, (Ptr{Nothing}, Ptr{Nothing}, Any), task.task, task.callbackfunc, task)
        if r != MSK_RES_OK.value
            throw(MosekError(r,getlasterrormsg(t)))
        end

        task
    end

    function Task(t::Ptr{Nothing},borrowed::Bool)
        task = new(msk_global_env,t,borrowed,nothing,nothing,nothing,nothing,nothing)
        task.callbackfunc = @cfunction(msk_info_callback_wrapper, Cint, (Ptr{Nothing}, Ptr{Nothing}, Int32, Ptr{Float64}, Ptr{Int32}, Ptr{Int64}))
        task.usercallbackfunc = nothing

        finalizer(deletetask,task)

        r = @msk_ccall(putcallbackfunc, Cint, (Ptr{Nothing}, Ptr{Nothing}, Any), task.task, task.callbackfunc, task)
        if r != MSK_RES_OK.value
            throw(MosekError(r,getlasterrormsg(t)))
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
    makeenv(func::Function)

Create a MOSEK environment, wither create direct for for use with  `do`-syntax.
"""
function makenv end

function makeenv()
    temp = Array{Ptr{Nothing}}(undef, 1)
    res = @msk_ccall(makeenv, Int32, (Ptr{Ptr{Nothing}}, Ptr{UInt8}), temp, C_NULL)
    if res != 0
        # TODO: Actually use result code
        error("MOSEK: Error creating environment")
    end
    Env(temp[1],nothing, nothing)
end

function makeenv(func::Function)
    temp = Array{Ptr{Nothing}}(undef,1)
    res = @msk_ccall(makeenv, Int32, (Ptr{Ptr{Nothing}}, Ptr{UInt8}), temp, C_NULL)
    if res != 0
        # TODO: Actually use result code
        error("MOSEK: Error creating environment")
    end
    env = Env(temp[1],nothing, nothing)

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
    maketask(;env::Env = msk_global_env, filename::String = "")
    
Create a task. If `filename` is not 0-length, initialize the task from this file.

    maketask(func :: Function;env::Env = msk_global_env, filename::String = "")

Create a task. If `filename` is not 0-length, initialize the task from
this file. The `func` parameter is a function 

```julia
func(task :: Task) :: Any
```

This can be used with the `do`-syntax:

```julia
maketask(env) do task
    readdata(task,"MyFile.task")
end
```
    maketask(task::Task)

Create a clone of `task`.
"""
function maketask end
function maketask(;env::Env = msk_global_env, filename::String = "")
    t = Task(env)
    if length(filename) > 0
        try
            readdata(t,filename)
        catch e
            deletetask(t)
            rethrow()
        end
    end
    t
end

function maketask(func :: Function;env::Env = msk_global_env, filename::String = "")
    t = Task(env)
    try
        if length(filename) > 0 readdata(t,filename) end
        func(t)
    finally
        deletetask(t)
    end
end

maketask(task::Task) = Task(task)

function maketask_ptr(t::Ptr{Nothing},borrowed::Bool)
    Task(t,borrowed)
end


"""
    deletetask(t::Task)

Destroy the task object.
"""
function deletetask(t::Task)
    if t.task != C_NULL
        if ! t.borrowed
            temp = Array{Ptr{Nothing}}(undef,1)
            temp[1] = t.task
            @msk_ccall(deletetask,Int32,(Ptr{Ptr{Nothing}},), temp)
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
        temp = Array{Ptr{Nothing}}(undef,1)
        temp[1] = e.env
        @msk_ccall(deleteenv,Int32,(Ptr{Ptr{Nothing}},), temp)
        e.env = C_NULL
    end
end

# function getlasterror(t::Task)
#     lasterrcode = Array{Cint}(undef,1)
#     lastmsglen = Array{Cint}(undef,1)

#     @msk_ccall(getlasterror,Cint,(Ptr{Nothing},Ptr{Cint},Cint,Ptr{Cint},Ptr{UInt8}),
#                t.task, lasterrcode, 0, lastmsglen, C_NULL)
#     lastmsg = Array{UInt8}(undef,lastmsglen[1])
#     @msk_ccall(getlasterror,Cint,(Ptr{Nothing},Ptr{Cint},Cint,Ptr{Cint},Ptr{UInt8}),
#                t.task, lasterrcode, lastmsglen[1], lastmsglen, lastmsg)
#     String(lastmsg[1:lastmsglen[1]-1])
# end

function getlasterrormsg(task::Task)
  lastrescode = Ref{Int32}()
  sizelastmsg = Ref{Int64}()
  if 0 != disable_sigint(()->ccall((:MSK_getlasterror64,libmosek),Int32,(Ptr{Nothing},Ref{Int32},Int64,Ref{Int64},Ptr{UInt8},),task.task,lastrescode,sizelastmsg,0,C_NULL))
    ""
  else
    lastmsg = Array{UInt8}(undef,sizelastmsg[]+1)
    if 0 != disable_sigint(()->ccall((:MSK_getlasterror64,libmosek),Int32,(Ptr{Nothing},Ref{Int32},Int64,Ref{Int64},Ptr{UInt8},),task.task,lastrescode,sizelastmsg,length(lastmasg),lastmsg))
      ""
    else
      lastmsg_len = findfirst(_c->_c==0,lastmsg)
      if lastmsg_len === nothing 
        String(lastmsg) 
      else 
        String(lastmsg[1:lastmsg_len-1]) 
      end
    end
  end   
end
function getlasterrormsg(task::Ptr{Nothing})
  lastrescode = Ref{Int32}()
  sizelastmsg = Ref{Int64}()
  if 0 != disable_sigint(()->ccall((:MSK_getlasterror64,libmosek),Int32,(Ptr{Nothing},Ref{Int32},Int64,Ref{Int64},Ptr{UInt8},),task,lastrescode,sizelastmsg,0,C_NULL))
    ""
  else
    lastmsg = Array{UInt8}(undef,sizelastmsg[]+1)
    if 0 != disable_sigint(()->ccall((:MSK_getlasterror64,libmosek),Int32,(Ptr{Nothing},Ref{Int32},Int64,Ref{Int64},Ptr{UInt8},),task,lastrescode,sizelastmsg,length(lastmasg),lastmsg))
      ""
    else
      lastmsg_len = findfirst(_c->_c==0,lastmsg)
      if lastmsg_len === nothing 
        String(lastmsg) 
      else 
        String(lastmsg[1:lastmsg_len-1]) 
      end
    end
  end   
end


using SparseArrays

include("msk_enums.jl")
include("msk_functions.jl")
include("msk_function_ext.jl")
include("msk_callback.jl")

## General Convex optimizer has been discontinued
#include("msk_geco.jl")

include("show.jl")
include("ext_functions.jl")

function getlasterrormsg(t::Ptr{Nothing})
    lastrescode_ = Ref{Int32}()
    __tmp_662 = Ref{Int64}()
    @MSK_getlasterror64(t,Ref{Int32}(),0,__tmp_662,C_NULL)
    __tmp_661 = __tmp_662[]
    sizelastmsg = Int64((__tmp_661 + Int64(1)))
    lastmsglen_ = Ref{Int64}()
    lastmsg_ = Array{UInt8}(undef,sizelastmsg)
    @MSK_getlasterror64(t,lastrescode_,sizelastmsg,lastmsglen_,lastmsg_)
    lastmsg_len = findfirst(_c->_c==0,lastmsg_)

    if lastmsg_len === nothing
        String(lastmsg_)
    else
        String(lastmsg_[1:lastmsg_len-1])
    end
end


#import MathProgBase
#struct MosekSolver <: MathProgBase.AbstractMathProgSolver
#    options
#end
#MosekSolver(;kwargs...) = MosekSolver(kwargs)
#export MosekSolver

#include("MosekSolverInterface.jl")


function writedatastreamcb(handle :: Ptr{Nothing},
                           src    :: Ptr{UInt8},
                           count  :: UInt64)
    let io = unsafe_pointer_to_objref(handle) :: IO
        write(io,unsafe_wrap(Vector{UInt8},src,count))
    end
    count
end

function writedatastream(task :: MSKtask, format :: Dataformat, compress :: Compresstype, stream :: IO)
    let handle = stream,
        cbfunc = @cfunction(writedatastreamcb, UInt64, (Ptr{Nothing},Ptr{UInt8},UInt64))
        @MSK_writedatahandle(task.task,cbfunc,handle,format.value,compress.value)
    end
end





end
