export
  putstreamfunc,
  putcallbackfunc

function msk_stream_callback_wrapper(userdata::Ptr{Nothing}, msg :: Ptr{UInt8})
  f = unsafe_pointer_to_objref(userdata) :: Function
  f(unsafe_string(msg))
  convert(Int32,0)::Int32
end

"""
    putstreamfunc(t::MSKtask, whichstream:: Streamtype, f :: Function)

Attach a log printing callback function to a task stream.

- `whichstream` Defines the stream to attach to. `MSK_STREAM_LOG` will catch the combined output from error and message streams.
- `f` is a function `function(msg::String)`
"""
function putstreamfunc(t::MSKtask, whichstream:: Streamtype, f :: Function)  
  cbfunc = @cfunction(msk_stream_callback_wrapper, Int32, (Ptr{Nothing},Ptr{UInt8}))
  r = @msk_ccall(linkfunctotaskstream,Int32,(Ptr{Nothing},Int32, Any, Ptr{Nothing}),t.task,whichstream.value,f,cbfunc)
  if r != MSK_RES_OK.value
    throw(MosekError(r,getlasterror(t)))
  end
  t.streamcallbackfunc = cbfunc
  t.userstreamcallbackfunc = f
  nothing
end


# MSKcallbackfunc :: ( MSKtask_t, void *, Cint, double *, int32 *, int64 * -> int32 )
# The callback function is called with following parameers:
#  - a 'where' identifying where in the solver we currently are
#  - a 'dinf' array carrying double information on the current state of the solution/solver
#  - a 'iinf' array carrying int32 information on the current state of the solution/solver
#  - a 'linf' array carrying int64 information on the current state of the solution/solver
function msk_info_callback_wrapper(t::Ptr{Nothing}, userdata::Ptr{Nothing}, where :: Int32, douinf :: Ptr{Float64}, intinf :: Ptr{Int32}, lintinf :: Ptr{Int64})
    r = if unsafe_load(cglobal(:jl_signal_pending, Cint),1) > Cint(0) # Not so good... jl_signal_pending is actually a sig_atomic_t...
            1
        else
            task = unsafe_pointer_to_objref(userdata) :: MSKtask
            if task.usercallbackfunc == nothing
                0
            else
                dinfa  = unsafe_wrap(Vector{Float64},douinf,(length(Dinfitem),),own=false)
                iinfa  = unsafe_wrap(Vector{Int32},intinf,(length(Iinfitem),),own=false)
                liinfa = unsafe_wrap(Vector{Int64},lintinf,(length(Liinfitem),),own=false)

                task.usercallbackfunc(Callbackcode(where), dinfa, iinfa, liinfa)
            end
        end
    
    convert(Cint,r)::Cint
end

# f :: where :: Cint, dinf :: Array{Float64,1}, iinf :: Array{Int32,1}, linf :: Array{Int64,1} -> Int32
# NOTE: On Win32 the callback function should be stdcall
"""
    putcallbackfunc(t::MSKtask, f::Function)

Set a callback. This will be called at various points during optimization. The `f` argument is a function of the form
```julia
    function(where::Callbackcode, dinf::Vector{Float64}, iinf::Vector{Int32}, liinf::Vector{Int64})
```
There are very strict rules about what you can call from the callback
function. Basically: Do not call anything using the same task as the
callback was called from. The one exception is `getxx()` that can be called when solving a mixed integer problem.

"""
function putcallbackfunc(t::MSKtask, f::Function)
    t.usercallbackfunc = f
    nothing
end

function clearcallbackfunc(t::MSKtask)
    t.usercallbackfunc = nothing
    nothing
end

