export
  putstreamfunc,
  putcallbackfunc

function msk_stream_callback_wrapper(userdata::Ptr{Nothing}, msg :: Ptr{UInt8})
  f = unsafe_pointer_to_objref(userdata) :: Function
  f(unsafe_string(msg))
  convert(Int32,0)::Int32
end

function putstreamfunc(t::MSKtask, whichstream:: Int32, f :: Function)  
  cbfunc = cfunction(msk_stream_callback_wrapper, Int32, (Ptr{Nothing},Ptr{UInt8}))
  r = @msk_ccall(linkfunctotaskstream,Int32,(Ptr{Nothing},Int32, Any, Ptr{Nothing}),t.task,whichstream,f,cbfunc)
  if r != MSK_RES_OK
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
                dinfa  = unsafe_wrap(Array{Float64,1},douinf,(MSK_DINF_END,),false)
                iinfa  = unsafe_wrap(Array{Int32,1},intinf,(MSK_IINF_END,),false)
                liinfa = unsafe_wrap(Array{Int64,1},lintinf,(MSK_LIINF_END,),false)

                task.usercallbackfunc(Int32(where), dinfa, iinfa, liinfa)
            end
        end
    
    convert(Cint,r)::Cint
end

# f :: where :: Cint, dinf :: Array{Float64,1}, iinf :: Array{Int32,1}, linf :: Array{Int64,1} -> Int32
# NOTE: On Win32 the callback function should be stdcall
function putcallbackfunc(t::MSKtask, f::Function)
    t.usercallbackfunc = f
    nothing
end

function clearcallbackfunc(t::MSKtask)
    t.usercallbackfunc = nothing
    nothing
end

