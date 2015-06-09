export
  putstreamfunc,
  putcallbackfunc

function msk_stream_callback_wrapper(userdata::Ptr{Void}, msg :: Ptr{Uint8})
  f = unsafe_pointer_to_objref(userdata) :: Function
  f (bytestring(msg))
  convert(Int32,0)::Int32
end

function putstreamfunc(t::MSKtask, whichstream:: Int32, f :: Function)
  cbfunc = cfunction(msk_stream_callback_wrapper, Int32, (Ptr{Void},Ptr{Uint8}))
  r = @msk_ccall(linkfunctotaskstream,Int32,(Ptr{Void},Int32, Any, Ptr{Void}),t.task,whichstream,f,cbfunc)
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
function msk_info_callback_wrapper(t::Ptr{Void}, userdata::Ptr{Void}, where :: Int32, douinf :: Ptr{Float64}, intinf :: Ptr{Int32}, lintinf :: Ptr{Int64})
  f      = unsafe_pointer_to_objref(userdata) :: Function
  dinfa  = pointer_to_array(douinf,(MSK_DINF_END,),false)
  iinfa  = pointer_to_array(intinf,(MSK_IINF_END,),false)
  liinfa = pointer_to_array(lintinf,(MSK_LIINF_END,),false)

  r = f(@compat(Int32(where)), dinfa, iinfa, liinfa)
  convert(Int32,r)::Int32
end

function msk_callback_wrapper(t::Ptr{Void}, userdata::Ptr{Void}, where :: Int32)
  f      = unsafe_pointer_to_objref(userdata) :: Function
  r = f(convert(Int32,where))
  convert(Int32,r)::Int32
end

# f :: where :: Cint, dinf :: Array{Float64,1}, iinf :: Array{Int32,1}, linf :: Array{Int64,1} -> Int32
# NOTE: On Win32 the callback function should be stdcall
function putcallbackfunc(t::MSKtask, f::Function)
  cbfunc = cfunction(msk_info_callback_wrapper, Int32, (Ptr{Void}, Ptr{Void}, Int32, Ptr{Float64}, Ptr{Int32}, Ptr{Int64}))
  #cbfunc = cfunction(msk_callback_wrapper, Int32, (Ptr{Void}, Ptr{Void}, Int32))

  r = @msk_ccall(putcallbackfunc, Int32, (Ptr{Void}, Ptr{Void}, Any), t.task, cbfunc, f)
  if r != MSK_RES_OK
    throw(MosekError(r,getlasterror(t)))
  end

  t.callbackfunc = cbfunc
  t.usercallbackfunc = f
  nothing
end

