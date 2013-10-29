  # MSKcallbackfunc :: ( MSKtask_t, void *, Cint, double *, int32 *, int64 * -> int32 )  
  # The callback function is called with following parameers:
  #  - a 'where' identifying where in the solver we currently are
  #  - a 'dinf' array carrying double information on the current state of the solution/solver
  #  - a 'iinf' array carrying int32 information on the current state of the solution/solver
  #  - a 'linf' array carrying int64 information on the current state of the solution/solver
  function msk_info_callback_wrapper(t::Ptr{Void}, userdata::Ptr{Void}, where :: Cint, douinf :: Ptr{Float64}, intinf :: Ptr{Int32}, lintinf :: Ptr{Int64})
    f      = unsafe_pointer_to_objref(userdata) :: Function
    dinfa  = pointer_to_array(douinf,(MSK_DINF_END,),false)
    iinfa  = pointer_to_array(intinf,(MSK_IINF_END,),false)
    liinfa = pointer_to_array(lintinf,(MSK_LIINF_END,),false)

    r = f(int(where), dinfa, iinfa, liinfa) 
    int(r)
  end

  # f :: where :: Cint, dinf :: Array{Float64,1}, iinf :: Array{Int32,1}, linf :: Array{Int64,1} -> Int32
  # NOTE: On Win32 the callback function should be stdcall
  function msk_putcallbackfunc(t::MSKtask, f::Function)
    mskcallback = cfunction(msk_info_callback_wrapper, Int32, (Ptr{Void}, Ptr{Void}, Cint, Ptr{Float64}, Ptr{Int32}, Ptr{Int64}))
    r = @msk_ccall(putcallbackfunc, Cint, (Ptr{Void}, Ptr{Void}, Any), task, mskcallback, f)
    @msk_handle_res(task,r)
    t.callback = mskcallback
    nothing
  end

  function msk_clearcallbackfuncI(t::MSKtask)
    r = @msk_ccall(putcallbackfnuc, Cint, (Ptr{Void}, Ptr{Void}, Ptr{Void}), task, C_NULL, C_NULL)
    @msk_handle_res(task,r)
    t.callback = mskcallback
    nothing
  end 




 
#  function msk_putnlfunc(t::MSKtask, nlgetsp :: Function, nlgetva :: Function)
#    msknlgetsp = cfunction(msk_nlgetsp_wrapper, Int32, (Ptr{Void}, Ptr{Int32}, Ptr{Int32}, Int32, Ptr{Bool}, Ptr{Int32}, Ptr{Int32}, Int32, Int32, Ptr{Int32}, Int32, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}))
#    msknlgetva = cfunction(msk_nlgetva_wrapper, Int32, (Ptr{Void}, Ptr{Float64}, Float64, Ptr{Float64}, Ptr{Float64}, Ptr{Int32}, Ptr{Int32}, Ptr{Float64}, Int32, Ptr{Int32}, Ptr{Float64}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Float}, Ptr{Float}, Int32, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Float64})
#    r = @msk_ccall(setcallbackfunc, Cint, (Ptr{Void}, Ptr{Void}, Any), model.ptr_model, grbcallback, (callback,model))
#    if ret != 0
#        throw(GurobiError(model.env, ret))
#    end
#    # we need to keep a reference to the callback function
#    # so that it isn't garbage collected
#    model.callback = callback
#    nothing 
#  end

# The nlgetva is in fact several functions in one, computing constraint 
# values, gradient of lagrangian and hessian. Furthermore, it has the quirk that
# the lengths of arrays cannot be validated when it is called. It is a mess.
# What we do here: We split it into: 
#  msk_nl_eval :: xx -> objval, conval
#  msk_nl_grdlag :: xx,yo,yc,subi -> grdlag
#  msk_nl_heslag :: xx,yo,yc,subi -> haslagsubi, heslagsubj, heslagval
#  msk_nl_grdobj :: xx -> grdobj

# nlgetsp is the same; it is several function in one. 
#  msk_nl_grdobjnnz    :: () -> nnz :: Int32 
#  msk_nl_grdobjsp     :: () -> grdobjsubj
#  msk_nl_coniisnl     :: i -> bool  
#  msk_nl_grdconinnz   :: i -> Int32  
#  msk_nl_grdconinsp   :: i -> grdconisuj  
#  msk_nl_laghesnnz    :: yo,ycsub -> Int32
#  msk_nl_laghessp     :: yo,ycsub -> laghessubi,laghessubj

function msk_nl_getsp_wrapper(nlhandle    Ptr{Void},   
                              numgrdobjnz Ptr{Int32},  
                              grdobjsub   Ptr{Int32},  
                              i           Int32,       
                              convali     Ptr{Bool},   
                              grdconinz   Ptr{Int32},  
                              grdconisub  Ptr{Int32},  
                              yo          Int32,       
                              numycnz     Int32,       
                              ycsub       Ptr{Int32},  
                              maxnumhesnz Int32,       
                              numhesnz    Ptr{Int32},  
                              hessubi     Ptr{Int32},  
                              hessubj     Ptr{Int32})  
  (nl_grdobjnnz,
   nl_grdobjsp,
   nl_coniisnl,
   nl_grdconinnz,
   nl_grdconisp,
   nl_laghesnnz,
   nl_laghessp) = unsafe_pointer_to_objref(nlhandle) :: (Function,Function,Function,Function,Function,Function,Function)
  
  numgrdobjnz = nl_grdobjnnz()

  if numgrdobjnz != C_NULL 
    unsafe_store!(numgrdobjnz, numgrdobjnz)
  end

  if  
  

end

#MSKint32t MSKnlgetspfunc ( 
#    MSKuserhandle_t       nlhandle, 
#    MSKint32t *           numgrdobjnz, 
#    MSKint32t *           grdobjsub, 
#    MSKint32t             i, 
#    MSKbooleant *         convali, 
#    MSKint32t *           grdconinz, 
#    MSKint32t *           grdconisub, 
#    MSKint32t             yo, 
#    MSKint32t             numycnz, 
#    MSKCONST MSKint32t *  ycsub, 
#    MSKint32t             maxnumhesnz, 
#    MSKint32t *           numhesnz, 
#    MSKint32t *           hessubi, 
#    MSKint32t *           hessubj);

#MSKint32t MSKnlgetvafunc ( 
#    MSKuserhandle_t       nlhandle, 
#    MSKCONST MSKrealt *   xx, 
#    MSKrealt              yo, 
#    MSKCONST MSKrealt *   yc, 
#    MSKrealt *            objval, 
#    MSKint32t *           numgrdobjnz, 
#    MSKint32t *           grdobjsub, 
#    MSKrealt *            grdobjval, 
#    MSKint32t             numi, 
#    MSKCONST MSKint32t *  subi, 
#    MSKrealt *            conval, 
#    MSKCONST MSKint32t *  grdconptrb, 
#    MSKCONST MSKint32t *  grdconptre, 
#    MSKCONST MSKint32t *  grdconsub, 
#    MSKrealt *            grdconval, 
#    MSKrealt *            grdlag, 
#    MSKint32t             maxnumhesnz, 
#    MSKint32t *           numhesnz, 
#    MSKint32t *           hessubi, 
#    MSKint32t *           hessubj, 
#    MSKrealt *            hesval);
