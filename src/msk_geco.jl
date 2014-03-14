# This file contains the implementation of the general convex interface

export
  putnlcallbacks,
  clearnlcallbacks

type MSKnlinfo
  numvar :: Int
  numcon :: Int

  nlgetva :: Ptr{Void}
  nlgetsp :: Ptr{Void}

  # these are stored with 0-base indexes
  grdobjsub  :: Array{Int32,1}
  grdconsub  :: Array{Int32,1}
  grdconptr  :: Array{Int32,1}

  hessubi    :: Array{Int32,1}
  hessubj    :: Array{Int32,1}


  evalobj  :: Function
  evalconi :: Function
  grdlag   :: Function
  heslag   :: Function
  grdobj   :: Function
  grdconi  :: Function
end

# This wraps the callback method returning the sparsity of the non-linear terms.
# Depending on which arguments are NULL it should return gradient and hessian sparsity.
# NOTE: This is not documented (at all!), but in reality it is only called in a few ways.
#       While the arguments allow fetching only a portion of the hessian pattern, we ignore this
#       and always return the whole thing. It doesn't really matter since the sizes are only
#       used for pre-allocating space in MOSEK - only a reasonable upper bound is really needed.
# NOTE: The 'hessian' here is the second derivative of the lagrangian
function msk_nl_getsp_wrapper(nlhandle::    Ptr{Void},
                              numgrdobjnz:: Ptr{Int32}, # number of nonzeros in gradient of objective
                              grdobjsub::   Ptr{Int32}, # subscripts of nonzeros in gradient of objective
                              i_::          Int32,      # constraint index
                              convali::     Ptr{Bool},  # 0/1 indicating whether constraint i is non-linear
                              grdconinz::   Ptr{Int32}, # number of nonzeros in gradient of constraint i
                              grdconisub::  Ptr{Int32}, # subscripts of nonzeros in gradient of constraint i
                              yo::          Int32,      # 0/1 include objective in computation of hessian pattern
                              numycnz::     Int32,      # number of constraints to include in computation of the hessian pattern
                              ycsub::       Ptr{Int32}, # indexes of constraints to include in computation of the hessian pattern
                              maxnumhesnz:: Int32,      # lengths of hessubi and hessubj
                              numhesnz_::   Ptr{Int32}, # number of hessian nonzeros
                              hessubi::     Ptr{Int32}, # column subscrips of hessian non-zeros
                              hessubj::     Ptr{Int32}) # row subscripts of hessian non-zeros
  nlinfo = unsafe_pointer_to_objref(nlhandle) :: MSKnlinfo

  i = i_+1

  grdobjlen = length(nlinfo.grdobjsub)

  if numgrdobjnz != C_NULL
    unsafe_store!(numgrdobjnz, convert(Int32,grdobjlen))
  end

  if grdobjsub != C_NULL
    if grdobjlen > 0
      grdobjsub_a = pointer_to_array(grdobjsub,(grdobjlen,))
      grdobjsub_a[1:grdobjlen] = nlinfo.grdobjsub
    end
  end

  if i <= nlinfo.numcon
    if convali != C_NULL
      if nlinfo.grdconptr[i+1] > nlinfo.grdconptr[i]
        unsafe_store!(convali, convert(Int32,1))
      else
        unsafe_store!(convali, convert(Int32,0))
      end
    end

    if grdconinz != C_NULL
      unsafe_store!(grdconinz, convert(Int32, nlinfo.grdconptr[i+1] - nlinfo.grdconptr[i]))
    end

    if grdconisub != C_NULL
      num = nlinfo.grdconptr[i+1] - nlinfo.grdconptr[i]
      if num > 0
        grdconisub_a = pointer_to_array(grdconisub,num)        
        grdconisub_a[1:num] = nlinfo.grdconsub[nlinfo.grdconptr[i]+1:nlinfo.grdconptr[i+1]]
      end
    end
  end

  numhesnz = length(nlinfo.hessubi)

  if numhesnz_ != C_NULL
    unsafe_store!(numhesnz_, convert(Int32, numhesnz))
  end

  if hessubi != C_NULL && hessubj != C_NULL && maxnumhesnz >= numhesnz
    hessubi_a = pointer_to_array(hessubi,(numhesnz,))
    hessubj_a = pointer_to_array(hessubj,(numhesnz,))

    hessubi_a[1:numhesnz] = nlinfo.hessubi
    hessubj_a[1:numhesnz] = nlinfo.hessubj
  end

  return int32(0) :: Int32
end

function msk_nl_getva_wrapper(nlhandle    :: Ptr{Void},
                              xx_         :: Ptr{Float64}, # input
                              yo          :: Float64,
                              yc_         :: Ptr{Float64}, # input, length = numcon
                              objval      :: Ptr{Float64},
                              numgrdobjnz :: Ptr{Int32},
                              grdobjsub   :: Ptr{Int32},
                              grdobjval   :: Ptr{Float64},
                              numi_       :: Int32,
                              subi_       :: Ptr{Int32},   # input
                              conval      :: Ptr{Float64},
                              grdconptrb_ :: Ptr{Int32},   # input
                              grdconptre_ :: Ptr{Int32},   # input
                              grdconsub_  :: Ptr{Int32},   # input
                              grdconval   :: Ptr{Float64},
                              grdlag      :: Ptr{Float64},
                              maxnumhesnz :: Int32,
                              numhesnz    :: Ptr{Int32},
                              hessubi     :: Ptr{Int32},
                              hessubj     :: Ptr{Int32},
                              hesval      :: Ptr{Float64})
  nlinfo = unsafe_pointer_to_objref(nlhandle) :: MSKnlinfo

  numi = convert(Int,numi_)
  xx = pointer_to_array(xx_,(nlinfo.numvar,))
  yc = pointer_to_array(yc_,(nlinfo.numcon,))
  subi = pointer_to_array(subi_,(numi,)) .+ 1

  if objval != C_NULL
    unsafe_store!(objval, nlinfo.evalobj(xx))
  end

  ngrdobjnz = length(nlinfo.grdobjsub)

  if numgrdobjnz != C_NULL
    unsafe_store!(numgrdobjnz, convert(Int32,ngrdobjnz))
  end

  if grdobjsub != C_NULL && grdobjval != C_NULL
    grdobjval_a = pointer_to_array(grdobjval,(ngrdobjnz,))
    grdobjsub_a = pointer_to_array(grdobjsub,(ngrdobjnz,))

    grdobjsub = nlinfo.grdobjsub .+ 1
    nlinfo.grdobj(xx,
                  grdobjsub,
                  grdobjval_a)
    for i in 1:length(grdobjsub_a)      
      grdobjsub_a[i] = nlinfo.grdobjsub[i]
    end
  end

  if numi > 0 && conval != C_NULL
    conv   = pointer_to_array(conval,(numi,))
    for i=1:numi
      conv[i] = nlinfo.evalconi(xx,subi[i])
    end
  end

  if grdconval != C_NULL
    grdconptrb = pointer_to_array(grdconptrb_,(numi,))
    grdconptre = pointer_to_array(grdconptre_,(numi,))
    for i=1:numi
      ptrb = grdconptrb[i]
      n    = grdconptre[i] - grdconptrb[i]
      nlinfo.grdconi(xx,subi[i],
                     pointer_to_array(grdconsub_+ptrb*4,(n,)) .+ 1,
                     pointer_to_array(grdconval +ptrb*4,(n,)))
    end
  end

  if grdlag != C_NULL
    nlinfo.grdlag(xx,yo,yc,subi,pointer_to_array(grdlag,(nlinfo.numvar,)))
  end

  nhesnz = length(nlinfo.hessubi)

  if numhesnz != C_NULL
    unsafe_store!(numhesnz,convert(Int32,nhesnz))
  end

  if maxnumhesnz > 0 && hessubi != C_NULL && hessubj != C_NULL && hesval != C_NULL
    hessubi_a = pointer_to_array(hessubi,(nhesnz,))
    hessubj_a = pointer_to_array(hessubj,(nhesnz,))
    nlinfo.heslag(xx,yo,yc,subi,
                  hessubi_a,
                  hessubj_a,
                  pointer_to_array(hesval, (nhesnz,)))

    # Hum... Very strange! If I use "hessubi_a -= 1" it appears that the
    # underlying data of the array is not the native pointer anymore, but
    # if I do it this way it is.
    for i=1:length(hessubi_a)
      hessubi_a[i] = hessubi_a[i]-convert(Int32,1)
      hessubj_a[i] = hessubj_a[i]-convert(Int32,1)
    end

    hesval_a  = pointer_to_array(hesval,(nhesnz,))
  end
  return convert(Int32,0) :: Int32
end

function putnlcallbacks(task::MSKtask,
                        grdobjsub :: Array{Int,1},
                        grdconsub :: Array{Int,1},
                        grdconptr :: Array{Int,1},
                        hessubi   :: Array{Int,1},
                        hessubj   :: Array{Int,1},
                        evalobj   :: Function,
                        evalconi  :: Function,
                        grdlag    :: Function,
                        grdobj    :: Function,
                        grdconi   :: Function,
                        heslag    :: Function)
  nvar = convert(Int,getnumvar(task))
  ncon = convert(Int,getnumcon(task))

  print( ncon,length(grdconptr))
  if ( length(hessubi) !=  length(hessubj))
    error("Arrays hessubi and hessubj have mismatching lengths")
  end
  if ( length(grdconptr) != ncon+1 )
    error("Length of grdconptr should match number of constraints")
  end

  nlgetsp = cfunction(msk_nl_getsp_wrapper,
                      Int32,
                      (Ptr{Void},Ptr{Int32},Ptr{Int32},Int32,Ptr{Bool},Ptr{Int32},Ptr{Int32},Int32,Int32,Ptr{Int32},Int32,Ptr{Int32},Ptr{Int32},Ptr{Int32}))
  nlgetva = cfunction(msk_nl_getva_wrapper,
                      Int32,
                      ( Ptr{Void}, # nlhandle
                        Ptr{Float64},Float64,Ptr{Float64}, # xx,yo,yc
                        Ptr{Float64},Ptr{Int32},Ptr{Int32},Ptr{Float64}, # objval,numgrdobjnz,grdobjsub,grdobjval
                        Int32,Ptr{Int32},Ptr{Float64}, # numi,subi,conval
                        Ptr{Int32},Ptr{Int32},Ptr{Int32},Ptr{Float64},Ptr{Float64},# grdconptrb,grdconptre,grdconsub,grdconval,grdlag
                        Int32,Ptr{Int32},Ptr{Int32},Ptr{Int32},Ptr{Float64}))# maxnumhesnz, numhesnz.hessubi,hessubj,hesval

  nlinfo = MSKnlinfo(nvar,ncon,
                     nlgetsp,nlgetva,
                     int32(grdobjsub .- 1),
                     int32(grdconsub .- 1),
                     int32(grdconptr .- 1),
                     int32(hessubi .- 1),
                     int32(hessubj .- 1),
                     evalobj,evalconi, grdlag,heslag,grdobj,grdconi)

  @msk_ccall("putnlfunc",
             Int32, (Ptr{Void},Any,Ptr{Void},Ptr{Void}),
             task.task, nlinfo, nlgetsp, nlgetva)
  task.nlinfo = nlinfo
end

function clearnlcallbacks(task::MSKtask)
  @msk_ccall("putnlfunc",
             Int32, (Ptr{Void},Ptr{Void},Ptr{Void},Ptr{Void}),
             task.task, C_NULL,C_NULL,C_NULL,C_NULL)
  task.nlinfo = nothing
end
