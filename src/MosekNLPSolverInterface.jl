
# Implements MathProgBase nonlinear interface


type CallbackData
    d::AbstractNLPEvaluator
    numVar::Int
    numConstr::Int
    Ihess::Vector{Int32}
    Jhess::Vector{Int32}
    jac_colval::Vector{Int32} # Compressed sparse row form of Jacobian sparsity
    jac_rowstarts::Vector{Int32}
    jac_nzval::Vector{Float64} # storage for full jacobian (mosek may ask for a subset)
    jac_nnz_original::Int # nnz in NLP evaluator jacobian == length(Ijac)
    jac_idxmap::Vector{Int} # map from indices in NLP evaluator jacobian to mosek jacobian
    J_tmp::Vector{Float64} # storage for NLP evaluator jacobian
    g_tmp::Vector{Float64} # storage for constraint values
end

function msk_nl_getsp_wrapper_mpb(nlhandle::    Ptr{Void},
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
    cb = unsafe_pointer_to_objref(nlhandle)::CallbackData

    i = i_+1

    if numgrdobjnz != C_NULL
        unsafe_store!(numgrdobjnz, convert(Int32,cb.numVar))
    end

    if grdobjsub != C_NULL
        grdobjsub_a = pointer_to_array(grdobjsub,(cb.numVar,))
        for i in 1:cb.numVar
            grdobjsub_a[i] = i-1
        end
    end

    if i <= cb.numConstr
        con_nnz = cb.jac_rowstarts[i+1] - cb.jac_rowstarts[i]
        if convali != C_NULL
            # can we only say that constraint is linear if we've added the linear part separately?
            if con_nnz > 0
                unsafe_store!(convali, convert(Int32,1))
            else
                unsafe_store!(convali, convert(Int32,0))
            end
        end

        if grdconinz != C_NULL
            unsafe_store!(grdconinz, convert(Int32, con_nnz))
        end

        if grdconisub != C_NULL
            if con_nnz > 0
                grdconisub_a = pointer_to_array(grdconisub,(con_nnz,))        
                grdconisub_a[1:con_nnz] = cb.jac_colval[cb.jac_rowstarts[i]:(cb.jac_rowstarts[i+1]-1)] - 1
            end
        end
    end

    hess_nnz = length(cb.Ihess)

    if numhesnz_ != C_NULL
        unsafe_store!(numhesnz_, convert(Int32, hess_nnz))
    end

    if hessubi != C_NULL && hessubj != C_NULL && maxnumhesnz >= hess_nnz
        hessubi_a = pointer_to_array(hessubi,(hess_nnz,))
        hessubj_a = pointer_to_array(hessubj,(hess_nnz,))

        for i in 1:hess_nnz
            hessubi_a[i] = cb.Ihess[i] - 1
            hessubj_a[i] = cb.Jhess[i] - 1
        end
    end

    return @compat(Int32(0))::Int32
end

function msk_nl_getva_wrapper_mpb(nlhandle    :: Ptr{Void},
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
                                  grdconval_  :: Ptr{Float64},
                                  grdlag      :: Ptr{Float64},
                                  maxnumhesnz :: Int32,
                                  numhesnz    :: Ptr{Int32},
                                  hessubi     :: Ptr{Int32},
                                  hessubj     :: Ptr{Int32},
                                  hesval      :: Ptr{Float64})
    cb = unsafe_pointer_to_objref(nlhandle)::CallbackData

    numi = convert(Int,numi_)
    xx = pointer_to_array(xx_,(cb.numVar,))
    yc = pointer_to_array(yc_,(cb.numConstr,))
    subi = pointer_to_array(subi_,(numi,))

    if objval != C_NULL
        unsafe_store!(objval, eval_f(cb.d, xx))
    end

    if numgrdobjnz != C_NULL
        unsafe_store!(numgrdobjnz, convert(Int32,cb.numVar))
    end

    if grdobjsub != C_NULL && grdobjval != C_NULL
        grdobjval_a = pointer_to_array(grdobjval,(cb.numVar,))
        grdobjsub_a = pointer_to_array(grdobjsub,(cb.numVar,))

        eval_grad_f(cb.d, grdobjval_a, xx)
        
        for i in 1:cb.numVar      
            grdobjsub_a[i] = i-1
        end
    end

    if numi > 0 && conval != C_NULL
        conv   = pointer_to_array(conval,(numi,))
        eval_g(cb.d, cb.g_tmp, xx)
        for i=1:numi
            conv[i] = cb.g_tmp[subi[i]+1]
        end
    end

    # do we need to compute Jacobian?
    if grdconval_ != C_NULL || grdlag != C_NULL
        eval_jac_g(cb.d, cb.J_tmp, xx)
        fill!(cb.jac_nzval, 0.0)
        for i in 1:cb.jac_nnz_original
            cb.jac_nzval[cb.jac_idxmap[i]] += cb.J_tmp[i]
        end
    end

    if grdconval_ != C_NULL
        grdconptrb = pointer_to_array(grdconptrb_,(numi,))
        grdconptre = pointer_to_array(grdconptre_,(numi,))
        grdconsub = pointer_to_array(grdconsub_,(grdconptre[numi],))
        grdconval = pointer_to_array(grdconval_,(grdconptre[numi],))
        nz = 1
        for i=1:numi
            ptrb = grdconptrb[i]
            n    = grdconptre[i] - grdconptrb[i]
            con = subi[i]+1
            # make sure input matches the sparsity pattern we provided
            @assert n == (cb.jac_rowstarts[con+1]-cb.jac_rowstarts[con])
            for k in 1:n
                @assert (grdconsub[grdconptrb[i]+k]+1) ==  cb.jac_colval[cb.jac_rowstarts[con]+k-1]
            end
        end

        # extract subset of Jacobian that mosek asked for
        for i in 1:numi
            con = subi[i]+1
            grdconval[(grdconptrb[i]+1):grdconptre[i]] = sub(cb.jac_nzval, cb.jac_rowstarts[con]:(cb.jac_rowstarts[con+1]-1))
        end
    end
    

    if grdlag != C_NULL
        # could use eval_jac_prod_t here, but just do a sparse matvec instead
        grdlag_a = pointer_to_array(grdlag,(cb.numVar,)) 
        eval_grad_f(cb.d, grdlag_a, xx)
        scale!(grdlag_a, yo)
        Jmat = SparseMatrixCSC(cb.numVar,cb.numConstr,cb.jac_rowstarts,cb.jac_colval,cb.jac_nzval)
        A_mul_B!(-1.0, Jmat, yc, 1.0, grdlag_a)
    end

    nhesnz = length(cb.Ihess)

    if numhesnz != C_NULL
        unsafe_store!(numhesnz,convert(Int32,nhesnz))
    end

    if maxnumhesnz > 0 && hessubi != C_NULL && hessubj != C_NULL && hesval != C_NULL
        hessubi_a = pointer_to_array(hessubi,(nhesnz,))
        hessubj_a = pointer_to_array(hessubj,(nhesnz,))
        
        scale!(yc,-1)
        eval_hesslag(cb.d,pointer_to_array(hesval,(nhesnz,)),xx,yo,yc)
        scale!(yc,-1)

        for i=1:length(hessubi_a)
            hessubi_a[i] = cb.Ihess[i]-1
            hessubj_a[i] = cb.Jhess[i]-1
        end
    end
    return convert(Int32,0)::Int32
end

function loadnonlinearproblem!(m::MosekMathProgModel, numVar::Integer, numConstr::Integer, l, u, lb, ub, sense::Symbol, d::AbstractNLPEvaluator)

    if !(numVar == length(l) == length(u)) ||
       !(numConstr == length(lb) == length(ub))
        throw(MosekMathProgModelError("Inconsistent data dimensions"))
    end

    Mosek.deletetask(m.task)
    m.task = maketask(Mosek.msk_global_env)
    loadoptions!(m)

    appendvars(m.task, numVar)
    appendcons(m.task, numConstr)

    m.numvar = numVar
    m.numcon = numConstr
    m.varmap = Int32[1:m.numvar;]
    m.conmap = Int32[1:m.numcon;]
    m.conslack = zeros(Int32,m.numvar)

    # input bounds
    putvarboundslice(m.task, 1, numVar+1, Int32[ MSK_BK_RA for i=1:numVar], float(l), float(u))
    putconboundslice(m.task, 1, numConstr+1, Int32[ MSK_BK_RA for i=1:numConstr], float(lb), float(ub))
    
    setsense!(m, sense)

    initialize(d, [:Grad, :Jac, :Hess])
    Ijac, Jjac = jac_structure(d)
    Ihess, Jhess = hesslag_structure(d)

    # Mosek needs jacobian sparsity in compressed sparse row format, which means that we need to 
    # reformat the indices while maintaining a map
    jac_nnz_original = length(Ijac)
    mergedindices = zeros(Int, jac_nnz_original)
    mergednnz = [0]
    mergedmap = zeros(Int, jac_nnz_original)

    # this is called for duplicate elements
    # assumption: at least one of the indices hasn't been combined before
    function combine(idx1,idx2)
        if mergedmap[idx1] == 0 && mergedmap[idx2] != 0
            mergednnz[1] += 1
            mergedmap[idx1] = idx2
            mergedindices[mergednnz[1]] = idx1
            return idx2
        else
            @assert mergedmap[idx2] == 0
            mergednnz[1] += 1
            mergedmap[idx2] = idx1
            mergedindices[mergednnz[1]] = idx2
            return idx1
        end
    end



    J = sparse(Jjac, Ijac, [i for i in 1:length(Ijac)], numVar, numConstr, combine)

    # map from original index into output index
    idxmap = zeros(Int, jac_nnz_original)
    for row in 1:numConstr
        for pos in J.colptr[row]:(J.colptr[row+1]-1)
            col = J.rowval[pos]
            origidx = J.nzval[pos] # this is the original index of this element
            idxmap[origidx] = pos
        end
    end
    for k in 1:mergednnz[1]
        origidx = mergedindices[k]
        mergedwith = mergedmap[origidx]
        @assert idxmap[origidx] == 0
        @assert idxmap[mergedwith] != 0
        idxmap[origidx] = idxmap[mergedwith]
    end
    for i in 1:jac_nnz_original
        @assert idxmap[i] != 0
    end

    cb = CallbackData(d, numVar, numConstr, Ihess, Jhess, J.rowval, J.colptr, zeros(J.colptr[numConstr+1]-1),
          jac_nnz_original, idxmap, zeros(jac_nnz_original), zeros(numConstr))

    nlgetsp = cfunction(msk_nl_getsp_wrapper_mpb,
                      Int32,
                      (Ptr{Void},Ptr{Int32},Ptr{Int32},Int32,Ptr{Bool},Ptr{Int32},Ptr{Int32},Int32,Int32,Ptr{Int32},Int32,Ptr{Int32},Ptr{Int32},Ptr{Int32}))
    nlgetva = cfunction(msk_nl_getva_wrapper_mpb,
                      Int32,
                      ( Ptr{Void}, # nlhandle
                        Ptr{Float64},Float64,Ptr{Float64}, # xx,yo,yc
                        Ptr{Float64},Ptr{Int32},Ptr{Int32},Ptr{Float64}, # objval,numgrdobjnz,grdobjsub,grdobjval
                        Int32,Ptr{Int32},Ptr{Float64}, # numi,subi,conval
                        Ptr{Int32},Ptr{Int32},Ptr{Int32},Ptr{Float64},Ptr{Float64},# grdconptrb,grdconptre,grdconsub,grdconval,grdlag
                        Int32,Ptr{Int32},Ptr{Int32},Ptr{Int32},Ptr{Float64}))# maxnumhesnz, numhesnz.hessubi,hessubj,hesval
    @Mosek.msk_ccall("putnlfunc",
             Int32, (Ptr{Void},Any,Ptr{Void},Ptr{Void}),
             m.task.task, cb, nlgetsp, nlgetva)

    #@Mosek.msk_ccall(linkfiletotaskstream, Int32,
    # (Ptr{Void},Int32,Ptr{Cchar}), m.task.task,MSK_STREAM_LOG ,"moseklog.txt")
end
