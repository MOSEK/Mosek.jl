MathOptInterface.candelete(m::MosekModel, ref::MathOptInterface.ConstraintReference) = isvalid(m,ref)
isvalid(m::MosekModel, ref::MathOptInterface.ConstraintReference) = allocated(m.c_block,BlockId(ref.value))

MathOptInterface.delete!(m::MosekModel, ref::MathOptInterface.ConstraintReference) = throw(MethodError())











function MathOptInterface.addconstraint!{D <: MathOptInterface.AbstractSet}(m   :: MosekModel,
                                                                            axb :: MathOptInterface.ScalarAffineFunction{Float64},
                                                                            dom :: D)
    N = MathOptInterface.dimension(dom)
    conid,conidxs = addlhsblock!(m, N, Int[1], axb.variables, axb.coefficients)
    append!(m.c_block_coneidx,0)
    if length(m.c_block_constant) < length(m.c_block)
        append!(m.c_block_constant,zeros(Float64,length(m.c_block) - length(m.c_block_constant)))
    end
    m.c_block_constant[conidxs] = axb.constant
    append!(m.c_block_slack,0)
    
    addbound!(m,conidxs,Float64[axb.constant],dom)
    MathOptInterface.ConstraintReference{S}(conid)
end

function MathOptInterface.addconstraint!{D <: MathOptInterface.AbstractSet}(m   :: MosekModel,
                                                                            axb :: MathOptInterface.VectorAffineFunction{Float64},
                                                                            dom :: D)
    conid,conidxs = addlhsblock!(m, MathOptInterface.dimension(dom), Int[1], axb.variables, axb.coefficients)
    append!(m.c_block_coneidx,0)
    if length(m.c_block_constant) < length(m.c_block)
        append!(m.c_block_constant,zeros(Float64,length(m.c_block) - length(m.c_block_constant)))
    end
    m.c_block_constant[conidxs] = axb.constant
    append!(m.c_block_slack,0)

    addbound!(m,conidxs,constant,dom)
    MathOptInterface.ConstraintReference{S}(conid)
end

function MathOptInterface.addconstraint!{D <: MathOptInterface.AbstractSet}(m   :: MosekModel,
                                                                            axb :: MathOptInterface.ScalarQuadraticFunction{Float64},
                                                                            dom :: D)
    conid,conidxs = addlhsblock!(m, MathOptInterface.dimension(dom), Int[1],
                                 axb.affine_variables, axb.affine_coefficients,
                                 ones(Int,length(axb.quadratic_coefficients)),axb.quadratic_rowvariables,axb.quadratic_colvariables,axb.quadratic_coefficients)

    append!(m.c_block_coneidx,0)
    if length(m.c_block_constant) < length(m.c_block)
        append!(m.c_block_constant,zeros(Float64,length(m.c_block) - length(m.c_block_constant)))
    end
    m.c_block_constant[conidxs] = axb.constant
    append!(m.c_block_slack,0)

    addlinearbound!(m,conidxs,constant,dom)
    MathOptInterface.ConstraintReference{S}(conid)
end

function MathOptInterface.addconstraint!{D <: MathOptInterface.AbstractSet}(m   :: MosekModel,
                                                                            axb :: MathOptInterface.VectorQuadraticFunction{Float64},
                                                                            dom :: D)
    conid,conidxs = addlhsblock!(m, MathOptInterface.dimension(dom), Int[1],
                                 axb.affine_variables, axb.affine_coefficients,
                                 axb.quadratic_outputindex,axb.quadratic_rowvariables,axb.quadratic_colvariables,axb.quadratic_coefficients)
    m.c_block_constant[conidxs] = axb.constant
    append!(m.c_block_coneidx,0)
    if length(m.c_block_constant) < length(m.c_block)
        append!(m.c_block_constant,zeros(Float64,length(m.c_block) - length(m.c_block_constant)))
    end
    append!(m.c_block_slack,0)

    addopenbound!(m,conidxs,constant,dom)
    MathOptInterface.ConstraintReference{S}(conid)
end

################################################################################
# Variable constraints

function MathOptInterface.addconstraint!(m :: MosekModel, xs :: MathOptInterface.SingleVariable, dom :: MathOptInterface.Reals)

    subj = getindexes(m.x_block, xs.variable)
    if m.x_boundflags[subj] & (boundflag_lower | boundflag_upper) != 0
        error("Cannot add multiple upper or lower bounds")        
    end

    ensurefree(m.xc_block,1)
    xcid = newblock(m.x_block,1)

    conid,conidxs = addlhsblock!(m, N, Int[1], axb.variables, axb.coefficients)
    
    append!(m.c_block_coneidx,0)
    if length(m.c_block_constant) < length(m.c_block)
        append!(m.c_block_constant,zeros(Float64,length(m.c_block) - length(m.c_block_constant)))
    end
    m.c_block_constant[conidxs] = axb.constant
    append!(m.c_block_slack,0)
    
    addbound!(m,conidxs,Float64[axb.constant],dom)
    MathOptInterface.ConstraintReference{S}(conid)
end


# function MathOptInterface.addconstraint!(m   :: MosekModel,
#                                          xs  :: MathOptInterface.VectorOfVariables,
#                                          dom :: MathOptInterface.)
#     if MathOptInterface.dimension(dom) != 1
#         error("Invalid set size")
#     end

# ......    

#     conid,conidxs = addlhsblock!(m, N, Int[1], axb.variables, axb.coefficients)
    
#     append!(m.c_block_coneidx,0)
#     if length(m.c_block_constant) < length(m.c_block)
#         append!(m.c_block_constant,zeros(Float64,length(m.c_block) - length(m.c_block_constant)))
#     end
#     m.c_block_constant[conidxs] = axb.constant
#     append!(m.c_block_slack,0)
    
#     addbound!(m,conidxs,Float64[axb.constant],dom)
#     MathOptInterface.ConstraintReference{S}(conid)
# end


################################################################################
################################################################################


# This adds a set of rows to A
function addlhsblock!(m        :: MosekModel,
                      N        :: Int,
                      subi     :: Vector{Int},
                      varidxs  :: Vector{MathOptInterface.VariableReference},
                      cofs     :: Vector{Float64})
    allocateconstraints(m,N)
    
    subj = Array{Int}(length(varidxs))
    for i in 1:length(subj)
        getindexes(m.x_block,varidx[i],subj,i)
    end

    conid = newblock(m.c_block,1,N)
    conidxs = getindexes(m.c_block,conid)

    At = sparse(subi, subi, subj, getnumvar(m.task), N)
    putarowlist(m.task,conidxs,At)

    conid,conidxs
end    

# This adds a set of rows to A and corresponding Q entries
function addlhsblock!(m        :: MosekModel,
                      N        :: Int,
                      subi     :: Vector{Int},
                      varidxs  :: Vector{MathOptInterface.VariableReference},
                      cofs     :: Vector{Float64},
                      qconi    :: Vector{Int},
                      qvaridxi :: Vector{MathOptInterface.VariableReference},
                      qvaridxj :: Vector{MathOptInterface.VariableReference},
                      qcof     :: Vector{Float64})
    if len(qsubi) > 0
        if m.problemtype == problemtype_conic
            error("Cannot mix cones and quadratic terms")
        end
        m.problemtype = problemtype_quadratic
    end

    conid, conidxs = addlhsblock!(m,N,subi,varidxs,cofs)

    qsubi = Array{Int}(length(qvaridxi))
    qsubj = Array{Int}(length(qvaridxj))
    for i in 1:length(qsubi)
        getindexes(m.x_block,qvaridxi[i],qsubi,i)
        getindexes(m.x_block,qvaridxj[i],qsubj,i)
    end

    qconk = conidxs[qconi]
    perm = sort(1:length(qconk),by=i->qconk[i])

    i = 1
    while i < length(perm)
        i0 = i
        k = qconk[perm[i0]]
        while i < length(perm) && k == qconk[perm[i]] i += 1 end

        putqconk(m.task,k,qsubi[perm[i0:i-1]],qsubj[perm[i0:i-1]],qcof[perm[i0:i-1]])
    end
    
    conid,conidxs
end    

addopenbound!(m :: MosekModel, conidxs :: Vector{Int}, constant :: Vector{Float64}, dom :: MathOptInterface.Reals) = putconboundlist(m.task,conidxs,fill(MSK_BK_FR,MathOptInterface.dimension(dom)),-constant,-constant)
addopenbound!(m :: MosekModel, conidxs :: Vector{Int}, constant :: Vector{Float64}, dom :: MathOptInterface.Zeros) = putconboundlist(m.task,conidxs,fill(MSK_BK_FX,MathOptInterface.dimension(dim)),-constant,-constant)
addopenbound!(m :: MosekModel, conidxs :: Vector{Int}, constant :: Vector{Float64}, dom :: MathOptInterface.Nonnegatives) = putconboundlist(m.task,conidxs,fill(MSK_BK_LO,MathOptInterface.dimension(dim)),-constant,-constant)
addopenbound!(m :: MosekModel, conidxs :: Vector{Int}, constant :: Vector{Float64}, dom :: MathOptInterface.Nonpositives) = putconboundlist(m.task,conidxs,fill(MSK_BK_UP,MathOptInterface.dimension(dim)),-constant,-constant)
addopenbound!(m :: MosekModel, conidxs :: Vector{Int}, constant :: Vector{Float64}, dom :: MathOptInterface.GreaterThan{Float64}) = putconboundlist(m.task,conidxs,Float64[MSK_BK_LO],dom.lower-constant,dom.lower-constant)
addopenbound!(m :: MosekModel, conidxs :: Vector{Int}, constant :: Vector{Float64}, dom :: MathOptInterface.LessThan{Float64})    = putconboundlist(m.task,conidxs,Float64[MSK_BK_UP],dom.upper-constant,dom.upper-constant)

addbound!(m :: MosekModel, conid :: Int, conidxs :: Vector{Int}, constant :: Vector{Float64}, dom :: MathOptInterface.Reals)                = putopenbound!(m,conidxs,constant,dom)
addbound!(m :: MosekModel, conid :: Int, conidxs :: Vector{Int}, constant :: Vector{Float64}, dom :: MathOptInterface.Zeros)                = putopenbound!(m,conidxs,constant,dom)
addbound!(m :: MosekModel, conid :: Int, conidxs :: Vector{Int}, constant :: Vector{Float64}, dom :: MathOptInterface.Nonnegatives)         = putopenbound!(m,conidxs,constant,dom)
addbound!(m :: MosekModel, conid :: Int, conidxs :: Vector{Int}, constant :: Vector{Float64}, dom :: MathOptInterface.Nonpositives)         = putopenbound!(m,conidxs,constant,dom)
addbound!(m :: MosekModel, conid :: Int, conidxs :: Vector{Int}, constant :: Vector{Float64}, dom :: MathOptInterface.GreaterThan{Float64}) = putopenbound!(m,conidxs,constant,dom)
addbound!(m :: MosekModel, conid :: Int, conidxs :: Vector{Int}, constant :: Vector{Float64}, dom :: MathOptInterface.LessThan{Float64})    = putopenbound!(m,conidxs,constant,dom)
addbound!(m :: MosekModel, conid :: Int, conidxs :: Vector{Int}, constant :: Vector{Float64}, dom :: MathOptInterface.EqualTo{Float64})     = putconboundlist(m.task,conidxs,Float64[MSK_BK_FX],dom.value-constant,dom.value-constant)
addbound!(m :: MosekModel, conid :: Int, conidxs :: Vector{Int}, constant :: Vector{Float64}, dom :: MathOptInterface.Interval{Float64})    = putconboundlist(m.task,conidxs,Float64[MSK_BK_RA],dom.lower-constant,dom.upper-constant)


function add_slack!(m :: MosekModel, conidxs :: Vector{Int}, constant :: Vector{Float64}, N :: Int)
    ensurefree(m.x_block,N)
    varid = newblock(m.x_block,Mosek_VAR,N)

    numvar = getnumvar(m.task)
    if length(s.x_block) > numvar
        appendvars(m.task, length(m.x_block) - numvar)
    end

    subj = getindexes(m.x_block,varid)

    putaijlist(m.task,conidxs,subj,-ones(Float64,N))
    putvarboundlist(m.task,subj,fill(MSK_BK_FR,N),zeros(Float64,N),zeros(Float64,N))
    
    putconboundlist(m.task,conidxs,Float64[MSK_BK_FX],-constant,-constant)

    varid,subj
end

function addbound!(m :: MosekModel, conid :: Int, conidxs :: Vector{Int}, constant :: Vector{Float64}, dom :: MathOptInterface.SecondOrderCone)
    N = MathOptInterface.dimension(dom)

    varid,subj = add_slack!(m,conidxs,constant,N)
    
    appendcone(m.task,MSK_CT_QUAD,subj)
    numcone = getnumcone(m.task)
    m.c_block_coneidx[conid] = numcone
    m.c_block_slack[conid]   = varid
end

function addbound!(m :: MosekModel, conid :: Int, conidxs :: Vector{Int}, constant :: Vector{Float64}, dom :: MathOptInterface.RotatedSecondOrderCone)
    N = MathOptInterface.dimension(dom)

    varid,subj = add_slack!(m,conidxs,constant,N)
    
    appendcone(m.task,MSK_CT_RQUAD,subj)
    numcone = getnumcone(m.task)
    m.c_block_coneidx[conid] = numcone
    m.c_block_slack[conid]   = varid
end

function addbound!(m :: MosekModel, conid :: Int, conidxs :: Vector{Int}, constant :: Vector{Float64}, dom :: MathOptInterface.PositiveSemidefiniteConeTriangle)
    dim = MathOptInterface.dimension(dom)

    appendbarvars(m.task,Int32[dim])
    barvaridx = getnumbarvar(m.task)

    idx = 1
    for j in 1:dim
        for i in i : dim
            matrixid = appendsparsesymmat(m.task,Int32(dim), Int32[i], Int32[j], Float64[-1.0])
            putbaraij(m.task, conidxs[idx], numvarvar, Int[matrixid], Float64[1.0])
            idx += 1
        end
    end

    putconboundlist(m.task,conidxs,Float64[MSK_BK_FX],-constant,-constant)

    c_block_coneidx[conid] = -barvaridx
    c_block_coneidx[conid] = -barvaridx
end






# function MathOptInterface.addconstraint!(m   :: MosekModel,
#                                          axb :: MathOptInterface.ScalarAffineFunction{Float64},
#                                          dom :: MathOptInterface.AbstractSet)
#     N = 1
#     if N != MathOptInterface.dimension(dom)
#         error("Dimensions mismatch")
#     end

#     conid,conidxs = makeconstr(m,ones(Int,length(axb.coefficients)),axb.variables,axb.coefficients,N,1)
#     bfix = axb.constant
#     if     typeof(dom) == MathOptInterface.Zeros
#         putconboundlist(m.task,conidxs,fill(MSK_BK_FX,N),Float64[-bfix],Float64[-bfix])
#     elseif typeof(dom) == MathOptInterface.Reals
#         putconboundlist(m.task,conidxs,fill(MSK_BK_FR,N),Float64[-bfix],Float64[-bfix])
#     elseif typeof(dom) == MathOptInterface.Nonnegatives
#         putconboundlist(m.task,conidxs,fill(MSK_BK_LO,N),Float64[-bfix],Float64[-bfix])
#     elseif typeof(dom) == MathOptInterface.Nonpositives
#         putconboundlist(m.task,conidxs,fill(MSK_BK_UP,N),Float64[-bfix],Float64[-bfix])
#     elseif typeof(dom) == MathOptInterface.GreaterThan{Float64}
#         putconboundlist(m.task,conidxs,fill(MSK_BK_LO,N),Float64[dom.lower-bfix],Float64[0])
#     elseif typeof(dom) == MathOptInterface.LessThan{Float64}
#         putconboundlist(m.task,conidxs,fill(MSK_BK_UP,N),Float64[0],Float64[dom.upper-bfix])
#     elseif typeof(dom) == MathOptInterface.EqualTo{Float64}
#         putconboundlist(m.task,conidxs,fill(MSK_BK_FX,N),Float64[dom.value-bfix],Float64[dom.value-bfix])
#     elseif typeof(dom) == MathOptInterface.Interval{Float64}
#         putconboundlist(m.task,conidxs,fill(MSK_BK_RA,N),Float64[dom.lower-bfix],Float64[dom.upper-bfix])
#     else
#         error("Invalid domain set")
#     end

#     push!(m.c_block_type, mosek_block_type_noneg)
#     push!(m.c_block_coneidx, 0)
#     push!(m.c_block_slack,   0)
    
#     MathOptInterface.ConstraintReference{S}(conid)
# end

# function MathOptInterface.addconstraint!(m   :: MosekModel,
#                                          axb :: MathOptInterface.VectorAffineFunction{Float64},
#                                          dom :: MathOptInterface.Abstract)
#     N = MathOptInterface.dimension(dom)

#     conid,conidxs = makeconstr(m,axb.outputindex,axb.variables,axb.coefficients,N)
#     bfix = axb.constant
#     if     typeof(dom) == MathOptInterface.Zeros
#         putconboundlist(m.task,conidxs,fill(MSK_BK_FX,N),-bfix,-bfix)
#     elseif typeof(dom) == MathOptInterface.Reals
#         putconboundlist(m.task,conidxs,fill(MSK_BK_FR,N),-bfix,-bfix)
#     elseif typeof(dom) == MathOptInterface.Nonnegatives
#         putconboundlist(m.task,conidxs,fill(MSK_BK_LO,N),-bfix,-bfix)
#     elseif typeof(dom) == MathOptInterface.Nonpositives
#         putconboundlist(m.task,conidxs,fill(MSK_BK_UP,N),-bfix,-bfix)
#     elseif typeof(dom) == MathOptInterface.GreaterThan{Float64}
#         putconboundlist(m.task,conidxs,fill(MSK_BK_LO,N),Float64[dom.lower-bfix[1]],Float64[0])
#     elseif typeof(dom) == MathOptInterface.LessThan{Float64}
#         putconboundlist(m.task,conidxs,fill(MSK_BK_UP,N),Float64[0],Float64[dom.upper-bfix[1]])
#     elseif typeof(dom) == MathOptInterface.EqualTo{Float64}
#         putconboundlist(m.task,conidxs,fill(MSK_BK_FX,N),Float64[dom.value-bfix[1]],Float64[dom.value-bfix[1]])
#     elseif typeof(dom) == MathOptInterface.Interval{Float64}
#         putconboundlist(m.task,conidxs,fill(MSK_BK_RA,N),Float64[dom.lower-bfix[1]],Float64[dom.upper-bfix[1]])
#     elseif typeof(dom) in [MathOptInterface.SecondOrderCone{Float64},
                           
#         if 
#         putaijlist()
        
#         putconboundlist(m.task,conidxs,fill(MSK_BK_RA,N),Float64[dom.lower-bfix[1]],Float64[dom.upper-bfix[1]])
        
#     else
#         error("Invalid domain set")
#     end

    
    
#     putconboundlist(m.task,conidxs,fill(MSK_BK_LO,N),-b,b)

#     push!(m.c_block_type, mosek_block_type_noneg)
#     push!(m.c_block_coneidx, 0)
#     push!(m.c_block_slack,   0)
#     MathOptInterface.ConstraintReference{S}(conid)
# end

# function MathOptInterface.addconstraint!(m   :: MosekModel,
#                                          fun :: MathOptInterface.ScalarQuadraticFunction,
#                                          dom :: MathOptInterface.Abstract)
# end

# function MathOptInterface.addconstraint!(m   :: MosekModel,
#                                          fun :: MathOptInterface.VectorQuadraticFunction,
#                                          dom :: MathOptInterface.Abstract)
# end




# function MathOptInterface.addconstraint!(
#     m           :: MosekModel,
#     b           :: Vector{Float64},
#     a_constridx :: Vector{Int},
#     a_varidx    :: Vector{VariableReference},
#     a_coef      :: Vector{Float64},
#     S           :: MathOptInterface.NonNegative)

#     N = S.dim
#     conid,conidxs = makeconstr(m,a_constridx, a_varidx,a_coef,N)
#     putconboundlist(m.task,conidxs,fill(MSK_BK_LO,N),-b,b)

#     push!(m.c_block_type, mosek_block_type_noneg)
#     push!(m.c_block_coneidx, 0)
#     push!(m.c_block_slack,   0)
#     MathOptInterface.ConstraintReference{S}(conid)
# end

# function MathOptInterface.addconstraint!(
#     m           :: MosekModel,
#     b           :: Vector{Float64},
#     a_constridx :: Vector{Int},
#     a_varidx    :: Vector{VariableReference},
#     a_coef      :: Vector{Float64},
#     S           :: MathOptInterface.NonPositive)

#     N = S.dim
#     conid,conidxs = makeconstr(m,a_constridx, a_varidx,a_coef,N)
#     putconboundlist(m.task,conidxs,fill(MSK_BK_UP,N),b,-b)

#     push!(m.c_block_type, mosek_block_type_nopos)
#     push!(m.c_block_coneidx, 0)
#     push!(m.c_block_slack,   0)
#     MathOptInterface.ConstraintReference{S}(conid)
# end

# function MathOptInterface.addconstraint!(
#     m           :: MosekModel,
#     b           :: Vector{Float64},
#     a_constridx :: Vector{Int},
#     a_varidx    :: Vector{VariableReference},
#     a_coef      :: Vector{Float64},
#     S           :: MathOptInterface.Zero)

#     N = S.dim
#     conid,conidxs = makeconstr(m,a_constridx, a_varidx,a_coef,N)
#     putconboundlist(m.task,conidxs,fill(MSK_BK_FX,N),-b,-b)

#     push!(m.c_block_type, mosek_block_type_nopos)
#     push!(m.c_block_coneidx, 0)
#     push!(m.c_block_slack,   0)
#     MathOptInterface.ConstraintReference{S}(conid)
# end

# function MathOptInterface.addconstraint!(
#     m           :: MosekModel,
#     b           :: Vector{Float64},
#     a_constridx :: Vector{Int},
#     a_varidx    :: Vector{VariableReference},
#     a_coef      :: Vector{Float64},
#     S           :: MathOptInterface.Interval)

#     N = length(S.lower)
#     conidxs = makeconstr(m,a_constridx, a_varidx,a_coef,N)

#     bl = S.lower - b
#     bu = S.upper - b
#     putconboundlist(m.task,conidxs,fill(MSK_BK_RA,N),bl,bu)

#     push!(m.c_block_type, mosek_block_type_nopos)
#     push!(m.c_block_coneidx, 0)
#     push!(m.c_block_slack,   0)
#     MathOptInterface.ConstraintReference{S}(conid)
# end

# function MathOptInterface.addconstraint!(
#     m           :: MosekModel,
#     varidx      :: VariableReference,
#     S           :: MathOptInterface.Integers)

# #     N = S.dim
# #     MathOptInterface.ConstraintReference{S}(varidx)
# # end









# function allocateconstraints(m           :: MosekModel,
#                              N           :: Int)
#     numcon = getnumcon(m.task)
#     ensurefree(m.c_block,N)
#     if length(s.c_block) > numcon
#         appendcons(s.task, length(s.c_block) - numcon)
#     end
# end




# function makeconstr(m           :: MosekModel,
#                     a_constridx :: Vector{Int},
#                     a_varidx    :: Vector{MathOptInterface.VariableReference},
#                     a_coef      :: Vector{Float64},
#                     N           :: Int)

#     allocateconstraints(m,N)
    
#     varidxs = Array{Int}(length(a_varidx))
#     for i in 1:length(a_varidx)
#         getindexes(m.c_block,a_varidx[i],varidxs,i)
#     end

#     conid = newblock(m.c_block,1,N)
#     conidxs = getindexes(m.c_block,conid)

#     At = sparse(idxs, a_constridx, a_coef, getnumvar(m.task), N)    
#     putarowlist(m.task,conidxs,At)

#     conid,conidxs
# end
