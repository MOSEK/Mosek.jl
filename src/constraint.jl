MathOptInterface.candelete(m::MosekModel, ref::MathOptInterface.ConstraintReference) = isvalid(m,ref)
MathOptInterface.delete!(  m::MosekModel, ref::MathOptInterface.ConstraintReference) = throw(MethodError())

function MathOptInterface.addconstraint!{D <: MathOptInterface.AbstractScalarSet}(m   :: MosekModel,
                                                                                  axb :: MathOptInterface.ScalarAffineFunction{Float64},
                                                                                  dom :: D)
    N = 1
    conid = allocateconstraints(m,N)
    addlhsblock!(m, conid, fill(1,length(axb.variables)),axb.variables, axb.coefficients)
    
    if length(m.c_constant) < length(m.c_block)
        append!(m.c_constant,zeros(Float64,length(m.c_block) - length(m.c_constant)))
    end
    conidxs = getindexes(m.c_block,conid)
    m.c_constant[conidxs] = axb.constant
    
    addbound!(m,conid,conidxs,Float64[axb.constant],dom)

    
    conref = MathOptInterface.ConstraintReference{MathOptInterface.ScalarAffineFunction{Float64},D}(UInt64(conid) << 1)
    select(m.constrmap,MathOptInterface.ScalarAffineFunction{Float64},D)[conref.value] = conid
    
    conref
end

function MathOptInterface.addconstraint!{D <: MathOptInterface.AbstractVectorSet}(m   :: MosekModel,
                                                                                  axb :: MathOptInterface.VectorAffineFunction{Float64},
                                                                                  dom :: D)
    N = MathOptInterface.dimension(dom)
    conid = allocateconstraints(m,N)
    addlhsblock!(m,
                 conid,
                 axb.outputindex,
                 axb.variables,
                 axb.coefficients)
    conidxs = getindexes(m.c_block,conid)
    m.c_constant[conidxs] = axb.constant

    addbound!(m,conid,conidxs,axb.constant,dom)
    conref = MathOptInterface.ConstraintReference{MathOptInterface.VectorAffineFunction{Float64},D}(UInt64(conid) << 1)
    select(m.constrmap,MathOptInterface.VectorAffineFunction{Float64},D)[conref.value] = conid    
    
    conref
end

function MathOptInterface.addconstraint!(m   :: MosekModel,
                                         axb :: MathOptInterface.VectorAffineFunction{Float64},
                                         dom :: MathOptInterface.PositiveSemidefiniteConeTriangle)
    N = MathOptInterface.dimension(dom)
    M = (N+1)*N >> 1
    conid = allocateconstraints(m,M)
    addlhsblock!(m,
                 conid,
                 axb.outputindex,
                 axb.variables,
                 axb.coefficients)
    conidxs = getindexes(m.c_block,conid)
    m.c_constant[conidxs] = axb.constant

    addbound!(m,conid,conidxs,axb.constant,dom)
    conref = MathOptInterface.ConstraintReference{MathOptInterface.VectorAffineFunction{Float64},MathOptInterface.PositiveSemidefiniteConeTriangle}(UInt64(conid) << 1)
    select(m.constrmap,MathOptInterface.VectorAffineFunction{Float64},MathOptInterface.PositiveSemidefiniteConeTriangle)[conref.value] = conid
    conref
end

################################################################################
# Variable constraints

# We allow following. Each variable can have
# - at most most upper and one lower bound
# - belong to at most one non-semidefinite cone
# - any number of semidefinite cones, which are implemented as ordinary constraints
# This is when things get a bit funky; By default a variable has no
# bounds, i.e. "free". Adding a `GreaterThan` or `Nonnegatives`
# constraint causes it to have a defined lower bound but no upper
# bound, allowing a `LessThan` or ``Nonpositives` constraint to be
# added later. Adding a `Interval` constraint defines both upper and
# lower bounds. Adding a `Reals` constraint will effectively be the
# same as an interval ]-infty;infty[, in that it will define both
# upper and lower bounds, and not allow those to be set afterwards.

domain_type_mask(dom :: MathOptInterface.Reals)        = (boundflag_lower | boundflag_upper)
domain_type_mask(dom :: MathOptInterface.Interval)     = (boundflag_lower | boundflag_upper)
domain_type_mask(dom :: MathOptInterface.EqualTo)      = (boundflag_lower | boundflag_upper)
domain_type_mask(dom :: MathOptInterface.GreaterThan)  = boundflag_lower
domain_type_mask(dom :: MathOptInterface.Nonnegatives) = boundflag_lower
domain_type_mask(dom :: MathOptInterface.LessThan)     = boundflag_upper
domain_type_mask(dom :: MathOptInterface.Nonpositives) = boundflag_upper

domain_type_mask(dom :: MathOptInterface.SecondOrderCone)        = boundflag_cone
domain_type_mask(dom :: MathOptInterface.RotatedSecondOrderCone) = boundflag_cone
domain_type_mask(dom :: MathOptInterface.ExponentialCone)        = boundflag_cone
domain_type_mask(dom :: MathOptInterface.PowerCone)              = boundflag_cone

domain_type_mask(dom :: MathOptInterface.PositiveSemidefiniteConeTriangle) = 0
domain_type_mask(dom :: MathOptInterface.PositiveSemidefiniteConeScaled)   = 0

domain_type_mask(dom :: MathOptInterface.Integer) = boundflag_int
domain_type_mask(dom :: MathOptInterface.ZeroOne) = (boundflag_int | boundflag_upper | boundflag_lower)

is_positivesemidefinite_set(dom :: Union{MathOptInterface.PositiveSemidefiniteConeTriangle,MathOptInterface.PositiveSemidefiniteConeScaled}) = true
is_positivesemidefinite_set(dom) = false
is_conic_set(dom :: Union{MathOptInterface.SecondOrderCone,
                          MathOptInterface.RotatedSecondOrderCone,
                          MathOptInterface.PowerCone,
                          MathOptInterface.DualPowerCone,
                          MathOptInterface.ExponentialCone,
                          MathOptInterface.DualExponentialCone}) = true
is_conic_set(dom) = false


addvarconstr(m :: MosekModel, subj :: Vector{Int}, dom :: MathOptInterface.Reals) = nothing
function addvarconstr(m :: MosekModel, subj :: Vector{Int}, dom :: MathOptInterface.Zeros)
    bnd = zeros(Float64,length(subj))
    putvarboundslice(m.task,subj,fill(MSK_BK_FX,length(subj)),bnd,bnd)
end
addvarconstr(m :: MosekModel, subj :: Vector{Int}, dom :: MathOptInterface.Interval) = putvarbound(m.task,subj[1],MSK_BK_RA,dom.lower, dom.upper)
addvarconstr(m :: MosekModel, subj :: Vector{Int}, dom :: MathOptInterface.EqualTo) = putvarbound(m.task,subj[1],MSK_BK_FX,dom.value, dom.value)

function addvarconstr(m :: MosekModel, subj :: Vector{Int}, dom :: MathOptInterface.LessThan)
    bk,lo,up = getvarbound(m.task,subj[1])
    bk = if (bk == MSK_BK_FR) MSK_BK_UP else MSK_BK_RA end
    putvarbound(m.task,Int32(subj[1]),bk,lo,dom.upper)
end

function addvarconstr(m :: MosekModel, subj :: Vector{Int}, dom :: MathOptInterface.GreaterThan)
    bk,lo,up = getvarbound(m.task,subj[1])
    bk = if (bk == MSK_BK_FR) MSK_BK_LO else MSK_BK_RA end
    putvarbound(m.task,Int32(subj[1]),bk,dom.lower,up)
end

function addvarconstr(m :: MosekModel, subj :: Vector{Int}, dom :: MathOptInterface.Nonnegatives)
    bkx = zeros(Int32,length(subj))
    blx = zeros(Float64,length(subj))
    bux = zeros(Float64,length(subj))
    for (i,j) in enumerate(subj)
        bk,lo,up = getvarbound(m.task,j)
        bkx[i] = if (bk == MSK_BK_FR) MSK_BK_RA else MSK_BK_LO end
        bux[i] = up
    end
    putvarboundlist(m.task,subj,bkx,blx,bux)
end

function addvarconstr(m :: MosekModel, subj :: Vector{Int}, dom :: MathOptInterface.Nonpositives)
    bkx = zeros(Int32,length(subj))
    blx = zeros(Float64,length(subj))
    bux = zeros(Float64,length(subj))
    for (i,j) in enumerate(subj)
        bk,lo,up = getvarbound(m.task,j)
        bkx[i] = if (bk == MSK_BK_FR) MSK_BK_RA else MSK_BK_UP end
        blx[i] = lo
    end
    putvarboundlist(m.task,subj,bkx,blx,bux)
end

addvarconstr(m :: MosekModel, subj :: Vector{Int}, dom :: MathOptInterface.SecondOrderCone) = appendcone(m.task,MSK_CT_QUAD, 0.0, subj)
addvarconstr(m :: MosekModel, subj :: Vector{Int}, dom :: MathOptInterface.RotatedSecondOrderCone) = appendcone(m.task,MSK_CT_RQUAD, 0.0, subj)


function MathOptInterface.addconstraint!{D <: MathOptInterface.AbstractSet}(m :: MosekModel, xs :: MathOptInterface.SingleVariable, dom :: D)    
    subj = getindexes(m.x_block, ref2id(xs.variable))
    
    mask = domain_type_mask(dom)
    if mask & m.x_boundflags[subj[1]] != 0
        error("Cannot multiple bound sets of the same type to a variable")
    end

    if is_positivesemidefinite_set(dom)
        error("Invalid size for positive semidefnite set")
    elseif is_conic_set(dom)
        error("Invalid size for conic set")
    else
        id = allocatevarconstraints(m,1)

        xc_sub = getindexes(m.xc_block,id)
        m.xc_bounds[xc_sub[1]] = mask
        m.xc_idxs[xc_sub[1]]   = subj[1]

        addvarconstr(m,subj,dom)

        m.x_boundflags[subj[1]] .|= mask

        #id2cref{MathOptInterface.SingleVariable,D}(id)        
        conref = MathOptInterface.ConstraintReference{MathOptInterface.SingleVariable,D}(UInt64(id) << 1)
        select(m.constrmap,MathOptInterface.SingleVariable,D)[conref.value] = id
        conref
    end
end

function MathOptInterface.addconstraint!(m :: MosekModel, xs :: MathOptInterface.VectorOfVariables, dom :: MathOptInterface.PositiveSemidefiniteConeTriangle)
    subj = Vector{Int}(length(xs.variables))
    for i in 1:length(subj)
        getindexes(m.x_block, ref2id(xs.variables[i]),subj,i)
    end
    
    mask = domain_type_mask(dom)
    if any(mask .& m.x_boundflags[subj[1]] .> 0)
        error("Cannot multiple bound sets of the same type to a variable")
    end

    N = MathOptInterface.dimension(dom)

    if N < 3
        error("Invalid dimension for semidefinite constraint")
    end

    NN = (N*(N+1)) >> 1

    if length(subj) != NN
        error("Mismatching variable length for semidefinite constraint")
    end
    
    id = allocateconstraints(m,NN)

    appendbarvars(m.task,Int32(N))
    barvaridx = getnumbarvar(m.task)
    
    subi = getindexes(m.c_block,id)
    m.c_block_slack[subi] = -barvaridx

    subii32 = convert(Vector{Int32},subi)
    putaijlist(m.task,
               subii32,
               convert(Vector{Int32},subj),
               ones(Float64,NN))
    putconboundlist(m.task,subii32,fill(MSK_BK_FX,NN),zeros(Float64,NN),zeros(Float64,NN))
    idx = 1
    for j in 1:N
        for i in 1:N
            symmatidx = appendsparsesymmat(m.task,N,Int32[i],Int32[j],Float64[-1.0])
            putbaraij(m.task,subii32[idx],barvaridx,[symmatidx],Float64[1.0])
            idx += 1
        end
    end

    # HACK: We need to return a negative to indicate that this is
    # not, in fact, a real variable constraint, but rather a real
    # constraint, but to resolve return value at compile time we
    # need to disguise it as a variable constraint.
    #id2cref{MathOptInterface.VectorOfVariables,D}(-id)
    conref = MathOptInterface.ConstraintReference{MathOptInterface.VectorOfVariables,MathOptInterface.PositiveSemidefiniteConeTriangle}((UInt64(-id) << 1) | 1)
    select(m.constrmap,MathOptInterface.VectorOfVariables,MathOptInterface.PositiveSemidefiniteConeTriangle)[conref.value] = id
    
    conref 
end

    
function MathOptInterface.addconstraint!{D <: MathOptInterface.AbstractSet}(m :: MosekModel, xs :: MathOptInterface.VectorOfVariables, dom :: D)
    subj = Vector{Int}(length(xs.variables))
    for i in 1:length(subj)
        getindexes(m.x_block, ref2id(xs.variables[i]),subj,i)
    end
    
    mask = domain_type_mask(dom)
    if any(mask .& m.x_boundflags[subj[1]] .> 0)
        error("Cannot multiple bound sets of the same type to a variable")
    end

    N = MathOptInterface.dimension(dom)
    allocated = ensurefree(m.xc_block,N)
    id = newblock(m.xc_block,N)
    if allocated > 0
        append!(m.xc_bounds, zeros(Int,allocated))
        append!(m.xc_idxs,   zeros(Int,allocated))
    end

    xc_sub = getindexes(m.xc_block,id)
    m.xc_bounds[xc_sub] = mask
    m.xc_idxs[xc_sub] = subj

    addvarconstr(m,subj,dom)

    m.x_boundflags[subj] .|= mask
    
    conref = MathOptInterface.ConstraintReference{MathOptInterface.VectorOfVariables,D}(UInt64(id) << 1)
    select(m.constrmap,MathOptInterface.VectorOfVariables,D)[conref.value] = id
    conref
end

################################################################################
################################################################################


# Put the linear left-hand side
function addlhsblock!(m        :: MosekModel,
                      conid    :: Int,
                      conidxs  :: Vector{Int},
                      varidxs  :: Vector{MathOptInterface.VariableReference},
                      cofs     :: Vector{Float64})
    consubi = getindexes(m.c_block,conid)
    subj = Array{Int}(length(varidxs))
    for i in 1:length(subj)
        getindexes(m.x_block,ref2id(varidxs[i]),subj,i)
    end

    N = length(consubi)    

    At = sparse(subj, conidxs, cofs, getnumvar(m.task), N)
    putarowlist(m.task,convert(Vector{Int32},consubi),At)
end


addbound!(m :: MosekModel, conid :: Int, conidxs :: Vector{Int}, constant :: Vector{Float64}, dom :: MathOptInterface.Reals)                = putconboundlist(m.task,convert(Vector{Int32},conidxs),fill(MSK_BK_FR,MathOptInterface.dimension(dom)),-constant,-constant)
addbound!(m :: MosekModel, conid :: Int, conidxs :: Vector{Int}, constant :: Vector{Float64}, dom :: MathOptInterface.Zeros)                = putconboundlist(m.task,convert(Vector{Int32},conidxs),fill(MSK_BK_FX,MathOptInterface.dimension(dom)),-constant,-constant)
addbound!(m :: MosekModel, conid :: Int, conidxs :: Vector{Int}, constant :: Vector{Float64}, dom :: MathOptInterface.Nonnegatives)         = putconboundlist(m.task,convert(Vector{Int32},conidxs),fill(MSK_BK_LO,MathOptInterface.dimension(dom)),-constant,-constant)
addbound!(m :: MosekModel, conid :: Int, conidxs :: Vector{Int}, constant :: Vector{Float64}, dom :: MathOptInterface.Nonpositives)         = putconboundlist(m.task,convert(Vector{Int32},conidxs),fill(MSK_BK_UP,MathOptInterface.dimension(dom)),-constant,-constant)
addbound!(m :: MosekModel, conid :: Int, conidxs :: Vector{Int}, constant :: Vector{Float64}, dom :: MathOptInterface.GreaterThan{Float64}) = putconbound(m.task,Int32(conidxs[1]),MSK_BK_LO,dom.lower-constant[1],dom.lower-constant[1])
addbound!(m :: MosekModel, conid :: Int, conidxs :: Vector{Int}, constant :: Vector{Float64}, dom :: MathOptInterface.LessThan{Float64})    = putconbound(m.task,Int32(conidxs[1]),MSK_BK_UP,dom.upper-constant[1],dom.upper-constant[1])
addbound!(m :: MosekModel, conid :: Int, conidxs :: Vector{Int}, constant :: Vector{Float64}, dom :: MathOptInterface.EqualTo{Float64})     = putconbound(m.task,Int32(conidxs[1]),MSK_BK_FX,dom.value-constant[1],dom.value-constant[1])
addbound!(m :: MosekModel, conid :: Int, conidxs :: Vector{Int}, constant :: Vector{Float64}, dom :: MathOptInterface.Interval{Float64})    = putconbound(m.task,Int32(conidxs[1]),MSK_BK_RA,dom.lower-constant[1],dom.upper-constant[1])

function add_slack!(m :: MosekModel, conidxs :: Vector{Int}, constant :: Vector{Float64}, N :: Int)
    ensurefree(m.x_block,N)
    varid = newblock(m.x_block,N)

    numvar = getnumvar(m.task)
    if length(m.x_block) > numvar
        appendvars(m.task, length(m.x_block) - numvar)
    end

    subj = getindexes(m.x_block,varid)

    putaijlist(m.task,conidxs,subj,-ones(Float64,N))
    putvarboundlist(m.task,subj,fill(MSK_BK_FR,N),zeros(Float64,N),zeros(Float64,N))
    
    putconboundlist(m.task,convert(Vector{Int32},conidxs),fill(MSK_BK_FX,N),-constant,-constant)

    varid,subj
end

function addbound!(m :: MosekModel, conid :: Int, conidxs :: Vector{Int}, constant :: Vector{Float64}, dom :: MathOptInterface.SecondOrderCone)
    N = MathOptInterface.dimension(dom)

    varid,subj = add_slack!(m,conidxs,constant,N)
    
    appendcone(m.task,MSK_CT_QUAD,0.0,subj)
    numcone = getnumcone(m.task)
    m.c_block_slack[conid]   = varid
end

function addbound!(m :: MosekModel, conid :: Int, conidxs :: Vector{Int}, constant :: Vector{Float64}, dom :: MathOptInterface.RotatedSecondOrderCone)
    N = MathOptInterface.dimension(dom)

    varid,subj = add_slack!(m,conidxs,constant,N)
    
    appendcone(m.task,MSK_CT_RQUAD,0.0,subj)
    numcone = getnumcone(m.task)
    m.c_block_slack[conid]   = varid
end

function addbound!(m :: MosekModel, conid :: Int, conidxs :: Vector{Int}, constant :: Vector{Float64}, dom :: MathOptInterface.PositiveSemidefiniteConeTriangle)
    dim = MathOptInterface.dimension(dom)

    appendbarvars(m.task,Int32[dim])
    barvaridx = getnumbarvar(m.task)

    idx = 1
    for j in 1:dim
        for i in j:dim
            matrixid = appendsparsesymmat(m.task,Int32(dim), Int32[i], Int32[j], Float64[-1.0])
            putbaraij(m.task, conidxs[idx], barvaridx, Int[matrixid], Float64[1.0])
            idx += 1
        end
    end

    putconboundlist(m.task,convert(Vector{Int32},conidxs),fill(MSK_BK_FX,length(constant)),-constant,-constant)

    m.c_block_slack[conid] = -barvaridx
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









function allocateconstraints(m :: MosekModel,
                             N :: Int)
    numcon = getnumcon(m.task)
    alloced = ensurefree(m.c_block,N)
    if alloced > 0
        appendcons(m.task, alloced)
        append!(m.c_constant, zeros(Float64,alloced))
        append!(m.c_block_slack, zeros(Float64,alloced))
    end
    newblock(m.c_block,N)
end


function allocatevarconstraints(m :: MosekModel,
                                N :: Int)
    alloced = ensurefree(m.xc_block,N)
    if alloced > 0
        append!(m.xc_bounds,zeros(Float64,alloced))
        append!(m.xc_idxs, zeros(Float64,alloced))
    end
    newblock(m.xc_block,N)
end

function allocatevariable(m :: MosekModel,N :: Int)
    numvar = getnumvar(m.task)
    alloced = ensurefree(m.x_block,N)
    if alloced > 0
        appendvars(m.task, length(m.x_block) - numvar)
        append!(m.x_boundflags, zeros(Int,length(m.x_block) - numvar))
        append!(m.x_numxc, zeros(Int,length(m.x_block) - numvar))
    end
    newblock(m.x_block,N)
end

isvalid(m::MosekModel, ref::MathOptInterface.ConstraintReference) = allocated(m.c_block,BlockId(ref.value))
