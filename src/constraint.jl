const VectorCone = Union{MathOptInterface.SecondOrderCone,
                         MathOptInterface.RotatedSecondOrderCone,
                         MathOptInterface.PowerCone,
                         MathOptInterface.DualPowerCone,
                         MathOptInterface.ExponentialCone,
                         MathOptInterface.DualExponentialCone}
const PositiveSemidefiniteCone = Union{MathOptInterface.PositiveSemidefiniteConeTriangle,
                                       MathOptInterface.PositiveSemidefiniteConeScaled}
const LinearFunction = Union{MathOptInterface.SingleVariable,
                             MathOptInterface.VectorOfVariables,
                             MathOptInterface.ScalarAffineFunction,
                             MathOptInterface.VectorAffineFunction}
const AffineFunction = Union{MathOptInterface.ScalarAffineFunction,
                             MathOptInterface.VectorAffineFunction}
const VariableFunction = Union{MathOptInterface.ScalarAffineFunction,
                               MathOptInterface.VectorAffineFunction}

const ScalarLinearDomain = Union{MathOptInterface.LessThan{Float64},
                                 MathOptInterface.GreaterThan{Float64},
                                 MathOptInterface.EqualTo{Float64},
                                 MathOptInterface.Interval{Float64} }
const VectorLinearDomain = Union{MathOptInterface.Nonpositives,
                                 MathOptInterface.Nonnegatives,
                                 MathOptInterface.Reals,
                                 MathOptInterface.Zeros}
const LinearDomain = Union{ScalarLinearDomain,VectorLinearDomain}
################################################################################
# ADD CONSTRAINT ###############################################################
################################################################################

function MathOptInterface.addconstraint!(
    m   :: MosekModel,
    axb :: MathOptInterface.ScalarAffineFunction{Float64},
    dom :: D) where {D <: MathOptInterface.AbstractScalarSet}
    
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

function MathOptInterface.addconstraint!(
    m   :: MosekModel,
    axb :: MathOptInterface.VectorAffineFunction{Float64},
    dom :: D) where { D <: MathOptInterface.AbstractVectorSet }
    
    N = MathOptInterface.dimension(dom)
    conid = allocateconstraints(m,N)
    addlhsblock!(m,
                 conid,
                 axb.outputindex,
                 axb.variables,
                 axb.coefficients)
    conidxs = getindexes(m.c_block,conid)
    m.c_constant[conidxs] = axb.constant

    m.c_block_slack[conid]
    
    addbound!(m,conid,conidxs,axb.constant,dom)
    conref = MathOptInterface.ConstraintReference{MathOptInterface.VectorAffineFunction{Float64},D}(UInt64(conid) << 1)
    select(m.constrmap,MathOptInterface.VectorAffineFunction{Float64},D)[conref.value] = conid
    
    conref
end

dim(dom :: MathOptInterface.PositiveSemidefiniteConeTriangle) = floor(Int,-.5 + sqrt(.25+2*MathOptInterface.dimension(dom)))
dim(dom :: MathOptInterface.PositiveSemidefiniteConeScaled) = floor(Int,-.5 + sqrt(.25+2*MathOptInterface.dimension(dom)))

function MathOptInterface.addconstraint!(m   :: MosekModel,
                                         axb :: MathOptInterface.VectorAffineFunction{Float64},
                                         dom :: PSDCone) where { PSDCone <: PositiveSemidefiniteCone }
    M = MathOptInterface.dimension(dom)
    N = dim(dom)
    
    conid = allocateconstraints(m,M)
    addlhsblock!(m,
                 conid,
                 axb.outputindex,
                 axb.variables,
                 axb.coefficients)
    conidxs = getindexes(m.c_block,conid)
    m.c_constant[conidxs] = axb.constant

    addbound!(m,conid,conidxs,axb.constant,dom)
    conref = MathOptInterface.ConstraintReference{MathOptInterface.VectorAffineFunction{Float64},PSDCone}(UInt64(conid) << 1)
    select(m.constrmap,MathOptInterface.VectorAffineFunction{Float64},PSDCone)[conref] = conid
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
is_conic_set(dom :: VectorCone) = true
is_conic_set(dom) = false

addvarconstr(m :: MosekModel, subj :: Vector{Int}, dom :: MathOptInterface.Reals) = nothing
function addvarconstr(m :: MosekModel, subj :: Vector{Int}, dom :: MathOptInterface.Zeros)
    bnd = zeros(Float64,length(subj))
    putvarboundslice(m.task,subj,fill(MSK_BK_FX,length(subj)),bnd,bnd)
end
addvarconstr(m :: MosekModel, subj :: Vector{Int}, dom :: MathOptInterface.Interval) = putvarbound(m.task,subj[1],MSK_BK_RA,dom.lower, dom.upper)
addvarconstr(m :: MosekModel, subj :: Vector{Int}, dom :: MathOptInterface.EqualTo) = putvarbound(m.task,subj[1],MSK_BK_FX,dom.value, dom.value)
addvarconstr(m :: MosekModel, subj :: Vector{Int}, dom :: MathOptInterface.Integer) = putvartype(m.task,subj[1],MSK_VAR_TYPE_INT)
function addvarconstr(m :: MosekModel, subj :: Vector{Int}, dom :: MathOptInterface.ZeroOne)
    putvartype(m.task,subj[1],MSK_VAR_TYPE_INT)
    putvarbound(m.task,subj[1],MSK_BK_RA,0.0,1.0)
end

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


abstractset2ct(dom::MathOptInterface.SecondOrderCone) = MSK_CT_QUAD
abstractset2ct(dom::MathOptInterface.RotatedSecondOrderCone) = MSK_CT_RQUAD
function MathOptInterface.addconstraint!(
    m   :: MosekModel,
    xs  :: MathOptInterface.SingleVariable,
    dom :: D) where {D <: MathOptInterface.AbstractScalarSet}
    
    subj = getindexes(m.x_block, ref2id(xs.variable))
    
    mask = domain_type_mask(dom)
    if mask & m.x_boundflags[subj[1]] != 0
        error("Cannot put multiple bound sets of the same type on a variable")
    end

    xcid = allocatevarconstraints(m,1)

    xc_sub = getindexes(m.xc_block,xcid)[1]
    
    m.xc_bounds[xcid]  = mask
    m.xc_idxs[xc_sub] = subj[1]

    addvarconstr(m,subj,dom)

    m.x_boundflags[subj[1]] .|= mask

    conref = MathOptInterface.ConstraintReference{MathOptInterface.SingleVariable,D}(UInt64(xcid) << 1)
    
    select(m.constrmap,MathOptInterface.SingleVariable,D)[conref.value] = xcid
    conref
end

function MathOptInterface.addconstraint!(m   :: MosekModel,
                                         xs  :: MathOptInterface.VectorOfVariables,
                                         dom :: D) where { D <: PositiveSemidefiniteCone }
    subj = Vector{Int}(length(xs.variables))
    for i in 1:length(subj)
        getindexes(m.x_block, ref2id(xs.variables[i]),subj,i)
    end
    
    mask = domain_type_mask(dom)
    if any(mask .& m.x_boundflags[subj[1]] .> 0)
        error("Cannot put multiple bound sets of the same type to a variable")
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

    addbound!(m,id,subi,zeros(Float64,NN),dom)
    
    #putconboundlist(m.task,subii32,fill(MSK_BK_FX,NN),zeros(Float64,NN),zeros(Float64,NN))
    #idx = 1
    #for j in 1:N
    #    for i in 1:N
    #        symmatidx = appendsparsesymmat(m.task,N,Int32[i],Int32[j],Float64[-1.0])
    #        putbaraij(m.task,subii32[idx],barvaridx,[symmatidx],Float64[1.0])
    #        idx += 1
    #    end
    #end

    # HACK: We need to return a negative to indicate that this is
    # not, in fact, a real variable constraint, but rather a real
    # constraint, but to resolve return value at compile time we
    # need to disguise it as a variable constraint.
    #id2cref{MathOptInterface.VectorOfVariables,D}(-id)
    conref = MathOptInterface.ConstraintReference{MathOptInterface.VectorOfVariables,MathOptInterface.PositiveSemidefiniteConeTriangle}((UInt64(-id) << 1) | 1)
    select(m.constrmap,MathOptInterface.VectorOfVariables,MathOptInterface.PositiveSemidefiniteConeTriangle)[conref.value] = id
    
    conref 
end



function aux_setvardom(m :: MosekModel,
                       xcid :: Int,
                       subj :: Vector{Int},
                       dom :: D) where { D <: VectorCone }
    appendcone(m.task,abstractset2ct(dom),  0.0, subj)
    coneidx = getnumcone(m.task)    
    m.conecounter += 1
    putconename(m.task,coneidx,"$(m.conecounter)")
    m.xc_coneid[xcid] = m.conecounter
end
aux_setvardom(m :: MosekModel, xcid :: Int, subj :: Vector{Int},dom :: D) where {D <: MathOptInterface.AbstractSet} = addvarconstr(m,subj,dom)

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
    xcid = allocatevarconstraint(m,N)
    xc_sub = getindexes(m.xc_block,xcid)
    
    m.xc_bounds[xcid]  = mask
    m.xc_idxs[xc_sub] = subj

    aux_setvardom(m,xcid,subj,dom)

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

function addbound!(m :: MosekModel, conid :: Int, conidxs :: Vector{Int}, constant :: Vector{Float64}, dom :: D) where { D <: VectorCone }
    N = MathOptInterface.dimension(dom)
    nalloc = ensurefree(m.x_block,N)
    
    varid = newblock(m.x_block,N)
    numvar = getnumvar(m.task)
    if length(m.x_block) > numvar
        appendvars(m.task, length(m.x_block) - numvar)
    end
    subj = getindexes(m.x_block,varid)

    putaijlist(m.task,conidxs,subj,-ones(Float64,N))
    putvarboundlist(m.task,subj,fill(MSK_BK_FR,N),zeros(Float64,N),zeros(Float64,N))
    putconboundlist(m.task,convert(Vector{Int32},conidxs),fill(MSK_BK_FX,N),-constant,-constant)

    m.c_block_slack[conid] = varid

    appendcone(m.task,abstractset2ct(dom),0.0,subj)
    coneidx = getnumcone(m.task)
    m.conecounter += 1
    putconename(m.task,coneidx,"$(m.conecounter)")
    m.c_coneid[conid] = m.conecounter
end

function addbound!(m :: MosekModel, conid :: Int, conidxs :: Vector{Int}, constant :: Vector{Float64}, dom :: MathOptInterface.PositiveSemidefiniteConeTriangle)
    dim = dim(dom)
    n = MathOptInterface.dimension(dom)
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

function addbound!(m :: MosekModel, conid :: Int, conidxs :: Vector{Int}, constant :: Vector{Float64}, dom :: MathOptInterface.PositiveSemidefiniteConeScaled)
    dim = dim(dom)
    n = MathOptInterface.dimension(dom)

    appendbarvars(m.task,Int32[dim])
    barvaridx = getnumbarvar(m.task)

    idx = 1
    for j in 1:dim
        for i in j:dim
            if i != j
                matrixid = appendsparsesymmat(m.task,Int32(dim), Int32[i], Int32[j], Float64[- sqrt(2.0)])
            else
                matrixid = appendsparsesymmat(m.task,Int32(dim), Int32[i], Int32[j], Float64[-1.0])
            end
            putbaraij(m.task, conidxs[idx], barvaridx, Int[matrixid], Float64[1.0])
            idx += 1
        end
    end

    putconboundlist(m.task,convert(Vector{Int32},conidxs),fill(MSK_BK_FX,length(constant)),-constant,-constant)

    m.c_block_slack[conid] = -barvaridx
end



################################################################################
##  MODIFY #####################################################################
################################################################################

MathOptInterface.canmodifyconstraint(m::MosekModel, c::MathOptInterface.ConstraintReference{F,D}, dom) where { F <: LinearFunction, D <: VectorCone } = false
MathOptInterface.canmodifyconstraint(m::MosekModel, c::MathOptInterface.ConstraintReference{F,D}, dom) where { F <: LinearFunction, D <: PositiveSemidefiniteCone } = false
MathOptInterface.canmodifyconstraint(m::MosekModel, c::MathOptInterface.ConstraintReference{F,D}, dom::D) where { F <: LinearFunction, D <: LinearDomain } = haskey(select(m.constrmap,F,D),c.value)
MathOptInterface.canmodifyconstraint(m::MosekModel, c::MathOptInterface.ConstraintReference{F,D}, func::F) where { F <: Union{MathOptInterface.SingleVariable,MathOptInterface.ScalarAffineFunction}, D <: ScalarLinearDomain } = haskey(select(m.constrmap,F,D),c.value)
MathOptInterface.canmodifyconstraint(m::MosekModel, c::MathOptInterface.ConstraintReference{F,D}, change::MathOptInterface.AbstractFunctionModification) where { F <: AffineFunction, D <: MathOptInterface.AbstractSet } = haskey(select(m.constrmap,F,D),c.value)

chgbound(bl::Float64,bu::Float64,k::Float64,dom :: MathOptInterface.LessThan{Float64})    = bl,dom.upper-k
chgbound(bl::Float64,bu::Float64,k::Float64,dom :: MathOptInterface.GreaterThan{Float64}) = dom.lower-k,bu
chgbound(bl::Float64,bu::Float64,k::Float64,dom :: MathOptInterface.EqualTo{Float64})     = dom.value-k,dom.value-k
chgbound(bl::Float64,bu::Float64,k::Float64,dom :: MathOptInterface.Interval{Float64})    = dom.lower-k,dom.upper-k

function MathOptInterface.modifyconstraint!(m::MosekModel,
                                            xcref::MathOptInterface.ConstraintReference{F,D},
                                            dom::D) where { F    <: MathOptInterface.SingleVariable,
                                                            D    <: ScalarLinearDomain }
    xcid = ref2id(xcref)
    j = m.xc_idxs[getindexes(m.xc_block,xcid)[1]]
    bk,bl,bu = getvarbound(m.task,j)
    bl,bu = chgbound(bl,bu,0.0,dom)
    putvarbound(m.task,j,bk,bl,bu)    
end


function MathOptInterface.modifyconstraint!(m::MosekModel,
                                            cref::MathOptInterface.ConstraintReference{F,D},
                                            dom::D) where { F    <: MathOptInterface.ScalarAffineFunction,
                                                            D    <: ScalarLinearDomain }
    cid = ref2id(cref)
    i = getindexes(m.c_block,cid)[1] # since we are in a scalar domain
    bk,bl,bu = getconbound(m.task,i)
    bl,bu = chgbound(bl,bu,0.0,dom)
    putconbound(m.task,i,bk,bl,bu)
end

function MathOptInterface.modifyconstraint!(m   ::MosekModel,
                                            c   ::MathOptInterface.ConstraintReference{MathOptInterface.ScalarAffineFunction{Float64},D},
                                            func::MathOptInterface.ScalarConstantChange{Float64}) where {D <: MathOptInterface.AbstractSet}

    cid = ref2id(c)

    i = getindexes(m.c_block, cid)[1]
    bk,bl,bu = getconbound(m.task,i)
    bl += m.c_constant[i] - func.new_constant
    bu += m.c_constant[i] - func.new_constant
    m.c_constant[i] = func.new_constant
    putconbound(m.task,i,bk,bl,bu)
end

function MathOptInterface.modifyconstraint!(m   ::MosekModel,
                                            c   ::MathOptInterface.ConstraintReference{MathOptInterface.ScalarAffineFunction{Float64},D},
                                            func::MathOptInterface.ScalarCoefficientChange{Float64}) where {D <: MathOptInterface.AbstractSet}
    cid = ref2id(c)
    xid = ref2id(func.variable)

    i = getindexes(m.c_block,cid)[1]
    j = getindexes(m.x_block,xid)[1]

    putaij(m.task,i,j,func.new_coefficient)
end

function MathOptInterface.modifyconstraint!(m::MosekModel,
                                            c::MathOptInterface.ConstraintReference{MathOptInterface.VectorAffineFunction{Float64},D},
                                            func::MathOptInterface.VectorConstantChange{Float64}) where {D <: MathOptInterface.AbstractSet}
    cid = ref2id(c)
    assert(cid > 0)

    subi = getindexes(m.c_block, cid)
    bk = Vector{Int32}(length(subi))
    bl = Vector{Float64}(length(subi))
    bu = Vector{Float64}(length(subi))

    bk,bl,bu = getconboundlist(m.task,convert(Vector{Int32},subi))
    bl += m.c_constant[subi] - func.new_constant
    bu += m.c_constant[subi] - func.new_constant
    m.c_constant[subi] = func.new_constant
    
    putconboundlist(m.task,convert(Vector{Int32},subi),bk,bl,bu)
end

function MathOptInterface.modifyconstraint!(m::MosekModel,
                                            c::MathOptInterface.ConstraintReference{MathOptInterface.VectorAffineFunction{Float64},D},
                                            func::MathOptInterface.MultirowChange{Float64}) where {D <: MathOptInterface.AbstractSet}
    cid = ref2id(c)
    assert(cid > 0)

    subi = getindexes(m.c_block, cid)[func.rows]
    xid = ref2id(func.variable)
    j = getindexes(m.x_block,xid)[1]

    putaijlist(m.task,convert(Vector{Int32},subi),fill(j,length(subi)),fund.new_coefficients)
end



MathOptInterface.cantransformconstraint(m::MosekModel, c::MathOptInterface.ConstraintReference{F,D1}, newdom::D2) where {F <: VariableFunction, D1, D2 } = false
MathOptInterface.cantransformconstraint(m::MosekModel, c::MathOptInterface.ConstraintReference{MathOptInterface.VectorAffineFunction,D1}, newdom::D2) where {D1 <: VectorLinearDomain, D2 <: VectorLinearDomain} = false
function MathOptInterface.cantransformconstraint(m::MosekModel,
                                                 cref::MathOptInterface.ConstraintReference{MathOptInterface.ScalarAffineFunction{Float64},D1},
                                                 newdom::D2) where {D1 <: ScalarLinearDomain,
                                                                    D2 <: ScalarLinearDomain}
    haskey(select(m.constrmap,MathOptInterface.ScalarAffineFunction{Float64},D1),cref.value)
end

function MathOptInterface.transformconstraint!(m::MosekModel,
                                               cref::MathOptInterface.ConstraintReference{F,D},
                                               newdom::D) where {F  <: MathOptInterface.AbstractFunction, D <: MathOptInterface.AbstractSet}
    MathOptInterface.modifyconstraint(m,cref,newdom)
    cref
end

function MathOptInterface.transformconstraint!(m::MosekModel,
                                               cref::MathOptInterface.ConstraintReference{MathOptInterface.ScalarAffineFunction{Float64},D1},
                                               newdom::D2) where {D1 <: ScalarLinearDomain,
                                                                  D2 <: ScalarLinearDomain}
    const F = MathOptInterface.ScalarAffineFunction{Float64}
    
    cid = ref2id(cref)

    subi = getindexes(m.c_block,cid)

    addbound!(m,cid,subi,m.c_constant[subi], newdom)
    
    newcref = MathOptInterface.ConstraintReference{F,D2}(UInt64(cid) << 1)    
    delete!(select(m.constrmap,F,D1), cref.value)
    select(m.constrmap,F, D2)[newcref.value] = cid
    newcref
end

function MathOptInterface.transformconstraint!(m::MosekModel,
                                               cref::MathOptInterface.ConstraintReference{MathOptInterface.VectorAffineFunction{Float64},D1},
                                               newdom::D2) where {D1 <: VectorLinearDomain,
                                                                  D2 <: VectorLinearDomain}
    const F = MathOptInterface.VectorAffineFunction{Float64}
    
    cid = ref2id(cref)

    subi = getindexes(m.c_block,cid)

    addbound!(m,cid,subi,m.c_constant[subi], newdom)
    
    newcref = MathOptInterface.ConstraintReference{F,D2}(UInt64(cid) << 1)    
    delete!(select(m.constrmap,F,D1), cref.value)
    select(m.constrmap,F,D2)[newcref.value] = cid
    newcref
end

#MathOptInterface.transformconstraint!(m::MosekModel, c::MathOptInterface.ConstraintReference{F,D1}, newdom::D2) where {F <: MathOptInterface.VectorAffineFunction , D1 <: VectorLinearDomain, D2 <: VectorLinearDomain} = false

################################################################################
##  DELETE #####################################################################
################################################################################

MathOptInterface.candelete(
    m   ::MosekModel,
    cref::MathOptInterface.ConstraintReference{F,D}) where {F <: Union{MathOptInterface.ScalarAffineFunction,
                                                                       MathOptInterface.VectorAffineFunction,
                                                                       MathOptInterface.SingleVariable,
                                                                       MathOptInterface.VectorOfVariables},
                                                            D <: Union{MathOptInterface.LessThan,
                                                                       MathOptInterface.GreaterThan,
                                                                       MathOptInterface.EqualTo,
                                                                       MathOptInterface.Interval,
                                                                       MathOptInterface.Zeros,
                                                                       MathOptInterface.Nonpositives,
                                                                       MathOptInterface.Nonnegatives,
                                                                       MathOptInterface.Reals}} = MathOptInterface.isvalid(m,cref)

MathOptInterface.candelete(
    m   ::MosekModel,
    cref::MathOptInterface.ConstraintReference{F,D}) where {F <: MathOptInterface.AbstractFunction,
                                                            D <: Union{MathOptInterface.SecondOrderCone,
                                                                       MathOptInterface.RotatedSecondOrderCone,
                                                                       MathOptInterface.ExponentialCone,
                                                                       MathOptInterface.DualExponentialCone,
                                                                       MathOptInterface.PowerCone,
                                                                       MathOptInterface.DualPowerCone}} = false

MathOptInterface.candelete(
    m   ::MosekModel,
    cref::MathOptInterface.ConstraintReference{F,D}) where {F <: Union{MathOptInterface.SingleVariable,
                                                                       MathOptInterface.VectorOfVariables},
                                                            D <: Union{MathOptInterface.SecondOrderCone,
                                                                       MathOptInterface.RotatedSecondOrderCone,
                                                                       MathOptInterface.ExponentialCone,
                                                                       MathOptInterface.DualExponentialCone,
                                                                       MathOptInterface.PowerCone,
                                                                       MathOptInterface.DualPowerCone}} = true


function Base.delete!(
    m::MosekModel,
    cref::MathOptInterface.ConstraintReference{F,D}) where {F <: Union{MathOptInterface.ScalarAffineFunction,
                                                                       MathOptInterface.VectorAffineFunction},
                                                            D <: Union{MathOptInterface.LessThan,
                                                                       MathOptInterface.GreaterThan,
                                                                       MathOptInterface.EqualTo,
                                                                       MathOptInterface.Interval,
                                                                       MathOptInterface.Zeros,
                                                                       MathOptInterface.Nonpositives,
                                                                       MathOptInterface.Nonnegatives,
                                                                       MathOptInterface.Reals}}

    delete!(select(m.constrmap, F, D), cref.value)
    
    cid = ref2id(cref)
    subi = getindexes(m.c_block,cid)

    n = length(subi)
    subi_i32 = convert(Vector{Int32},subi)
    ptr = fill(Int64(0),n)
    putarowlist(m.task,subi_i32,ptr,ptr,Int32[],Float64[])
    b = fill(0.0,n)
    putconboundlist(m.task,subi_i32,fill(MSK_BK_FX,n),b,b)

    m.c_constant[subi] = 0.0
    deleteblock(m.c_block,cid)
end

function Base.delete!(
    m::MosekModel,
    cref::MathOptInterface.ConstraintReference{F,D}) where {F <: Union{MathOptInterface.VectorAffineFunction},
                                                            D <: Union{MathOptInterface.SecondOrderCone,
                                                                       MathOptInterface.RotatedSecondOrderCone,
                                                                       MathOptInterface.ExponentialCone,
                                                                       MathOptInterface.DualExponentialCone,
                                                                       MathOptInterface.PowerCone,
                                                                       MathOptInterface.DualPowerCone}}
    
    delete!(select(m.constrmap, F, D), cref.value)
    xcid = ref2id(cref)
    sub = getindexes(m.xc_block,xcid)

    subj = [ getindexes(m.x_block,i)[1] for i in sub ]
    N = length(subj)
    m.x_boundflags[subj] .&= ~m.xc_bounds[xcid]    
    asgn,coneidx = getconenameindex(m.task,"$(m.xc_coneid[xcid])")
    m.xc_coneid[xcid] = 0
    removecone(m.task,coneidx)

    m.x_numxc[subj]  -= 1
    m.xc_idxs[sub]    = 0
    m.xc_bounds[xcid] = 0
    
    deleteblock(m.xc_block,xcid)
end


function Base.delete!(
    m::MosekModel,
    cref::MathOptInterface.ConstraintReference{F,D}) where {F <: Union{MathOptInterface.SingleVariable,
                                                                       MathOptInterface.VectorOfVariables},
                                                            D <: Union{MathOptInterface.LessThan,
                                                                       MathOptInterface.GreaterThan,
                                                                       MathOptInterface.EqualTo,
                                                                       MathOptInterface.Interval,
                                                                       MathOptInterface.Zeros,
                                                                       MathOptInterface.Nonpositives,
                                                                       MathOptInterface.Nonnegatives,
                                                                       MathOptInterface.Reals,
                                                                       MathOptInterface.ZeroOne,
                                                                       MathOptInterface.Integer}}
    delete!(select(m.constrmap, F, D), cref.value)
    
    xcid = ref2id(cref)
    sub = getindexes(m.xc_block,xcid)

    subj = [ getindexes(m.x_block,i)[1] for i in sub ]
    N = length(subj)

    m.x_boundflags[subj] .&= ~m.xc_bounds[xcid]
    if m.xc_bounds[xcid] & boundflag_int != 0
        for i in 1:length(subj)
            putvartype(m.task,subj[i],MSK_VAR_TYPE_CONT)
        end
    end
    
    if m.xc_bounds[xcid] & boundflag_lower != 0 && m.xc_bounds[xcid] & boundflag_upper != 0
        bnd = fill(0.0, length(N))
        putvarboundlist(m.task,convert(Vector{Int32},subj),fill(MSK_BK_FR,N),bnd,bnd)
    elseif m.xc_bounds[xcid] & boundflag_lower != 0
        bnd = fill(0.0, length(N))
        bk,bl,bu = getvarboundlist(m.task, convert(Vector{Int32},subj))

        for i in 1:N
            if MSK_BK_RA == bk[i] || MSK_BK_FX == bk[i]
                bk[i] = MSK_BK_UP
            else
                bk[i] = MSK_BK_FR
            end
        end

        putvarboundlist(m.task,convert(Vector{Int32},subj),bk,bl,bu)
    elseif m.xc_bounds[xcid] & boundflag_upper != 0
        bnd = fill(0.0, length(N))
        bk,bl,bu = getvarboundlist(m.task, subj)

        for i in 1:N
            if MSK_BK_RA == bk[i] || MSK_BK_FX == bk[i]
                bk[i] = MSK_BK_LO
            else
                bk[i] = MSK_BK_FR
            end
        end
        
        putvarboundlist(m.task,convert(Vector{Int32},subj),bk,bl,bu)
    else
        assert(false)
        # should not happen
    end

    m.x_numxc[subj] -= 1
    m.xc_idxs[sub] = 0
    m.xc_bounds[xcid] = 0
    
    deleteblock(m.xc_block,xcid)
end






################################################################################
################################################################################
################################################################################

function allocateconstraints(m :: MosekModel,
                             N :: Int)
    numcon = getnumcon(m.task)
    alloced = ensurefree(m.c_block,N)
    id = newblock(m.c_block,N)

    M = numblocks(m.c_block) - length(m.c_block_slack)
    if alloced > 0
        appendcons(m.task, alloced)
        append!(m.c_constant, zeros(Float64,alloced))
    end
    if M > 0
        append!(m.c_block_slack, zeros(Float64,M))
        append!(m.c_coneid, zeros(Float64,M))
    end
    id
end


function allocatevarconstraints(m :: MosekModel,
                                N :: Int)
    nalloc = ensurefree(m.xc_block,N)
    id = newblock(m.xc_block,N)

    M = numblocks(m.xc_block) - length(m.xc_bounds)
    if M > 0
        append!(m.xc_bounds,zeros(Float64,M))
        append!(m.xc_coneid,zeros(Float64,M))
    end
    if nalloc > 0
        append!(m.xc_idxs, zeros(Float64,nalloc))
    end

    id
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

MathOptInterface.isvalid(m::MosekModel, ref::MathOptInterface.ConstraintReference{F,D}) where { F,D } = haskey(select(m.constrmap,F,D),ref.value)
MathOptInterface.isvalid(m::MosekModel, ref::MathOptInterface.VariableReference) = allocated(m.x_block,ref2id(ref))


function getvarboundlist(t::Mosek.Task, subj :: Vector{Int32})
    n = length(subj)
    bk = Vector{Boundkey}(n)
    bl = Vector{Float64}(n)
    bu = Vector{Float64}(n)
    for i in 1:n
        bki,bli,bui = getvarbound(t,subj[i])
        bk[i] = bki
        bl[i] = bli
        bu[i] = bui
    end
    bk,bl,bu
end

function getconboundlist(t::Mosek.Task, subj :: Vector{Int32})
    n = length(subj)
    bk = Vector{Boundkey}(n)
    bl = Vector{Float64}(n)
    bu = Vector{Float64}(n)
    for i in 1:n
        bki,bli,bui = getconbound(t,subj[i])
        bk[i] = bki
        bl[i] = bli
        bu[i] = bui
    end
    bk,bl,bu
end
