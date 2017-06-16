include("LinkedInts.lj")

mosek_block_type_unallocated = 0
mosek_block_type_zero   = 1
mosek_block_type_nonneg = 2
mosek_block_type_nonpos = 3
mosek_block_type_range  = 4
mosek_block_type_qcone  = 5
mosek_block_type_rqcone = 6
mosek_block_type_psd    = 7
mosek_block_type_integer = 8

"""
    MosekMathProgModel <: MathProgBase.AbstractMathModel

Linear variables and constraint can be deleted. MOSEK does not support
deleting PSD variables.
"""
mutable struct MosekMathProgModel  <: MathProgBase.AbstractMathModel
    task :: MSKtask
    
    x_block :: LinkedInts
    """
    Counts for each variable the number of integer constraints that
    are imposed on that variable. Zero means continuous.
    """
    x_block_integer :: Vector{Int}

    ###########################

    c_block :: LinkedInts
    """
    Type of the constraint block.
    """
    c_block_type    :: Array{Int,1}
    c_block_coneidx :: Array{Int,1}
    """
    Each element is either 
    - 0, meaning: no slack, when domain is defined directly as a bound,
    - a `x_block` reference, e.g. for qcones, or
    - a PSD variable index
    """
    c_block_slack   :: Array{Int,1}
    
    ###########################
    trm :: Int32
end
Mosek_VAR   = 1
Mosek_SLACK = 2


function MathProgBase.freemodel!(m::MosekMathProgModel)
    deletetask(m.task)
end

function MathProgBase.optimize!(m::MosekMathProgModel)
    m.trm = MSK_optimize(m.task)
end

#function MathProgBase.loadproblem!(...)
#end

#function MathProgBase.setparameters! ()
#end

MathProgBase.candelete(m::MosekMathProgModel,ref::VariableReference) =
    isvalid(m,ref)
MathProgBase.isvalid(m::MosekMathProgModel, ref::VariableReference) = 
    allocated(m.x_block,BlockId(ref.value))

function MathProgBase.addvariables!(m::MosekMathProgModel, N :: Int)
    ensurefree(m.x_block,N)
    r = VariableReference[ VariableReference(newblock(m.x_block,Mosek_VAR,1).id) for i in 1:N ]
    ids = BlockId[ newblock(m.x_block,Mosek_VAR,1) for i in 1:N ]
    if length(s.x_block) > numvar
        appendvars(s.task, length(s.x_block) - numvar)
    end

    # clear the variables
    idxs = Array{Int,1}(N)
    for i in 1:N
        getindexes(s.x_block,ids[i],idxs,i)
    end
    
    bld = Array{Float64}(N)
    putvarboundlist(m.task,
                    convert(Array{Int32,1}, getindexes(s.x_block, id)),
                    Int32[MSK_BK_FR for i in 1:N],
                    bnd,bld)
    
    VariableReference[VariableReference(id.id) for id in ids]
end

function MathProgBase.addvariable!(m::MosekMathProgModel)
    N = 1
    ensurefree(m.x_block,N)
    id = newblock(m.x_block,Mosek_VAR,N)
    numvar = getnumvar(s.task)
    if length(s.x_block) > numvar
        appendvars(s.task, length(s.x_block) - numvar)
    end
    
    bld = Array{Float64}(N)
    putvarboundlist(m.task,
                    convert(Array{Int32,1}, getindexes(s.x_block, id)),
                    Int32[MSK_BK_FR for i in 1:N],
                    bnd,bld)

    VariableReference(id.id)
end

function Base.delete!(m::MosekMathProgModel, refs::Vector{VariableReference})
    if ! all(r -> candelete(m,ref),refs)
        throw(CannotDelete())
    else
        ids = BlockId[ BlockId(ref.value) for ref in refs]
        sizes = Int[blocksize(m.x_block,id) for id in ids]
        N = sum(sizes)
        indexes = Array{Int}(N)
        offset = 1
        for i in 1:length(ids)
            getindexes(m.x_block,ids[i],indexes,offset)
            offset += sizes[i]
        end

        # clear all non-zeros in columns
        putacollist(m.task, 
                    indexes,
                    zeros{Int64}(N),
                    zeros{Int64}(N),
                    Int32[],
                    Float64[])
        # clear bounds
        bnd = Array{Float64,1}(N)
        putvarboundlist(m.task,
                        indexes,
                        Int32[MSK_BK_FR for i in 1:N],
                        bnd,bnd)

        for i in 1:length(ids)
            deleteblock(s.x_block,ids[i])
        end
    end
        
end

function Base.delete!(m::MosekMathProgModel, ref::VariableReference)
    if ! candelete(m,ref)
        throw(CannotDelete())
    else
        id = BlockId(ref.value)
        
        indexes = convert(Array{Int32,1},getindexes(s.x_block,id))
        N = blocksize(s.x_blocks,id)

        # clear all non-zeros in columns
        putacollist(m.task, 
                    indexes,
                    zeros{Int64}(N),
                    zeros{Int64}(N),
                    Int32[],
                    Float64[])
        # clear bounds
        bnd = Array{Float64,1}(N)
        putvarboundlist(m.task,
                        indexes,
                        Int32[MSK_BK_FR for i in 1:N],
                        bnd,bnd)

        deleteblock(s.x_block,id)
    end
end

candelete(m::MosekMathProgModel, ref::ConstraintReference) = isvalid(m,ref)
isvalid(m::MosekMathProgModel, ref::ConstraintReference) = allocated(m.c_block,BlockId(ref.value))

#Base.delete!(m::AbstractMathProgModel, ref::ConstraintReference) = throw(MethodError())

#function MathProgBase.addconstraint!(m::MosekMathProgModel, b :: Vector{Float64}, a_varidx, a_coef, Q_vari, Q_varj, Q_coef, S::MosekSet)::QuadraticConstraintReference{typeof(S)}
#end



function allocateconstraints(
    m           :: MosekMathProgModel,
    N           :: Int)
    numcon = getnumcon(m.task)
    ensurefree(m.c_block,N)
    if length(s.c_block) > numcon
        appendcons(s.task, length(s.c_block) - numcon)
    end
end




function makeconstr{T <: MathProgBase.AbstractSet}(
    m           :: MosekMathProgModel,
    a_constridx :: Vector{Int},
    a_varidx    :: Vector{VariableReference},
    a_coef      :: Vector{Float64},
    N           :: Int)

    allocateconstraints(m,N)
    
    varidxs = Array{Int}(length(a_varidx))
    for i in 1:length(a_varidx)
        getindexes(m.c_block,a_varidx[i],varidxs,i)
    end

    conid = newblock(m.c_block,1,N)
    conidxs = getindexes(m.c_block,conid)

    At = sparse(idxs, a_constridx, a_coef, getnumvar(m.task), N)    
    putarowlist(m.task,conidxs,At)

    conid,conidxs
end

function MathProgBase.addconstraint!(
    m           :: MosekMathProgModel,
    b           :: Vector{Float64},
    a_constridx :: Vector{Int},
    a_varidx    :: Vector{VariableReference},
    a_coef      :: Vector{Float64},
    S           :: MathProgBase.NonNegative)

    N = S.dim
    conid,conidxs = makeconstr(m,a_constridx, a_varidx,a_coef,N)
    putconboundlist(m.task,conidxs,fill(MSK_BK_LO,N),-b,b)

    push!(m.c_block_type, mosek_block_type_noneg)
    push!(m.c_block_coneidx, 0)
    push!(m.c_block_slack,   0)
    MathProgBase.ConstraintReference{S}(conid)
end

function MathProgBase.addconstraint!(
    m           :: MosekMathProgModel,
    b           :: Vector{Float64},
    a_constridx :: Vector{Int},
    a_varidx    :: Vector{VariableReference},
    a_coef      :: Vector{Float64},
    S           :: MathProgBase.NonPositive)

    N = S.dim
    conid,conidxs = makeconstr(m,a_constridx, a_varidx,a_coef,N)
    putconboundlist(m.task,conidxs,fill(MSK_BK_UP,N),b,-b)

    push!(m.c_block_type, mosek_block_type_nopos)
    push!(m.c_block_coneidx, 0)
    push!(m.c_block_slack,   0)
    MathProgBase.ConstraintReference{S}(conid)
end

function MathProgBase.addconstraint!(
    m           :: MosekMathProgModel,
    b           :: Vector{Float64},
    a_constridx :: Vector{Int},
    a_varidx    :: Vector{VariableReference},
    a_coef      :: Vector{Float64},
    S           :: MathProgBase.Zero)

    N = S.dim
    conid,conidxs = makeconstr(m,a_constridx, a_varidx,a_coef,N)
    putconboundlist(m.task,conidxs,fill(MSK_BK_FX,N),-b,-b)

    push!(m.c_block_type, mosek_block_type_nopos)
    push!(m.c_block_coneidx, 0)
    push!(m.c_block_slack,   0)
    MathProgBase.ConstraintReference{S}(conid)
end

function MathProgBase.addconstraint!(
    m           :: MosekMathProgModel,
    b           :: Vector{Float64},
    a_constridx :: Vector{Int},
    a_varidx    :: Vector{VariableReference},
    a_coef      :: Vector{Float64},
    S           :: MathProgBase.Interval)

    N = length(S.lower)
    conidxs = makeconstr(m,a_constridx, a_varidx,a_coef,N)

    bl = S.lower - b
    bu = S.upper - b
    putconboundlist(m.task,conidxs,fill(MSK_BK_RA,N),bl,bu)

    push!(m.c_block_type, mosek_block_type_nopos)
    push!(m.c_block_coneidx, 0)
    push!(m.c_block_slack,   0)
    MathProgBase.ConstraintReference{S}(conid)
end

function MathProgBase.addconstraint!(
    m           :: MosekMathProgModel,
    varidx      :: VariableReference,
    S           :: MathProgBase.Integers)

    N = S.dim
    MathProgBase.ConstraintReference{S}(varidx)
end

function MathProgBase.addconstraint!(m::MosekMathProgModel, varidx, S::AbstractSet)::VariablewiseConstraintReference{typeof(S)}
end


function MathProgBase.setobjective!(m::AbstractMathProgModel, N::Int, b, a_varidx, a_coef, Q_vari, Q_varj, Q_coef)
end
function MathProgBase.modifyobjective!(m::AbstractMathProgModel, i::Int, args...) end
function MathProgBase.modifyobjective!(m::AbstractMathProgModel, i::Int, b) end
function MathProgBase.modifyobjective!(m::AbstractMathProgModel, i::Int, a_varidx, a_coef) end
function MathProgBase.modifyobjective!(m::AbstractMathProgModel, i::Int, Q_vari, Q_varj, Q_coef) end
function MathProgBase.getobjective(m, i:Int) end



# COMMENTS:
#
# Missing: addcostraint!(m, varidxs,S)
# Treating integrality as a constraint gives me a headache.
