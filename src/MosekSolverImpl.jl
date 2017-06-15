include("LinkedInts.lj")

mosek_block_type_unallocated = 0
mosek_block_type_free   = 1
mosek_block_type_nonneg = 2
mosek_block_type_nonpos = 3
mosek_block_type_qcone  = 4
mosek_block_type_rqcone = 5
mosek_block_type_psd    = 6


"""
    MosekMathProgModel <: MathProgBase.AbstractMathModel

Linear variables and constraint can be deleted. MOSEK does not support
deleting PSD variables.
"""
mutable struct MosekMathProgModel  <: MathProgBase.AbstractMathModel
    task :: MSKtask
    
    x_block :: LinkedInts

    ###########################

    c_block :: LinkedInts
    """
    Type of the constraint block.
    """
    c_block_type    :: Array{Int,1}
    c_block_coneidx :: Array{Int,1}
    """
    Each element is either 0 (meaning: no slack), or an `BlockId` from
    `x_block`.
    """
    c_block_slack   :: Array{Int,1}
    
    ###########################
    trm :: Int32
end
Mosek_VAR = 1
Mosek_SLACK = 2




function MathProgBase.freemodel!(m::MosekMathProgModel)
    deletetask(m.task)
end

function MathProgBase.optimize!(m::MosekMathProgModel)
    m.trm = MSK_optimize(m.task)
end

function MathProgBase.loadproblem!(...)
end

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

function MathProgBase.addvariable!(m::MosekMathProgModel, N :: Int)
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
isvalid(m::AbstractMathProgModel, ref::ConstraintReference) = 
    allocated(m.c_block,BlockId(ref.value))

Base.delete!(m::AbstractMathProgModel, ref::ConstraintReference) = throw(MethodError())
function addconstraint!(m::MosekMathProgModel, b, a_constridx, a_varidx, a_coef, Q_constridx, Q_vari, Q_varj, Q_coef, S::AbstractSet)::ConstraintReference{typeof(S)}
end
function addconstraint!(m::MosekMathProgModel, b, a_varidx, a_coef, Q_vari, Q_varj, Q_coef, S::AbstractSet)::ConstraintReference{typeof(S)}
end
function addconstraint!(m::MosekMathProgModel, varidx, S::AbstractSet)::ConstraintReference{typeof(S)}
end
