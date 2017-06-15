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
    
    x_block_ptr :: Array{Int,1}
    x_block_num :: Array{Int,1}

    x_block_numfree :: Int
    """
    Linked list of linear native variable vector indexes. Effectively,
    a permutation of the variable vector.
    """
    x_block_next :: Array{Int,1}
    x_block_prev :: Array{Int,1}
    """
    Index into `x_block_prev` of the first used variable.
    """
    x_block_used_ptr :: Int
    """
    Index into `x_block_prev` of the first free variable.
    """
    x_block_free_ptr :: Int

    ###########################

    c_block_ptrb :: Array{Int,1}
    """
    Type of the constraint block.
    """
    c_block_type :: Array{Int,1}
    """
    Each entry is either 0 or the index of the cone of the
    corresponding block.
    """
    c_block_coneidx :: Array{Int,1}
    c_block_prev :: Array{Int,1}
    c_block_used_ptr :: Array{Int,1}
    c_block_free_ptr :: Array{Int,1}
    
    ###########################
    trm :: Int32
end

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

candelete(m::MosekMathProgModel,ref::VariableReference) =
    isvalid(m,ref) &&
    m.x_block_type != mosek_block_type_psd

isvalid(m::MosekMathProgModel, ref::VariableReference) = 
    ref.value > 0 &&
    ref.value < length(m.x_block_ptrb) && 
    m.x_block_ptr > 0



#internal
function allocate_variables(m::MosekMathProgModel, N :: Int)
    numvar = getnumvar(m.task)
    appendvars(m.task,N)
    
    curnum = length(m.x_block_next)
    m.x_block_next[curnum] = curnum+1

    next = Array{Int}(N)
    prev = Array{Int}(N)
    
    next[1:N-1] = curnum:curnum+N-1
    next[N]     = 0
    prev[2:N]   = curnum+1:curnum+N
    prev[1]     = m.x_block_free_ptr

    append!(m.x_block_next, next)
    append!(m.x_block_ptrb, prev)
    if prev[1] > 0
        m.x_block_next[prev[1]] = curnum+1
    end

    m.x_block_free_ptr = curnum+N
end

function addvariables!(m::MosekMathProgModel, N :: Int)
    if m.x_block_numfree < N 
        allocate_variables(m,N-m.x_block_numfree)
    end
    
    ptre = m.x_block_free_ptr
    ptrb = ptre
    for _ in 1:N-1
        ptrb = m.x_block_prev[ptre]
    end

    m.x_block_free_ptr = m.x_block_prev[ptrb]
    if m.x_block_prev[ptrb] > 0
        m.x_block_next[m.x_block_prev[ptrb]] = 0
    end

    m.x_block_prev[ptrb] = m.x_block_used_ptr
    m.x_block_used_ptr = ptre

    m.x_block_numfree -= N

    ptr = ptrb
    num = length(m.x_block_ptr)
    for i in 1:N
        push!(m.x_block_ptr,ptr)
        push!(m.x_block_num,1)
        ptr = m.x_block_next[ptr]
    end

    VariableRef[ VariableReference(i) for i in num+1:num+N]
end

function addvariable!(m::MosekMathProgModel, N :: Int)
    if m.x_block_numfree < N 
        allocate_variables(m,N-m.x_block_numfree)
    end
    
    ptre = m.x_block_free_ptr
    ptrb = ptre
    for _ in 1:N-1
        ptrb = m.x_block_prev[ptre]
    end

    m.x_block_free_ptr = m.x_block_prev[ptrb]
    if m.x_block_prev[ptrb] > 0
        m.x_block_next[m.x_block_prev[ptrb]] = 0
    end

    m.x_block_prev[ptrb] = m.x_block_used_ptr
    m.x_block_used_ptr = ptre

    m.x_block_numfree -= N

    push!(m.x_block_ptr,ptrb)
    push!(m.x_block_num,N)

    VariableReference(length(m.x_block_ptr))
end

function Base.delete!(m::MosekMathProgModel, ref::VariableReference)
    if ! candelete(m,ref)
        throw(CannotDelete())
    else
        idx = ref.value
        ptrb = m.x_block_ptr[idx]
        num  = m.x_block_num[idx]

        
    end
end






candelete(m::MosekMathProgModel, ref::ConstraintReference) = isvalid(m,ref)
isvalid(m::AbstractMathProgModel, ref::ConstraintReference) = 
    ref.value > 0 &&
    ref.value < length(m.c_block_ptrb) &&
    m.c_block_valid[ref.value]

Base.delete!(m::AbstractMathProgModel, ref::ConstraintReference) = throw(MethodError())
function addconstraint!(m::MosekMathProgModel, b, a_constridx, a_varidx, a_coef, Q_constridx, Q_vari, Q_varj, Q_coef, S::AbstractSet)::ConstraintReference{typeof(S)}
end
function addconstraint!(m::MosekMathProgModel, b, a_varidx, a_coef, Q_vari, Q_varj, Q_coef, S::AbstractSet)::ConstraintReference{typeof(S)}
end
function addconstraint!(m::MosekMathProgModel, varidx, S::AbstractSet)::ConstraintReference{typeof(S)}
end
