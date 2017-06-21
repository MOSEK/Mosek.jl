
MathProgBase.candelete(m::MosekModel,ref::MathProgBase.VariableReference) =
    isvalid(m,ref)
MathProgBase.isvalid(m::MosekModel, ref::MathProgBase.VariableReference) = 
    allocated(m.x_block,BlockId(ref.value))

function MathProgBase.addvariables!(m::MosekModel, N :: Int)
    ensurefree(m.x_block,N)
    r = MathProgBase.VariableReference[ MathProgBase.VariableReference(newblock(m.x_block,Mosek_VAR,1).id) for i in 1:N ]
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
    
    MathProgBase.VariableReference[MathProgBase.VariableReference(id.id) for id in ids]
end

function MathProgBase.addvariable!(m::MosekModel)
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

    MathProgBase.VariableReference(id.id)
end

function MathProgBase.addvariable!(m::MosekModel,
                                   cref::Vector{Union{
                                   MathProgBase.AffineConstraintReference{MathProgBase.NonPositive},
                                   MathProgBase.AffineConstraintRef{MathProgBase.NonNegative},
                                   MathProgBase.AffineConstraintRef{MathProgBase.Zero},
                                   MathProgBase.AffineConstraintRef{MathProgBase.Interval}
                                   }},
                                   coefs :: Vector{Float64})
    N = 1
    ensurefree(m.x_block,N)
    id = newblock(m.x_block,Mosek_VAR,N)
    numvar = getnumvar(s.task)
    if length(s.x_block) > numvar
        appendvars(s.task, length(s.x_block) - numvar)
    end

    sumlen = foldl((v,i) -> v+blocksize(m.c_block,BlockId(i)),0,cref)
    if sumlen != length(coefs)
        error("Lengths of cref and coefs do not match")
    end

    idxs = Array{Int}(length(coef))
    offset = 1
    for i in 1:length(cref)
        getindexes(s.c_block, BlockId(cref[i].value), offset, idxs)
        offset += blocksize(s.c_block, BlockId(cref[i].value))
    end
    
    bld = Array{Float64}(N)
    putvarboundlist(m.task,
                    convert(Array{Int32,1}, getindexes(s.x_block, id)),
                    Int32[MSK_BK_FR for i in 1:N],
                    bnd,bld)

    MathProgBase.VariableReference(id.id)
end


function Base.delete!(m::MosekModel, refs::Vector{MathProgBase.VariableReference})
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

function Base.delete!(m::MosekModel, ref::MathProgBase.VariableReference)
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
