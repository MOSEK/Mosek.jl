
MathOptInterface.candelete(m::MosekModel,ref::MathOptInterface.VariableReference) = isvalid(m,ref) && m.x_numxc[ref2id(ref)] == 0
isvalid(m::MosekModel, ref::MathOptInterface.VariableReference) = allocated(m.x_block,ref2id(ref))

function MathOptInterface.addvariables!(m::MosekModel, N :: Int)
    ids = [ allocatevariable(m,1) for i in 1:N ]

    m.publicnumvar += N
    
    idxs = Vector{Int}(N)
    for i in 1:N
        getindexes(m.x_block,ids[i],idxs,i)
    end
    
    bnd = zeros(Float64,N)
    putvarboundlist(m.task,
                    convert(Vector{Int32}, idxs),
                    fill(MSK_BK_FR,N),
                    bnd,bnd)

    # if DEBUG
    if DEBUG
        for i in idxs
            putvarname(m.task,Int32(i),"x$i")
        end
    end

    [ id2vref(id) for id in ids]
end

function MathOptInterface.addvariable!(m::MosekModel)
    N = 1
    id = allocatevariable(m,1)
    m.publicnumvar += N
    bnd = Vector{Float64}(N)
    subj = convert(Vector{Int32}, getindexes(m.x_block, id))
    putvarboundlist(m.task,
                    subj,
                    fill(MSK_BK_FR,N),
                    bnd,bnd)
    # if DEBUG
    if DEBUG
        for i in subj
            putvarname(m.task,Int32(i),"x$i")
        end
    end

    id2vref(id)
end


function Base.delete!(m::MosekModel, refs::Vector{MathOptInterface.VariableReference})
    ids = Int[ ref2id(ref) for ref in refs ]

    if ! all(id -> m.x_numxc[id] == 0, idxs)
        error("Cannot delete a variable while a bound constraint is defined on it")
    elseif ! all(r -> MathOptInterface.candelete(m,ref),refs)
        throw(CannotDelete())
    else
        sizes = Int[blocksize(m.x_block,id) for id in ids]
        N = sum(sizes)
        m.publicnumvar -= length(refs)
        indexes = Array{Int}(N)
        offset = 1
        for i in 1:length(ids)
            getindexes(m.x_block,ids[i],indexes,offset)
            offset += sizes[i]
        end

        # clear all non-zeros in columns
        putacollist(m.task, 
                    indexes,
                    zeros(Int64,N),
                    zeros(Int64,N),
                    Int32[],
                    Float64[])        
        putclist(m.task,indexes,zeros(Int64,N))
        # clear bounds
        bnd = zeros(Float64,N)
        putvarboundlist(m.task,
                        indexes,
                        fill(MSK_BK_FX,N),
                        bnd,bnd)

        if DEBUG
            for i in idxs
                putvarname(m.task,Int32(idxs),"deleted$i")
            end
        end
        
        for i in 1:length(ids)
            deleteblock(s.x_block,ids[i])
        end
    end
end

function Base.delete!(m::MosekModel, ref::MathOptInterface.VariableReference)
    if m.x_numxc[ref2id(ref)] != 0
        error("Cannot delete a variable while a bound constraint is defined on it")
    elseif ! MathOptInterface.candelete(m,ref)
        throw(CannotDelete())
    else
        id = ref2id(ref)

        m.publicnumvar -= 1
        
        indexes = convert(Array{Int32,1},getindexes(m.x_block,id))
        N = blocksize(m.x_block,id)

        # clear all non-zeros in columns
        for i in indexes
            putcj(m.task,i,0.0)
        end
        putacollist(m.task, 
                    indexes,
                    zeros(Int64,N),
                    zeros(Int64,N),
                    Int32[],
                    Float64[])
        # clear bounds
        bnd = zeros(Float64,N)
        putvarboundlist(m.task,
                        indexes,
                        fill(MSK_BK_FX,N),
                        bnd,bnd)
        if DEBUG
            for i in indexes
                putvarname(m.task,Int32(i),"deleted$i")
            end
        end

        deleteblock(m.x_block,id)
    end
end



###############################################################################
## ATTRIBUTES
