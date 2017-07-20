
MathOptInterface.candelete(m::MosekModel,ref::MathOptInterface.VariableReference) = isvalid(m,ref)
isvalid(m::MosekModel, ref::MathOptInterface.VariableReference) = allocated(m.x_block,Int(ref.value))

function MathOptInterface.addvariables!(m::MosekModel, N :: Int)
    ensurefree(m.x_block,N)    
    ids = [ newblock(m.x_block,1) for i in 1:N ]
    numvar = getnumvar(m.task)
    if length(m.x_block) > numvar
        appendvars(m.task, length(m.x_block) - numvar)
        append!(m.x_boundflags, zeros(Int,length(m.x_block) - numvar))
    end
        
    # clear the variables
    idxs = Vector{Int}(N)
    for i in 1:N
        getindexes(m.x_block,ids[i],idxs,i)
    end
    
    bnd = zeros(Float64,N)
    putvarboundlist(m.task,
                    convert(Vector{Int32}, idxs),
                    fill(MSK_BK_FR,N),
                    bnd,bnd)
    
    MathOptInterface.VariableReference[MathOptInterface.VariableReference(id) for id in ids]
end

function MathOptInterface.addvariable!(m::MosekModel)
    N = 1
    ensurefree(m.x_block,N)
    id = newblock(m.x_block,N)
    numvar = getnumvar(m.task)
    if length(m.x_block) > numvar
        appendvars(m.task, length(m.x_block) - numvar)
        append!(m.x_boundflags, zeros(Int,length(m.x_block) - numvar))
    end
    
    bnd = Vector{Float64}(N)
    putvarboundlist(m.task,
                    convert(Vector{Int32}, getindexes(m.x_block, id)),
                    fill(MSK_BK_FR,N),
                    bnd,bnd)

    MathOptInterface.VariableReference(id)
end





function Base.delete!(m::MosekModel, refs::Vector{MathOptInterface.VariableReference})
    assert(0)
    if ! all(r -> MathOptInterface.candelete(m,ref),refs)
        throw(CannotDelete())
    else
        ids = Int[ ref.value for Int(ref) in refs ]
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
                    zeros(Int64,N),
                    zeros(Int64,N),
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

function Base.delete!(m::MosekModel, ref::MathOptInterface.VariableReference)
    if ! MathOptInterface.candelete(m,ref)
        throw(CannotDelete())
    else
        id = Int(ref.value)
        
        indexes = convert(Array{Int32,1},getindexes(m.x_block,id))
        N = blocksize(m.x_block,id)

        # clear all non-zeros in columns
        putacollist(m.task, 
                    indexes,
                    zeros(Int64,N),
                    zeros(Int64,N),
                    Int32[],
                    Float64[])
        # clear bounds
        bnd = Array{Float64,1}(N)
        putvarboundlist(m.task,
                        indexes,
                        fill(MSK_BK_FR,N),
                        bnd,bnd)

        deleteblock(m.x_block,id)
    end
end
