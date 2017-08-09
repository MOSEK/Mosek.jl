
mutable struct LinkedInts
    next :: Vector{Int}
    prev :: Vector{Int}

    free_ptr :: Int
    free_cap :: Int
    root     :: Int

    block :: Vector{Int}
    size  :: Vector{Int}
end

LinkedInts(capacity=128) =
    LinkedInts(Int[],Int[],
               0,0,0,
               Int[],
               Int[])

allocated(s::LinkedInts, id :: Int) = id > 0 && id <= length(s.block) && s.block[id] > 0
blocksize(s::LinkedInts, id :: Int) = s.size[id]
Base.length(s::LinkedInts) = length(s.next)
numblocks(s::LinkedInts) = length(s.block)

function Base.string(l::LinkedInts)
    r = String[]
    push!(r,"LinkedInts(\n")
    push!(r,@sprintf("  Number of blocks: %d\n", length(l.block)))
    push!(r,"  Blocks:\n")

    for i in 1:length(l.block)
        if l.block[i] > 0
            idxs = getindexes(l,i)            
            push!(r,@sprintf("    #%d: %s\n",i,string(idxs)))
        end
    end
    p = l.free_ptr
    freelst = Int[]
    while p > 0
        push!(freelst,p)
        p = l.prev[p]
    end
    push!(r,@sprintf("  Free: %s\n",string(freelst)))
    
    push!(r,")")
    join(r)
end

"""
    ensurefree(s::LinkedInts, N :: Int)

Ensure that there are at least `N` elements free, and allocate as necessary.
"""
function ensurefree(s::LinkedInts, N :: Int)
    if s.free_cap < N
        num = N - s.free_cap

        cap = length(s.next)
        first = cap+1
        last  = cap+num
        
        append!(s.next,Int[i+1 for i in first:last])
        append!(s.prev,Int[i-1 for i in first:last])

        s.next[last] = 0
        s.prev[first] = s.free_ptr
        if s.prev[first] > 0
            s.next[s.prev[first]] = first
        end
        s.free_ptr = last
        s.free_cap += num

        num
    else
        0
    end
end

"""
    newblock(s::LinkedInts, N :: Int)

Add a new block in list `idx`
"""
function newblock(s::LinkedInts, N :: Int) :: Int
    assert(N>0)
    ensurefree(s,N)
    # remove from free list
    ptre = s.free_ptr
    ptrb = ptre
    for i = 1:N-1
        ptrb = s.prev[ptrb]
    end

    prev = s.prev[ptrb]

    if prev > 0
        s.next[prev] = 0
    end

    s.free_ptr = s.prev[ptrb]
    s.free_cap -= N
    
    # insert into list `idx`
    s.prev[ptrb] = s.root
    if s.root > 0
        s.next[s.root] = ptrb
    end
    s.root = ptre
    push!(s.block,ptrb)
    push!(s.size,N)

    id = length(s.block)
    
    id
end

"""
Move a block to the free list.
"""
function deleteblock(s::LinkedInts, id :: Int)
    if s.size[id] > 0
        ptrb = s.block[id]
        N = s.size[id]
        ptre = ptrb
        for i in 2:N
            ptre = s.next[ptre]
        end
        prev = s.prev[ptrb]
        next = s.next[ptre]
        
        # remove from list and clear the block id
        if s.root == ptre s.root = prev end
        if prev > 0 s.next[prev] = next end
        if next > 0 s.prev[next] = prev end
            
        s.size[id]  = 0
        s.block[id] = 0
        
        # add to free list        
        if s.free_ptr > 0
            s.next[s.free_ptr] = ptrb
        end
        s.prev[ptrb] = s.free_ptr
        s.free_ptr = ptre
        s.next[ptre] = 0

        s.free_cap += N
    end
end

"""
    getindexes(s::LinkedInts, id :: Int)


"""
function getindexes(s::LinkedInts, id :: Int)
    N = s.size[id]
    r = Array{Int}(N)
    p = s.block[id]
    for i in 1:N
        r[i] = p
        p = s.next[p]
    end
    r
end


function getindexes(s::LinkedInts, id :: Int, target :: Array{Int,1}, offset :: Int)
    N = s.size[id]
    p = s.block[id]
    for i in 1:N
        target[i+offset-1] = p
        p = s.next[p]
    end
    N
end

function getoneindex(s::LinkedInts, id :: Int)
    N = s.size[id]
    if N < 1
        error("No values at id")
    end

    s.block[i]
end


"""
Get a list if the currently free elements.
"""
function getfreeindexes(s::LinkedInts)    
    N = s.free_cap
    r = Array{Int}(N)
    ptr = s.free_ptr
    for i in 1:N
        r[N-i+1] = ptr
        ptr  = s.prev[ptr]
    end
    r
end



"""
Get a list if the currently free elements.
""" 
function getusedindexes(s::LinkedInts)    
    N = length(s.next) - s.free_cap
    r = Array{Int}(N)
    ptr = s.root
    for i in 1:N
        r[N-i+1] = ptr
        ptr  = s.prev[ptr]
    end
    r
end



"""
Check consistency of the internal structures.
"""
function checkconsistency(s::LinkedInts) :: Bool
    if length(s.prev) != length(s.next)
        return false
    end
    
    N = length(s.prev)

    all(i -> s.prev[i] == 0 || s.next[s.prev[i]] == i, 1:N) &&
    all(i -> s.next[i] == 0 || s.prev[s.next[i]] == i, 1:N)
end
