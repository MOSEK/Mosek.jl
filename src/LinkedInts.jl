
struct BlockId id :: Int end
mutable struct LinkedInts
    next :: Array{Int,1}
    prev :: Array{Int,1}

    free_ptr :: Int
    free_cap :: Int
    roots :: Array{Int,1}
    cap   :: Array{Int,1}

    block :: Array{Int,1}
    size  :: Array{Int,1}
    owner :: Array{Int,1}

    function LinkedInts(nlists :: Int) 
        new(Int[],Int[],0,0,zeros{Int}(nlists),zeros{Int}(nlists),Int[],Int[],Int[])
    end
end




allocated(s::LinkedInts, id :: BlockId) = id.id > 0 && id.id <= length(s.block) && s.block[id.id] > 0
blocksize(s::LinkedInts, id :: BlockId) = s.size[id.id]
Base.length(s::LinkedInts) = length(s.next)
"""
    allocate(s::LinkedInts, N :: Int, idx :: Int)

Allocate `N` item and place it in list `idx`.

Return the index if the first allocated element.
"""
function ensurefree(s::LinkedInts, N :: Int)
    if s.free_cap < N
        num = N - s.free_cap

        cap = length(s.next)
        append!(s.next,Int[i+1 for i in cap+1:cap+num])
        append!(s.prev,Int[i-1 for i in cap+1:cap+num])

        s.next[cap+num] = 0
        s.prev[cap+1] = s.roots[idx]
        if s.prev[cap+1] > 0
            s.next[s.prev[cap+1]] = cap+1
        end
        s.free_ptr = cap+num
        s.free_cap += num
    end
end

"""
Add a new block in list `idx`
"""
function newblock(s::LinkedInts, idx :: Int, N :: Int)
    ensurefree(N)
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
    
    # insert into list `idx`
    s.prev[ptrb] = s.roots[idx]
    if s.roots[idx] > 0
        s.next[s.roots[idx]] = ptrb
    end
    s.roots[idx] = ptre
    push!(s.block,ptrb)
    push!(s.size,N)
    push!(s.owner,idx)

    id = length(s.block)
    
    BlockId(id)
end

"""
Move a block to the free list.
"""
function deleteblock(s::LinkedInts, id_ :: BlockId)
    id = id_.id
    if s.size[id] > 0
        idx = s.owner[id]
        ptrb = s.block[id]
        N = s.size[id]
        ptre = ptrb
        for i in 2:N
            ptre = s.next[ptre]
        end
        prev = s.prev[ptrb]
        next = s.next[ptre]
        
        # remove from list and clear the block id
        if prev s.next[prev] = next end
        if next s.prev[next] = prev end
        s.cap[idx] -= N
        s.size[id] = 0
        s.block = 0
        
        # add to free list
        if s.free_ptr > 0
            s.prev[ptrb] = s.free_ptr
        end
        s.free_ptr = ptre
        s.next[ptre] = 0

        s.free_cap += N
    end
end

"""
    get(s::LinkedInts, ptrb :: Int, N :: Int)

Get `N` integers starting at `ptrb`. It is not checked if the list
contains the requested number of items.

Return an array of integers.
"""
function getindexes(s::LinkedInts, id :: BlockId)
    N = s.size[id.id]
    r = Array{Int}(N)
    p = s.block[id.id]
    for i in 1:N
        r[i] = p
        p = s.next[p]
    end
    r
end

function getindexes(s::LinkedInts, id :: BlockId, target :: Array{Int,1}, offset :: Int)
    N = s.size[id.id]
    p = s.block[id.id]
    for i in 1:N
        target[i+offset] = p
        p = s.next[p]
    end
end
