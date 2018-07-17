module MosekBlocks
import Mosek

export Model,variable!

mutable struct ManagedIndexes
    next       :: Vector{Int}
    prev       :: Vector{Int}
    first_used :: Int
    last_used  :: Int
    first_free :: Int
    last_free  :: Int

    nfree      :: Int
end

ManagedIndexes() = ManagedIndexes( Int[], Int[], 0, 0, 0, 0, 0 )

Base.length(mi :: ManagedIndexes) = length(mi.next)

function ensure!(self :: ManagedIndexes, num :: Int)
    if self.nfree < num
        addnum = num-self.nfree
        oldnum = length(self.next)
        newnum = oldnum + addnum
        resize!(self.next,newnum)
        resize!(self.prev,newnum)
        
        for i in oldnum+1:newnum
            self.next[i] = i+1
            self.prev[i] = i-1
        end
        first = oldnum+1
        last  = newnum
        
        self.next[last] = 0
        self.prev[first] = self.last_free
        if self.last_free == 0
            self.first_free = first
            self.last_free = last
        else
            self.next[self.last_free] = first
            self.last_free = last
        end
        self.nfree += addnum
    end
end

function alloc!(self :: ManagedIndexes, num :: Int)
    ensure!(self,num)
    first = self.first_free
    last = first
    for i in 1:num-1
        last = self.next[last]
    end
    
    if self.next[last] > 0
        self.prev[self.next[last]] = 0
    end

    self.first_free = self.next[last]
    if self.first_free == 0
        self.last_free = 0
    end

    self.next[last] = 0
    self.prev[first] = self.last_used
    if self.last_used == 0
        self.first_used = first
    else
        self.next[self.last_used] = first
    end
    self.last_used = last
    self.nfree -= num
    
    first
end

function free(mi :: ManagedIndexes, index :: Int, num :: Int)
    first = index
    last = index
    for i in 1:num-1
        last = mi.next[last]
    end

    if mi.prev[first] == 0
        if mi.next[last] == 0
            mi.first_used = 0
            mi.last_used = 0
        else
            mi.first_used = mi.next[last]
            mi.prev[mi.first_used] = 0
        end
    else
        if mi.next[last] == 0
            mi.last_used = mi.prev[first]
            mi.next[mi.last_used] = 0
        else
            mi.prev[mi.next[last]] = mi.prev[first]
            mi.next[mi.prev[first]] = mi.next[last]
        end
    end

    if mi.first_free == 0
        mi.first_free = first
        mi.last_free = last
    else
        mi.next[last] = mi.first_free
        mi.prev[mi.first_free] = last
        mi.first_free = first
    end

    mi.nfree += num
end

function get_indexes(self::ManagedIndexes, head::Int, num::Int, res::Vector{Int})
    ii = head
    res[1] = ii
    for i in 2:num
        ii = self.next[ii]
        res[i] = ii
    end
end
#----------------------------
mutable struct ManagedBlocks
    ii :: ManagedIndexes
    blocklist :: ManagedIndexes

    blocksize :: Vector{Int}
    blockhead :: Vector{Int}
end

ManagedBlocks() = ManagedBlocks(ManagedIndexes(),ManagedIndexes(),Int[],Int[])

function alloc!(self :: ManagedBlocks, num :: Int)
    p       = alloc!(self.ii,num)
    blockid = alloc!(self.blocklist,1)

    numblocks = length(self.blocklist)
    if length(self.blocksize) < numblocks
        resize!(self.blocksize,numblocks)
        resize!(self.blockhead,numblocks)
    end
    self.blocksize[blockid] = num
    self.blockhead[blockid] = p

    blockid
end

function free!(self :: ManagedBlocks, blockid :: Int)
    free!(self.ii, self.blockhead[blockid], self.blocksize[blockid])
end

Base.length(self::ManagedBlocks) = length(self.ii)

num_blocks(self::ManagedBlocks) = length(self.blocklist)

function get_indexes(self::ManagedBlocks, blockid :: Int)
    num = self.blocksize[blockid]
    head = self.blockhead[blockid]
    res = zeros(Int,num)
    get_indexes(self.ii, head, num, res)
    res
end


type Model
    task :: Mosek.MSKtask

    varblocks :: ManagedBlocks
    conblocks :: ManagedBlocks

    conbfix   :: Vector{Float64}
    conslack  :: Vector{Int}
end

Model() = Model(Mosek.maketask(),
                ManagedBlocks(),
                ManagedBlocks(),
                Float64[],
                Float64[] )

type Variable{N}
    shape   :: NTuple{N,Int}
    indexes :: Vector{Int}
end

type Constraint{N}
    shape :: NTuple{N,Int}
    indexes :: Vector{Int}
end

type Expr{N}
    shape :: NTuple{N,Int}
    ptr   :: Vector{Int}
    nidx  :: Vector{Int}
    cof   :: Vector{Float64}
    bfix  :: Vector{Float64}
end

function alloc_var!(m::Model, num :: Int)
    numvar = Mosek.getnumvar(m.task)
    blockid = alloc!(m.varblocks,num)
    if length(m.varblocks) < numvar
        nadd = numvar-length(m.varblocks)
        Mosek.appendvars(m.task,nadd)
        Mosek.putvarboundslice(m.task,numvar+1,numvar+nadd+1,
                               fill(Mosek.MSK_BK_FR,nadd),
                               zeros(Float64,nadd),
                               zeros(Float64,nadd))
    end
    blockid
end

function alloc_con!(m::Model, num :: Int)
    numcon = Mosek.getnumcon(m.task)
    blockid = alloc!(m.varblocks,num)
    if length(m.conblocks) < numvar
        nadd = numcon-length(m.conblocks)
        newnumcon = numcon+nadd
        Mosek.appendcons(m.task,nadd)
        Mosek.putconboundslice(m.task,numcon+1,numcon+nadd+1,
                               fill(Mosek.MSK_BK_FR,nadd),
                               zeros(Float64,nadd),
                               zeros(Float64,nadd))
        resize!(m.conbfix,newnumcon)
        resize!(m.conslack,newnumcon)

        m.conbfix[numcon+1:]  = 0.0
        m.conslack[numcon+1:] = 0
    end
    blockid
end

function variable!{N}(m::Model, shape :: NTuple{N,Int})
    num = prod(shape)
    blockid = alloc_var!(m,num)
    idxs = get_indexes(m.varblocks, blockid)
    Variable( shape, idxs )
end

variable!(m::Model, shape :: Int) = variable!(m,(shape,))

function constraint!(m:Model, shape :: NTuple{N,Int}, expr :: Expr, dom ::D)
    where N, D <: AbstractDomain

    
end





###############
# TEST

m = MosekBlocks.Model()

v1 = MosekBlocks.variable!(m, (2,2))
v2 = MosekBlocks.variable!(m,2)
