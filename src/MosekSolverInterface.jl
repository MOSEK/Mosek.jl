module MosekMathProgSolverInterface
using ..Mosek
using Compat

export MosekSolver


# Known issues:
#  - SOCP and QP cannot be mixed, but this is not checked (an error from mosek will be produced, though)
#  - Adding a conic quadratic constraint will add an empty constraint to ensure that the number of values
#    in constraint solution is as expected. The actual constraint solution value is bogus.
#  - Adding rotated conic quadratic constraints will result in a constraint being added, but the constraint soloution
#    for this is pointless. Also, a variable is added, but this is filtered out in the results.
#  - Loading an SOCP problem file will cause some funky problems as information on extra variables etc. is lost.
#  - Dual information is currently useless.
#
#  - The concept of dual values is a bit shaky. Specifically; for a variable x there is a dual for the upper bound,
#    one for the lower bound and one for the conic "bound". The dual value reported will be (slx-sux+snx).

importall MathProgBase.SolverInterface

const MosekMathProgModel_LINR = 0
const MosekMathProgModel_QOQP = 1
const MosekMathProgModel_SOCP = 2
const MosekMathProgModel_SDP  = 3

const MosekMathProgModel_VTCNT = @compat(Int8(0))
const MosekMathProgModel_VTINT = @compat(Int8(1))
const MosekMathProgModel_VTBIN = @compat(Int8(2))

type MosekMathProgModel <: AbstractMathProgModel
  task :: Mosek.MSKtask
  probtype :: Int

  # numvar 
  #   Number of elements used in varmap,barvarij
  numvar   :: Int
  # varmap 
  #   Maps model variables to MOSEK variables.

  #   - varmap[i] > 0: it refers to MOSEK linear variable index (varmap[i]),
  #     barvarij[i] is 0.
  #   - varmap[i] < 0: it refers to MOSEK SDP variable index (-varmap[i]), 
  #     and barvarij[i] is the linear linear index into the
  #     column-oriented lower triangular part of the SDP variable.
  varmap   :: Array{Int32,1} 
  barvarij :: Array{Int64,1}
  
  # numbarvar 
  #   Number of used elements in barvarmap.
  numbarvar :: Int32
  # barvarmap 
  #   Maps Model PSD variable indexes into MOSEK barvar indexes. When
  #   using the Semidefinite interface, these are the semidefinite
  #   variables thar can be accessed.  Semidefinite variables added
  #   thorugh loadconicproblem!() are not mapped here.
  barvarmap :: Array{Int32,1}

  # binvarflag 
  #   Defines per variable if it is binary
  binvarflag :: Array{Bool,1}

  # numcon 
  #   Number of elements used in conmap
  numcon   :: Int32
  # conmap 
  #   Maps Model constraints to native MOSEK constraints. Auxiliary
  #   constraints are not mapped through. Positive values map to the
  #   corresponding constraint in MOSEK, while 0 indicated a
  #   'NULL-constraint', i.e. A placeholder constraint that exists
  #   from the user's point of view, but doesn't map to a MOSEK
  #   constraint (used specifically as a place-holder constraint for
  #   conic-quadratic constraints).
  conmap   :: Array{Int32,1}
  # conslack 
  #   Maps Model constraint to corresponding slack variable
  # 
  #   - conslack[i] == 0: No slack (for linear constraints)
  #   - conslack[i] > 0: Constraint is conic quadratic and corresponds to MOSEK variable conslack[i]
  #   - conslack[i] < 0: Constraint is PSD conic and corresponds to
  #
  #   MOSEK variable -coneslack[i]. barconij[i] is the linear index
  #   into the column-oriented lower triangular part of the SDP
  #   variable.
  conslack :: Array{Int32,1} 
  barconij :: Array{Int64,1}

  # quadratic constraints
  numqcon :: Int32
  qconmap :: Array{Int32,1}

  # options Options from MosekSolver
  options
end



immutable MosekSolver <: AbstractMathProgSolver
  options
end
MosekSolver(;kwargs...) = MosekSolver(kwargs)

type MosekMathProgModelError <: Exception
  msg :: String
end

function loadoptions!(m::MosekMathProgModel)
  # write to console by default
  printstream(msg::String) = print(msg)

  putstreamfunc(m.task,MSK_STREAM_LOG,printstream)
  for (o,val) in m.options
      if isa(val, Integer)
          parname = "MSK_IPAR_$o"
          putnaintparam(m.task, parname, val)
      elseif isa(val, FloatingPoint)
          parname = "MSK_DPAR_$o"
          putnadouparam(m.task, parname, val)
      elseif isa(val, String)
          parname = "MSK_SPAR_$o"
          putnastrparam(m.task, parname, val)
      else
          error("Value $val for parameter $o has unrecognized type")
      end
  end
end


function model(s::MosekSolver)
  task = maketask(Mosek.msk_global_env)

  m = MosekMathProgModel(task,
                         MosekMathProgModel_LINR, # pt
                         0,       # numvar
                         Int32[], # varmap
                         Int64[], # barvarij
                         0,       # numbarvar
                         Int32[], # barvarmap
                         Bool[],  # binvarflag
                         0,       # numcon
                         Int32[], # conmap
                         Int32[], # conslack
                         Int64[], # barconij
                         0,       # numqcon
                         Int64[], # qconmap
                         s.options)
  loadoptions!(m)
  return m
end

function getBoundsKey(lb, ub)
    ret = convert(Int32,0)
    if lb == -Inf && ub == Inf
        ret = MSK_BK_FR
    elseif lb == ub
        ret = MSK_BK_FX
    elseif ub  == Inf
        ret = MSK_BK_LO
    elseif lb == -Inf
        ret = MSK_BK_UP
    else
        ret = MSK_BK_RA
    end
    return convert(Int32, ret) #just to be safe
end

# NOTE: This method will load data into an existing task, but
# it will not necessarily reset all exitsing data, depending on the
# file read (e.g. reading an MPS will not reset parameters)
# - auxilary variables and cones will not be correctly mapped.
# - SDP variables are not mapped to the linear variable
# 
function loadproblem!(m::        MosekMathProgModel, 
                      filename:: String)    
    readdata(m.task, filename)

    m.numvar     = getnumvar(m.task)
    m.varmap     = Int32[1:m.numvar]
    m.barvarij   = zeros(Int64,m.numvar)
    m.binvarflag = fill(false,length(varmap)) 
    m.numbarvar  = getnumbarvar(m.task)
    m.barvarmap  = Int32[1:m.numbarvar]
    
    let numqonz = getnumqobjnz(m.task),
        numqconknz = Int32[ getnumqconknz(m.task,i) for i in 1:getnumcon(m.task) ]
        
        if numqonz + sum(numqconknz) == 0
            m.numcon     = getnumcon(m.task)
            m.conmap     = Int32[1:m.numcon]
            m.conslack   = zeros(Int32,m.numcon)
            m.barconij   = zeros(Int64,m.numcon)
            m.numqcon    = 0
            m.qconmap    = Int32[]

            m.probtype =
                if m.numbarvar > 0
                    MosekMathProgModel_SDP
                elseif getnumcone(m.task) > 0
                    MosekMathProgModel_SOCP
                else
                    MosekMathProgModel_LINR
                end
        else # quadratic problem
            m.numcon  = count(nqnz -> nqnz > 0, numqconknz)
            m.conmap   = find(nqnz -> nqnz == 0, numqconknz)
            m.conslack = zeros(Int32,m.numcon)
            m.barconij = zeros(Int64,m.numcon)
            m.numqcon = getnumcon(m.task) - m.numcon
            m.qconmap  = find(nqnz -> nqnz >  0, numqconknz)

            m.probtype   = MosekMathProgModel_QOQP
        end
    end
end

function writeproblem(m:: MosekMathProgModel, filename:: String)
  writedata(m.task, filename)
end


function loadproblem!(m::     MosekMathProgModel,
                      A::     SparseMatrixCSC{Float64,Int},
                      collb:: Array{Float64,1},
                      colub:: Array{Float64,1},
                      obj::   Array{Float64,1},
                      rowlb:: Array{Float64,1},
                      rowub:: Array{Float64,1},
                      sense:: Symbol)
  putmaxnumvar(m.task,0)
  putmaxnumcon(m.task,0)
  putmaxnumcone(m.task,0)
  putmaxnumbarvar(m.task,0)
  putmaxnumqnz(m.task,0)
  
  nrows,ncols = size(A)  
  if ncols != length(collb) ||
     ncols != length(colub) ||
     ncols != size(obj,1)   ||
     nrows != length(rowlb) ||
     nrows != length(rowub) ||
     ncols != length(obj)

    throw(MosekMathProgModelError("Inconsistent data dimensions"))
  end

  appendvars(m.task,ncols)
  appendcons(m.task,nrows)

  m.numvar     = ncols
  m.varmap     = Int32[1:m.numvar;]
  m.barvarij   = zeros(Int64,m.numvar)
  m.binvarflag = fill(false,m.numvar) 
  m.numcon     = nrows
  m.conmap     = Int32[1:m.numcon;]
  m.conslack   = zeros(Int32,m.numcon)
  m.barconij   = zeros(Int64,m.numcon)
  
  m.numqcon    = 0
  m.qconmap    = Int32[]

  # input coefficients
  putclist(m.task, Int32[1:ncols;], obj)
  putacolslice(m.task, 1, ncols+1, A.colptr[1:ncols], A.colptr[2:ncols+1], A.rowval, A.nzval)
  setsense!(m, sense)

  # input bounds
  putvarboundslice(m.task, 1, ncols+1, Int32[ MSK_BK_RA for i=1:ncols ], collb, colub)
  putconboundslice(m.task, 1, nrows+1, Int32[ MSK_BK_RA for i=1:nrows ], rowlb, rowub)
  nothing
end

function loadproblem!(m::     MosekMathProgModel,
                      A,
                      collb,
                      colub,
                      obj,
                      rowlb,
                      rowub,
                      sense :: Symbol)

  loadproblem!(m,
               convert(SparseMatrixCSC{Float64,Int},sparse(A)),
               convert(Array{Float64,1},collb),
               convert(Array{Float64,1},colub),
               convert(Array{Float64,1},obj),
               convert(Array{Float64,1},rowlb),
               convert(Array{Float64,1},rowub),
               sense)
end

#internal 
function upgradeProbType(m::MosekMathProgModel, pt :: Int)
    if     pt == MosekMathProgModel_QOQP
        if     m.probtype < pt
            m.probtype = pt
        elseif m.probtype > pt
            error("Incompatible problem type")
        end
    elseif pt == MosekMathProgModel_LINR
        # do nothing
    elseif pt == MosekMathProgModel_SOCP
        if     m.probtype == MosekMathProgModel_QOQP
            error("Incompatible problem type")
        elseif m.probtype < pt
            m.probtype = pt
        end
    elseif pt == MosekMathProgModel_SDP
        if     m.probtype == MosekMathProgModel_QOQP
            error("Incompatible problem type")
        elseif m.probtype < pt
            m.probtype = pt
        end
    else
        error("Invalid problem type")
    end
end


#internal
function complbk(bk,bl)
  if bl > -Inf
    if bk in [ MSK_BK_UP, MSK_BK_RA, MSK_BK_FX ]
      MSK_BK_RA
    else
      MSK_BK_LO
    end
  else
    if bk in [ MSK_BK_UP, MSK_BK_RA, MSK_BK_FX ]
      MSK_BK_UP
    else
      MSK_BK_FR
    end
  end
end

#internal
function compubk(bk,bu)
  if bu < Inf
    if bk in [ MSK_BK_LO, MSK_BK_RA, MSK_BK_FX ]
      MSK_BK_RA
    else
      MSK_BK_UP
    end
  else
    if bk in [ MSK_BK_LO, MSK_BK_RA, MSK_BK_FX ]
      MSK_BK_LO
    else
      MSK_BK_FR
    end
  end
end

#internal
function getLB(m::MosekMathProgModel,bk::Array{Int32,1},b::Array{Float64,1})
    map(1:m.numvar) do k
        if   (m.varmap[k] > 0) 
            eltbk = bk[m.varmap[k]]
            if eltbk in [ MSK_BK_FR, MSK_BK_UP ] -Inf
            else                                 b[m.varmap[k]]
            end
        else # SDP
            -Inf 
        end
    end
end

#internal
function getUB(m::MosekMathProgModel,bk::Array{Int32,1},b::Array{Float64,1})
    map(1:m.numvar) do k
        if   (m.varmap[k] > 0) 
            eltbk = bk[m.varmap[k]]
            if eltbk in [ MSK_BK_FR, MSK_BK_LO ] Inf
            else                                 b[m.varmap[k]]
            end
        else # SDP
            -Inf 
        end
    end
end

function getvarLB(m::MosekMathProgModel)
    bk,bl,bu = getvarboundslice(m.task,1,getnumvar(m.task))
    getLB(m,bk,bl)
end

function getvarUB(m::MosekMathProgModel)
    bk,bl,bu = getvarboundslice(m.task,1,getnumvar(m.task)+1)
    getUB(m,bk,bu)
end

function getconstrLB(m::MosekMathProgModel)
    bk,bl,bu = getconboundslice(m.task,1,getnumcon(m.task)+1)
    getLB(m,bk,bl)
end

function getconstrUB(m::MosekMathProgModel)
    bk,bl,bu = getconboundslice(m.task,1,getnumcon(m.task)+1)
    getLB(m,bk,bu)
end

# Note: Bounds on SDP variables are simply ignored without warning
function setvarLB!(m::MosekMathProgModel, collb)
    if m.numvar != length(collb)
        throw(MosekMathProgModelError("Bound vector has wrong size"))
    end    
    
    const idxs = find(i -> i > 0, m.varmap[1:m.numvar])
    local bnd  = collb[idxs]
    local subj = m.varmap[idxs]
    
    for i in 1:length(idxs)
        if m.binvarflag[idxs[i]]
            bnd[i] = max(bnd[i], 0.0)
        end
    end
    bk,bl,bu = getvarboundslice(m.task,1,getnumvar(m.task)+1)
    
    newbk = Int32[ complbk(bk[m.varmap[idxs[i]]],bnd[i]) for i=1:length(idxs) ]

    putvarboundlist(m.task, m.varmap[idxs], newbk, bnd, bu)
end

function setvarUB!(m::MosekMathProgModel, colub)
    if m.numvar != length(colub)
        throw(MosekMathProgModelError("Bound vector has wrong size"))
    end

    const idxs = find(i -> m.varmap[i] > 0, 1:m.numvar)
    local bnd  = colub[idxs]
    local subj = m.varmap[idxs]
    
    for i in 1:length(idxs)
        if m.binvarflag[idxs[i]]
            bnd[i] = min(bnd[i], 1.0)
        end
    end
    bk,bl,bu = getvarboundslice(m.task,1,getnumvar(m.task)+1)
    
    newbk = Int32[ compubk(bk[m.varmap[idxs[i]]],bnd[i]) for i=1:length(idxs) ]

    putvarboundlist(m.task, m.varmap[1:m.numvar], newbk, bl, bnd)
end


# Note: Bounds on any conic constraints are ignored without warning
function setconstrLB!(m::MosekMathProgModel, rowlb)
    if m.numcon != length(rowlb)
        throw(MosekMathProgModelError("Bound vector has wrong size"))
    end
    
    const idxs = find(i -> m.conmap[i] > 0 && m.conslack[i] == 0, 1:m.numcon)
    local bnd  = rowlb[idxs]
    local subj = m.conmap[idxs]
    
    bk,bl,bu = getconboundslice(m.task,1,getnumcon(m.task)+1)
    
    newbk = Int32[ complbk(bk[m.conmap[idxs[i]]],bnd[i]) for i=1:length(idxs) ]
    newbl = rowlb[m.conmap[idxs]] 
    putconboundlist(m.task, m.conmap[idxs], newbk, newbl, bu)
end


function setconstrUB!(m::MosekMathProgModel, rowub)
    if m.numcon != length(rowub)
        throw(MosekMathProgModelError("Bound vector has wrong size"))
    end
    
    const idxs = find(i -> m.conmap[i] > 0 && m.conslack[i] == 0, 1:m.numcon)
    local bnd  = rowub[idxs]
    local subj = m.conmap[idxs]

    bk,bl,bu = getconboundslice(m.task,1,getnumcon(m.task)+1)

    newbk = Int32[ compubk(bk[m.conmap[idxs[i]]],bnd[i]) for i=1:length(idxs) ]
    newbu = rowub[m.conmap[idxs]] 
    putconboundlist(m.task, m.conmap[idxs], newbk, bl, newbu)
    #println("putconboundlist: ",idxs,",",m.conmap[idxs])
    #println("\t",newbk)
    #println("\t",bl)
    #println("\t",newbu)
    #println(MSK_BK_LO," ",MSK_BK_UP," ",MSK_BK_FR," ",MSK_BK_RA)
end

# WARNING: Only gets non-PSD parts. Coefficients for PSD variables are returned as 0
function getobj(m::MosekMathProgModel)
    c = getc(m.task)
    local res = zeros(Float64,m.numvar)
    if all(i -> i > 0, m.varmap[1:m.numvar])
        res[:] = c[m.varmap[1:numvar]]
    else        
        const lidxs = find(i -> i > 0,m.varmap[1:m.numvar])
        res[idxs] = c[m.varmap[lidxs]]
        const baridxs = find(i -> i < 0,m.varmap[1:m.numvar])
    end
    res
end

# WARNING: Only sets non-PSD parts. Coefficients for PSD variables are not modified
function setobj!(m::MosekMathProgModel, obj::Array{Float64,1})
    if size(obj,1) != m.numvar
        throw(MosekMathProgModelError("Objective vector has wrong size"))
    end
    
    const idxs = find(i -> i > 0, m.varmap[1:m.numvar])
    const c = obj[idxs]
    putclist(m.task,m.varmap[idxs],c)
    nothing
end

function setobj!(m::MosekMathProgModel, obj)
  setobj(m,dense(float(obj)))
end

#internal
function ensureVarMapSize(m::MosekMathProgModel, numvar::Int)
    if (length(m.varmap) < numvar)
        newsz = max(1024,numvar,2*length(m.varmap))
        oldsz = length(m.varmap)
        resize!(m.varmap, newsz)
        resize!(m.barvarij, newsz)
        resize!(m.binvarflag, newsz)
        for i in (oldsz+1):newsz
            m.varmap[i] = 0
            m.barvarij[i] = 0
            m.binvarflag[i] = false
        end
    end
end

#internal
function ensureConMapSize(m::MosekMathProgModel, numcon::Int)
    if (length(m.conmap) < numcon)
        newsz = max(1024,numcon,2*length(m.conmap))
        oldsz = length(m.conmap)
        resize!(m.conmap,newsz)
        resize!(m.conslack,newsz)
        resize!(m.barconij,newsz)
        for i in (oldsz+1):newsz
            m.conmap[i] = 0
            m.conslack[i] = 0
            m.barconij[i] = 0
        end
    end
end

#internal
function ensureQuadConMapSize(m::MosekMathProgModel, numqcon::Int)
    if (length(m.qconmap) < numqcon)
        newsz = max(1024,numqcon,2*length(m.qconmap))
        resize!(m.qconmap,newsz)
    end
end


#internal
function ensureBarvarMapSize(m::MosekMathProgModel, numbarvar::Int)
    if (length(m.barvarmap) < numbarvar)
        newsz = max(1024,numbarvar,2*length(m.barvarmap))
        oldsz = length(m.barvarmap)
        resize!(m.barvarmap,newsz)
        for i in (oldsz+1):newsz
            m.barvarmap[i] = 0
        end
    end
end

#internal
function addUserVar(m::MosekMathProgModel, natidx::Int32)
    ensureVarMapSize(m,m.numvar+1)
    m.numvar += 1
    m.varmap[m.numvar]     = natidx
    m.barvarij[m.numvar]   = 0
    m.binvarflag[m.numvar] = false
    
    return m.numvar
end

#internal
# NOTE: When adding a barvar through this function it *will*not* get
# added to the variable vector - it is only accessable throught the low-level interface.
function addUserBarvar(m::MosekMathProgModel, natidx::Int32)
    ensureBarvarMapSize(m,m.numbarvar+1)
    m.numbarvar += 1

    barvaridx = m.numbarvar
    m.barvarmap[barvaridx] = natidx
    
    return barvaridx
end

#internal
function addUserCon(m::MosekMathProgModel, natidx::Int32)
    ensureConMapSize(m,m.numcon+1)
    m.numcon += 1
    m.conmap[m.numcon]   = natidx
    m.conslack[m.numcon] = 0
    m.barconij[m.numcon] = 0
    return m.numcon
end

function addUserQuadCon(m::MosekMathProgModel, natidx::Int32)
    ensureQuadConMapSize(m,m.numcon+1)
    m.numqcon += 1
    m.qconmap[m.numqcon] = natidx
    return m.numqcon
end


function addvar!(m::MosekMathProgModel, rowidx, rowcoef, collb, colub, objcoef)
    appendvars(m.task,1)
    varidx = getnumvar(m.task)
    bk = getBoundsKey(collb, colub)
    putvarbound(m.task,varidx,bk,collb,colub)
    putcj(m.task,varidx,objcoef)

    nzidxs = find(v -> v < 0.0 || v > 0.0, rowcoef);
    putaijlist(m.task,nzidxs,[ varidx for i in 1:size(nzidxs,1)],rowcoef[nzidxs])

    return addUserVar(m,varidx)
end

# NOTE: We silently ignore any SDP column
function addconstr!(m::MosekMathProgModel, colidx, colcoef, lb, ub)
    appendcons(m.task,1)
    constridx = getnumcon(m.task)

    local idxs = find(i -> m.varmap[i] > 0, colidx)

    putarow(m.task,constridx,m.varmap[colidx[idxs]],colcoef[idxs])

    bk = getBoundsKey(lb, ub)
    putconbound(m.task,constridx,bk,lb,ub)

    return addUserCon(m,constridx)
end

updatemodel!(m::MosekMathProgModel) = nothing

function setsense!(m::MosekMathProgModel,sense::Symbol)
  if     sense == :Min
    putobjsense(m.task, MSK_OBJECTIVE_SENSE_MINIMIZE)
  elseif sense == :Max
    putobjsense(m.task, MSK_OBJECTIVE_SENSE_MAXIMIZE)
  else
   throw(MosekMathProgModelError("Invalid objective sense"))
  end
  nothing
end

function getsense(m::MosekMathProgModel)
  sense = getobjsense(m.task)
  if sense == MSK_OBJECTIVE_SENSE_MINIMIZE
    :Min
  elseif sense == MSK_OBJECTIVE_SENSE_MAXIMIZE
    :Max
  else
    None
  end
end

numvar(       m::MosekMathProgModel) = m.numvar
numconstr(    m::MosekMathProgModel) = m.numcon+m.numqcon
numlinconstr( m::MosekMathProgModel) = m.numcon
numquadconstr(m::MosekMathProgModel) = m.numqcon
optimize!(m::MosekMathProgModel) = 
    begin
        optimize(m.task)
    end
# function optimize!(m::MosekMathProgModel)  optimize(m.task); writedata(m.task,"mskprob.opf") end

function getsoldef(m::MosekMathProgModel)
  for sol in (MSK_SOL_ITG, MSK_SOL_BAS, MSK_SOL_ITR)
    if solutiondef(m.task,sol) && getsolsta(m.task,sol) != MSK_SOL_STA_UNKNOWN
      return sol
    end
  end
  return -1
end


# NOTE: What are the 'legal' values to return?
# Another NOTE: status seems to mash together the problem status,
# the solution status and the solver status. I'll try to cope
# in some sensible manner.
function status(m::MosekMathProgModel)
  soldef = getsoldef(m)
  if soldef < 0 return :Unknown end
  prosta = getprosta(m.task,soldef)
  solsta = getsolsta(m.task,soldef)

  if     solsta == MSK_SOL_STA_UNKNOWN
    :Unknown
  elseif solsta == MSK_SOL_STA_DUAL_FEAS ||
         solsta == MSK_SOL_STA_PRIM_FEAS ||
         solsta == MSK_SOL_STA_NEAR_PRIM_FEAS ||
         solsta == MSK_SOL_STA_NEAR_DUAL_FEAS ||
         solsta == MSK_SOL_STA_PRIM_AND_DUAL_FEAS ||
         solsta == MSK_SOL_STA_NEAR_PRIM_AND_DUAL_FEAS
    :Unknown
  elseif solsta == MSK_SOL_STA_DUAL_INFEAS_CER ||
         solsta == MSK_SOL_STA_NEAR_DUAL_INFEAS_CER
    :Unbounded
  elseif solsta == MSK_SOL_STA_PRIM_INFEAS_CER ||
         solsta == MSK_SOL_STA_NEAR_PRIM_INFEAS_CER
    :Infeasible
  elseif solsta == MSK_SOL_STA_OPTIMAL ||
         solsta == MSK_SOL_STA_NEAR_OPTIMAL ||
         solsta == MSK_SOL_STA_INTEGER_OPTIMAL ||
         solsta == MSK_SOL_STA_NEAR_INTEGER_OPTIMAL
    :Optimal
  else
    error("Internal value error")
  end
end

function getobjval(m::MosekMathProgModel)
  soldef = getsoldef(m)
  if soldef < 0 return NaN end
  getprimalobj(m.task,soldef)
end

# NOTE: I am not entirely sure how to implement this... If the solution status
# is feasible for an integer problem, then the objective value is the best
# known bound.
getobjbound(m::MosekMathProgModel) = getdouinf(m.task,MSK_DINF_MIO_OBJ_BOUND)

# NOTE: We may want different warmstart values interior point and simplex
function setwarmstart!(m::MosekMathProgModel, xx::Array{Float64,1})
    idxs = find(i -> i > 0, m.varmap)
    let subj = m.varmap[idxs]
        for whichsol in [ MSK_SOL_BAS,MSK_SOL_ITR,MSK_SOL_ITG ]
            let newxx = 
                if solutiondef(m.task,whichsol)
                    getxx(m.task,whichsol)
                else
                    zeros(Float64,getnumvar(m.task))
                end
                newxx[subj] = xx[idxs]
                putxx(m.task,whichsol,xx)
            end
        end
    end
end

function getsolution(m::MosekMathProgModel)
  soldef = getsoldef(m)
  if soldef < 0 
      throw(MosekMathProgModelError("No solution available"))
  end
  solsta = getsolsta(m.task,soldef)
  if solsta in [ MSK_SOL_STA_OPTIMAL, MSK_SOL_STA_PRIM_FEAS, MSK_SOL_STA_PRIM_AND_DUAL_FEAS, MSK_SOL_STA_NEAR_OPTIMAL, MSK_SOL_STA_NEAR_PRIM_FEAS, MSK_SOL_STA_NEAR_PRIM_AND_DUAL_FEAS, MSK_SOL_STA_INTEGER_OPTIMAL, MSK_SOL_STA_NEAR_INTEGER_OPTIMAL ]
      if m.numbarvar > 0
          const xx = getxx(m.task,soldef)
          const barx      = [ getbarxj(m.task,soldef,j) for j in 1:getnumbarvar(m.task) ]
          const barvardim = [ getdimbarvarj(m.task,j)    for j in 1:getnumbarvar(m.task) ]
          
          map(1:m.numvar) do k
              const j = m.varmap[k]
              if   (j > 0) xx[j]
              else         barx[-j][m.barvarij[k]]
              end
          end
      else
          getxx(m.task,soldef)[m.varmap[1:m.numvar]]
      end
  else
    throw(MosekMathProgModelError("No solution available"))
  end
end

function getconstrsolution(m::MosekMathProgModel)
  soldef = getsoldef(m)
  if soldef < 0 throw(MosekMathProgModelError("No solution available")) end
  solsta = getsolsta(m.task,soldef)
  if solsta in [ MSK_SOL_STA_OPTIMAL, MSK_SOL_STA_PRIM_FEAS, MSK_SOL_STA_PRIM_AND_DUAL_FEAS, MSK_SOL_STA_NEAR_OPTIMAL, MSK_SOL_STA_NEAR_PRIM_FEAS, MSK_SOL_STA_NEAR_PRIM_AND_DUAL_FEAS, MSK_SOL_STA_INTEGER_OPTIMAL, MSK_SOL_STA_NEAR_INTEGER_OPTIMAL ]
      xc = getxc(m.task,soldef)
      xx = getxx(m.task,soldef)

      if m.numbarvar > 0
          const barx = [ getbarxj(m.task,soldef,j)   for j in 1:getnumbarvar(m.task) ]
          const barvardim = [ getdimbarvar(m.task,j) for j in 1:getnumbarvar(m.task) ]

          map(1:m.numcon) do k
              const j = m.conslack[k]
              if     (j == 0) xc[m.conmap[k]]
              elseif (j >  0) xx[j] # conic slack
              else            barx[-j][m.barconij[k]]
              end
          end
      else
          map(1:m.numcon) do k
              const j = m.conslack[k]
              if     (j == 0) xc[m.conmap[k]]
              else            xx[j] # conic slack
              end
          end
      end
  else
    throw(MosekMathProgModelError("No solution available"))
  end
end


# internal
function getvardual(m::MosekMathProgModel,soldef::Int32)
    sux = getsux(m.task,soldef)
    slx = getslx(m.task,soldef)
    snx = if (soldef == MSK_SOL_ITR) getsnx(m.task,soldef) else zeros(Float64,length(slx)) end

    s = slx-sux+snx

    if m.numbarvar > 0
        const bars = [ getbarsj(m.task,soldef,j) for j in 1:getnumbarvar(m.task) ]
        const barvardim = Int32[ getdimbarvarj(m.task,j) for j in 1:getnumbarvar(m.task) ]

        map(1:m.numvar) do k
            const j = m.varmap[k]

            if   (j > 0) s[j]
            else         bars[-j][m.barvarij[k]]
            end
        end
    else
        s[m.varmap[1:m.numvar]]
    end
end

#internal
function getcondual(m::MosekMathProgModel,soldef::Int32)
    let snx  = if (soldef == MSK_SOL_ITR) getsnx(m.task,soldef) else zeros(Float64,getnumvar(m.task)) end,
        y    = gety(m.task,soldef),
        bars = map(j -> getbarsj(m.task,soldef,j), 1:getnumbarvar(m.task)),
        bkc,blc,buc = getconboundslice(m.task,1,getnumcon(m.task)+1)

        map(1:m.numcon) do k
            const j = m.conslack[k]
            const i = m.conmap[k]
            if     (j == 0)
                y[i]
            elseif (j >  0)
                snx[j]
            else
                bars[-j][m.barconij[k]]
            end
        end
    end
end


function getreducedcosts(m::MosekMathProgModel)
  soldef = getsoldef(m)
  if soldef < 0 throw(MosekMathProgModelError("No solution available")) end
  solsta = getsolsta(m.task,soldef)

  if solsta in [MSK_SOL_STA_OPTIMAL, 
                MSK_SOL_STA_DUAL_FEAS, 
                MSK_SOL_STA_PRIM_AND_DUAL_FEAS, 
                MSK_SOL_STA_NEAR_OPTIMAL, 
                MSK_SOL_STA_NEAR_DUAL_FEAS,
                MSK_SOL_STA_NEAR_PRIM_AND_DUAL_FEAS ]
      if (getsense(m) == :Min) 
          getvardual(m,soldef)
      else
          getvardual(m,soldef)
      end
  else
      throw(MosekMathProgModelError("No solution available"))
  end
end


function getconstrduals(m::MosekMathProgModel)
    soldef = getsoldef(m)
    if soldef < 0 throw(MosekMathProgModelError("No solution available")) end
    solsta = getsolsta(m.task,soldef)
    if solsta in [MSK_SOL_STA_OPTIMAL, 
                  MSK_SOL_STA_DUAL_FEAS, 
                  MSK_SOL_STA_PRIM_AND_DUAL_FEAS, 
                  MSK_SOL_STA_NEAR_OPTIMAL, 
                  MSK_SOL_STA_NEAR_DUAL_FEAS, 
                  MSK_SOL_STA_NEAR_PRIM_AND_DUAL_FEAS, 
                  MSK_SOL_STA_PRIM_INFEAS_CER]
        getcondual(m,soldef)
    else
        throw(MosekMathProgModelError("No solution available"))
    end
end



function getinfeasibilityray(m::MosekMathProgModel)
    soldef = getsoldef(m)
    if soldef < 0 throw(MosekMathProgModelError("No solution available")) end
    solsta = getsolsta(m.task,soldef)
    if solsta in [ MSK_SOL_STA_PRIM_INFEAS_CER, MSK_SOL_STA_NEAR_PRIM_INFEAS_CER ]
        -getvardual(m,soldef)
    else
        throw(MosekMathProgModelError("No solution available"))
    end
end

function getunboundedray(m::MosekMathProgModel)
  soldef = getsoldef(m)
  if soldef < 0 throw(MosekMathProgModelError("No solution available")) end
  solsta = getsolsta(m.task,soldef)
  if solsta in [ MSK_SOL_STA_DUAL_INFEAS_CER, MSK_SOL_STA_NEAR_DUAL_INFEAS_CER ]
      const barx = [ getbarxj(m.task,soldef,j) for j in 1:getnumbarvar(m.task) ]
      const xx   = getxx(m.task,soldef)

      map(1:m.numvar) do k
          const j = m.varmap[k]
          if   (j > 0) xx[j] # conic slack
          else         barx[-j][m.barvarij[k]]
          end
      end
  else
      throw(MosekMathProgModelError("No solution available"))
  end
end

getrawsolver(m::MosekMathProgModel) = m.task


function getvartype(m::MosekMathProgModel)
    numvar = getnumvar(m.task)
    vartypes = Array(Int32,numvar)
    getvartypelist(m.task,Int32[1:numvar],vartypes)

    return [(if (vartypes == MSK_VAR_TYPE_CONT)
               :Cont
             else
               bk,bl,bu = getvarbound(m.task,j)
               if (bk == MSK_BK_RA && bl == 0.0 && bu - 1.0 == 0.0) # hmm... not sure '==' here is ok.
                 :Bin
               else
                 :Int
               end
             end
             )
            for j in m.varmap ]
end

# NOTE: Var types for semidefinite variables are ignored.
function setvartype!(m::MosekMathProgModel, intvarflag :: Array{Symbol,1})
  n = min(length(intvarflag),m.numvar)
  all(x->in(x,[:Cont,:Int,:Bin]), intvarflag) || error("Invalid variable type present")
  #any(k -> intvarflag[k] in [:Bin,:Int] && m.varmap[k] < 0, 1:n) || error("Invalid variable type for semidefinite variable element")

  m.binvarflag[1:n] = map(f -> f == :Bin, intvarflag[1:n])

  idxs = find(i -> m.varmap[i] > 0,1:n)
  putvartypelist(m.task,m.varmap[idxs],Int32[ (if (c == :Cont) MSK_VAR_TYPE_CONT else MSK_VAR_TYPE_INT end) for c in intvarflag[idxs] ])

  for k in find(f -> f == :Bin, intvarflag)
      local j = m.varmap[k]

      bk,bl,bu = getvarbound(m.task,j)
      bl =
          if (bk in [ MSK_BK_LO, MSK_BK_RA ]) max(bl,0.0)
          else                                0.0
          end
      bu =
          if (bk in [ MSK_BK_UP, MSK_BK_RA ]) min(bu,1.0)
          elseif (bk == MSK_BK_FX)            0.0
          else                                1.0
      end
      putvarbound(m.task,j,MSK_BK_RA,bl,bu)
  end
end

function getintvarflag(m::MosekMathProgModel)
    local vt = getvartypelist(m.task,[1:getnumvar(m.task)])
    map(1:m.numvar) do k
        j = m.varmap[k]
        if     j < 0                      :Cont
        elseif m.binvarflag[k]            :Bin
        elseif vt[j] == MSK_VAR_TYPE_CONT :Int
        end
    end
end

include("MosekLowLevelQCQO.jl")
include("MosekLowLevelSDP.jl")
include("MosekNLPSolverInterface.jl")
include("MosekConicInterface.jl")

end
