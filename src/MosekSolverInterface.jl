module MosekMathProgSolverInterface
using ..Mosek

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

type MosekMathProgModel <: AbstractMathProgModel
  task :: Mosek.MSKtask
  probtype :: Int

  numvar :: Int64
  varmap :: Array{Int32,1} # map user defined variables to MOSEK variables
  vtypemap :: Array{Symbol,1}

  numcon :: Int64
  conmap :: Array{Int32,1} # map user defined constraints to MOSEK variables
  options
end

immutable MosekSolver <: AbstractMathProgSolver
  options
end
MosekSolver(;kwargs...) = MosekSolver(kwargs)

type MosekMathProgModelError
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
                         MosekMathProgModel_LINR,
                         0,
                         Array(Int32,1024),
                         Symbol[],
                         0,
                         Array(Int32,1024),
                         s.options)
  loadoptions!(m)
  return m
end

# NOTE: This method will load data into an existing task, but
# it will not necessarily reset all exitsing data, depending on the
# file read (e.g. reading an MPS will not reset parameters)
# Also, auxilary variables and cones will not be correctly mapped.
function loadproblem!(m:: MosekMathProgModel, filename:: String)
  readdata(m.task, filename)
  m.numvar = getnumvar(m.task)
  m.numcon = getnumcon(m.task)
  m.varmap = Int32[1:m.numvar]
  m.conmap = Int32[1:m.numcon]
  vtype = getvartypelist(m.task, m.varmap)
  m.vtypemap = map(vtype) do c
    if c == MSK_VAR_TYPE_CONT
      return :Cont
    elseif c == MSK_VAR_TYPE_INT
      return :Int
    else
      error("Unrecognized variable type")
    end
  end
end

function writeproblem(m:: MosekMathProgModel, filename:: String)
  writedata(m.task, filename)
end


function loadproblem!(m::     MosekMathProgModel,
                      A::     SparseMatrixCSC,
                      collb:: Array{Float64,1},
                      colub:: Array{Float64,1},
                      obj::   SparseMatrixCSC,
                      rowlb:: Array{Float64,1},
                      rowub:: Array{Float64,1},
                      sense:: Symbol)
  Mosek.deletetask(m.task)
  m.task = maketask(Mosek.msk_global_env)
  loadoptions!(m)

  nrows,ncols = size(A)
  if ncols != length(collb) ||
     ncols != length(colub) ||
     ncols != size(obj,1)   ||
     nrows != length(rowlb) ||
     nrows != length(rowub) ||
     size(obj,2) != 1
    throw(MosekMathProgModelError("Inconsistent data dimensions"))
  end

  appendvars(m.task, ncols)
  appendcons(m.task, nrows)

  m.numvar = ncols
  m.numcon = nrows
  m.varmap = Int32[1:m.numvar]
  m.conmap = Int32[1:m.numcon]
  m.vtypemap = fill(:Cont,ncols)

  # input coefficients
  putclist(m.task, obj.rowval, obj.nzval)
  putacolslice(m.task, 1, ncols+1, A.colptr[1:ncols], A.colptr[2:ncols+1], A.rowval, A.nzval)
  setsense!(m, sense)

  # input bounds
  putvarboundslice(m.task, 1, ncols+1, Int32[ MSK_BK_RA for i=1:ncols], collb, colub)
  putconboundslice(m.task, 1, nrows+1, Int32[ MSK_BK_RA for i=1:nrows], rowlb, rowub)
  nothing
end

function loadproblem! (m::     MosekMathProgModel,
                       A,
                       collb,
                       colub,
                       obj,
                       rowlb,
                       rowub,
                       sense)
  loadproblem!(m,sparse(float(A)),float(collb),float(colub),sparse(float(obj)),float(rowlb),float(rowub),sense)
end


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

function getvarLB(m::MosekMathProgModel)
  bkx,blx,bux = getvarboundslice(m.task,1,getnumvar(m.task)+1)
  Float64[ (if bkx[i] in [ MSK_BK_FR, MSK_BK_UP ] -Inf else blx[i] end) for i=m.varmap[1:m.numvar] ]
end

function getvarUB(m::MosekMathProgModel)
  bkx,blx,bux = getvarboundslice(m.task,1,getnumvar(m.task)+1)
  Float64[ (if bkx[i] in [ MSK_BK_FR, MSK_BK_LO ]  Inf else bux[i] end) for i=m.varmap[1:m.numvar] ]
end

function getconstrLB(m::MosekMathProgModel)
  bkc,blc,buc = getconboundslice(m.task,1,getnumcon(m.task)+1)
  Float64[ (if bkc[i] in [ MSK_BK_FR, MSK_BK_UP ] -Inf else blc[i] end) for i=m.conmap[1:m.numcon] ]
end

function getconstrUB(m::MosekMathProgModel)
  bkc,blc,buc = getconboundslice(m.task,1,getnumcon(m.task)+1)
  Float64[ (if bkc[i] in [ MSK_BK_FR, MSK_BK_LO ]  Inf else buc[i] end) for i=m.conmap[1:m.numcon] ]
end

function setvarLB!(m::MosekMathProgModel, collb)
  if m.numvar != length(collb)
    throw(MosekMathProgModelError("Bound vector has wrong size"))
  end

  bnd = copy(collb)
  for i in 1:length(collb)
    if m.vtypemap[i] == :Bin
      bnd[i] = min(collb[i], 1.0)
   end
  end
  bk,bl,bu = getvarboundslice(m.task,1,getnumvar(m.task)+1)

  newbk = [ complbk(bk[i],bnd[i]) for i=m.varmap[1:m.numvar] ]
  newbu = [ bu[i] for i=m.varmap[1:m.numvar] ]

  putvarboundlist(m.task, m.varmap[1:m.numvar], bk, bnd, bu)
end

function setvarUB!(m::MosekMathProgModel, colub)
  if m.numvar != length(colub)
    throw(MosekMathProgModelError("Bound vector has wrong size"))
  end

  bnd = copy(colub)
  for i in 1:length(colub)
    if m.vtypemap[i] == :Bin
      bnd[i] = max(colub[i], 0.0)
   end
  end
  bk,bl,bu = getvarboundslice(m.task,1,getnumvar(m.task)+1)

  newbk = [ compubk(bk[i],bnd[i]) for i=m.varmap[1:m.numvar] ]
  newbl = [ bl[i] for i=m.varmap[1:m.numvar] ]

  putvarboundlist(m.task, m.varmap[1:m.numvar], bk, bl, bnd)
end


function setconstrLB!(m::MosekMathProgModel, rowlb)
  if m.numcon != length(rowlb)
    throw(MosekMathProgModelError("Bound vector has wrong size"))
  end

  bk,bl,bu = getconboundslice(m.task,1,getnumcon(m.task)+1)

  newbk = [ complbk(bk[i],rowlb[i]) for i=m.conmap[1:m.numcon] ]
  newbu = [ bu[i] for i=m.conmap[1:m.numcon] ]
  putconboundlist(m.task, m.conmap[1:m.numcon], bk, rowlb, bu)
end


function setconstrUB!(m::MosekMathProgModel, rowub)
  if m.numcon != length(rowub)
    throw(MosekMathProgModelError("Bound vector has wrong size"))
  end

  bk,bl,bu = getconboundslice(m.task,1,getnumcon(m.task)+1)

  newbk = [ compubk(bk[i],rowub[i]) for i=m.conmap[1:m.numcon] ]
  newbl = [ bl[i] for i=m.conmap[1:m.numcon] ]
  putconboundlist(m.task, m.conmap[1:m.numcon],bk,bl,rowub)
end

function getobj(m::MosekMathProgModel)
  c = getc(m.task)
  [ c[i] for i=m.varmap[1:m.numvar] ]
end

function setobj!(m::MosekMathProgModel, obj::Array{Float64,1})
  if size(obj,1) != m.numvar
    throw(MosekMathProgModelError("Objective vector has wrong size"))
  end

  putclist(m.task,m.varmap[1:m.numvar],obj)
end

function setobj!(m::MosekMathProgModel, obj)
  setobj(m,dense(float(obj)))
end

function ensureVarMapSize(m::MosekMathProgModel, numvar::Int32)
  if (length(m.varmap) < numvar)
    newvarmaplen = max(numvar,2*length(m.varmap))
    newvarmap = Array(Int32,newvarmaplen)
    newvarmap[1:m.numvar] = m.varmap[1:m.numvar]
    m.varmap = newvarmap
  end
end

function ensureConMapSize(m::MosekMathProgModel, numcon::Int32)
  if (length(m.conmap) < numcon)
    newconmaplen = max(numcon,2*length(m.conmap))
    newconmap = Array(Int32,newconmaplen)
    newconmap[1:m.numcon] = m.conmap[1:m.numcon]
    m.conmap = newconmap
  end
end

function addUserVar(m::MosekMathProgModel, natidx::Int32)
  ensureVarMapSize(m,convert(Int32,m.numvar+1))
  m.numvar += 1
  m.varmap[m.numvar] = natidx
  return m.numvar
end

function addUserCon(m::MosekMathProgModel, natidx::Int32)
  ensureConMapSize(m,convert(Int32,m.numcon+1))
  m.numcon += 1
  m.conmap[m.numcon] = natidx
  return m.numcon
end

function addvar!(m::MosekMathProgModel, rowidx, rowcoef, collb, colub, objcoef)
    appendvars(m.task,1)
    varidx = getnumvar(m.task)
    bk = getBoundsKey(collb, colub)
    putvarbound(m.task,varidx,bk,collb,colub)
    putcj(m.task,varidx,objcoef)

    return addUserVar(m,varidx)
end

function addconstr!(m::MosekMathProgModel, colidx, colcoef, lb, ub)
    appendcons(m.task,1)
    constridx = getnumcon(m.task)
    putarow(m.task,constridx,[ m.varmap[i] for i in colidx ],colcoef)
    bk = getBoundsKey(lb, ub)
    putconbound(m.task,constridx,bk,lb,ub)

    return addUserCon(m,constridx)
end

updatemodel!(m::MosekMathProgModel) = nothing

function setsense!(m::MosekMathProgModel,sense)
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

numvar(m::MosekMathProgModel) = m.numvar
numconstr(m::MosekMathProgModel) = m.numcon
optimize!(m::MosekMathProgModel) = optimize(m.task)
# function optimize!(m::MosekMathProgModel)  optimize(m.task); writedata(m.task,"mskprob.opf") end


function getsoldef(m::MosekMathProgModel)
  if solutiondef(m.task,MSK_SOL_ITG) MSK_SOL_ITG
  elseif solutiondef(m.task,MSK_SOL_BAS) MSK_SOL_BAS
  elseif solutiondef(m.task,MSK_SOL_ITR) MSK_SOL_ITR
  else -1
  end
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

function getsolution(m::MosekMathProgModel)
  soldef = getsoldef(m)
  if soldef < 0 throw(MosekMathProgModelError("No solution available"))
  end
  solsta = getsolsta(m.task,soldef)
  if solsta in [ MSK_SOL_STA_OPTIMAL, MSK_SOL_STA_PRIM_FEAS, MSK_SOL_STA_PRIM_AND_DUAL_FEAS, MSK_SOL_STA_NEAR_OPTIMAL, MSK_SOL_STA_NEAR_PRIM_FEAS, MSK_SOL_STA_NEAR_PRIM_AND_DUAL_FEAS, MSK_SOL_STA_INTEGER_OPTIMAL, MSK_SOL_STA_NEAR_INTEGER_OPTIMAL ]
    xx = getxx(m.task,soldef)
    Float64[ xx[i] for i=m.varmap[1:m.numvar] ]
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
    Float64[ xc[i] for i=m.conmap[1:m.numcon] ]
  else
    throw(MosekMathProgModelError("No solution available"))
  end
end

function getreducedcosts(m::MosekMathProgModel)
  soldef = getsoldef(m)
  if soldef < 0 throw(MosekMathProgModelError("No solution available")) end
  solsta = getsolsta(m.task,soldef)

  if solsta in [ MSK_SOL_STA_OPTIMAL, MSK_SOL_STA_DUAL_FEAS, MSK_SOL_STA_PRIM_AND_DUAL_FEAS, MSK_SOL_STA_NEAR_OPTIMAL, MSK_SOL_STA_NEAR_DUAL_FEAS, MSK_SOL_STA_NEAR_PRIM_AND_DUAL_FEAS ]
    sux = getsux(m.task,soldef)
    slx = getslx(m.task,soldef)
    snx = if (soldef == MSK_SOL_ITR) getsnx(m.task,soldef) else zeros(Float64,length(slx)) end

    if getsense(m) == :Min
      Float64[   sux[i] - slx[i] + snx[i]  for i=m.varmap[1:m.numvar] ]
    else
      Float64[ -(sux[i] - slx[i] + snx[i]) for i=m.varmap[1:m.numvar] ]
    end
  else
    throw(MosekMathProgModelError("No solution available"))
  end
end

function getconstrduals(m::MosekMathProgModel)
  soldef = getsoldef(m)
  if soldef < 0 throw(MosekMathProgModelError("No solution available")) end
  solsta = getsolsta(m.task,soldef)
  if solsta in [ MSK_SOL_STA_OPTIMAL, MSK_SOL_STA_DUAL_FEAS, MSK_SOL_STA_PRIM_AND_DUAL_FEAS, MSK_SOL_STA_NEAR_OPTIMAL, MSK_SOL_STA_NEAR_DUAL_FEAS, MSK_SOL_STA_NEAR_PRIM_AND_DUAL_FEAS, MSK_SOL_STA_PRIM_INFEAS_CER ]
    y = gety(m.task,soldef)
    Float64[ y[i] for i=m.conmap[1:m.numcon] ]
  else
    throw(MosekMathProgModelError("No solution available"))
  end
end

function getinfeasibilityray(m::MosekMathProgModel)
  soldef = getsoldef(m)
  if soldef < 0 throw(MosekMathProgModelError("No solution available")) end
  solsta = getsolsta(m.task,soldef)
  if solsta in [ MSK_SOL_STA_PRIM_INFEAS_CER, MSK_SOL_STA_NEAR_PRIM_INFEAS_CER ]
    sux = getsux(m.task,soldef)
    slx = getslx(m.task,soldef)
    snx = if (soldef == MSK_SOL_ITR) getsnx(m.task,soldef) else zeros(Float64,length(slx)) end

    Float64[ sux[i] - slx[i] + snx[i] for i=m.varmap[1:m.numvar] ]
  else
    throw(MosekMathProgModelError("No solution available"))
  end
end

function getunboundedray(m::MosekMathProgModel)
  soldef = getsoldef(m)
  if soldef < 0 throw(MosekMathProgModelError("No solution available")) end
  solsta = getsolsta(m.task,soldef)
  if solsta in [ MSK_SOL_STA_DUAL_INFEAS_CER, MSK_SOL_STA_NEAR_DUAL_INFEAS_CER ]
    xx = getxx(m.task,soldef)
    Float64[ xx[i] for i=m.varmap[1:m.numvar] ]
  else
    throw(MosekMathProgModelError("No solution available"))
  end
end

getrawsolver(m::MosekMathProgModel) = m.task

function setvartype!(m::MosekMathProgModel, vartype :: Array{Symbol,1})
  numvar = getnumvar(m.task)
  n = min(length(vartype),numvar)
  all(x->in(x,[:Cont,:Int,:Bin]), vartype) || error("Invalid variable type present")
  m.vtypemap = copy(vartype)
  putvartypelist(m.task,m.varmap[1:m.numvar],Int32[ (if (c == :Cont) MSK_VAR_TYPE_CONT else MSK_VAR_TYPE_INT end) for c in vartype ])
end

getvartype(m::MosekMathProgModel) = copy(m.vtypemap)

# QCQO interface, so far only non-conic.

function setquadobj!(m::MosekMathProgModel, rowidx,colidx,quadval)
  qosubi = [ m.varmap[i] for i=rowidx ]
  qosubj = [ m.varmap[i] for i=colidx ]
  qoval  = copy(quadval)
  for i=1:length(rowidx)
    if qosubj[i] > qosubi[i]
      cj = qosubj[i]
      qosubj[i] = qosubi[i]
      qosubi[i] = cj
    end
  end
  putqobj(m.task,qosubi,qosubj,convert(Array{Float64},qoval))
end


# Note:
#  If the quadratic terms define a quadratic cone, the linear terms, sense and rhs are ignored.
function addquadconstr!(m::MosekMathProgModel, linearidx, linearval, quadrowidx, quadcolidx, quadval, sense, rhs)
  subj = Int32[ m.varmap[i] for i=linearidx ]
  valj = linearval
  qcksubi = Int32[ m.varmap[i] for i=quadrowidx ]
  qcksubj = Int32[ m.varmap[i] for i=quadcolidx ]
  qckval  = quadval

  # detect SOCP form

  ct,x =
    begin
      let num_posonediag = 0,
          offdiag_idx    = 0,
          negdiag_idx    = 0
        for i=1:length(qcksubi)
          if qcksubi[i] == qcksubj[i]
            if abs(qckval[i]-1.0) < 1e-12
              num_posonediag += 1
            elseif abs(qckval[i]+1.0) < 1e-12
              negdiag_idx = i
            end
          elseif qcksubi[i] != qcksubj[i]
            if abs(qckval[i]+1.0) < 1e-12
              offdiag_idx = i
            end
          end
        end

        if num_posonediag == length(qcksubj)-1 && negdiag_idx > 0
          x = Array(Int64,length(qcksubj))
          x[1] = qcksubj[negdiag_idx]
          for i=1:negdiag_idx-1 x[i+1] = qcksubj[i] end
          for i=negdiag_idx+1:length(qcksubj) x[i] = qcksubj[i] end

          MSK_CT_QUAD, x
        elseif num_posonediag == length(qcksubj)-1 && offdiag_idx > 0
          x = Array(Int64,length(qcksubj)+1)
          x[1] = qcksubi[offdiag_idx]
          x[2] = qcksubj[offdiag_idx]
          for i=1:offdiag_idx-1 x[i+2] = qcksubj[i] end
          for i=offdiag_idx+1:length(qcksubj) x[i+1] = qcksubj[i] end

          MSK_CT_RQUAD, x
        else
          -1,()
        end
      end
    end

  if     ct == MSK_CT_QUAD || ct == MSK_CT_RQUAD
    if     m.probtype == MosekMathProgModel_QOQP
      throw(MosekMathProgModelError("Cannot mix conic and quadratic terms"))
    elseif m.probtype == MosekMathProgModel_LINR
      m.probtype = MosekMathProgModel_SOCP
    end
    # SOCP and SDP can be mixed, SDP includes SOCP

    n = length(x)

    nvar  = getnumvar(m.task)+1
    ncon  = getnumcon(m.task)+1

    appendvars(m.task,n) # create aux variable z
    appendcons(m.task,n)

    # z in R^n, free
    z = nvar
    putvarboundslice(m.task, nvar,nvar+n, Int32[ MSK_BK_FR for i=1:n ] , zeros(n),zeros(n))

    cof = Array(Float64,2,n)
    cof[1,:] =  1.0
    cof[2,:] = -1.0
    if ct == MSK_CT_RQUAD
      cof[1,1] =  0.5;
    end

    subj = Array(Int32,2,n)
    subj[1,:] = x
    subj[2,:] = z:(z+n-1)

    ptrb = [1:n]*2 .- 1
    ptre = [1:n]*2 .+ 1

    # 0.5 x_1 - z_1 = 0
    #     x_i - z_i = 0, i=2..n
    putarowslice(m.task, ncon, ncon+n, ptrb, ptre, subj[:], cof[:] )
    putconboundslice(m.task, ncon, ncon+n, Int32[ MSK_BK_FX for i=1:n ], zeros(n), zeros(n) )

    appendcone(m.task, ct, 0.0, [z:z+n-1])

    # we add a dummy constraint to make sure that there is a place-holder for the constarint. The value is always 0.
    appendcons(m.task,1)
    dummycon = getnumcon(m.task)
    putconbound(m.task,dummycon, MSK_BK_FX, 0.0, 0.0)

    addUserCon(m,dummycon)
  else
    if     m.probtype == MosekMathProgModel_SOCP || m.probtype == MosekMathProgModel_SDP
      throw(MosekMathProgModelError("Cannot mix conic and quadratic terms"))
    elseif m.probtype == MosekMathProgModel_LINR
      m.probtype = MosekMathProgModel_QOQP
    end

    for i=1:length(quadrowidx)
      if qcksubj[i] > qcksubi[i]
        cj = qcksubj[i]
        qcksubj[i] = qcksubi[i]
        qcksubi[i] = cj
      elseif qcksubj[i] == qcksubi[i]
        qckval[i] = qckval[i] * 2
      end
    end

    k = getnumcon(m.task)+1
    appendcons(m.task,1)
    considx = getnumcon(m.task)

    putarow(m.task, k, convert(Array{Float64},subj), convert(Array{Float64},valj))
    putqconk(m.task,k, qcksubi,qcksubj,qckval)

    if sense == '<'
      putconbound(m.task,k,MSK_BK_UP, -Inf,convert(Float64,rhs))
    elseif sense == '>'
      putconbound(m.task,k,MSK_BK_LO, convert(Float64,rhs),Inf)
    else
      putconbound(m.task,k,MSK_BK_FR, -Inf,Inf)
    end

    addUserCon(m,considx)
  end
end

#####
# SDP
#####
function sparseToSparseTriple(mat::SparseMatrixCSC)
    if issym(mat)
        nnz = convert(Int64, (countnz(mat)+countnz(diag(mat))) / 2)
    elseif istriu(mat)
        nnz = nfilled(mat)
    else
        error("Matrix must be symmetric or upper triangular")
    end
    II = Array(Cint, nnz)
    JJ = Array(Cint, nnz)
    VV = Array(Cdouble, nnz)
    m, n = size(mat)
    k = 0
    colptr::Vector{Int64} = mat.colptr
    nzval::Vector{Float64} = mat.nzval

    for i = 1:n
        qi = convert(Cint, i)
        for j = colptr[i]:(colptr[i+1]-1)
            qj = convert(Cint, mat.rowval[j])
            if qi <= qj && nzval[j] != 0.0
                k += 1
                II[k] = qj
                JJ[k] = qi
                VV[k] = nzval[j]
            end
        end
    end
    return II,JJ,VV
end

function denseToSparseTriple(mat::Matrix)
    if issym(mat)
        nnz = convert(Int64, (countnz(mat)+countnz(diag(mat))) / 2)
        II = Array(Int32, nnz)
        JJ = Array(Int32, nnz)
        VV = Array(Float64, nnz)
        m, n = size(mat)
        cnt = 1
        for j in 1:m # get LOWER TRIANGULAR
            for i in j:n
                if mat[i,j] != 0.0
                    II[cnt] = i
                    JJ[cnt] = j
                    VV[cnt] = mat[i,j]
                    cnt += 1
                end
            end
        end
    elseif istriu(mat)
        nnz = nfilled(mat)
        II = Array(Int32, nnz)
        JJ = Array(Int32, nnz)
        VV = Array(Float64, nnz)
        m, n = size(mat)
        cnt = 1
        for i in 1:m # UPPER TRIANGULAR -> LOWER TRIANGULAR
            for j in i:n
                if mat[i,j] != 0.0
                    II[cnt] = j
                    JJ[cnt] = i
                    VV[cnt] = mat[i,j]
                    cnt += 1
                end
            end
        end
    # elseif istril(mat)
    #     nnz = nfilled(mat)
    #     II = Array(Int32, nnz)
    #     JJ = Array(Int32, nnz)
    #     VV = Array(Float64, nnz)
    #     m, n = size(mat)
    #     cnt = 1
    #     for j in 1:m # get LOWER TRIANGULAR
    #         for i in j:n
    #             if mat[i,j] != 0.0
    #                 II[cnt] = i
    #                 JJ[cnt] = j
    #                 VV[cnt] = mat[i,j]
    #                 cnt += 1
    #             end
    #         end
    #     end
    else
        error("Matrix must be symmetric or upper triangular")
    end
    return II,JJ,VV
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

function addsdpvar!(m::MosekMathProgModel, dim)
    appendbarvars(m.task, Cint[dim])
    return convert(Int64, getnumbarvar(m.task))
end

function addsdpmatrix!(m::MosekMathProgModel, mat)
    if isa(mat, Matrix)
        II,JJ,VV = denseToSparseTriple(mat)
    elseif isa(mat, SparseMatrixCSC)
        II,JJ,VV = sparseToSparseTriple(mat)
    else
        II,JJ,VV = sparseToSparseTriple(sparse(mat))
    end
    idx = Mosek.appendsparsesymmat(m.task, size(mat,1), II, JJ, VV)
    return convert(Int64, idx)
end

function addsdpconstr!(m::MosekMathProgModel, matvaridx, matcoefidx, scalidx, scalcoef, lb, ub)
    m.probtype = MosekMathProgModel_SDP
    appendcons(m.task,1)
    constridx = getnumcon(m.task)
    for i in 1:length(matvaridx)
        putbaraij(m.task, constridx, matvaridx[i], [matcoefidx[i]], [1])
    end
    putarow(m.task,constridx,scalidx,scalcoef)
    bk = getBoundsKey(lb, ub)
    putconbound(m.task,constridx,bk,lb,ub)

    userconidx = addUserCon(m,constridx)

    return convert(Int64, userconidx)
end

function setsdpobj!(m::MosekMathProgModel, matvaridx, matcoefidx)
    for (it,varidx) in enumerate(matvaridx)
        putbarcj(m.task, varidx, [matcoefidx[it]], [1])
    end
end

function getsdpsolution(m::MosekMathProgModel, idx)
    V = getbarxj(m.task, MSK_SOL_ITR, idx)
    n = convert(Int64, sqrt(8*length(V)+1)/2-1/2 )
    cnt = 0
    A = Array(Float64,n,n)
    for j in 1:n
      cnt += 1
      A[j,j] = V[cnt]
      for i in (j+1):n
        cnt += 1
        A[i,j] = V[cnt]
        A[j,i] = V[cnt]
      end
    end
    return A
end

getsdpdual(m::MosekMathProgModel) = getconstrduals(m)

include("MosekNLPSolverInterface.jl")

end
