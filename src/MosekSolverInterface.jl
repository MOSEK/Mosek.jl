export 
  MosekSolver,
  loadproblem,
  writeproblem,
  getvarLB,
  setvarLB,
  getvarLB,
  setvarLB,
  getconstrLB,
  setconstrLB,
  getconstrUB,
  setconstrUB,
  getobj,
  setobj,
  addvar,
  addconstr,
  updatemodel,
  setsense,
  getsense,
  numvar,
  numconstr,
  optimize,
  status,
  getobjval,
  getobjbound,
  getsolution,
  getconstrsolution,
  getreducedcosts,
  getconstrduals,
  getrawsolver,
  setvartype

require(joinpath(Pkg.dir("MathProgBase"),"src","LinprogSolverInterface.jl"))
import LinprogSolverInterface.LinprogSolver

type MosekSolver <: LinprogSolver
  task :: MSKtask
end

type MosekSolverError
  msg :: String
end

# NOTE: This method will load data into an existing task, but
# it will not necessarily reset all exiting data, depending on the 
# file read (e.g. reading an MPS will not reset parameters)
function loadproblem(m:: MosekSolver, filename:: String)
  readdata(m.task, filename)
end

function writeproblem(m:: MosekSolver, filename:: String)
  writedata(m.task, filename)
end


function loadproblem( m::     MosekSolver, 
                      A::     SparseMatrixCSC, 
                      collb:: Array{Float64,1},
                      colub:: Array{Float64,1},
                      obj::   SparseMatrixCSC, 
                      rowlb:: Array{Float64,1}, 
                      rowub:: Array{Float64,1})
  deletetask(m.task)
  m.task = maketask(msk_global_env)
  
  nrows,ncols = size(A) 
  if ncols != length(collb) ||
     ncols != length(colub) ||
     ncols != size(obj,1)   ||
     nrows != length(rowlb) ||
     nrows != length(rowub) ||
     size(obj,2) != 1
    throw(MosekSolverError("Inconsistent data dimensions"))
  end

  appendvars(task, ncols)
  appendcons(task, nrows)
  
  # input coefficients
  putclist(m.task, obj.rowval, obj.nzval)
  putacolslice(m.task, 1, nrows+1, A.colptr[1:ncols], A.colptr[2:ncols+1], A.rowval, A.nzval)

  # input bounds
  putvarboundslice(m.task, 1, ncols+1, [ MSK_BK_RA for i=1:ncols], collb, colub)
  putconboundslice(m.task, 1, nrows+1, [ MSK_BK_RA for i=1:nrows], rowlb, rowub)
  nothing
end

function loadproblem ( m::     MosekSolver,
                       A,
                       collb,
                       colub,
                       obj,
                       rowlb,
                       rowub)
  loadproblem(m,sparse(A),collb,colub,sparse(obj'),rowlb,rowub)
end

function getvarLB(m::MosekSolver)
  bkx,blx,bux = getvarboundslice(m.task,1,getnumvar(m))  
  [ (if bkx[i] in [ MSK_BK_FR, MSK_BK_UP ] -Inf else blx[i] end) for i=1:length(bkx) ] :: Array{Float64,1}
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
  if bl < Inf
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

function setvarLB(m::MosekSolver, collb)
  nvars = getnumvar(m.task)
  if nvars != length(collb)
    throw(MosekSolverError("Bound vector has wrong size"))
  end
  bk,bl,bu = getvarboundslice(m.task,1,nvars+1)

  bk = [ complbk(bk[i],collb[i]) for i=1:nvars ]
  putvarbound(m.task, bk, collb, bu)
end

function getvarUB(m::MosekSolver) 
  bkx,blx,bux = getvarboundslice(m.task,1,getnumvar(m))  
  [ (if bkx[i] in [ MSK_BK_FR, MSK_BK_LO ] Inf else bux[i] end) for i=1:length(bkx) ] :: Array{Float64,1} 
end

function setvarUB(m::MosekSolver, colub) 
  nvars = getnumvar(m.task)
  if nvars != length(colub)
    throw(MosekSolverError("Bound vector has wrong size"))
  end
  bk,bl,bu = getvarboundslice(m.task,1,nvars+1)

  bk = [ compubk(bk[i],colub[i]) for i=1:nvars ]
  putvarbound(m.task, bk, bl, colub)
end

function getconstrLB(m::MosekSolver) 
  bkc,blc,buc = getconboundslice(m.task,1,getnumvar(m))  
  [ (if bkc[i] in [ MSK_BK_FR, MSK_BK_UP ] -Inf else blc[i] end) for i=1:length(bkc) ] :: Array{Float64,1}
end

function setconstrLB(m::MosekSolver, rowlb) 
  nvars = getnumcon(m.task)
  if ncons != length(rowlb)
    throw(MosekSolverError("Bound vector has wrong size"))
  end
  bk,bl,bu = getconboundslice(m.task,1,ncons+1)

  bk = [ complbk(bk[i],rowlb[i]) for i=1:ncons ]
  putconbound(m.task, bk, rowlb, bu)
end

function getconstrUB(m::MosekSolver) 
  bkc,blc,buc = getconboundslice(m.task,1,getnumvar(m))  
  [ (if bkc[i] in [ MSK_BK_FR, MSK_BK_LO ] Inf else buc[i] end) for i=1:length(bkc) ] :: Array{Float64,1} 
end

function setconstrUB(m::MosekSolver, rowub) 
  nvars = getnumcon(m.task)
  if ncons != length(rowub)
    throw(MosekSolverError("Bound vector has wrong size"))
  end
  bk,bl,bu = getconboundslice(m.task,1,ncons+1)

  bk = [ compubk(bk[i],rowub[i]) for i=1:ncons ]
  putconbound(m.task, bk, rowub, bu)
end

function getobj(m::MosekSolver) 
  getc(m.task)
end

function setobj(m::MosekSolver, obj::Array{Float64,1})     
  nvars = getnumvar(m.task)
  if size(obj,1) != nvars
    throw(MosekSolverError("Objective vector has wrong size"))
  end

  putclist(m.task,[1:numvar+1],obj)
end

function setobj(m::MosekSolver, obj)
  setobj(m,dense(float(obj)))
end

# addvar(m::MosekSolver, rowidx, rowcoef, collb, colub, objcoef)
# addconstr(m::MosekSolver, colidx, colcoef, rowlb, rowub) 

updatemodel(m::MosekSolver) = nothing

function setsense(m::MosekSolver,sense) 
  if     sense == :Min
    putobjsense(m.task, MSK_OBJECTIVE_SENSE_MINIMIZE)
  elseif sense == :Max
    putobjsense(m.task, MSK_OBJECTIVE_SENSE_MAXIMIZE)
  else
   throw(MosekSolverError("Invalid objective sense"))
  end
  nothing  
end

function getsense(m::MosekSolver) 
  sense = getobjsense(m.task)
  if sense == MSK_OBJECTIVE_SENSE_MINIMIZE
    :Min
  elseif sense == MSK_OBJECTIVE_SENSE_MAXIMIZE
    :Max
  else
    None
  end
end

numvar(m::MosekSolver) = getnumvar(m.task)
numconstr(m::MosekSolver) = getnumcon(m.task)
optimize(m::MosekSolver) = optimize(m.task)




function getsoldef(m::MosekSolver)
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
function status(m::MosekSolver) 
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

function getobjval(m::MosekSolver) 
  soldef = getsoldef(m)
  if soldef < 0 return NaN end

  getprimalobj(m,soldef)
end

# getobjbound(m::MosekSolver)

function getsolution(m::MosekSolver) 
  soldef = getsoldef(m)
  if soldef < 0 throw(MosekSolverError("No solution available")) end
  getxx(m,soldef)
end

function getconstrsolution(m::MosekSolver) 
  soldef = getsoldef(m)
  if soldef < 0 throw(MosekSolverError("No solution available")) end
  getxc(m,soldef)
end

# NOTE: is the dual in your model slx-sux or sux-slx?
function getreducedcosts(m::MosekSolver) 
  soldef = getsoldef(m)
  if soldef < 0 throw(MosekSolverError("No solution available")) end
  getslx(m,soldef) - getsux(m,soldef)
end

function getconstrduals(m::MosekSolver)
  soldef = getsoldef(m)
  if soldef < 0 throw(MosekSolverError("No solution available")) end
  gety(m,soldef)
end

getrawsolver(m::MosekSolver) = m.task

## NOTE: What does this do?!
#function setvartype(m::MosekSolver, vartype) 

# getvartype(m::MosekSolver) 


