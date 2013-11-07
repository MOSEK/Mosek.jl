module MosekMathProgSolverInterface
using ..Mosek

export MosekSolver


require(joinpath(Pkg.dir("MathProgBase"),"src","MathProgSolverInterface.jl"))
importall MathProgSolverInterface

type MosekMathProgModel <: AbstractMathProgModel 
  task :: Mosek.MSKtask
end

immutable MosekSolver <: AbstractMathProgSolver
  options
end
MosekSolver(;kwargs...) = MosekSolver(kwargs)

type MosekMathProgModelError
  msg :: String
end

function model(s::MosekSolver)
  # TODO: process solver options
  task = maketask(Mosek.msk_global_env)
  return MosekMathProgModel(task)
end

# NOTE: This method will load data into an existing task, but
# it will not necessarily reset all exiting data, depending on the 
# file read (e.g. reading an MPS will not reset parameters)
function loadproblem!(m:: MosekMathProgModel, filename:: String)
  readdata(m.task, filename)
end

function writeproblem(m:: MosekMathProgModel, filename:: String)
  writedata(m.task, filename)
end


function loadproblem!( m::     MosekMathProgModel, 
                      A::     SparseMatrixCSC, 
                      collb:: Array{Float64,1},
                      colub:: Array{Float64,1},
                      obj::   SparseMatrixCSC, 
                      rowlb:: Array{Float64,1}, 
                      rowub:: Array{Float64,1},
                      sense:: Symbol)
  Mosek.deletetask(m.task)
  m.task = maketask(Mosek.msk_global_env)
  
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
  
  # input coefficients
  putclist(m.task, obj.rowval, obj.nzval)
  putacolslice(m.task, 1, nrows+1, A.colptr[1:ncols], A.colptr[2:ncols+1], A.rowval, A.nzval)
  setsense!(m, sense)

  # input bounds
  putvarboundslice(m.task, 1, ncols+1, Int32[ MSK_BK_RA for i=1:ncols], collb, colub)
  putconboundslice(m.task, 1, nrows+1, Int32[ MSK_BK_RA for i=1:nrows], rowlb, rowub)
  nothing
end

function loadproblem! ( m::     MosekMathProgModel,
                       A,
                       collb,
                       colub,
                       obj,
                       rowlb,
                       rowub,
                       sense)
  loadproblem!(m,sparse(float(A)),float(collb),float(colub),sparse(float(obj)),float(rowlb),float(rowub),sense)
end

function getvarLB(m::MosekMathProgModel)
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

function setvarLB!(m::MosekMathProgModel, collb)
  nvars = getnumvar(m.task)
  if nvars != length(collb)
    throw(MosekMathProgModelError("Bound vector has wrong size"))
  end
  bk,bl,bu = getvarboundslice(m.task,1,nvars+1)

  bk = [ complbk(bk[i],collb[i]) for i=1:nvars ]
  putvarbound(m.task, bk, collb, bu)
end

function getvarUB(m::MosekMathProgModel) 
  bkx,blx,bux = getvarboundslice(m.task,1,getnumvar(m))  
  [ (if bkx[i] in [ MSK_BK_FR, MSK_BK_LO ] Inf else bux[i] end) for i=1:length(bkx) ] :: Array{Float64,1} 
end

function setvarUB!(m::MosekMathProgModel, colub) 
  nvars = getnumvar(m.task)
  if nvars != length(colub)
    throw(MosekMathProgModelError("Bound vector has wrong size"))
  end
  bk,bl,bu = getvarboundslice(m.task,1,nvars+1)

  bk = [ compubk(bk[i],colub[i]) for i=1:nvars ]
  putvarbound(m.task, bk, bl, colub)
end

function getconstrLB(m::MosekMathProgModel) 
  bkc,blc,buc = getconboundslice(m.task,1,getnumvar(m))  
  [ (if bkc[i] in [ MSK_BK_FR, MSK_BK_UP ] -Inf else blc[i] end) for i=1:length(bkc) ] :: Array{Float64,1}
end

function setconstrLB!(m::MosekMathProgModel, rowlb) 
  nvars = getnumcon(m.task)
  if ncons != length(rowlb)
    throw(MosekMathProgModelError("Bound vector has wrong size"))
  end
  bk,bl,bu = getconboundslice(m.task,1,ncons+1)

  bk = [ complbk(bk[i],rowlb[i]) for i=1:ncons ]
  putconbound(m.task, bk, rowlb, bu)
end

function getconstrUB(m::MosekMathProgModel) 
  bkc,blc,buc = getconboundslice(m.task,1,getnumvar(m))  
  [ (if bkc[i] in [ MSK_BK_FR, MSK_BK_LO ] Inf else buc[i] end) for i=1:length(bkc) ] :: Array{Float64,1} 
end

function setconstrUB!(m::MosekMathProgModel, rowub) 
  nvars = getnumcon(m.task)
  if ncons != length(rowub)
    throw(MosekMathProgModelError("Bound vector has wrong size"))
  end
  bk,bl,bu = getconboundslice(m.task,1,ncons+1)

  bk = [ compubk(bk[i],rowub[i]) for i=1:ncons ]
  putconbound(m.task, bk, rowub, bu)
end

function getobj(m::MosekMathProgModel) 
  getc(m.task)
end

function setobj!(m::MosekMathProgModel, obj::Array{Float64,1})     
  nvars = getnumvar(m.task)
  if size(obj,1) != nvars
    throw(MosekMathProgModelError("Objective vector has wrong size"))
  end

  putclist(m.task,[1:numvar+1],obj)
end

function setobj!(m::MosekMathProgModel, obj)
  setobj(m,dense(float(obj)))
end

# addvar(m::MosekMathProgModel, rowidx, rowcoef, collb, colub, objcoef)
# addconstr(m::MosekMathProgModel, colidx, colcoef, rowlb, rowub) 

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

numvar(m::MosekMathProgModel) = getnumvar(m.task)
numconstr(m::MosekMathProgModel) = getnumcon(m.task)
optimize!(m::MosekMathProgModel) = optimize(m.task)




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

# getobjbound(m::MosekMathProgModel)

function getsolution(m::MosekMathProgModel) 
  soldef = getsoldef(m)
  if soldef < 0 throw(MosekMathProgModelError("No solution available")) end
  getxx(m.task,soldef)
end

function getconstrsolution(m::MosekMathProgModel) 
  soldef = getsoldef(m)
  if soldef < 0 throw(MosekMathProgModelError("No solution available")) end
  getxc(m.task,soldef)
end

# NOTE: is the dual in your model slx-sux or sux-slx?
function getreducedcosts(m::MosekMathProgModel) 
  soldef = getsoldef(m)
  if soldef < 0 throw(MosekMathProgModelError("No solution available")) end
  getslx(m.task,soldef) - getsux(m.task,soldef)
end

function getconstrduals(m::MosekMathProgModel)
  soldef = getsoldef(m)
  if soldef < 0 throw(MosekMathProgModelError("No solution available")) end
  gety(m.task,soldef)
end

getrawsolver(m::MosekMathProgModel) = m.task

## NOTE: What does this do?!
#function setvartype(m::MosekMathProgModel, vartype) 

# getvartype(m::MosekMathProgModel) 

end


