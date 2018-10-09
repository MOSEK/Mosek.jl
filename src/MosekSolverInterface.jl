module MosekMathProgSolverInterface
import ..Mosek

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

import MathProgBase

status(t::Mosek.Task) = status(t,Mosek.MSK_RES_OK)
function status(t::Mosek.Task, r::Mosek.Rescode)
    if  r == Mosek.MSK_RES_TRM_MAX_ITERATIONS ||
        r == Mosek.MSK_RES_TRM_MAX_NUM_SETBACKS ||
        r == Mosek.MSK_RES_TRM_MAX_TIME ||
        r == Mosek.MSK_RES_TRM_MIO_NEAR_ABS_GAP ||
        r == Mosek.MSK_RES_TRM_MIO_NEAR_REL_GAP ||
        r == Mosek.MSK_RES_TRM_MIO_NUM_BRANCHES ||
        r == Mosek.MSK_RES_TRM_MIO_NUM_RELAXS ||
        r == Mosek.MSK_RES_TRM_NUM_MAX_NUM_INT_SOLUTIONS ||
        r == Mosek.MSK_RES_TRM_OBJECTIVE_RANGE
        :UserLimit
    #elseif r == Mosek.MSK_RES_TRM_STALL
    #    # Yes, I know this is not a standard code, but it doesn't fit any other code
    #    :Stall
    elseif r == Mosek.MSK_RES_TRM_USER_CALLBACK
        :UserBreak
    elseif r.value >= 1000 && r.value < 10000 # awful hack! Means: MSK_RES_ERR_*. r < 1000 is MSK_RES_WRN_*, r >= 10000 means MSK_RES_TRM_*
        :Error
    elseif ! ( (Mosek.solutiondef(t,Mosek.MSK_SOL_ITG) && Mosek.getsolsta(t,Mosek.MSK_SOL_ITG) != Mosek.MSK_SOL_STA_UNKNOWN) || 
               (Mosek.solutiondef(t,Mosek.MSK_SOL_BAS) && Mosek.getsolsta(t,Mosek.MSK_SOL_BAS) != Mosek.MSK_SOL_STA_UNKNOWN) ||
               (Mosek.solutiondef(t,Mosek.MSK_SOL_ITR) && Mosek.getsolsta(t,Mosek.MSK_SOL_ITR) != Mosek.MSK_SOL_STA_UNKNOWN) )
        :Unknown
    else
        sol =
            if     Mosek.solutiondef(t,Mosek.MSK_SOL_ITG) && Mosek.getsolsta(t,Mosek.MSK_SOL_ITG) != Mosek.MSK_SOL_STA_UNKNOWN
                Mosek.MSK_SOL_ITG
            elseif Mosek.solutiondef(t,Mosek.MSK_SOL_BAS) && Mosek.getsolsta(t,Mosek.MSK_SOL_BAS) != Mosek.MSK_SOL_STA_UNKNOWN
                Mosek.MSK_SOL_BAS
            else
                Mosek.MSK_SOL_ITR
            end

        prosta = Mosek.getprosta(t,sol)
        solsta = Mosek.getsolsta(t,sol)

        if  solsta == Mosek.MSK_SOL_STA_PRIM_AND_DUAL_FEAS
            #:Suboptimal # Not a standard code
            :Unknown
        elseif solsta == Mosek.MSK_SOL_STA_PRIM_FEAS
            if Mosek.getnumintvar(t) > 0
                #:Suboptimal # Not a standard code
                :Unknown
            else
                :Unknown
            end
        elseif solsta == Mosek.MSK_SOL_STA_DUAL_ILLPOSED_CER ||
            solsta == Mosek.MSK_SOL_STA_PRIM_ILLPOSED_CER
            :DualityFailure
        elseif solsta == Mosek.MSK_SOL_STA_DUAL_FEAS
            :Unknown
        elseif solsta == Mosek.MSK_SOL_STA_DUAL_INFEAS_CER
            :Unbounded
        elseif solsta == Mosek.MSK_SOL_STA_PRIM_INFEAS_CER
            :Infeasible
        elseif solsta == Mosek.MSK_SOL_STA_OPTIMAL ||
            solsta == Mosek.MSK_SOL_STA_INTEGER_OPTIMAL
            :Optimal
        else
            :Unknown
        end
    end
end

function getsense(t::Mosek.Task)
    sense = Mosek.getobjsense(t)
    if sense == Mosek.MSK_OBJECTIVE_SENSE_MAXIMIZE
        :Max
    elseif sense == Mosek.MSK_OBJECTIVE_SENSE_MINIMIZE
        :Min
    else
        :None
    end
end

function setsense!(t::Mosek.Task,sense)
    if     sense == :Max
        Mosek.putobjsense(t,Mosek.MSK_OBJECTIVE_SENSE_MAXIMIZE)
    elseif sense == :Min
        Mosek.putobjsense(t,Mosek.MSK_OBJECTIVE_SENSE_MINIMIZE)
    end
end


function getobjgap(t::Mosek.Task)
    sol = getsoldef(t)
    if sol == MSK_SOL_ITG
        Mosek.getdouinf(m.task,MSK_DINF_MIO_OBJ_REL_GAP)
    else
        0.0
    end
end

getobjval(t::Mosek.Task) = Mosek.getprimalobj(t,getsoldef(t))

function makebounds(bl_ :: Vector{Float64},
                    bu_ :: Vector{Float64}) :: Tuple{Vector{Mosek.Boundkey},Vector{Float64},Vector{Float64}}
    bk = Vector{Mosek.Boundkey}(undef, length(bl_))
    bl = Vector{Float64}(undef, length(bl_))
    bu = Vector{Float64}(undef, length(bl_))

    for i in 1:length(bl_)
        if bl_[i] > -Inf
            if bu_[i] < Inf
                if abs(bu_[i]-bl_[i]) < 1e-8
                    bk[i] = Mosek.MSK_BK_FX
                    bl[i] = bl_[i]
                    bu[i] = bl_[i]
                else
                    bk[i] = Mosek.MSK_BK_RA
                    bl[i] = bl_[i]
                    bu[i] = bu_[i]
                end
            else # bu_[i] == Inf
                bk[i] = Mosek.MSK_BK_LO
                bl[i] = bl_[i]
                bu[i] = Inf
            end
        else # bl_[i] == -Inf
            if bu_[i] < Inf
                bk[i] = Mosek.MSK_BK_UP
                bl[i] = -Inf
                bu[i] = bu_[i]
            else
                bk[i] = Mosek.MSK_BK_FR
                bl[i] = -Inf
                bu[i] = Inf
            end
        end
    end

    bk,bl,bu
end











#mutable struct MosekMathProgModel <: MathProgBase.SolverInterface.AbstractMathProgModel
mutable struct MosekMathProgModel <: MathProgBase.AbstractMathProgModel
  task :: Mosek.Task
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
  varmap   :: Vector{Int32}
  barvarij :: Vector{Int64}

  # numbarvar
  #   Number of used elements in barvarmap.
  numbarvar :: Int32
  # barvarmap
  #   Maps Model PSD variable indexes into MOSEK barvar indexes. When
  #   using the Semidefinite interface, these are the semidefinite
  #   variables thar can be accessed.  Semidefinite variables added
  #   thorugh loadconicproblem!() are not mapped here.
  barvarmap :: Vector{Int32}

  # binvarflag
  #   Defines per variable if it is binary
  binvarflag :: Vector{Bool}

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
  conmap   :: Vector{Int32}
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
  conslack :: Vector{Int32}
  barconij :: Vector{Int64}

  # quadratic constraints
  numqcon :: Int32
  qconmap :: Vector{Int32}

  # options Options from MosekSolver
  options
end


mutable struct MosekMathProgModelError <: Exception
  msg :: AbstractString
end

function loadoptions_internal!(t::Mosek.Task, options)
    # write to console by default
    printstream(msg::AbstractString) = print(msg)

    be_quiet = false
    for (option,val) in options
        parname = string(option)
        if parname == "QUIET"
            be_quiet = convert(Bool,val)
        elseif startswith(parname, "MSK_IPAR_")
            Mosek.putnaintparam(t, parname, convert(Integer, val))
        elseif startswith(parname, "MSK_DPAR_")
            Mosek.putnadouparam(t, parname, convert(AbstractFloat, val))
        elseif startswith(parname, "MSK_SPAR_")
            Mosek.putnastrparam(t, parname, convert(AbstractString, val))
        elseif isa(val, Integer)
            parname = "MSK_IPAR_$parname"
            Mosek.putnaintparam(t, parname, val)
        elseif isa(val, AbstractFloat)
            parname = "MSK_DPAR_$parname"
            Mosek.putnadouparam(t, parname, val)
        elseif isa(val, AbstractString)
            parname = "MSK_SPAR_$parname"
            Mosek.putnastrparam(t, parname, val)
        else
            error("Value $val for parameter $option has unrecognized type")
        end
    end
    if ! be_quiet
        Mosek.putstreamfunc(t,Mosek.MSK_STREAM_LOG,printstream)
    end
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

function getsoldef(t::Mosek.Task)
    for sol in (Mosek.MSK_SOL_ITG, Mosek.MSK_SOL_BAS, Mosek.MSK_SOL_ITR)
        if Mosek.solutiondef(t,sol) && Mosek.getsolsta(t,sol) != Mosek.MSK_SOL_STA_UNKNOWN
            return sol
        end
    end
    throw(MosekMathProgModelError("No solution available"))
end

#include("MosekNLPSolverInterface.jl")
include("MosekLPQCQPInterface.jl")
include("MosekConicInterface.jl")

end

