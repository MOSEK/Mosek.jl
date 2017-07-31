include("LinkedInts.jl")

mosek_block_type_unallocated = 0
mosek_block_type_zero   = 1
mosek_block_type_nonneg = 2
mosek_block_type_nonpos = 3
mosek_block_type_range  = 4
mosek_block_type_qcone  = 5
mosek_block_type_rqcone = 6
mosek_block_type_psd    = 7
mosek_block_type_integer = 8

problemtype_linear    = 0
problemtype_conic     = 1
problemtype_quadratic = 2

import MathOptInterface

immutable MosekSolver <: MathOptInterface.AbstractSolver
  options
end

MosekSolver(;kwargs...) = MosekSolver(kwargs)


boundflag_lower = 0x1
boundflag_upper = 0x2
boundflag_cone  = 0x4
boundflag_int   = 0x8
    
"""
    MosekModel <: MathOptInterface.AbstractModel

Linear variables and constraint can be deleted. For some reason MOSEK
does not support deleting PSD variables.

Note also that adding variables and constraints will permanently add
some (currently between 1 and 3) Int64s that a `delete!` will not
remove. This ensures that References (Variable and constraint) that
are deleted are thereafter invalid. 
"""
mutable struct MosekModel  <: MathOptInterface.AbstractSolverInstance
    task :: MSKtask

    problemtype :: Int

    """
    The total length of `x_block` matches the number of variables in
    the underlying task, and the number of blocks corresponds to the
    number variables allocated in the Model.
    """
    x_block      :: LinkedInts
    
    """
    One entry per scalar variable in the task indicating which bound types are defined
    """
    x_boundflags :: Vector{Int}
    
    """
    One entry per scalar variable in the task defining the number of
    variable constraints imposed on that variable. It is not allowed
    to delete a variable without deleting the constraints on it first
    (to get around the problem of deleting roots in conic constraints).
    """
    x_numxc :: Vector{Int}

    """
    One entry per variable-constraint
    """
    xc_block     :: LinkedInts
    
    """
    One entry per variable-constraint block indicating which bound
    types it defines. The values are binary ORed `boundflag_...` values.
    """
    xc_bounds    :: Vector{Int} # ORed boundflag values

    """
    One entry per scalar variable-constraint, defining which native
    variables the bound block covers.
    """
    xc_idxs      :: Vector{Int}

    ###########################
    """
    One scalar entry per constraint in the underlying task. One block
    per constraint allocated in the Model.
    """
    c_block :: LinkedInts

    """
    One entry per allocated scalar constraint. Defines the fixed term
    on the left-hand for each scalar constraint.
    """
    c_constant :: Vector{Float64}

    """
    One entry per allocated scalar constraint.
    Each element is either
    - 0: meaning: no slack, when domain is defined directly as a bound,
    - positive: a `x_block` reference, e.g. for qcones, or
    - negative: Negated index of a PSD variable in the underlying task
    """
    c_block_slack   :: Vector{Int}

    ###########################
    trm :: Int32
end
Mosek_VAR   = 1
Mosek_SLACK = 2


function MathOptInterface.SolverInstance(solver::MosekSolver)
    MosekModel(maketask(),# task
               problemtype_linear, # problemtype
               LinkedInts(),# c_block
               Int[], # x_boundflags
               Int[], # x_boundflags
               LinkedInts(), # xc_block
               Int[], # xc_bounds
               Int[], # xc_idxs
               LinkedInts(), # c_block
               Float64[], # c_constant
               Int[], # c_block_slack
               Mosek.MSK_RES_OK) # trm
end

function MathOptInterface.free!(m::MosekModel)
    deletetask(m.task)
end

function MathOptInterface.optimize!(m::MosekModel)
    m.trm = optimize(m.task)
end

function MathOptInterface.writeproblem(m::MosekModel, filename :: String)
    writedata(m.task,filename)
end

# For linear objectives we accept:
# EITER affine left-hand side and ranged, unbounded, half-open, fixed (equality), PSD or SOC domains
# OR affine and quadratic left-hand side, and ranged, unbounded, half-open, fixed (equality) domains (quadratic constraints must be unbounded or half-open)
#
# For non-quadratic problems we allow binary and integer variables (but not constraints)
function MathOptInterface.supportsproblem(m::MosekModel, objective_type::MathOptInterface.ScalarAffineFunction, constraint_types::Vector) :: Bool
    if any( (fun,dom) => (typeof(fun) == MathOptInterface.ScalarQuadraticFunction ||
                          typeof(fun) == MathOptInterface.VectorQuadraticFunction))
        return false
    end
    

    for (fun,dom) in constraint_types
        if  typeof(fun) in [MathOptInterface.ScalarAffineFunction,
                            MathOptInterface.ScalarVariablewiseFunction,
                            MathOptInterface.VectorAffineFunction] &&
                                typeof(dom) in [MathOptInterface.Zeros,
                                                MathOptInterface.Reals,
                                                MathOptInterface.Nonnegatives,
                                                MathOptInterface.Nonpositives,
                                                MathOptInterface.GreaterThan,
                                                MathOptInterface.LessThan,
                                                MathOptInterface.EqualTo,
                                                MathOptInterface.Interval,
                                                MathOptInterface.SecondOrderCone,
                                                MathOptInterface.RotatedSecondOrderCone,
                                                MathOptInterface.PositiveSemidefiniteConeTriangle ]
            # ok
        elseif typeof(dom) in [MathOptInterface.ZeroOne,
                               MathOptInterface.Integer] &&
                                   typeof(fun) in [MathOptInterface.ScalarVariablewiseFunction,
                                                   MathOptInterface.VectorVariablewiseFunction]
            # ok
        else
            return false
        end
    end

    return true
end

MathOptInterface.supportsproblem(m::MosekModel, objective_type::MathOptInterface.ScalarQuadraticFunction, constraint_types::Vector) = false


ref2id(ref :: MathOptInterface.VariableReference) :: Int =
    if ref.value & 1 == 0
        Int(ref.value >> 1)
    else
        - Int(ref.value >> 1)
    end

ref2id(ref :: MathOptInterface.ConstraintReference) :: Int =
    if ref.value & 1 == 0
        Int(ref.value >> 1)
    else
        - Int(ref.value >> 1)
    end

id2vref(id :: Int) :: MathOptInterface.VariableReference =
    if id < 0
        MathOptInterface.VariableReference((UInt64(-id) << 1) | 1)
    else
        MathOptInterface.VariableReference(UInt64(id) << 1)
    end

#id2cref{F,S}(id :: Int) :: MathOptInterface.ConstraintReference{F,S} =
#    if id < 0
#        MathOptInterface.ConstraintReference{F,S}((UInt64(-id) << 1) | 1)
#    else
#        MathOptInterface.ConstraintReference{F,S}(UInt64(id) << 1)
#    end


include("variable.jl")
include("constraint.jl")

#function MathOptInterface.setobjective!(m::MosekModel, N::Int, b, a_varidx, a_coef, Q_vari, Q_varj, Q_coef)
#end
#function MathOptInterface.modifyobjective!(m::MosekModel, i::Int, args...) end
#function MathOptInterface.modifyobjective!(m::MosekModel, i::Int, b) end
#function MathOptInterface.modifyobjective!(m::MosekModel, i::Int, a_varidx, a_coef) end
#function MathOptInterface.modifyobjective!(m::MosekModel, i::Int, Q_vari, Q_varj, Q_coef) end
#function MathOptInterface.getobjectiveaffine(m::MosekModel) end
#function MathOptInterface.getobjectiveconstant(m::MosekModel) end
#function MathOptInterface.modifyobjective!(m::MosekModel,i::Int, b) end





export MosekSolver, MosekModel
