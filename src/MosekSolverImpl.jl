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

    x_block      :: LinkedInts
    xc_block     :: LinkedInts
    x_boundflags :: Vector{Int}
    
    """
    Counts for each variable the number of integer constraints that
    are imposed on that variable. Zero means continuous.
    """
    x_block_integer :: Vector{Int}

    ###########################
    c_block :: LinkedInts

    c_constant :: Vector{Float64}

    """
    Native index of the cone to which the constraint block belongs.
    """
    c_block_coneidx :: Vector{Int}
    """
    Each element is either
    - 0, meaning: no slack, when domain is defined directly as a bound,
    - a `x_block` reference (positive number), e.g. for qcones, or
    - a PSD variable index (negative number)
    """
    c_block_slack   :: Vector{Int}

    ###########################
    trm :: Int32
end
Mosek_VAR   = 1
Mosek_SLACK = 2


function MathOptInterface.SolverInstance(solver::MosekSolver)
    MosekModel(maketask(),
               problemtype_linear,
               LinkedInts(),
               LinkedInts(),
               Int[],
               Int[],
               LinkedInts(),
               Int[],
               Int[],
               Int[],
               Mosek.MSK_RES_OK)
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
    isquad == any( (fun,dom) => (typeof(fun) == MathOptInterface.ScalarQuadraticFunction ||
                                 typeof(fun) == MathOptInterface.VectorQuadraticFunction))

    if isquad
        for (fun,dom) in constraint_types
            if  typeof(fun) in [MathOptInterface.ScalarQuadraticFunction,
                                MathOptInterface.VectorQuadraticFunction ] &&
                typeof(dom) in [MathOptInterface.Reals,
                                MathOptInterface.Nonnegatives,
                                MathOptInterface.Nonpositives,
                                MathOptInterface.GreaterThan,
                                MathOptInterface.LessThan]
                # ok
            elseif typeof(fun) in [MathOptInterface.ScalarAffineFunction,
                                   MathOptInterface.ScalarVariablewiseFunction,
                                   MathOptInterface.VectorAffineFunction] &&
                typeof(dom) in [MathOptInterface.Zeros,
                                MathOptInterface.Reals,
                                MathOptInterface.Nonnegatives,
                                MathOptInterface.Nonpositives,
                                MathOptInterface.GreaterThan,
                                MathOptInterface.LessThan,
                                MathOptInterface.EqualTo]
                # ok
            else
                return false
            end           
        end
    else # ! isquad
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
    end
    return true
end

# For affine+quadratic objective, we accept linear and convec quadratic constraints
function MathOptInterface.supportsproblem(m::MosekModel, objective_type::MathOptInterface.ScalarQuadraticFunction, constraint_types::Vector) :: Bool
    for (fun,dom) in constraint_types
        if typeof(dom) in [MathOptInterface.Zeros,
                           MathOptInterface.Reals,
                           MathOptInterface.Nonnegatives,
                           MathOptInterface.Nonpositives,
                           MathOptInterface.GreaterThan,
                           MathOptInterface.LessThan,
                           MathOptInterface.EqualTo] &&
            typeof(fun) in [MathOptInterface.ScalarAffineFunction,
                            MathOptInterface.ScalarVariablewiseFunction,
                            MathOptInterface.VectorAffineFunction,
                            MathOptInterface.ScalarQuadraticFunction,
                            MathOptInterface.VectorQuadraticFunction]
            # ok
        elseif typeof(dom) in [MathOptInterface.Reals,
                               MathOptInterface.Nonnegatives,
                               MathOptInterface.Nonpositives,
                               MathOptInterface.GreaterThan,
                               MathOptInterface.LessThan] &&
            typeof(fun) in [MathOptInterface.ScalarQuadraticFunction,
                            MathOptInterface.VectorQuadraticFunction]
            # ok
        else
            return false
        end
    end
    
    return true
end




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
