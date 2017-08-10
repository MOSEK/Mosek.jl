include("LinkedInts.jl")

const DEBUG = true

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
boundflag_all   = 0x0f


# Mapping of all constraint types to its index
struct ConstraintMap
    x_lessthan           :: Dict{UInt64,Int} # (SingleVariable,LessThan) -> constraint number
    x_greaterthan        :: Dict{UInt64,Int} # (SingleVariable,GreaterThan) -> constraint number
    x_equalto            :: Dict{UInt64,Int} # (SingleVariable,EqualTo) -> constraint number
    x_interval           :: Dict{UInt64,Int} # (SingleVariable,Interval) -> constraint number
    x_nonpositives       :: Dict{UInt64,Int} # (SingleVariable,Nonpositives) -> constraint number
    x_nonnegatives       :: Dict{UInt64,Int} # (SingleVariable,Nonnegatives) -> constraint number
    x_binary             :: Dict{UInt64,Int} # (SingleVariable,ZeroOne) -> constraint number
    x_integer            :: Dict{UInt64,Int} # (SingleVariable,Integer) -> constraint number

    xs_nonpositives      :: Dict{UInt64,Int} # (VectorOfVariables,Nonpositives) -> constraint number
    xs_nonnegatives      :: Dict{UInt64,Int} # (VectorOfVariables,Nonnegatives) -> constraint number
    xs_zeros             :: Dict{UInt64,Int} # (VectorOfVariables,Zeros) -> constraint number
    xs_reals             :: Dict{UInt64,Int} # (VectorOfVariables,Reals) -> constraint number
    xs_qcone             :: Dict{UInt64,Int} # (VectorOfVariables,SecondOrderCone) -> constraint number
    xs_rqcone            :: Dict{UInt64,Int} # (VectorOfVariables,RotatedSecondOrderCone) -> constraint number
    xs_psdconetriangle   :: Dict{UInt64,Int} # (VectorOfVariables,PositiveSemidefiniteConeTriangle) -> constraint number

    axb_lessthan         :: Dict{UInt64,Int} # (ScalarAffineFunction,LessThan) -> constraint number
    axb_greaterthan      :: Dict{UInt64,Int} # (ScalarAffineFunction,GreaterThan) -> constraint number
    axb_equalto          :: Dict{UInt64,Int} # (ScalarAffineFunction,EqualTo) -> constraint number
    axb_interval         :: Dict{UInt64,Int} # (ScalarAffineFunction,Interval) -> constraint number
    axb_nonpositives     :: Dict{UInt64,Int} # (ScalarAffineFunction,Nonpositives) -> constraint number
    axb_nonnegatives     :: Dict{UInt64,Int} # (ScalarAffineFunction,Nonnegatives) -> constraint number
    axb_binary           :: Dict{UInt64,Int} # (ScalarAffineFunction,ZeroOne) -> constraint number
    axb_integer          :: Dict{UInt64,Int} # (ScalarAffineFunction,Integer) -> constraint number

    axbs_nonpositives    :: Dict{UInt64,Int} # (VectorAffineFunction,Nonpositives) -> constraint number
    axbs_nonnegatives    :: Dict{UInt64,Int} # (VectorAffineFunction,Nonnegatives) -> constraint number
    axbs_zeros           :: Dict{UInt64,Int} # (VectorAffineFunction,Zeros) -> constraint number
    axbs_reals           :: Dict{UInt64,Int} # (VectorAffineFunction,Reals) -> constraint number
    axbs_qcone           :: Dict{UInt64,Int} # (VectorAffineFunction,SecondOrderCone) -> constraint number
    axbs_rqcone          :: Dict{UInt64,Int} # (VectorAffineFunction,RotatedSecondOrderCone) -> constraint number
    axbs_psdconetriangle :: Dict{UInt64,Int} # (VectorAffineFunction,PositiveSemidefiniteConeTriangle) -> constraint number
end

ConstraintMap() = ConstraintMap([Dict{UInt64,Int}() for i in 1:30]...)
select(cm::ConstraintMap,::Type{MathOptInterface.SingleVariable},      ::Type{MathOptInterface.LessThan{Float64}}) =                               cm.x_lessthan
select(cm::ConstraintMap,::Type{MathOptInterface.SingleVariable},      ::Type{MathOptInterface.GreaterThan{Float64}}) =                            cm.x_greaterthan
select(cm::ConstraintMap,::Type{MathOptInterface.SingleVariable},      ::Type{MathOptInterface.EqualTo{Float64}}) =                                cm.x_equalto
select(cm::ConstraintMap,::Type{MathOptInterface.SingleVariable},      ::Type{MathOptInterface.Interval{Float64}}) =                               cm.x_interval
select(cm::ConstraintMap,::Type{MathOptInterface.SingleVariable},      ::Type{MathOptInterface.Nonpositives}) =                           cm.x_nonpositives
select(cm::ConstraintMap,::Type{MathOptInterface.SingleVariable},      ::Type{MathOptInterface.Nonnegatives}) =                           cm.x_nonnegatives
select(cm::ConstraintMap,::Type{MathOptInterface.SingleVariable},      ::Type{MathOptInterface.ZeroOne}) =                           cm.x_binary
select(cm::ConstraintMap,::Type{MathOptInterface.SingleVariable},      ::Type{MathOptInterface.Integer}) =                           cm.x_integer
select(cm::ConstraintMap,::Type{MathOptInterface.VectorOfVariables},   ::Type{MathOptInterface.Nonpositives}) =                        cm.xs_nonpositives
select(cm::ConstraintMap,::Type{MathOptInterface.VectorOfVariables},   ::Type{MathOptInterface.Nonnegatives}) =                        cm.xs_nonnegatives
select(cm::ConstraintMap,::Type{MathOptInterface.VectorOfVariables},   ::Type{MathOptInterface.Zeros}) =                        cm.xs_zeros
select(cm::ConstraintMap,::Type{MathOptInterface.VectorOfVariables},   ::Type{MathOptInterface.Reals}) =                        cm.xs_reals
select(cm::ConstraintMap,::Type{MathOptInterface.VectorOfVariables},   ::Type{MathOptInterface.SecondOrderCone}) =                     cm.xs_qcone
select(cm::ConstraintMap,::Type{MathOptInterface.VectorOfVariables},   ::Type{MathOptInterface.RotatedSecondOrderCone}) =              cm.xs_rqcone
select(cm::ConstraintMap,::Type{MathOptInterface.VectorOfVariables},   ::Type{MathOptInterface.PositiveSemidefiniteConeTriangle}) =    cm.xs_psdconetriangle
select(cm::ConstraintMap,::Type{MathOptInterface.ScalarAffineFunction{Float64}},::Type{MathOptInterface.LessThan{Float64}}) =                         cm.axb_lessthan
select(cm::ConstraintMap,::Type{MathOptInterface.ScalarAffineFunction{Float64}},::Type{MathOptInterface.GreaterThan{Float64}}) =                      cm.axb_greaterthan
select(cm::ConstraintMap,::Type{MathOptInterface.ScalarAffineFunction{Float64}},::Type{MathOptInterface.EqualTo{Float64}}) =                          cm.axb_equalto
select(cm::ConstraintMap,::Type{MathOptInterface.ScalarAffineFunction{Float64}},::Type{MathOptInterface.Interval{Float64}}) =                         cm.axb_interval
select(cm::ConstraintMap,::Type{MathOptInterface.ScalarAffineFunction{Float64}},::Type{MathOptInterface.ZeroOne}) =                         cm.axb_binary
select(cm::ConstraintMap,::Type{MathOptInterface.ScalarAffineFunction{Float64}},::Type{MathOptInterface.Integer}) =                         cm.axb_integer
select(cm::ConstraintMap,::Type{MathOptInterface.ScalarAffineFunction{Float64}},::Type{MathOptInterface.Nonpositives}) =                     cm.axb_nonpositives
select(cm::ConstraintMap,::Type{MathOptInterface.ScalarAffineFunction{Float64}},::Type{MathOptInterface.Nonnegatives}) =                     cm.axb_nonnegatives
select(cm::ConstraintMap,::Type{MathOptInterface.VectorAffineFunction{Float64}},::Type{MathOptInterface.Nonpositives}) =                     cm.axbs_nonpositives
select(cm::ConstraintMap,::Type{MathOptInterface.VectorAffineFunction{Float64}},::Type{MathOptInterface.Nonnegatives}) =                     cm.axbs_nonnegatives
select(cm::ConstraintMap,::Type{MathOptInterface.VectorAffineFunction{Float64}},::Type{MathOptInterface.Zeros}) =                     cm.axbs_zeros
select(cm::ConstraintMap,::Type{MathOptInterface.VectorAffineFunction{Float64}},::Type{MathOptInterface.Reals}) =                     cm.axbs_reals
select(cm::ConstraintMap,::Type{MathOptInterface.VectorAffineFunction{Float64}},::Type{MathOptInterface.SecondOrderCone}) =                  cm.axbs_qcone
select(cm::ConstraintMap,::Type{MathOptInterface.VectorAffineFunction{Float64}},::Type{MathOptInterface.RotatedSecondOrderCone}) =           cm.axbs_rqcone
select(cm::ConstraintMap,::Type{MathOptInterface.VectorAffineFunction{Float64}},::Type{MathOptInterface.PositiveSemidefiniteConeTriangle}) = cm.axbs_psdconetriangle

Base.getindex{F,D}(cm::ConstraintMap,r :: MathOptInterface.ConstraintReference{F,D}) = select(cm,F,D)[r.value]


struct MosekSolution
    whichsol :: Int32
    solsta   :: Int32
    prosta   :: Int32

    xxstatus :: Vector{Int32}
    xx       :: Vector{Float64}
    slx      :: Vector{Float64}
    sux      :: Vector{Float64}
    snx      :: Vector{Float64}

    cstatus  :: Vector{Int32}
    xc       :: Vector{Float64}
    slc      :: Vector{Float64}
    suc      :: Vector{Float64}
    y        :: Vector{Float64}
end

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


    """
    Number of variables explicitly created by user
    """
    publicnumvar :: Int

    """
    """
    constrmap :: ConstraintMap

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
    xc_bounds    :: Vector{Int}

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
    solutions :: Vector{MosekSolution}

end


function MathOptInterface.SolverInstance(solver::MosekSolver)
    MosekModel(maketask(),# task
               0, # public numvar
               ConstraintMap(), # public constraints
               LinkedInts(),# c_block
               Int[], # x_boundflags
               Int[], # x_boundflags
               LinkedInts(), # xc_block
               Int[], # xc_bounds
               Int[], # xc_idxs
               LinkedInts(), # c_block
               Float64[], # c_constant
               Int[], # c_block_slack
               Mosek.MSK_RES_OK,
               MosekSolution[]) # trm
end

function MathOptInterface.free!(m::MosekModel)
    deletetask(m.task)
end

function MathOptInterface.optimize!(m::MosekModel)
    m.trm = optimize(m.task)
    m.solutions = MosekSolution[]
    if solutiondef(m.task,MSK_SOL_ITG)

        push!(m.solutions,
              MosekSolution(MSK_SOL_ITG,
                            getsolsta(m.task,MSK_SOL_ITG),
                            getprosta(m.task,MSK_SOL_ITG),
                            getskx(m.task,MSK_SOL_ITG),
                            getxx(m.task,MSK_SOL_ITG),
                            Float64[],
                            Float64[],
                            Float64[],
                            getskc(m.task,MSK_SOL_ITG),
                            getxc(m.task,MSK_SOL_ITG),
                            Float64[],
                            Float64[],
                            Float64[]))
    end
    if solutiondef(m.task,MSK_SOL_BAS)
        push!(m.solutions,
              MosekSolution(MSK_SOL_BAS,
                            getsolsta(m.task,MSK_SOL_BAS),
                            getprosta(m.task,MSK_SOL_BAS),
                            getskx(m.task,MSK_SOL_BAS),
                            getxx(m.task,MSK_SOL_BAS),
                            getslx(m.task,MSK_SOL_BAS),
                            getsux(m.task,MSK_SOL_BAS),
                            Float64[],
                            getskc(m.task,MSK_SOL_BAS),
                            getxc(m.task,MSK_SOL_BAS),
                            getslc(m.task,MSK_SOL_BAS),
                            getsuc(m.task,MSK_SOL_BAS),
                            gety(m.task,MSK_SOL_BAS)))
    end
    if solutiondef(m.task,MSK_SOL_ITR)
        push!(m.solutions,
              MosekSolution(MSK_SOL_ITR,
                            getsolsta(m.task,MSK_SOL_ITR),
                            getprosta(m.task,MSK_SOL_ITR),
                            getskx(m.task,MSK_SOL_ITR),
                            getxx(m.task,MSK_SOL_ITR),
                            getslx(m.task,MSK_SOL_ITR),
                            getsux(m.task,MSK_SOL_ITR),
                            getsnx(m.task,MSK_SOL_ITR),
                            getskc(m.task,MSK_SOL_ITR),
                            getxc(m.task,MSK_SOL_ITR),
                            getslc(m.task,MSK_SOL_ITR),
                            getsuc(m.task,MSK_SOL_ITR),
                            gety(m.task,MSK_SOL_ITR)))
    end
end

function MathOptInterface.writeproblem(m::MosekModel, filename :: String)
    putintparam(m.task,MSK_IPAR_OPF_WRITE_SOLUTIONS, MSK_ON)
    writedata(m.task,filename)
end

# For linear objectives we accept:
# EITER affine left-hand side and ranged, unbounded, half-open, fixed (equality), PSD or SOC domains
# OR affine and quadratic left-hand side, and ranged, unbounded, half-open, fixed (equality) domains (quadratic constraints must be unbounded or half-open)
#
# For non-quadratic problems we allow binary and integer variables (but not constraints)
function supportsconstraints(m::MosekSolver, constraint_types) :: Bool
    for (fun,dom) in constraint_types
        if  fun in [MathOptInterface.ScalarAffineFunction{Float64},
                    MathOptInterface.SingleVariable,
                    MathOptInterface.VectorAffineFunction{Float64}] &&
            dom in [MathOptInterface.Zeros,
                    MathOptInterface.Reals,
                    MathOptInterface.Nonnegatives,
                    MathOptInterface.Nonpositives,
                    MathOptInterface.GreaterThan{Float64},
                    MathOptInterface.LessThan{Float64},
                    MathOptInterface.EqualTo{Float64},
                    MathOptInterface.Interval{Float64},
                    MathOptInterface.SecondOrderCone,
                    MathOptInterface.RotatedSecondOrderCone,
                    MathOptInterface.PositiveSemidefiniteConeTriangle ]
            # ok
        elseif dom in [MathOptInterface.ZeroOne,
                       MathOptInterface.Integer] &&
                           fun in [MathOptInterface.SingleVariable,
                                   MathOptInterface.VectorOfVariables]
            # ok
        else
            return false
        end
    end
    true
end


MathOptInterface.supportsproblem(m::MosekSolver, ::Type{MathOptInterface.SingleVariable},                constraint_types) :: Bool = supportsconstraints(m,constraint_types)
MathOptInterface.supportsproblem(m::MosekSolver, ::Type{MathOptInterface.ScalarAffineFunction{Float64}}, constraint_types) :: Bool = supportsconstraints(m,constraint_types)
MathOptInterface.supportsproblem{F}(m::MosekSolver, ::Type{F}, constraint_types) :: Bool = false

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

include("objective.jl")
include("variable.jl")
include("constraint.jl")
include("attributes.jl")

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
