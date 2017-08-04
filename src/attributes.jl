

#### solver attributes
MathOptInterface.getattribute(m::Union{MosekSolver,MosekModel},::MathOptInterface.SupportsDuals) = true
MathOptInterface.getattribute(m::Union{MosekSolver,MosekModel},::MathOptInterface.SupportsAddConstraintAfterSolve) = true
MathOptInterface.getattribute(m::Union{MosekSolver,MosekModel},::MathOptInterface.SupportsAddVariableAfterSolve) = true
MathOptInterface.getattribute(m::Union{MosekSolver,MosekModel},::MathOptInterface.SupportsDeleteConstraint) = true
MathOptInterface.getattribute(m::Union{MosekSolver,MosekModel},::MathOptInterface.SupportsDeleteVariable) = true
MathOptInterface.getattribute(m::Union{MosekSolver,MosekModel},::MathOptInterface.SupportsConicThroughQuadratic) = false # though actually the solver does

#### objective
function MathOptInterface.getattribute(m::MosekModel,attr::MathOptInterface.ObjectiveValue)
    solitems = Int32[]
    if solutiondef(m.task,MSK_SOL_ITG) append!(solitems,getdouinf(m.task,MSK_SOL_ITG)) end
    if solutiondef(m.task,MSK_SOL_BAS) append!(solitems,getdouinf(m.task,MSK_SOL_BAS)) end
    if solutiondef(m.task,MSK_SOL_ITR) append!(solitems,getdouinf(m.task,MSK_SOL_ITR)) end

    getprimalobj(m.task,solitems[attr.resultindex])
end

MathOptInterface.cangetattribute(m::MosekSolver,attr::MathOptInterface.ObjectiveValue) = true
function MathOptInterface.cangetattribute(m::MosekModel,attr::MathOptInterface.ObjectiveValue)
    num = 0
    if solutiondef(m.task,MSK_SOL_ITG) num += 1 end
    if solutiondef(m.task,MSK_SOL_BAS) num += 1 end
    if solutiondef(m.task,MSK_SOL_ITR) num += 1 end


    attr.resultindex > 0 && attr.resultindex <= num
end


MathOptInterface.cangetattribute(m::MosekSolver,attr::MathOptInterface.ObjectiveBound) = true
MathOptInterface.cangetattribute(m::MosekModel,attr::MathOptInterface.ObjectiveBound) = true
MathOptInterface.getattribute(m::MosekModel,attr::MathOptInterface.ObjectiveBound) = getdouinf(m.task,MSK_DINF_MIO_OBJ_BOUND)

MathOptInterface.cangetattribute(m::MosekSolver,attr::MathOptInterface.RelativeGap) = true
MathOptInterface.cangetattribute(m::MosekModel,attr::MathOptInterface.RelativeGap) = true
MathOptInterface.getattribute(m::MosekModel,attr::MathOptInterface.RelativeGap) = getdouinf(m.task,MSK_DINF_MIO_REL_GAP)

MathOptInterface.cangetattribute(m::MosekSolver,attr::MathOptInterface.SolveTime) = true
MathOptInterface.cangetattribute(m::MosekModel,attr::MathOptInterface.SolveTime) = true
MathOptInterface.getattribute(m::MosekModel,attr::MathOptInterface.SolveTime) = getdouinf(m.task,MSK_DINF_OPTIMIZER_TIME)

MathOptInterface.cangetattribute(m::MosekSolver,attr::MathOptInterface.ObjectiveSense) = true
MathOptInterface.cangetattribute(m::MosekModel,attr::MathOptInterface.ObjectiveSense)  = true
MathOptInterface.cansetattribute(m::MosekSolver,attr::MathOptInterface.ObjectiveSense) = true
MathOptInterface.cansetattribute(m::MosekModel,attr::MathOptInterface.ObjectiveSense)  = true

# NOTE: The MOSEK interface currently only supports Min and Max
function MathOptInterface.getattribute(m::MosekModel,attr::MathOptInterface.ObjectiveSense)
    sense = getobjsense(m.task)
    if sense == MSK_OBJECTIVE_SENSE_MINIMIZE
        MathOptInterface.MinSense
    else
        MathOptInterface.MaxSense
    end
end

function MathOptInterface.setattribute!(m::MosekModel,attr::MathOptInterface.ObjectiveSense, sense::MathOptInterface.OptimizationSense)
    if sense == MathOptInterface.MinSense
        setobjsense(m.task,MSK_OBJECTIVE_SENSE_MINIMIZE)
    elseif sense == MathOptInterface.MaxSense
        setobjsense(m.task,MSK_OBJECTIVE_SENSE_MAXIMIZE)
    else
        error("Sense '$sense' is not supported")
    end
end

#### Solver/Solution information

MathOptInterface.cangetattribute(m::MosekSolver,attr::MathOptInterface.SimplexIterations) = true
MathOptInterface.cangetattribute(m::MosekModel,attr::MathOptInterface.SimplexIterations) = true
function MathOptInterface.getattribute(m::MosekModel,attr::MathOptInterface.SimplexIterations)
    miosimiter = getliintinf(m.task,MSK_LIINF_MIO_SIMPLEX_ITER)
    if miosimiter > 0
        Int(miosimiter)
    else
        Int(getiintinf(m.task,MSK_IINF_SIM_PRIMAL_ITER) + getliintinf(m.task,MSK_IINF_SIM_DUAL_ITER))
    end
end

MathOptInterface.cangetattribute(m::MosekSolver,attr::MathOptInterface.BarrierIterations) = true
MathOptInterface.cangetattribute(m::MosekModel,attr::MathOptInterface.BarrierIterations) = true
function MathOptInterface.getattribute(m::MosekModel,attr::MathOptInterface.BarrierIterations)
    miosimiter = getliintinf(m.task,MSK_LIINF_MIO_INTPNT_ITER)
    if miosimiter > 0
        Int(miosimiter)
    else
        Int(getiintinf(m.task,MSK_IINF_INTPNT_ITER))
    end
end

MathOptInterface.cangetattribute(m::MosekSolver,attr::MathOptInterface.NodeCount) = true
MathOptInterface.cangetattribute(m::MosekModel,attr::MathOptInterface.NodeCount) = true
function MathOptInterface.getattribute(m::MosekModel,attr::MathOptInterface.NodeCount)
        Int(getiintinf(m.task,MSK_IINF_NUM_BRANCH))
end

MathOptInterface.cangetattribute(m::MosekSolver,attr::MathOptInterface.RawSolver) = true
MathOptInterface.cangetattribute(m::MosekModel,attr::MathOptInterface.RawSolver) = true
MathOptInterface.getattribute(m::MosekModel,attr::MathOptInterface.RawSolver) = m.task

MathOptInterface.cangetattribute(m::MosekSolver,attr::MathOptInterface.ResultCount) = true
MathOptInterface.cangetattribute(m::MosekModel,attr::MathOptInterface.ResultCount) = true
function MathOptInterface.getattribute(m::MosekModel,attr::MathOptInterface.ResultCount)
    num = 0
    if solutiondef(m.task,MSK_SOL_ITG) num += 1 end
    if solutiondef(m.task,MSK_SOL_BAS) num += 1 end
    if solutiondef(m.task,MSK_SOL_ITR) num += 1 end

    return num
end
    
#### Problem information

MathOptInterface.cangetattribute(m::MosekSolver,attr::MathOptInterface.NumberOfVariables) = true
MathOptInterface.cangetattribute(m::MosekModel,attr::MathOptInterface.NumberOfVariables) = true
MathOptInterface.getattribute(m::MosekModel,attr::MathOptInterface.NumberOfVariables) = m.publicnumvar

MathOptInterface.cangetattribute(m::MosekSolver,attr::MathOptInterface.NumberOfConstraints) = true
MathOptInterface.cangetattribute{F,D}(m::MosekModel,attr::MathOptInterface.NumberOfConstraints{F,D}) = true
MathOptInterface.getattribute{F,D}(m::MosekModel,attr::MathOptInterface.NumberOfConstraints{F,D}) = length(select(m.constrmap,F,D))

#MathOptInterface.cangetattribute(m::MosekSolver,attr::MathOptInterface.ListOfVariableReferences) = false
#MathOptInterface.cangetattribute(m::MosekSolver,attr::MathOptInterface.ListOfConstraints) = false

MathOptInterface.cangetattribute{F,D}(m::MosekSolver,attr::MathOptInterface.ListOfConstraintReferences{F,D}) = true
MathOptInterface.getattribute{F,D}(m::MosekSolver,attr::MathOptInterface.ListOfConstraintReferences{F,D}) = keys(select(m.constrmap,F,D))


#### Warm start values

MathOptInterface.cansetattribute(m::MosekSolver,attr::MathOptInterface.VariablePrimalStart) = true
MathOptInterface.cansetattribute(m::MosekModel,attr::MathOptInterface.VariablePrimalStart) = true 
function MathOptInterface.setattribute!(m::MosekModel,attr::MathOptInterface.VariablePrimalStart, v :: MathOptInterface.VariableReference, val::Float64)
    subj = getindexes(m.x_block,ref2id(v))

    xx = Float64[val]
    for sol in [ MSK_SOL_BAS, MSK_SOL_ITG ]
        putxxslice(m.task,sol,Int(subj),Int(subj+1),xx)
    end
end

function MathOptInterface.setattribute!(m::MosekModel,attr::MathOptInterface.VariablePrimalStart, vs::Vector{MathOptInterface.VariableReference}, vals::Vector{Float64})
    subj = Array{Int}(length(vs))
    for i in 1:length(subj)
        getindexes(m.x_block,ref2id(vs[i]),subj,i)
    end
    
    for sol in [ MSK_SOL_BAS, MSK_SOL_ITG ]
        if solutiondef(m.task,sol)
            xx = getxx(m.task,sol)
            xx[subj] = vals
            putxx(m.task,sol,xx)
        else
            xx = zeros(Float64,getnumvar(m.task))
            xx[subj] = vals
            putxx(m.task,sol,xx)
        end
    end
end

MathOptInterface.cansetattribute(m::MosekSolver,attr::MathOptInterface.ConstraintPrimalStart) = false # not sure what exactly this would be...
MathOptInterface.cansetattribute(m::MosekModel,attr::MathOptInterface.ConstraintPrimalStart) = false

MathOptInterface.cansetattribute(m::MosekSolver,attr::MathOptInterface.ConstraintDualStart) = false # for now
MathOptInterface.cansetattribute(m::MosekModel,attr::MathOptInterface.ConstraintDualStart) = false

# function MathOptInterface.setattribute!(m::MosekModel,attr::MathOptInterface.ConstraintDualStart, vs::Vector{MathOptInterface.ConstraintReference}, vals::Vector{Float64})
#     subj = Array{Int}(length(vs))
#     for i in 1:length(subj)
#         getindexes(m.x_block,ref2id(vs[i]),subj,i)
#     end
    
#     for sol in [ MSK_SOL_BAS, MSK_SOL_ITG ]
#         if solutiondef(m.task,sol)
#             xx = getxx(m.task,sol)
#             xx[subj] = vals
#             putxx(m.task,sol,xx)
#         else
#             xx = zeros(Float64,getnumvar(m.task))
#             xx[subj] = vals
#             putxx(m.task,sol,xx)
#         end
#     end
# end



#### Variable solution values

function MathOptInterface.cangetattribute(m::MosekModel,attr::MathOptInterface.VariablePrimal)
    num = 0
    if solutiondef(m.task,MSK_SOL_ITG) num += 1 end
    if solutiondef(m.task,MSK_SOL_BAS) num += 1 end
    if solutiondef(m.task,MSK_SOL_ITR) num += 1 end

    println("can get ... $(attr.N > 0 && attr.N <= num)")
    attr.N > 0 && attr.N <= num
end

MathOptInterface.cangetattribute(m::MosekSolver,attr::MathOptInterface.VariablePrimal) = true 

MathOptInterface.cangetattribute(m::MosekModel,attr::MathOptInterface.VariablePrimal,vs::Vector{MathOptInterface.VariableReference}) = MathOptInterface.cangetattribute(m,attr)
MathOptInterface.cangetattribute(m::MosekModel,attr::MathOptInterface.VariablePrimal,v::MathOptInterface.VariableReference) = MathOptInterface.cangetattribute(m,attr)


function MathOptInterface.getattribute!(output::Vector{Float64},m::MosekModel,attr::MathOptInterface.VariablePrimal, vs::Vector{MathOptInterface.VariableReference})
    println("get! VariablePrimal... ")
    solitems = Float64[]
    if solutiondef(m.task,MSK_SOL_ITG) append!(solitems,MSK_SOL_ITG) end
    if solutiondef(m.task,MSK_SOL_BAS) append!(solitems,MSK_SOL_BAS) end
    if solutiondef(m.task,MSK_SOL_ITR) append!(solitems,MSK_SOL_ITR) end

    subj = Array{Int}(length(vs))
    for i in 1:length(subj)
        getindexes(m.x_block,ref2id(vs[i]),subj,i)
    end
    
    xx = getxx(m.task, solitems[attr.N])
    println("xx = $xx")
    output[1:length(output)] = xx[subj]
end

function MathOptInterface.getattribute(m::MosekModel,attr::MathOptInterface.VariablePrimal, vs::Vector{MathOptInterface.VariableReference})
    println("get VariablePrimal... ")
    output = Vector{Int64}(length(vs))
    MathOptInterface.getattribute!(output,m,attr,vs)
    output
end

function MathOptInterface.getattribute(m::MosekModel,attr::MathOptInterface.VariablePrimal, vref::MathOptInterface.VariableReference)
    solitems = Float64[]
    if solutiondef(m.task,MSK_SOL_ITG) append!(solitems,MSK_SOL_ITG) end
    if solutiondef(m.task,MSK_SOL_BAS) append!(solitems,MSK_SOL_BAS) end
    if solutiondef(m.task,MSK_SOL_ITR) append!(solitems,MSK_SOL_ITR) end

    subj = getindexes(m.x_block,ref2id(vref))
    xx = getxxslice(m.task, subj[1],subj[1]+1,solitems[attr.N])
    xx[1]
end







MathOptInterface.cangetattribute(m::MosekSolver,attr::MathOptInterface.VariableBasisStatus) = false # is a solution selector missing?
MathOptInterface.cangetattribute(m::MosekModel,attr::MathOptInterface.VariableBasisStatus)  = false

#### Constraint solution values

function getsolcode(m, N::Int)
    solitems = Float64[]
    if solutiondef(m.task,MSK_SOL_ITG) append!(solitems,MSK_SOL_ITG) end
    if solutiondef(m.task,MSK_SOL_BAS) append!(solitems,MSK_SOL_BAS) end
    if solutiondef(m.task,MSK_SOL_ITR) append!(solitems,MSK_SOL_ITR) end
    return solitems[N]
end

MathOptInterface.cangetattribute(m::MosekSolver,attr::MathOptInterface.ConstraintPrimal) = true
function MathOptInterface.cangetattribute(m::MosekModel,attr::MathOptInterface.ConstraintPrimal)
    num = 0
    if solutiondef(m.task,MSK_SOL_ITG) num += 1 end
    if solutiondef(m.task,MSK_SOL_BAS) num += 1 end
    if solutiondef(m.task,MSK_SOL_ITR) num += 1 end

    attr.N > 0 && attr.N <= num
end

function MathOptInterface.getattribute!{D}(
    output::Vector{Float64},
    m     ::MosekModel,
    attr  ::MathOptInterface.ConstraintPrimal,
    cref  ::MathOptInterface.ConstraintReference{MathOptInterface.SingleVariable,D})

    whichsol = getsolcode(m,attr.N)
    
    conid = ref2id(cref)
    idxs  = getindexes(m.xc_block,conid)
    subj  = m.xc_idxs[idxs[1]]
    
    xx = getxxslice(m.task,whichsol, Int(subj),Int(subj+1))
    xx[1]
end

# Semidefinite domain for a variable 
function MathOptInterface.getattribute!(
    output::Vector{Float64},
    m     ::MosekModel,
    attr  ::MathOptInterface.ConstraintPrimal,
    cref  ::MathOptInterface.ConstraintReference{MathOptInterface.VectorOfVariables,MathOptInterface.PositiveSemidefiniteConeTriangle})

    whichsol = getsolcode(m,attr.N)
    cid = ref2id(cref)
    assert(cid < 0)

    # It is in fact a real constraint and cid is the id of an ordinary constraint
    barvaridx = - m.c_slack[-cid]
    getbarxj(m.task,whichsol,barvaridx)
end

# Any other domain for variable vector
function MathOptInterface.getattribute!{D}(
    output::Vector{Float64},
    m     ::MosekModel,
    attr  ::MathOptInterface.ConstraintPrimal,
    cref  ::MathOptInterface.ConstraintReference{MathOptInterface.VectorOfVariables,D})

    whichsol = getsolcode(m,attr.N)
    xcid = ref2id(cref)
    assert(xcid > 0)
    
    mask = m.xc_bounds[xcid]
    idxs = getindexes(m.xc_block,xcid)
    subj = m.xc_idxs[idxs]

    if     mask & (boundflag_lower | boundflag_upper) != 0 # ranged,
        slx = getslx(m.task,whichsol)
        sux = getslx(m.task,whichsol)

        slx[subj] - sux[subj]
    elseif mask & boundflag_lower != 0
        slx = getslx(m.task,whichsol)
        slx[subj]
    elseif mask & boundflag_upper != 0
        sux = getslx(m.task,whichsol)
        -sux[subj]
    elseif mask & boundflag_cone != 0
        sux = getslx(m.task,whichsol)
        snx[subj]
    else
        error("No primal value for this type of constraint") # e.g. integer constraints
    end
end

solsize{F <: MathOptInterface.AbstractScalarFunction,D}(m::MosekModel, cref :: MathOptInterface.ConstraintReference{F,D}) = 1
function solsize{D}(m::MosekModel, cref :: MathOptInterface.ConstraintReference{MathOptInterface.VectorOfVariables,D})
    cid = ref2id(cref)
    if cid < 0
        blocksize(m.c_block,-cid)
    else
        blocksize(m.xc_block,cid)
    end
end

function solsize{D}(m::MosekModel, cref :: MathOptInterface.ConstraintReference{MathOptInterface.VectorAffineFunction,D})
    cid = ref2id(cref)
    blocksize(m.c_block,cid)
end



function MathOptInterface.getattribute{F,D}(m::MosekModel,attr::MathOptInterface.ConstraintPrimal, cref :: MathOptInterface.ConstraintReference{F,D})
    cid = ref2id(cref)
    output = Vector{Int64}(solsize(m,cref))
    MathOptInterface.getattribute!(output,attr,vs)
    output
end




#### Status codes
MathOptInterface.cangetattribute(m::MosekModel,attr::MathOptInterface.TerminationStatus) = true
function MathOptInterface.getattribute(m::MosekModel,attr::MathOptInterface.TerminationStatus)
    if     m.trm == MSK_RES_OK
        MathOptInterface.Success
    elseif m.trm == MSK_RES_TRM_MAX_ITERATIONS
        MathOptInterface.IterationLimit
    elseif m.trm == MSK_RES_TRM_MAX_TIME
        MathOptInterface.TimeLimit
    elseif m.trm == MSK_RES_TRM_OBJECTIVE_RANGE
        MathOptInterface.ObjectiveLimit
    elseif m.trm == MSK_RES_TRM_MIO_NEAR_REL_GAP
        MathOptInterface.AlmostSuccess
    elseif m.trm == MSK_RES_TRM_MIO_NEAR_ABS_GAP
        MathOptInterface.AlmostSuccess
    elseif m.trm == MSK_RES_TRM_MIO_NUM_RELAXS
        MathOptInterface.OtherLimit
    elseif m.trm == MSK_RES_TRM_MIO_NUM_BRANCHES
        MathOptInterface.NodeLimit
    elseif m.trm == MSK_RES_TRM_NUM_MAX_NUM_INT_SOLUTIONS
        MathOptInterface.SolutionLimit
    elseif m.trm == MSK_RES_TRM_STALL
        MathOptInterface.SlowProgress
    elseif m.trm == MSK_RES_TRM_USER_CALLBACK
        MathOptInterface.Interrupted
    elseif m.trm == MSK_RES_TRM_MAX_NUM_SETBACKS
        MathOptInterface.OtherLimit
    elseif m.trm == MSK_RES_TRM_NUMERICAL_PROBLEM
        MathOptInterface.SlowProgress
    elseif m.trm == MSK_RES_TRM_INTERNAL
        MathOptInterface.OtherError
    elseif m.trm == MSK_RES_TRM_INTERNAL_STOP
        MathOptInterface.OtherError
    else
        MathOptInterface.OtherError
    end
end

MathOptInterface.cangetattribute(m::MosekModel,attr::MathOptInterface.PrimalStatus) = true
function MathOptInterface.getattribute(m::MosekModel,attr::MathOptInterface.PrimalStatus)
    solitems = Int32[]
    if solutiondef(m.task,MSK_SOL_ITG) append!(solitems,MSK_SOL_ITG) end
    if solutiondef(m.task,MSK_SOL_BAS) append!(solitems,MSK_SOL_BAS) end
    if solutiondef(m.task,MSK_SOL_ITR) append!(solitems,MSK_SOL_ITR) end

    solsta = getsolsta(m.task,solitems[attr.N])

    if     solsta == MSK_SOL_STA_UNKNOWN
        MathOptInterface.UnknownResultStatus
    elseif solsta == MSK_SOL_STA_OPTIMAL
        MathOptInterface.FeasiblePoint
    elseif solsta == MSK_SOL_STA_PRIM_FEAS
        MathOptInterface.FeasiblePoint
    elseif solsta == MSK_SOL_STA_DUAL_FEAS
        MathOptInterface.UnknownResultStatus
    elseif solsta == MSK_SOL_STA_PRIM_AND_DUAL_FEAS
        MathOptInterface.FeasiblePoint
    elseif solsta == MSK_SOL_STA_NEAR_OPTIMAL
        MathOptInterface.NearlyFeasiblePoint
    elseif solsta == MSK_SOL_STA_NEAR_PRIM_FEAS
        MathOptInterface.NearlyFeasiblePoint
    elseif solsta == MSK_SOL_STA_NEAR_DUAL_FEAS
        MathOptInterface.UnknownResultStatus
    elseif solsta == MSK_SOL_STA_NEAR_PRIM_AND_DUAL_FEAS
        MathOptInterface.NearFeasiblePoint
    elseif solsta == MSK_SOL_STA_PRIM_INFEAS_CER
        MathOptInterface.UnknownResultStatus
    elseif solsta == MSK_SOL_STA_DUAL_INFEAS_CER
        MathOptInterface.InfeasibilityCertificate
    elseif solsta == MSK_SOL_STA_NEAR_PRIM_INFEAS_CER
        MathOptInterface.UnknownResultStatus
    elseif solsta == MSK_SOL_STA_NEAR_DUAL_INFEAS_CER
        MathOptInterface.NearlyInfeasibilityCertificate
    elseif solsta == MSK_SOL_STA_PRIM_ILLPOSED_CER
        MathOptInterface.UnknownResultStatus
    elseif solsta == MSK_SOL_STA_DUAL_ILLPOSED_CER
        MathOptInterface.ReductionCertificate
    elseif solsta == MSK_SOL_STA_INTEGER_OPTIMAL
        MathOptInterface.FeasiblePoint
    elseif solsta == MSK_SOL_STA_NEAR_INTEGER_OPTIMAL
        MathOptInterface.NearlyFeasiblePoint
    else
        MathOptInterface.UnknownResultStatus
    end
end

MathOptInterface.cangetattribute(m::MosekModel,attr::MathOptInterface.DualStatus) = true
function MathOptInterface.getattribute(m::MosekModel,attr::MathOptInterface.DualStatus)
    solitems = Int32[]
    if solutiondef(m.task,MSK_SOL_ITG) append!(solitems,MSK_SOL_ITG) end
    if solutiondef(m.task,MSK_SOL_BAS) append!(solitems,MSK_SOL_BAS) end
    if solutiondef(m.task,MSK_SOL_ITR) append!(solitems,MSK_SOL_ITR) end

    solsta = getsolsta(m.task,solitems[attr.N])

    if     solsta == MSK_SOL_STA_UNKNOWN
        MathOptInterface.UnknownResultStatus
    elseif solsta == MSK_SOL_STA_OPTIMAL
        MathOptInterface.FeasiblePoint
    elseif solsta == MSK_SOL_STA_PRIM_FEAS
        MathOptInterface.UnknownResultStatus
    elseif solsta == MSK_SOL_STA_DUAL_FEAS
        MathOptInterface.FeasiblePoint
    elseif solsta == MSK_SOL_STA_PRIM_AND_DUAL_FEAS
        MathOptInterface.FeasiblePoint
    elseif solsta == MSK_SOL_STA_NEAR_OPTIMAL
        MathOptInterface.NearlyFeasiblePoint
    elseif solsta == MSK_SOL_STA_NEAR_PRIM_FEAS
        MathOptInterface.UnknownResultStatus
    elseif solsta == MSK_SOL_STA_NEAR_DUAL_FEAS
        MathOptInterface.NearlyFeasiblePoint
    elseif solsta == MSK_SOL_STA_NEAR_PRIM_AND_DUAL_FEAS
        MathOptInterface.NearFeasiblePoint
    elseif solsta == MSK_SOL_STA_PRIM_INFEAS_CER
        MathOptInterface.InfeasibilityCertificate
    elseif solsta == MSK_SOL_STA_DUAL_INFEAS_CER
        MathOptInterface.UnknownResultStatus
    elseif solsta == MSK_SOL_STA_NEAR_PRIM_INFEAS_CER
        MathOptInterface.NearlyInfeasibilityCertificate
    elseif solsta == MSK_SOL_STA_NEAR_DUAL_INFEAS_CER
        MathOptInterface.UnknownResultStatus
    elseif solsta == MSK_SOL_STA_PRIM_ILLPOSED_CER
        MathOptInterface.ReductionCertificate
    elseif solsta == MSK_SOL_STA_DUAL_ILLPOSED_CER
        MathOptInterface.UnknownResultStatus
    elseif solsta == MSK_SOL_STA_INTEGER_OPTIMAL
        MathOptInterface.FeasiblePoint
    elseif solsta == MSK_SOL_STA_NEAR_INTEGER_OPTIMAL
        MathOptInterface.NearlyFeasiblePoint
    else
        MathOptInterface.UnknownResultStatus
    end
end
