
#### solver attributes
MathOptInterface.getattribute(m::Union{MosekSolver,MosekModel},::MathOptInterface.SupportsDuals) = true
MathOptInterface.getattribute(m::Union{MosekSolver,MosekModel},::MathOptInterface.SupportsAddConstraintAfterSolve) = true
MathOptInterface.getattribute(m::Union{MosekSolver,MosekModel},::MathOptInterface.SupportsAddVariableAfterSolve) = true
MathOptInterface.getattribute(m::Union{MosekSolver,MosekModel},::MathOptInterface.SupportsDeleteConstraint) = true
MathOptInterface.getattribute(m::Union{MosekSolver,MosekModel},::MathOptInterface.SupportsDeleteVariable) = true
MathOptInterface.getattribute(m::Union{MosekSolver,MosekModel},::MathOptInterface.SupportsConicThroughQuadratic) = false # though actually the solver does

#### objective
MathOptInterface.getattribute(m::MosekModel,attr::MathOptInterface.ObjectiveValue) = getprimalobj(m.task,m.solutions[attr.resultindex].whichsol)
MathOptInterface.cangetattribute(m::MosekSolver,attr::MathOptInterface.ObjectiveValue) = true
MathOptInterface.cangetattribute(m::MosekModel,attr::MathOptInterface.ObjectiveValue) = attr.resultindex > 0 && attr.resultindex <= length(m.solutions)

MathOptInterface.cangetattribute(m::MosekModel,attr::MathOptInterface.ObjectiveBound) = getintinf(m.task,MSK_IINF_MIO_NUM_RELAX) > 0
MathOptInterface.getattribute(m::MosekModel,attr::MathOptInterface.ObjectiveBound) = getdouinf(m.task,MSK_DINF_MIO_OBJ_BOUND)

MathOptInterface.cangetattribute(m::MosekModel,attr::MathOptInterface.RelativeGap) = getdouinf(m.task,MSK_DINF_MIO_OBJ_REL_GAP) > -1.0
MathOptInterface.getattribute(m::MosekModel,attr::MathOptInterface.RelativeGap) = getdouinf(m.task,MSK_DINF_MIO_OBJ_REL_GAP)

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
    miosimiter = getlintinf(m.task,MSK_LIINF_MIO_SIMPLEX_ITER)
    if miosimiter > 0
        Int(miosimiter)
    else
        Int(getintinf(m.task,MSK_IINF_SIM_PRIMAL_ITER) + getintinf(m.task,MSK_IINF_SIM_DUAL_ITER))
    end
end

MathOptInterface.cangetattribute(m::MosekSolver,attr::MathOptInterface.BarrierIterations) = true
MathOptInterface.cangetattribute(m::MosekModel,attr::MathOptInterface.BarrierIterations) = true
function MathOptInterface.getattribute(m::MosekModel,attr::MathOptInterface.BarrierIterations)
    miosimiter = getlintinf(m.task,MSK_LIINF_MIO_INTPNT_ITER)
    if miosimiter > 0
        Int(miosimiter)
    else
        Int(getintinf(m.task,MSK_IINF_INTPNT_ITER))
    end
end

MathOptInterface.cangetattribute(m::MosekSolver,attr::MathOptInterface.NodeCount) = true
MathOptInterface.cangetattribute(m::MosekModel,attr::MathOptInterface.NodeCount) = true
function MathOptInterface.getattribute(m::MosekModel,attr::MathOptInterface.NodeCount)
        Int(getintinf(m.task,MSK_IINF_MIO_NUM_BRANCH))
end

MathOptInterface.cangetattribute(m::MosekSolver,attr::MathOptInterface.RawSolver) = true
MathOptInterface.cangetattribute(m::MosekModel,attr::MathOptInterface.RawSolver) = true
MathOptInterface.getattribute(m::MosekModel,attr::MathOptInterface.RawSolver) = m.task

MathOptInterface.cangetattribute(m::MosekSolver,attr::MathOptInterface.ResultCount) = true
MathOptInterface.cangetattribute(m::MosekModel,attr::MathOptInterface.ResultCount) = true
MathOptInterface.getattribute(m::MosekModel,attr::MathOptInterface.ResultCount) = length(m.solutions)
    
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

MathOptInterface.cangetattribute(m::MosekModel,attr::MathOptInterface.VariablePrimal) = attr.N > 0 && attr.N <= length(m.solutions)
MathOptInterface.cangetattribute(m::MosekSolver,attr::MathOptInterface.VariablePrimal) = true 

MathOptInterface.cangetattribute(m::MosekModel,attr::MathOptInterface.VariablePrimal,vs::Vector{MathOptInterface.VariableReference}) = MathOptInterface.cangetattribute(m,attr)
MathOptInterface.cangetattribute(m::MosekModel,attr::MathOptInterface.VariablePrimal,v::MathOptInterface.VariableReference) = MathOptInterface.cangetattribute(m,attr)


function MathOptInterface.getattribute!(output::Vector{Float64},m::MosekModel,attr::MathOptInterface.VariablePrimal, vs::Vector{MathOptInterface.VariableReference})
    subj = Array{Int}(length(vs))
    for i in 1:length(subj)
        getindexes(m.x_block,ref2id(vs[i]),subj,i)
    end
    
    output[1:length(output)] = m.solutions[attr.N].xx[subj]
end

function MathOptInterface.getattribute(m::MosekModel,attr::MathOptInterface.VariablePrimal, vs::Vector{MathOptInterface.VariableReference})
    output = Vector{Float64}(length(vs))
    MathOptInterface.getattribute!(output,m,attr,vs)
    output
end

function MathOptInterface.getattribute(m::MosekModel,attr::MathOptInterface.VariablePrimal, vref::MathOptInterface.VariableReference)
    subj = getindexes(m.x_block,ref2id(vref))[1]
    m.solutions[attr.N].xx[subj]
end






MathOptInterface.cangetattribute(m::MosekSolver,attr::MathOptInterface.VariableBasisStatus) = false # is a solution selector missing?
MathOptInterface.cangetattribute(m::MosekModel,attr::MathOptInterface.VariableBasisStatus)  = false

#### Constraint solution values
function cangetattribute(m::MosekModel,attr::MathOptInterface.ConstraintPrimal)
    num = 0
    if solutiondef(m.task,MSK_SOL_ITG) num += 1 end
    if solutiondef(m.task,MSK_SOL_BAS) num += 1 end
    if solutiondef(m.task,MSK_SOL_ITR) num += 1 end
    attr.N > 0 && attr.N <= num
end

MathOptInterface.cangetattribute(m::MosekModel,attr::MathOptInterface.ConstraintPrimal, c) = cangetattribute(m,attr)
MathOptInterface.cangetattribute{F,D}(m::MosekModel,attr::MathOptInterface.ConstraintPrimal, cref::MathOptInterface.ConstraintReference{F,D}) = cangetattribute(m,attr)
MathOptInterface.cangetattribute(m::MosekModel,attr::MathOptInterface.ConstraintPrimal, cref::MathOptInterface.ConstraintReference{MathOptInterface.ScalarAffineFunction{Float64},MathOptInterface.LessThan{Float64}}) = cangetattribute(m,attr)

function MathOptInterface.getattribute{D}(
    m     ::MosekModel,
    attr  ::MathOptInterface.ConstraintPrimal,
    cref  ::MathOptInterface.ConstraintReference{MathOptInterface.SingleVariable,D})
    
    conid = ref2id(cref)
    idxs  = getindexes(m.xc_block,conid)
    subj  = m.xc_idxs[idxs[1]]
    
    m.solutions[attr.N].xx[subj]
end

# Semidefinite domain for a variable 
function MathOptInterface.getattribute!(
    output::Vector{Float64},
    m     ::MosekModel,
    attr  ::MathOptInterface.ConstraintPrimal,
    cref  ::MathOptInterface.ConstraintReference{MathOptInterface.VectorOfVariables,MathOptInterface.PositiveSemidefiniteConeTriangle})

    cid = ref2id(cref)
    assert(cid < 0)

    # It is in fact a real constraint and cid is the id of an ordinary constraint
    barvaridx = - m.c_slack[-cid]    
    output[1:length(output)] = getbarxj(m.task,m.solutions[attr.N].whichsol,barvaridx)
end

# Any other domain for variable vector
function MathOptInterface.getattribute!{D}(
    output::Vector{Float64},
    m     ::MosekModel,
    attr  ::MathOptInterface.ConstraintPrimal,
    cref  ::MathOptInterface.ConstraintReference{MathOptInterface.VectorOfVariables,D})

    xcid = ref2id(cref)
    assert(xcid > 0)
    
    mask = m.xc_bounds[xcid]
    idxs = getindexes(m.xc_block,xcid)
    subj = m.xc_idxs[idxs]

    output[1:length(output)] = m.solutions[attr.N].xx[subj]
end

function MathOptInterface.getattribute{D}(m     ::MosekModel,
                                          attr  ::MathOptInterface.ConstraintPrimal,
                                          cref  ::MathOptInterface.ConstraintReference{MathOptInterface.ScalarAffineFunction{Float64},D})
    cid = ref2id(cref)
    subi = getindexes(m.c_block,cid)[1]
    m.solutions[attr.N].xc[subi]
end



function MathOptInterface.getattribute!{D}(
    output::Vector{Float64},
    m     ::MosekModel,
    attr  ::MathOptInterface.ConstraintPrimal,
    cref  ::MathOptInterface.ConstraintReference{MathOptInterface.VectorAffineFunction{Float64},D})

    cid = ref2id(cref)
    subi = getindexes(m.c_block,cid)

    if     m.c_block_slack[cid] == 0 # no slack
        output[1:length(output)] = m.solutions[attr.N].xc[subi] + m.c_constant[subi]
    elseif m.c_block_slack[cid] >  0 # qcone slack
        xid = m.c_block_slack[cid]
        xsubj = getindexes(m.x_block, xid)
        output[1:length(output)] = m.solutions[attr.N].xx[xsubj]  + m.c_constant[subi]
    elseif m.c_block_slack[cid]  # psd slack
        xid = - m.c_block_slack[cid]
        output[1:length(output)] = getbarxj(m.task,m.solutions[attr.N].whichsol,Int32(xid)) + m.c_constant[subi]
    end
end













MathOptInterface.cangetattribute(m::MosekModel,attr::MathOptInterface.ConstraintDual) = attr.N > 0 && attr.N <= length(m.solutions)
MathOptInterface.cangetattribute(m::MosekModel,attr::MathOptInterface.ConstraintDual, crefs::Vector{MathOptInterface.ConstraintReference}) = MathOptInterface.cangetattribute(m,attr)
MathOptInterface.cangetattribute(m::MosekModel,
                                 attr::MathOptInterface.ConstraintDual,
                                 cref::MathOptInterface.ConstraintReference{F,D}) where {F,D} = MathOptInterface.cangetattribute(m,attr)


function MathOptInterface.getattribute(
    m     ::MosekModel,
    attr  ::MathOptInterface.ConstraintDual,
    cref  ::MathOptInterface.ConstraintReference{MathOptInterface.SingleVariable,D}) where { D <: MathOptInterface.AbstractSet }


    xcid = ref2id(cref)
    idxs = getindexes(m.xc_block,xcid) # variable ids

    assert(blocksize(m.xc_block,xcid) > 0)

    subj  = getindexes(m.x_block, m.xc_idxs[idxs][1])[1]
    if (getobjsense(m.task) == MSK_OBJECTIVE_SENSE_MINIMIZE)
        if     m.xc_bounds[xcid] & boundflag_lower != 0 && m.xc_bounds[xcid] & boundflag_upper != 0
            m.solutions[attr.N].slx[subj] - m.solutions[attr.N].sux[subj]
        elseif (m.xc_bounds[xcid] & boundflag_lower) != 0
            m.solutions[attr.N].slx[subj]
        elseif (m.xc_bounds[xcid] & boundflag_upper) != 0 
            - m.solutions[attr.N].sux[subj]
        elseif (m.xc_bounds[xcid] & boundflag_cone) != 0
            m.solutions[attr.N].snx[subj]
        else
            error("Dual value available for this constraint")
        end
    else
        if     m.xc_bounds[xcid] & boundflag_lower != 0 && m.xc_bounds[xcid] & boundflag_upper != 0
            m.solutions[attr.N].sux[subj] - m.solutions[attr.N].slx[subj]
        elseif (m.xc_bounds[xcid] & boundflag_lower) != 0
            - m.solutions[attr.N].slx[subj]
        elseif (m.xc_bounds[xcid] & boundflag_upper) != 0 
            m.solutions[attr.N].sux[subj]
        elseif (m.xc_bounds[xcid] & boundflag_cone) != 0
            - m.solutions[attr.N].snx[subj]
        else
            error("Dual value available for this constraint")
        end
    end
end

# Semidefinite domain for a variable 
function MathOptInterface.getattribute!(
    output::Vector{Float64},
    m     ::MosekModel,
    attr  ::MathOptInterface.ConstraintDual,
    cref  ::MathOptInterface.ConstraintReference{MathOptInterface.VectorOfVariables,MathOptInterface.PositiveSemidefiniteConeTriangle})

    whichsol = getsolcode(m,attr.N)
    cid = ref2id(cref)
    assert(cid < 0)

    # It is in fact a real constraint and cid is the id of an ordinary constraint
    barvaridx = - m.c_slack[-cid]
    if (getobjsense(m.task) == MSK_OBJECTIVE_SENSE_MINIMIZE)
        output[1:length(output)] = getbarsj(m.task,whichsol,barvaridx)
    else
        output[1:length(output)] = - getbarsj(m.task,whichsol,barvaridx)
    end
end

# Any other domain for variable vector
function MathOptInterface.getattribute!{D}(
    output::Vector{Float64},
    m     ::MosekModel,
    attr  ::MathOptInterface.ConstraintDual,
    cref  ::MathOptInterface.ConstraintReference{MathOptInterface.VectorOfVariables,D})

    xcid = ref2id(cref)
    assert(xcid > 0)
    
    mask = m.xc_bounds[xcid]
    idxs = getindexes(m.xc_block,xcid)
    subj = m.xc_idxs[idxs]

    if (getobjsense(m.task) == MSK_OBJECTIVE_SENSE_MINIMIZE)
        if     m.xc_bounds[conid] & boundflag_lower != 0 && m.xc_bounds[conid] & boundflag_upper != 0
            output[1:length(output)] = m.solutions[attr.N].slx[subj] - m.solutions[attr.N].sux[subj]
        elseif (m.xc_bounds[conid] & boundflag_lower) != 0
            output[1:length(output)] = m.solutions[attr.N].slx[subj]
        elseif (m.xc_bounds[conid] & boundflag_upper) != 0 
            output[1:length(output)] = - m.solutions[attr.N].sux[subj]
        elseif (m.xc_bounds[conid] & boundflag_cone) != 0
            output[1:length(output)] = m.solutions[attr.N].snx[subj]
        else
            error("Dual value available for this constraint")
        end
    else
        if     m.xc_bounds[conid] & boundflag_lower != 0 && m.xc_bounds[conid] & boundflag_upper != 0
            output[1:length(output)] = m.solutions[attr.N].sux[subj] - m.solutions[attr.N].slx[subj]
        elseif (m.xc_bounds[conid] & boundflag_lower) != 0
            output[1:length(output)] = - m.solutions[attr.N].slx[subj]
        elseif (m.xc_bounds[conid] & boundflag_upper) != 0 
            output[1:length(output)] = m.solutions[attr.N].sux[subj]
        elseif (m.xc_bounds[conid] & boundflag_cone) != 0
            output[1:length(output)] = - m.solutions[attr.N].snx[subj]
        else
            error("Dual value available for this constraint")
        end
    end
end

function MathOptInterface.getattribute{D}(m     ::MosekModel,
                                          attr  ::MathOptInterface.ConstraintDual,
                                          cref  ::MathOptInterface.ConstraintReference{MathOptInterface.ScalarAffineFunction{Float64},D})
    cid = ref2id(cref)
    subi = getindexes(m.c_block,cid)[1]

    if getobjsense(m.task) == MSK_OBJECTIVE_SENSE_MINIMIZE
        m.solutions[attr.N].y[subi]
    else
        - m.solutions[attr.N].y[subi]
    end
end


function MathOptInterface.getattribute(
    m     ::MosekModel,
    attr  ::MathOptInterface.ConstraintDual,
    cref  ::MathOptInterface.ConstraintReference{MathOptInterface.VectorAffineFunction{Float64},D}) where { D <: MathOptInterface.AbstractSet }
    
    n = blocksize(m.c_block,ref2id(cref))
    res = Vector{Float64}(n)
    MathOptInterface.getattribute!(res,m,attr,cref)
    res
end
function MathOptInterface.getattribute!(
    output::Vector{Float64},
    m     ::MosekModel,
    attr  ::MathOptInterface.ConstraintDual,
    cref  ::MathOptInterface.ConstraintReference{MathOptInterface.VectorAffineFunction{Float64},D}) where { D <: MathOptInterface.AbstractSet }

    cid = ref2id(cref)
    subi = getindexes(m.c_block,cid)

    if getobjsense(m.task) == MSK_OBJECTIVE_SENSE_MINIMIZE
        if     m.c_block_slack[cid] == 0 # no slack
            output[1:length(output)] = m.solutions[attr.N].y[subi]
        elseif m.c_block_slack[cid] >  0 # qcone slack
            xid = m.c_block_slack[cid]
            subj = getindexes(m.x_block, xid)
            output[1:length(output)] = m.solutions[attr.N].snx[xsubj]
        elseif m.c_block_slack[cid]  # psd slack
            xid = - m.c_block_slack[cid]
            output[1:length(output)] = getbarsj(m.task,m.solutions[attr.N].whichsol,Int32(xid))
        end
    else
        if     m.c_block_slack[cid] == 0 # no slack
            output[1:length(output)] = - m.solutions[attr.N].y[subi]
        elseif m.c_block_slack[cid] >  0 # qcone slack
            xid = m.c_block_slack[cid]
            subj = getindexes(m.x_block, xid)
            output[1:length(output)] = - m.solutions[attr.N].snx[subj]
        elseif m.c_block_slack[cid]  # psd slack
            xid = - m.c_block_slack[cid]
            output[1:length(output)] = - getbarsj(m.task,m.solutions[attr.N].whichsol,Int32(xid))
        end
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

function MathOptInterface.getattribute{F <: MathOptInterface.AbstractScalarFunction,D}(m::MosekModel,attr::MathOptInterface.ConstraintPrimal, cref :: MathOptInterface.ConstraintReference{F,D})
    cid = ref2id(cref)
    output = Vector{Float64}(solsize(m,cref))
    MathOptInterface.getattribute!(output,m,attr,cref)
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

function MathOptInterface.cangetattribute(m::MosekModel,attr::MathOptInterface.PrimalStatus)
    if attr.N < 0 || attr.N > length(m.solutions)
        false
    else
        solsta = m.solutions[attr.N].solsta
        if     solsta == MSK_SOL_STA_UNKNOWN true
        elseif solsta == MSK_SOL_STA_OPTIMAL true
        elseif solsta == MSK_SOL_STA_PRIM_FEAS true
        elseif solsta == MSK_SOL_STA_DUAL_FEAS true
        elseif solsta == MSK_SOL_STA_PRIM_AND_DUAL_FEAS true
        elseif solsta == MSK_SOL_STA_NEAR_OPTIMAL true
        elseif solsta == MSK_SOL_STA_NEAR_PRIM_FEAS true
        elseif solsta == MSK_SOL_STA_NEAR_DUAL_FEAS true
        elseif solsta == MSK_SOL_STA_NEAR_PRIM_AND_DUAL_FEAS true
        elseif solsta == MSK_SOL_STA_PRIM_INFEAS_CER false
        elseif solsta == MSK_SOL_STA_DUAL_INFEAS_CER true
        elseif solsta == MSK_SOL_STA_NEAR_PRIM_INFEAS_CER false
        elseif solsta == MSK_SOL_STA_NEAR_DUAL_INFEAS_CER true
        elseif solsta == MSK_SOL_STA_PRIM_ILLPOSED_CER false
        elseif solsta == MSK_SOL_STA_DUAL_ILLPOSED_CER true
        elseif solsta == MSK_SOL_STA_INTEGER_OPTIMAL true
        elseif solsta == MSK_SOL_STA_NEAR_INTEGER_OPTIMAL true
        else false
        end
    end
end
function MathOptInterface.getattribute(m::MosekModel,attr::MathOptInterface.PrimalStatus)
    solsta = m.solutions[attr.N].solsta
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

function MathOptInterface.cangetattribute(m::MosekModel,attr::MathOptInterface.DualStatus)
    if attr.N < 0 || attr.N > length(m.solutions)
        false
    else
        solsta = m.solutions[attr.N].solsta
        if     solsta == MSK_SOL_STA_UNKNOWN true
        elseif solsta == MSK_SOL_STA_OPTIMAL true
        elseif solsta == MSK_SOL_STA_PRIM_FEAS true
        elseif solsta == MSK_SOL_STA_DUAL_FEAS true
        elseif solsta == MSK_SOL_STA_PRIM_AND_DUAL_FEAS true
        elseif solsta == MSK_SOL_STA_NEAR_OPTIMAL true
        elseif solsta == MSK_SOL_STA_NEAR_PRIM_FEAS true
        elseif solsta == MSK_SOL_STA_NEAR_DUAL_FEAS true
        elseif solsta == MSK_SOL_STA_NEAR_PRIM_AND_DUAL_FEAS true
        elseif solsta == MSK_SOL_STA_PRIM_INFEAS_CER true
        elseif solsta == MSK_SOL_STA_DUAL_INFEAS_CER false
        elseif solsta == MSK_SOL_STA_NEAR_PRIM_INFEAS_CER true
        elseif solsta == MSK_SOL_STA_NEAR_DUAL_INFEAS_CER false
        elseif solsta == MSK_SOL_STA_PRIM_ILLPOSED_CER true
        elseif solsta == MSK_SOL_STA_DUAL_ILLPOSED_CER false
        elseif solsta == MSK_SOL_STA_INTEGER_OPTIMAL false
        elseif solsta == MSK_SOL_STA_NEAR_INTEGER_OPTIMAL false
        else false
        end
    end
end
function MathOptInterface.getattribute(m::MosekModel,attr::MathOptInterface.DualStatus)
    solsta = m.solutions[attr.N].solsta
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
