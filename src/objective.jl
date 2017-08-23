



function MathOptInterface.setobjective!(m::MosekModel, sense :: MathOptInterface.OptimizationSense, func::MathOptInterface.SingleVariable)
    if sense == MathOptInterface.MinSense
        putobjsense(m.task,MSK_OBJECTIVE_SENSE_MINIMIZE)
    else
        putobjsense(m.task,MSK_OBJECTIVE_SENSE_MAXIMIZE)
    end

    numvar = getnumvar(m.task)
    c = zeros(Float64,numvar)
    vid = ref2id(func.variable)
    subj = getindexes(m.x_block,vid)

    c[subj[1]] = 1.0
    
    putclist(m.task,Int32[1:numvar...],c)
    putcfix(m.task,0.0)
end

function MathOptInterface.setobjective!(m::MosekModel, sense :: MathOptInterface.OptimizationSense, func::MathOptInterface.ScalarAffineFunction{Float64})
    if sense == MathOptInterface.MinSense
        putobjsense(m.task,MSK_OBJECTIVE_SENSE_MINIMIZE)
    else
        putobjsense(m.task,MSK_OBJECTIVE_SENSE_MAXIMIZE)
    end

    numvar = getnumvar(m.task)
    c = zeros(Float64,numvar)
    vids = [ ref2id(vid) for vid in func.variables ]    
    subj = Vector{Int}(length(vids))
    for i in 1:length(vids)
        getindexes(m.x_block,vids[i],subj,i)
    end

    c[subj] = func.coefficients
    
    putclist(m.task,Int32[1:numvar...],c)
    putcfix(m.task,func.constant)
end

MathOptInterface.canmodifyobjective(m::MosekModel, change :: MathOptInterface.ScalarConstantChange) = true
MathOptInterface.canmodifyobjective(m::MosekModel, change :: MathOptInterface.ScalarCoefficientChange) = true

function MathOptInterface.modifyobjective!(m::MosekModel, change :: MathOptInterface.ScalarConstantChange)
    putcfix(m.task,change.new_constant)
end

function MathOptInterface.modifyobjective!(m::MosekModel, change :: MathOptInterface.ScalarCoefficientChange)
    vid = ref2id(change.variable)
    subj = getindexes(m.x_block,vid)
    putcj(m.task,Int32(subj[1]),change.new_coefficient)
end
