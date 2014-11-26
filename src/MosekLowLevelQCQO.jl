

# Note: 
#  We do not covert quadratic objective to SOCP
setquadobj!(m::MosekMathProgModel, rowidx,colidx,quadval) =
    setquadobj!(m,
                convert(Array{Int32,1},rowidx),
                convert(Array{Int32,1},colidx),
                convert(Array{Float64,1},quadval))

function setquadobj!(m::MosekMathProgModel, 
                     rowidx::Array{Int32,1},
                     colidx::Array{Int32,1},
                     quadval::Array{Float64,1})
  upgradeProbType(m,MosekMathProgModel_QOQP)
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



addquadconstr!(m::MosekMathProgModel, linearidx, linearval, quadrowidx, quadcolidx, quadval, sense, rhs) =    
    addquadconstr!(m,
                   convert(Array{Int32,1},  linearidx),
                   convert(Array{Float64,1},linearval),
                   convert(Array{Int32,1},  quadrowidx),
                   convert(Array{Int32,1},  quadcolidx),
                   convert(Array{Float64,1},quadval),
                   sense,
                   convert(Float64,rhs))

# Note:
#  If the quadratic terms define a quadratic cone (only the special
#  forms are detected), the linear terms, sense and rhs are ignored.
#  
function addquadconstr!(m::MosekMathProgModel,
                        linearidx ::Array{Int32,1}, 
                        linearval ::Array{Float64,1}, 
                        quadrowidx::Array{Int32,1}, 
                        quadcolidx::Array{Int32,1}, 
                        quadval   ::Array{Float64,1}, 
                        sense, 
                        rhs       ::Float64)
  subj    = Int32[ m.varmap[i] for i=linearidx ]
  valj    = linearval
  qcksubi = Int32[ m.varmap[i] for i=quadrowidx ]
  qcksubj = Int32[ m.varmap[i] for i=quadcolidx ]
  qckval  = quadval

  # detect SOCP form
  ct,x =
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
          x = zeros(Int64,length(qcksubj))
          x[1] = qcksubj[negdiag_idx]
          if negdiag_idx > 1
              x[2:negdiag_idx] = qcksubj[1:negdiag_idx-1]
          end
          if negdiag_idx < length(qcksubj)
              x[negdiag_idx+1:length(qcksubj)] = qcksubj[negdiag_idx+1:length(qcksubj)]
          end

          MSK_CT_QUAD, x
        elseif num_posonediag == length(qcksubj)-1 && offdiag_idx > 0
          x = zeros(Int64,length(qcksubj)+1)            
          x[1] = qcksubi[offdiag_idx]
          x[2] = qcksubj[offdiag_idx]
          if offdiag_idx > 1
              x[3:offdiag_idx+1] = qcksubj[1:offdiag_idx-1]
          end
          if offdiag_idx < length(qcksubj)
              x[offdiag_idx+2:length(qcksubj)+1] = qcksubj[offdiag_idx+1:length(qcksubj)]
          end
            
          MSK_CT_RQUAD, x
        else
          -1,()
        end
      end

  if ct in [ MSK_CT_QUAD, MSK_CT_RQUAD ]
    upgradeProbType(m,MosekMathProgModel_SOCP)
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
    addUserCon(m,dummycon)

    m.conslack[dummycon] = 0
    putconbound(m.task,dummycon, MSK_BK_FX, 0.0, 0.0)

  else
    upgradeProbType(m,MosekMathProgModel_QOQP)

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

    putarow(m.task, k, subj, valj)
    putqconk(m.task,k, qcksubi,qcksubj,qckval)

    if sense == '<'
      putconbound(m.task,k,MSK_BK_UP, -Inf,rhs)
    elseif sense == '>'
      putconbound(m.task,k,MSK_BK_LO, rhs,Inf)
    else
      putconbound(m.task,k,MSK_BK_FR, -Inf,Inf)
    end

    addUserCon(m,considx)
    m.conslack[considx] = 0
  end
end
