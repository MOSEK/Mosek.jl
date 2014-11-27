

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
  qckval  = quadval

  # detect SOCP form
  ct,idxs =
      let num_posonediag = 0,
          offdiag_idx    = 0,
          negdiag_idx    = 0
        for i=1:length(quadrowidx)
          if quadrowidx[i] == quadcolidx[i]
            if abs(qckval[i]-1.0) < 1e-12
              num_posonediag += 1
            elseif abs(qckval[i]+1.0) < 1e-12
              negdiag_idx = i
            end
          elseif quadrowidx[i] != quadcolidx[i]
            if abs(qckval[i]+1.0) < 1e-12
              offdiag_idx = i
            end
          end
        end

        if num_posonediag == length(quadcolidx)-1 && negdiag_idx > 0
          idxs = zeros(Int,length(quadcolidx))
          idxs[1] = quadrowidx[negdiag_idx]

          for i=1:negdiag_idx-1                  idxs[i+1] = quadrowidx[i] end
          for i=negdiag_idx+1:length(quadcolidx) idxs[i]   = quadrowidx[i] end

          MSK_CT_QUAD, idxs
        elseif num_posonediag == length(quadcolidx)-1 && offdiag_idx > 0
          idxs = zeros(Int,length(quadcolidx)+1)
          idxs[1] = quadrowidx[offdiag_idx]
          idxs[2] = quadcolidx[offdiag_idx]
          for i=1:offdiag_idx-1                  idxs[i+2] = quadcolidx[i] end
          for i=offdiag_idx+1:length(quadcolidx) idxs[i+1] = quadcolidx[i] end
            
          MSK_CT_RQUAD, idxs
        else
          -1,()
        end
      end

  if ct in [ MSK_CT_QUAD, MSK_CT_RQUAD ]
      # SOCP and SDP can be mixed, SDP includes SOCP
      upgradeProbType(m,MosekMathProgModel_SOCP)
      
      x = m.varmap[idxs]
      n = length(x)
      numbarelm = count(xitem -> xitem < 0, x)
      
      
      nvar  = getnumvar(m.task)+1
      ncon  = getnumcon(m.task)+1
      appendvars(m.task,n) # create aux variable z - these will be put into a cone
      appendcons(m.task,n) 
      
      z     = nvar
      cof  = Array(Float64,n*2-numbarelm)
      subj = Array(Int32,  n*2-numbarelm)
      ptr  = Array(Int64,  n+1)

      if numbarelm == 0
          cof[1:2:2*n]  = 1.0
          cof[2:2:2*n]  = -1.0
          subj[1:2:2*n] = z:z+n-1
          subj[2:2:2*n] = x
          ptr[:]        = Int64[1:2:2*n+1]
      else
          let k = 1
              for i in 1:n
                  ptr[i] = k

                  cof[k]  = 1.0
                  subj[k] = z+i-0

                  if x[i] > 0
                      cof[k] = -1.0
                      subj[k] = x[i]
                      k += 1
                  end
              end
              ptr[length(ptr)] = n*2-numbarelm+1
          end
      end
      if ct == MSK_CT_RQUAD cof[1] = 2.0 end
          
      # z in R^n, free
      let b = zeros(Float64,n)
          putvarboundslice(m.task, nvar,nvar+n, Int32[ MSK_BK_FR for i=1:n ], b,b)
          putconboundslice(m.task, ncon,ncon+n, Int32[ MSK_BK_FX for i=1:n ], b,b)
      end

      putarowslice(m.task, ncon, ncon+n, ptr[1:n], ptr[2:n+1], subj, cof)
      appendcone(m.task, ct, 0.0, [z:z+n-1])
      
      # add bar non-zeros
      if numbarelm > 0
          baridxs = find(x -> x < 0, x)
          for (i,k) in enumerate(baridxs)
              j     = -x[k]
              d     = getdimbarvar(m.task,j)
              ii,jj = lintriltoij(m.barvarij[k],d)
              midx  = appendsparsesymmat(m.task,d,Int32[ii],Int32[jj], Float64[1.0])

              putbaraij(m.task, ncon+i-1, j, Int64[midx], Float64[1.0])
          end
      end

      # we add an empty placeholder quadratic constraint. The value is always 0, the dual is undefined.
      appendcons(m.task,1)
      dummycon = getnumcon(m.task)
      addUserQuadCon(m,dummycon)
      putconbound(m.task,dummycon, MSK_BK_FX, 0.0, 0.0)
  else      
    upgradeProbType(m,MosekMathProgModel_QOQP)
    
    qcksubi = m.varmap[quadrowidx]
    qcksubj = m.varmap[quadcolidx]

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

    addUserQuadCon(m,considx)
  end
end



function getquadconstrduals(m::MosekMathProgModel)
    soldef = getsoldef(m)
    if soldef < 0 throw(MosekMathProgModelError("No solution available")) end
    solsta = getsolsta(m.task,soldef)
    if solsta in [MSK_SOL_STA_OPTIMAL, 
                  MSK_SOL_STA_DUAL_FEAS, 
                  MSK_SOL_STA_PRIM_AND_DUAL_FEAS, 
                  MSK_SOL_STA_NEAR_OPTIMAL, 
                  MSK_SOL_STA_NEAR_DUAL_FEAS, 
                  MSK_SOL_STA_NEAR_PRIM_AND_DUAL_FEAS, 
                  MSK_SOL_STA_PRIM_INFEAS_CER]
        let y = gety(m.task,soldef)
            y[m.qconmap]
        end
    else
        throw(MosekMathProgModelError("No solution available"))
    end    
end

function getquadconstrsolution(m::MosekMathProgModel)
    soldef = getsoldef(m)
    if soldef < 0 throw(MosekMathProgModelError("No solution available")) end
    solsta = getsolsta(m.task,soldef)
    if solsta in [MSK_SOL_STA_OPTIMAL, 
                  MSK_SOL_STA_PRIM_FEAS, 
                  MSK_SOL_STA_PRIM_AND_DUAL_FEAS, 
                  MSK_SOL_STA_NEAR_OPTIMAL, 
                  MSK_SOL_STA_NEAR_PRIM_FEAS, 
                  MSK_SOL_STA_NEAR_PRIM_AND_DUAL_FEAS, 
                  MSK_SOL_STA_DUAL_INFEAS_CER]
        let xc = getxc(m.task,soldef)
            xc[m.qconmap]
        end
    else
        throw(MosekMathProgModelError("No solution available"))
    end    
end

