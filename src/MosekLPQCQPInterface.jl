
type MosekLinearQuadraticModel <: MathProgBase.AbstractLinearQuadraticModel
    task :: Mosek.MSKtask

    binvarflag :: Array{Bool,1}
    # NOTE: bkx/blx/bux are the bound on the variable in the
    # continuous problem. Setting a variable to :Bin, will not change
    # these. Setting :Cont or :Int on a :Bin variable will effectively
    # revert to the original bounds.
    numvar     :: Int
    numcon     :: Int
    bkx        :: Array{Int32,1}
    blx        :: Array{Float64,1}
    bux        :: Array{Float64,1}

    bkc        :: Array{Int32,1}
    blc        :: Array{Float64,1}
    buc        :: Array{Float64,1}
    lincon     :: Array{Int32,1}
    quadcon    :: Array{Int32,1}

    options
end

function makebounds(bl_ :: Array{Float64,1},
                    bu_ :: Array{Float64,1})
    bk = Array(Int32,length(bl_))
    bl = Array(Float64,length(bl_))
    bu = Array(Float64,length(bl_))

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

##############################################################
## Linear
#############################################################

function loadproblem!(m::     MosekMathProgModel,
                      A::     SparseMatrixCSC{Float64,Int},
                      collb:: Array{Float64,1},
                      colub:: Array{Float64,1},
                      obj::   Array{Float64,1},
                      rowlb:: Array{Float64,1},
                      rowub:: Array{Float64,1},
                      sense:: Symbol)
    putmaxnumvar(m.task,0)
    putmaxnumcon(m.task,0)
    putmaxnumcone(m.task,0)
    putmaxnumbarvar(m.task,0)
    putmaxnumqnz(m.task,0)

    nrows,ncols = size(A)
    if ncols != length(collb) ||
        ncols != length(colub) ||
        ncols != size(obj,1)   ||
        nrows != length(rowlb) ||
        nrows != length(rowub) ||
        ncols != length(obj)

        throw(MosekMathProgModelError("Inconsistent data dimensions"))
    end

    appendvars(m.task,ncols)
    appendcons(m.task,nrows)

    m.binvarflag = fill(false,m.numvar)

    (m.bkx,m.blx,m.bux) = makebounds(collb,colub)

    (m.bkc,m.blc,m.buc) = makebounds(rowlb,rowub)

    m.numvar = length(m.bkx)
    m.numcon = length(m.bkc)
    m.lincon = Int32[1:nrows;]
    m.quadcon = Array(Int32,0)

    # input coefficients
    Mosek.putclist(m.task, Int32[1:ncols;], obj)
    Mosek.putacolslice(m.task, 1, ncols+1, A.colptr[1:ncols], A.colptr[2:ncols+1], A.rowval, A.nzval)
    setsense!(m, sense)

    # input bounds
    Mosek.putvarboundslice(m.task, 1, ncols+1, m.bkx, m.blx, m.bux)
    Mosek.putconboundslice(m.task, 1, nrows+1, m.bkc, m.blc, m.buc)

    nothing
end

function loadproblem!(m::        MosekLinearQuadraticModel,
                      filename:: AbstractString)
    tmptask = Mosek.maketask()
    try
        readdata(tmptask, filename)

        if  Mosek.getnumcone(tmptask) > 0 ||
            Mosek.getnumbarvar(tmptask) > 0
            throw(MosekMathProgModelError("Not a linear/quadratic model"))
        end


        numcon = Mosek.getnumcon(tmptask)
        numvar = Mosek.getnumvar(tmptask)

        bkx,blx,bux = Mosek.getvarboundslice(tmptask,1,numvar+1)
        bkc,blc,buc = Mosek.getconboundslice(tmptask,1,numcon+1)

        vts = Mosek.getvartypelist(tmptask,Int32[1:numvar])
        binflags = Bool[ (vts[i] == MSK_VARIABLE_TYPE_INT &&
                          bkx[i] == MSK_BK_RA &&
                          abs(blx[i]) < 1e-8 &&
                          abs(bux[i]-1.0) < 1e-8)
                        for i in 1:numvar ]

        lincon = find(i -> Mosek.getnumqconknz(m.task,i) == 0, 1:numcon)
        quadcon = find(i -> Mosek.getnumqconknz(m.task,i) > 0, 1:numcon)

        Mosek.deletetask(m.task)

        m.task = tmptask
        m.binvaflags = binflags
        m.bkx = bkx
        m.blx = blx
        m.bux = bux
        m.bkc = bkc
        m.blc = blc
        m.buc = buc
        m.numvar = lenght(bkx)
        m.numcon = length(bkc)
        m.lincon = lincon
        m.quadcon = quadcon
    catch
        Mosek.deletetask(tmptask)
        rethrow()
    end
end

function loadoptions!(m::MosekLinearQuadraticModel)
  loadoptions_internal!(m.task, m.options)
end

function writeproblem(m::MosekLinearQuadraticModel, filename::AbstractString)
    Mosek.writedata(m.task,filename)
end

function getvarLB(m::MosekLinearQuadraticModel)
    bk,bu,bl = Mosek.getvarboundslice(m.task,1,m.numvar+1)
    for i in 1:length(bk)
        if bk == MSK_BK_FR || bk == MSK_BK_UP
            bl[i] = -Inf
        end
    end
    bl
end

function getvarUB(m::MosekLinearQuadraticModel)
    bk,bu,bl = Mosek.getvarboundslice(m.task,1,m.numvar+1)
    for i in 1:length(bk)
        if bk == MSK_BK_FR || bk == MSK_BK_LO
            bu[i] = Inf
        end
    end
    bu
end

function getconstrLB(m::MosekLinearQuadraticModel)
    bk,bu,bl = Mosek.getconboundslice(m.task,1,m.numcon+1)
    for i in 1:length(bk)
        if bk == MSK_BK_FR || bk == MSK_BK_UP
            bl[i] = -Inf
        end
    end
    bl
end

function getconstrUB(m::MosekLinearQuadraticModel)
    bk,bu,bl = Mosek.getconboundslice(m.task,1,m.numcon+1)
    for i in 1:length(bk)
        if bk == MSK_BK_FR || bk == MSK_BK_LO
            bu[i] = Inf
        end
    end
    bu
end

function setvarLB!(m::MosekLinearQuadraticModel, bnd::Array{Float64,1})
    n = min(length(bnd),m.numvar)

    for i in 1:n
        if bnd[i] > -Inf
            if m.bux[i] < Inf
                if abs(bnd[i]-m.bux[i]) < 1e-8
                    m.bkx[i] = Mosek.MSK_BK_FX
                    m.blx[i] = m.bux[i]
                    m.bux[i] = m.bux[i]
                else
                    m.bkx[i] = Mosek.MSK_BK_RA
                    m.blx[i] = bnd[i]
                end
            else # bux[i] == Inf
                m.bkx[i] = Mosek.MSK_BK_LO
                m.blx[i] = bnd[i]
            end
        else # bnd[i] == -Inf
            if m.bux[i] < Inf
                m.bkx[i] = Mosek.MSK_BK_UP
            else
                m.bkx[i] = Mosek.MSK_BK_FR
            end
            m.blx[i] = -Inf
        end
    end

    Mosek.putvarboundslice(m.task,1,n+1,m.bkx,m.blx,m.bux)
    if any(m.binvarflags)
        idxs = convert(Array{Int32,1},find(v->v, m.binvarflags))
        bkx = Int32[ Mosek.MSK_BK_RA for i in 1:length(idxs)]
        blx = Float64[ max(m.blx[i],0.0) for i in idxs ]
        bux = Float64[ min(m.bux[i],1.0) for i in idxs ]

        Mosek.putvarboundlist(m.task,idxs, bkx,blx,bux)
    end

    nothing
end

function setvarUB!(m::MosekLinearQuadraticModel, bnd::Array{Float64,1})
    n = min(length(bnd),m.numvar)

    for i in 1:n
        if bnd[i] < Inf
            if m.blx[i] > -Inf
                if abs(bnd[i]-m.blx[i]) < 1e-8
                    m.bkx[i] = Mosek.MSK_BK_FX
                    m.blx[i] = m.blx[i]
                    m.bux[i] = m.blx[i]
                else
                    m.bkx[i] = Mosek.MSK_BK_RA
                    m.bux[i] = bnd[i]
                end
            else # blx[i] == -Inf
                m.bkx[i] = Mosek.MSK_BK_UP
                m.bux[i] = bnd[i]
            end
        else # bnd[i] == Inf
            if m.blx[i] > -Inf
                m.bkx[i] = Mosek.MSK_BK_LO
            else
                m.bkx[i] = Mosek.MSK_BK_FR
            end
            m.blx[i] = -Inf
        end
    end

    Mosek.putvarboundslice(m.task,1,n+1,m.bkx,m.blx,m.bux)
    if any(m.binvarflags)
        idxs = convert(Array{Int32,1},find(v->v, m.binvarflags))
        bkx = Int32[ Mosek.MSK_BK_RA for i in 1:length(idxs)]
        blx = Float64[ max(m.blx[i],0.0) for i in idxs ]
        bux = Float64[ min(m.bux[i],1.0) for i in idxs ]

        Mosek.putvarboundlist(m.task,idxs, bkx,blx,bux)
    end

    nothing
end


function setconstrLB!(m::MosekLinearQuadraticModel, bnd::Array{Float64,1})
    n = min(length(bnd),length(m.lincon))

    for i in 1:n
        if bnd[i] > -Inf
            if m.buc[i] < Inf
                if abs(bnd[i]-m.buc[i]) < 1e-8
                    m.bkc[i] = Mosek.MSK_BK_FX
                    m.blc[i] = m.buc[i]
                    m.buc[i] = m.buc[i]
                else
                    m.bkc[i] = Mosek.MSK_BK_RA
                    m.blc[i] = bnd[i]
                end
            else # buc[i] == Inf
                m.bkc[i] = Mosek.MSK_BK_LO
                m.blc[i] = bnd[i]
            end
        else # bnd[i] == -Inf
            if m.buc[i] < Inf
                m.bkc[i] = Mosek.MSK_BK_UP
            else
                m.bkc[i] = Mosek.MSK_BK_FR
            end
            m.blc[i] = -Inf
        end
    end

    Mosek.putvarboundslice(m.task,1,n+1,m.bkc,m.blc,m.buc)
    if any(m.binvarflags)
        idxs = convert(Array{Int32,1},find(v->v, m.binvarflags))
        bkc = Int32[ Mosek.MSK_BK_RA for i in 1:length(idxs)]
        blc = Float64[ max(m.blc[i],0.0) for i in idxs ]
        buc = Float64[ min(m.buc[i],1.0) for i in idxs ]

        Mosek.putvarboundlist(m.task,idxs, bkc,blc,buc)
    end

    nothing
end

function setconstUB!(m::MosekLinearQuadraticModel, bnd::Array{Float64,1})
    n = min(length(bnd),m.numcon)

    for i in 1:n
        if bnd[i] < Inf
            if m.blc[m.lincon[i]] > -Inf
                if abs(bnd[i]-m.blc[m.lincon[i]]) < 1e-8
                    m.bkc[m.lincon[i]] = Mosek.MSK_BK_FX
                    m.blc[m.lincon[i]] = m.blc[i]
                    m.buc[m.lincon[i]] = m.blc[i]
                else
                    m.bkc[m.lincon[i]] = Mosek.MSK_BK_RA
                    m.buc[m.lincon[i]] = bnd[i]
                end
            else # blc[i] == -Inf
                m.bkc[m.lincon[i]] = Mosek.MSK_BK_UP
                m.buc[m.lincon[i]] = bnd[i]
            end
        else # bnd[i] == Inf
            if m.blc[m.lincon[i]] > -Inf
                m.bkc[m.lincon[i]] = Mosek.MSK_BK_LO
            else
                m.bkc[m.lincon[i]] = Mosek.MSK_BK_FR
            end
            m.blc[m.lincon[i]] = -Inf
        end
    end

    Mosek.putvarboundlist(m.task,m.lincon[1:n],m.bkc,m.blc,m.buc)
    if any(m.binvarflags[m.lincon[1:n]])
        idxs = convert(Array{Int32,1},find(i->m.binvarflags[m.lincon[i]], 1:n))
        bkc = Int32[ Mosek.MSK_BK_RA for i in 1:length(idxs)]
        blc = Float64[ max(m.blc[m.lincon[i]],0.0) for i in idxs ]
        buc = Float64[ min(m.buc[m.lincon[i]],1.0) for i in idxs ]

        Mosek.putvarboundlist(m.task,idxs, bkc,blc,buc)
    end

    nothing
end


function getobj(m::MosekLinearQuadraticModel)
    Mosek.getcslice(m.task,1,m.numvar+1)
end

function setobj(m::MosekLinearQuadraticModel, c :: Array{Float64,1})
    n = min(length(c),m.numvar)
    Mosek.putclist(m.task,Int32[1:n;],c[1:n])
end


function getconstrmatrix(m::MosekLinearQuadraticModel)
    numnz = sum(Int[ Mosek.getarownumnz(m.task,i) for i in 1:m.numcon ])
    asubi = Array(Int32,numnz)
    asubj = Array(Int32,numnz)
    aval  = Array(Float64,numnz)

    let ptr = 1
        for i in 1:m.numcon
            subj,valj = Mosek.getarow(m.task,i)
            n = length(subj)
            asubi[ptr:ptr+n-1] = i
            asubj[ptr:ptr+n-1] = subj
            aval[ptr:ptr+n-1]  = valj
            ptr += n
        end
    end

    sparse(asubi,asubj,aval,length(m.lincon),m.numvar)
end

function addvar!(m::MosekLinearQuadraticModel,
                 bl  ::Float64,
                 bu  ::Float64,
                 c   ::Float64)
    m.numvar += 1
    if bl > -Inf
        if bu < Inf
            if abs(bl-bu) < 1e-8
                push!(m.bkx,Mosek.MSK_BK_FX)
                push!(m.blx,bl)
                push!(m.bux,bl)
            else
                push!(m.bkx,Mosek.MSK_BK_RA)
                push!(m.blx,bl)
                push!(m.bux,bu)
            end
        else
            push!(m.bkx,Mosek.MSK_BK_LO)
            push!(m.blx,bl)
            push!(m.bux,Inf)
        end
    else
        if bu < Inf
            push!(m.bkx,Mosek.MSK_BK_UP)
            push!(m.blx,-Inf)
            push!(m.bux,bu)
        else
            push!(m.bkx,Mosek.MSK_BK_FR)
            push!(m.blx,-Inf)
            push!(m.bux,Inf)
        end
    end

    Mosek.appendvars(m.task,1);
    Mosek.putvarbound(m.task,m.numvar,m.bkx[m.numvar],m.blx[m.numvar],m.bux[m.numvar])
    push!(m.binvarflag,false)
end

function addvar!(m::MosekLinearQuadraticModel,
                 subi::Array{Int32,1},
                 val ::Array{Float64,1},
                 bl  ::Float64,
                 bu  ::Float64,
                 c   ::Float64)
    addvar!(m,bl,bu,c)
    Mosek.putacol(m.task,m.numvar,subi,val)
end


function addconstr!(m::MosekLinearQuadraticModel,
                    subj::Array{Int32,1},
                    val ::Array{Float64,1},
                    bl  ::Float64,
                    bu  ::Float64)
    m.numcon += 1
    push!(m.lincon,m.numcon)
    if bl > -Inf
        if bu < Inf
            if abs(bl-bu) < 1e-8
                push!(m.bkx,Mosek.MSK_BK_FX)
                push!(m.blx,bl)
                push!(m.bux,bl)
            else
                push!(m.bkx,Mosek.MSK_BK_RA)
                push!(m.blx,bl)
                push!(m.bux,bu)
            end
        else
            push!(m.bkx,Mosek.MSK_BK_LO)
            push!(m.blx,bl)
            push!(m.bux,Inf)
        end
    else
        if bu < Inf
            push!(m.bkx,Mosek.MSK_BK_UP)
            push!(m.blx,-Inf)
            push!(m.bux,bu)
        else
            push!(m.bkx,Mosek.MSK_BK_FR)
            push!(m.blx,-Inf)
            push!(m.bux,Inf)
        end
    end
    Mosek.appendcons(m.task,1);
    Mosek.putconbound(m.task,m.numcon,m.bkx[m.numcon],m.blx[m.numcon],m.bux[m.numcon])
    Mosek.putarow(m.task,m.numcon,subj,val)
end

updatemodel!(m::MosekLinearQuadraticModel) = nothing

numlinconstr(m::MosekLinearQuadraticModel) = length(m.lincon)

function getconstrsolution(m::MosekLinearQuadraticModel)
    sol = getsoldef(m.task)
    if sol < 0
        throw(Mosek.MosekMathProgModelError("No solution available"))
    end

    Mosek.getxc(sol)[m.lincon]
end

function getreducedcosts(m::MosekLinearQuadraticModel)
    sol = getsoldef(m.task)
    if sol < 0 || sol == Mosek.MSK_SOL_ITG
        throw(Mosek.MosekMathProgModelError("Solution not available"))
    end
    solsta = Mosek.getsolsta(m.task,sol)

    if solsta in [MSK_SOL_STA_OPTIMAL,
                  MSK_SOL_STA_DUAL_FEAS,
                  MSK_SOL_STA_PRIM_AND_DUAL_FEAS,
                  MSK_SOL_STA_NEAR_OPTIMAL,
                  MSK_SOL_STA_NEAR_DUAL_FEAS,
                  MSK_SOL_STA_NEAR_PRIM_AND_DUAL_FEAS ]
        Mosek.getslx(m.task) - Mosek.getsux(m.task)
    else
        throw(MosekMathProgModelError("No solution available"))
    end
end

function getconstrduals(m::MosekLinearQuadraticModel)
    sol = getsoldef(m.task)
    if sol < 0 || sol == Mosek.MSK_SOL_ITG
        throw(Mosek.MosekMathProgModelError("Solution not available"))
    end
    solsta = Mosek.getsolsta(m.task,sol)

    if solsta in [MSK_SOL_STA_OPTIMAL,
                  MSK_SOL_STA_DUAL_FEAS,
                  MSK_SOL_STA_PRIM_AND_DUAL_FEAS,
                  MSK_SOL_STA_NEAR_OPTIMAL,
                  MSK_SOL_STA_NEAR_DUAL_FEAS,
                  MSK_SOL_STA_NEAR_PRIM_AND_DUAL_FEAS ]
        Mosek.gety(m.task)[m.lincon]
    else
        throw(MosekMathProgModelError("Solution not available"))
    end
end

function getbasis(m::MosekLinearQuadraticModel)
    if ! Mosek.solutiondef(m.task,Mosek.MSK_SOL_BAS)
        throw(Mosek.MosekMathProgModelError("Basis not available"))
    end
    sol = Mosek.MSK_SOL_BAS

    skx = Mosek.getskx(m.task,sol)
    skc = Mosek.getskc(m.task,sol)

    cbasis = [if     skx[i] == MSK_SK_BAS :Basic
              elseif skx[i] == MSK_SK_LO  :NonBasicAtLower
              elseif skx[i] == MSK_SK_UP  :NonBasicAtUpper
              elseif skx[i] == MSK_SK_FX  :NonBasicAtLower # or upper. Doesn't matter.
              else                        :SuperBasic
              for i in 1:m.numvar]
    rbasis = [if     skc[i] == MSK_SK_BAS :Basic
              elseif skc[i] == MSK_SK_LO  :NonBasicAtLower
              elseif skc[i] == MSK_SK_UP  :NonBasicAtUpper
              elseif skc[i] == MSK_SK_FX  :NonBasicAtLower # or upper. Doesn't matter.
              else                        :SuperBasic
              for i in m.lincon]

    cbasis,rbasis
end

function getunboundedray(m::MosekLinearQuadraticModel)
    soldef = getsoldef(m.task)
    if soldef < 0 throw(MosekMathProgModelError("No solution available")) end
    solsta = getsolsta(m.task,soldef)
    if ! solsta in [ MSK_SOL_STA_DUAL_INFEAS_CER, MSK_SOL_STA_NEAR_DUAL_INFEAS_CER ]
        throw(MosekMathProgModelError("No ray available"))
    else
        getxx(m.task,soldef)
    end
end

getrawsolver(m::MosekLinearQuadraticModel) = m.task

getsimplexiter(m::MosekLinearQuadraticModel) = Mosek.getintinf(m.task,Mosek.MSK_IINFITEM_SIM_PRIMAL_ITER)+Mosek.getintinf(m.task,Mosek.MSK_IINFITEM_SIM_DUAL_ITER)+Mosek.getintinf(m.task,Mosek.MSK_IINFITEM_SIM_PRIMAL_DUAL_ITER)

getbarrieriter(m::MosekLinearQuadraticModel) = Mosek.getintinf(m.task,Mosek.MSK_IINFITEM_SIM_INTPNT_ITER)

function setwarmstart!(m::MosekLinearQuadraticModel, v::Array{Float64,1})
    n = min(m.numvar,length(v))
    vals = Array(Float64, n)
    vals[:] = v[1:n]

    nanidxs = find(i -> isnan(vals[i]),vals)
    vals[nanidxs] = 0.0

    skx = Int32[ if isnan(vals[i]) Mosek.MSK_SK_UNK else Mosek.MSK_SK_BAS end for i in 1:n ]

    Mosek.putxxslice(m.task,Mosek.MSK_SOL_BAS,1,n+1,vals);
    Mosek.putskxslice(m.task,Mosek.MSK_SOL_BAS,1,n+1,skx);
end

##############################################################
## Integer Programming
##############################################################

getnodecount(m::MosekLinearQuadraticModel) = 0

##############################################################
## Quadratic
##############################################################

numquadconstr(m::MosekLinearQuadraticModel) = length(m.quadcon)

function setquadobj!{T}(m::MosekLinearQuadraticModel,Q::SparseMatrixCSC{Float64,T})
    n = Q.ptr[length(Q.ptr)]
    qosubi = convert(Array{Int32,1},Q.rowval)
    qosubj = Array(Int32,n)
    for i in 1:m.n
        qosubj[Q.ptr[i]:Q.ptr[i+1]] = i
    end
    qoval = Q.nzval

    setquadobj!(m,qosubi,qosubj,qoval)
end

function setquadobj!(m::MosekLinearQuadraticModel,Q::Array{Float64,2})
    m,n = size(Q)
    qosubi = reshape([i for i in 1:m, j in 1:n], m*n)
    qosubj = reshape([j for i in 1:m, j in 1:n], m*n)
    qoval  = reshape(Q, m*n)

    setquadobj!(m,qosubi,qosubj,qoval)
end

function setquadobj!(m::MosekLinearQuadraticModel,
                     subi  :: Array{Int32,1},
                     subj  :: Array{Int32,1},
                     valij :: Array{Float64,1})
    n = length(subi)
    let qsubi = subi[:],
        qsubj = subj[:],
        qval  = valij[:]
    for i in 1:n
        if qsubi[i] < qsubj[i]
            tmp = qsubi[i]
            qsubi[i] = qsubj[i]
            qsubj[i] = tmp
            qval[i] *= 0.5
        elseif qsubi[i] > qsubj[i]
            qval[i] *= 0.5
        end
    end
    Mosek.putqobj(m.task,subi,subj,valij)
end

function addquadconstr!(m      :: MosekLinearQuadraticModel,
                        subj   :: Array{Int32,1},
                        valj   :: Array{Float64,1},
                        qsubi  :: Array{Int32,1},
                        qsubj  :: Array{Int32,1},
                        qvalij :: Array{Float64,1},
                        sense  :: Char,
                        bnd    :: Float64)

    if sense == '<'
        push!(m.bkc,Mosek.MSK_BK_UP)
    elseif sense = '>'
        push!(m.bkc,Mosek.MSK_BK_LO)
    else
        throw(MosekMathProgSolverInterface.MosekMathProgModelError("Invalid sense"))
    end

    m.numcon += 1
    push!(m.quadcon,m.numcon)
    push!(m.blc,bnd)
    push!(m.buc,bnd)

    Mosek.appendcons(m.task,1)
    Mosek.putconbound(m.task,m.numcon,m.bkc[m.numcon],m.blc[m.numcon],m.buc[m.numcon])

    let qsubi = qsubi[:],
        qsubj = qsubj[:],
        qval  = qvalij[:]
        for i in 1:length(qsubi)
            if qsubi[i] < qsubj[i]
                t = qsubi[i]
                qsubi[i] = qsubj[i]
                qsubj[i] = t
                qval[i] *= 0.5
            elseif qsubi[i] > qsubj[i]
                qval[i] *= 0.5
            end
        end
        Mosek.putqconk(m.task,m.numcon,qsubi,qsubj,qval)
    end

    Mosek.putarow(m.task,m.numcon,subj,valj)
end

function getquadconstrsolution(m::MosekLinearQuadraticModel)
    sol = getsoldef(m.task)
    if sol < 0
        throw(Mosek.MosekMathProgModelError("No solution available"))
    end

    Mosek.getxc(sol)[m.quadcon]
end

function getquadconstrduals(m::MosekLinearQuadraticModel)
    sol = getsoldef(m.task)
    if sol < 0 || sol == Mosek.MSK_SOL_ITG
        throw(Mosek.MosekMathProgModelError("Solution not available"))
    end
    solsta = Mosek.getsolsta(m.task,sol)

    if solsta in [MSK_SOL_STA_OPTIMAL,
                  MSK_SOL_STA_DUAL_FEAS,
                  MSK_SOL_STA_PRIM_AND_DUAL_FEAS,
                  MSK_SOL_STA_NEAR_OPTIMAL,
                  MSK_SOL_STA_NEAR_DUAL_FEAS,
                  MSK_SOL_STA_NEAR_PRIM_AND_DUAL_FEAS ]
        Mosek.gety(m.task)[m.quadcon]
    else
        throw(MosekMathProgModelError("Solution not available"))
    end
end

function getquadinfeasibilityray(m::MosekLinearQuadraticModel)
    sol = getsoldef(m)
    if sol < 0 throw(MosekMathProgModelError("No solution available")) end
    solsta = getsolsta(m.task,sol)

    s = Mosek.getsux(m.task,sol) - Mosek.getslx(m.task,sol)
    if solsta in [ MSK_SOL_STA_PRIM_INFEAS_CER, MSK_SOL_STA_NEAR_PRIM_INFEAS_CER ]
        -s[m.quadcon]
    else
        throw(MosekMathProgModelError("No solution available"))
    end
end

function getquadconstrRHS(m::MosekLinearQuadraticModel)
    m.blc[m.quadcon]
end

function setquadconstrRHS!(m::MosekLinearQuadraticModel, bnd::Array{Float64,1})
    n = min(length(bnd),length(m.quadcon))
    m.blc[m.quadcon[1:n]] = bnd[1:n]

    Mosek.putconboundlist(m.task,m.quadcon[1:n],m.bkc,m.blc,m.buc)
end
