
mutable struct MosekLinearQuadraticModel <: MathProgBase.AbstractLinearQuadraticModel
    task :: Mosek.MSKtask

    binvarflags:: Vector{Bool}
    # NOTE: bkx/blx/bux are the bound on the variable in the
    # continuous problem. Setting a variable to :Bin, will not change
    # these. Setting :Cont or :Int on a :Bin variable will effectively
    # revert to the original bounds.
    numvar     :: Int
    numcon     :: Int
    bkx        :: Vector{Mosek.Boundkey}
    blx        :: Vector{Float64}
    bux        :: Vector{Float64}

    bkc        :: Vector{Mosek.Boundkey}
    blc        :: Vector{Float64}
    buc        :: Vector{Float64}
    lincon     :: Vector{Int32}
    quadcon    :: Vector{Int32}

    lasttrm    :: Mosek.Rescode

    options
end

type MosekNonlinearModel <: MathProgBase.AbstractLinearQuadraticModel
    m :: MosekLinearQuadraticModel
end

function MathProgBase.LinearQuadraticModel(s::MosekSolver) 
  r = MosekLinearQuadraticModel(Mosek.maketask(),
                            Array{Bool}(0),
                            0,
                            0,
                            Array{Int32}(0),
                            Array{Float64}(0),
                            Array{Float64}(0),

                            Array{Int32}(0),
                            Array{Float64}(0),
                            Array{Float64}(0),

                            Array{Int32}(0),
                            Array{Int32}(0),

                            Mosek.MSK_RES_OK,

                            s.options)
    loadoptions!(r)
    r
end
##############################################################
## Linear
##############################################################

function MathProgBase.loadproblem!{T1,T2,T3,T4,T5}(m::MosekLinearQuadraticModel,
                                                   A,
                                                   collb :: Array{T1,1},
                                                   colub :: Array{T2,1},
                                                   obj   :: Array{T3,1},
                                                   rowlb :: Array{T4,1},
                                                   rowub :: Array{T5,1},
                                                   sense:: Symbol)
    MathProgBase.loadproblem!(m,
                              convert(SparseMatrixCSC{Float64,Int},A),
                              convert(Array{Float64,1},collb),
                              convert(Array{Float64,1},colub),
                              convert(Array{Float64,1},obj),
                              convert(Array{Float64,1},rowlb),
                              convert(Array{Float64,1},rowub),
                              sense)
end

function MathProgBase.loadproblem!(m::MosekLinearQuadraticModel,
                                   A::     SparseMatrixCSC{Float64,Int},
                                   collb:: Array{Float64,1},
                                   colub:: Array{Float64,1},
                                   obj::   Array{Float64,1},
                                   rowlb:: Array{Float64,1},
                                   rowub:: Array{Float64,1},
                                   sense:: Symbol)
    Mosek.resizetask(m.task,0,0,0,0,0);

    nrows,ncols = size(A)
    if ncols != length(collb) ||
        ncols != length(colub) ||
        ncols != size(obj,1)   ||
        nrows != length(rowlb) ||
        nrows != length(rowub) ||
        ncols != length(obj)

        throw(MosekMathProgModelError("Inconsistent data dimensions"))
    end

    Mosek.appendvars(m.task,ncols)
    Mosek.appendcons(m.task,nrows)


    (m.bkx,m.blx,m.bux) = makebounds(collb,colub)

    (m.bkc,m.blc,m.buc) = makebounds(rowlb,rowub)

    m.numvar = length(m.bkx)
    m.numcon = length(m.bkc)
    m.lincon = Int32[1:nrows;]
    m.quadcon = Array{Int32}(0)
    m.binvarflags = fill(false,m.numvar)

    # input coefficients
    Mosek.putclist(m.task, Int32[1:ncols;], obj)
    Mosek.putacolslice(m.task, Int32(1), Int32(ncols+1), A.colptr[1:ncols], A.colptr[2:ncols+1], A.rowval, A.nzval)
    MathProgBase.setsense!(m, sense)

    # input bounds
    Mosek.putvarboundslice(m.task, Int32(1), Int32(ncols+1), m.bkx, m.blx, m.bux)
    Mosek.putconboundslice(m.task, Int32(1), Int32(nrows+1), m.bkc, m.blc, m.buc)
    m
end

function MathProgBase.loadproblem!(m::MosekLinearQuadraticModel,
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

        bkx,blx,bux = Mosek.getvarboundslice(tmptask,Int32(1),Int32(numvar+1))
        bkc,blc,buc = Mosek.getconboundslice(tmptask,Int32(1),Int32(numcon+1))

        vts = Mosek.getvartypelist(tmptask,Int32[1:numvar])
        binflags = Bool[ (vts[i] == Mosek.MSK_VARIABLE_TYPE_INT &&
                          bkx[i] == Mosek.MSK_BK_RA &&
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

    m
end

function loadoptions!(m::MosekLinearQuadraticModel)
  loadoptions_internal!(m.task, m.options)
end

function MathProgBase.writeproblem(m::MosekLinearQuadraticModel, filename::AbstractString)
    Mosek.putintparam(m.task,Mosek.MSK_IPAR_OPF_WRITE_SOLUTIONS, Mosek.MSK_ON)
    Mosek.writedata(m.task,filename)
end

function MathProgBase.getvarLB(m::MosekLinearQuadraticModel)
    bk,bl,bu = Mosek.getvarboundslice(m.task,Int32(1),Int32(m.numvar+1))
    for i in 1:length(bk)
        if bk[i] == Mosek.MSK_BK_FR || bk[i] == Mosek.MSK_BK_UP
            bl[i] = -Inf
        end
    end
    bl
end

function MathProgBase.getvarUB(m::MosekLinearQuadraticModel)
    bk,bl,bu = Mosek.getvarboundslice(m.task,Int32(1),Int32(m.numvar+1))
    for i in 1:length(bk)
        if bk[i] == Mosek.MSK_BK_FR || bk[i] == Mosek.MSK_BK_LO
            bu[i] = Inf
        end
    end
    bu
end

function MathProgBase.getconstrLB(m::MosekLinearQuadraticModel)
    bk,bl,bu = Mosek.getconboundslice(m.task,1,m.numcon+1)
    for i in 1:length(bk)
        if bk[i] == Mosek.MSK_BK_FR || bk[i] == Mosek.MSK_BK_UP
            bl[i] = -Inf
        end
    end
    bl
end

function MathProgBase.getconstrUB(m::MosekLinearQuadraticModel)
    bk,bl,bu = Mosek.getconboundslice(m.task,1,m.numcon+1)
    for i in 1:length(bk)
        if bk[i] == Mosek.MSK_BK_FR || bk[i] == Mosek.MSK_BK_LO
            bu[i] = Inf
        end
    end
    bu
end

MathProgBase.setvarLB!{T}(m::MosekLinearQuadraticModel, bnd::Array{T,1}) = MathProgBase.setvarLB!(m,convert(Array{Float64,1},bnd))
function MathProgBase.setvarLB!(m::MosekLinearQuadraticModel, bnd::Array{Float64,1})
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

MathProgBase.setvarUB!{T}(m::MosekLinearQuadraticModel, bnd::Array{T,1}) = MathProgBase.setvarUB!(m,convert(Array{Float64,1},bnd))
function MathProgBase.setvarUB!(m::MosekLinearQuadraticModel, bnd::Array{Float64,1})
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
            m.bux[i] = Inf
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


MathProgBase.setconstrLB!{T}(m::MosekLinearQuadraticModel, bnd::Array{T,1}) = MathProgBase.setconstrLB!(m,convert(Array{Float64,1},bnd))
function MathProgBase.setconstrLB!(m::MosekLinearQuadraticModel, bnd::Array{Float64,1})
    n = min(length(bnd),length(m.lincon))

    for i in 1:n
        bk = m.bkc[m.lincon[i]]
        if bnd[i] > -Inf
            if bk == Mosek.MSK_BK_UP || bk == Mosek.MSK_BK_RA || bk == Mosek.MSK_BK_FX
                if abs(bnd[i]-m.buc[m.lincon[i]]) < 1e-8
                    m.bkc[m.lincon[i]] = Mosek.MSK_BK_FX
                    m.blc[m.lincon[i]] = m.buc[m.lincon[i]]
                    m.buc[m.lincon[i]] = m.buc[m.lincon[i]]
                else
                    m.bkc[m.lincon[i]] = Mosek.MSK_BK_RA
                    m.blc[m.lincon[i]] = bnd[i]
                end
            else # buc[m.lincon[i]] == Inf
                m.bkc[m.lincon[i]] = Mosek.MSK_BK_LO
                m.blc[m.lincon[i]] = bnd[i]
            end
        else # bnd[m.lincon[i]] == -Inf
            if bk == Mosek.MSK_BK_UP || bk == Mosek.MSK_BK_RA || bk == Mosek.MSK_BK_FX
                m.bkc[m.lincon[i]] = Mosek.MSK_BK_UP
            else
                m.bkc[m.lincon[i]] = Mosek.MSK_BK_FR
            end
            m.blc[m.lincon[i]] = -Inf
        end
    end

    Mosek.putconboundslice(m.task,Int32(1),Int32(m.numcon+1),m.bkc,m.blc,m.buc)

    nothing
end

MathProgBase.setconstrUB!{T}(m::MosekLinearQuadraticModel, bnd::Array{T,1}) = MathProgBase.setconstrUB!(m,convert(Array{Float64,1},bnd))
function MathProgBase.setconstrUB!(m::MosekLinearQuadraticModel, bnd::Array{Float64,1})
    n = min(length(bnd),m.numcon)

    for i in 1:n
        bk = m.bkc[m.lincon[i]]
        if bnd[i] < Inf
            if bk == Mosek.MSK_BK_LO || bk == Mosek.MSK_BK_RA || bk == Mosek.MSK_BK_FX
                if abs(bnd[i]-m.blc[m.lincon[i]]) < 1e-8
                    m.bkc[m.lincon[i]] = Mosek.MSK_BK_FX
                    m.buc[m.lincon[i]] = m.blc[m.lincon[i]]
                else
                    m.bkc[m.lincon[i]] = Mosek.MSK_BK_RA
                    m.buc[m.lincon[i]] = bnd[i]
                end
            else # blc[i] == -Inf
                m.bkc[m.lincon[i]] = Mosek.MSK_BK_UP
                m.buc[m.lincon[i]] = bnd[i]
            end
        else # bnd[i] == Inf
            if bk == Mosek.MSK_BK_LO || bk == Mosek.MSK_BK_RA || bk == Mosek.MSK_BK_FX
                m.bkc[m.lincon[i]] = Mosek.MSK_BK_LO
            else
                m.bkc[m.lincon[i]] = Mosek.MSK_BK_FR
            end
            m.buc[m.lincon[i]] = Inf
        end
    end
    Mosek.putconboundslice(m.task,Int32(1),Int32(m.numcon+1),m.bkc,m.blc,m.buc)

    nothing
end


function MathProgBase.getobj(m::MosekLinearQuadraticModel)
    Mosek.getcslice(m.task,1,m.numvar+1)
end

MathProgBase.setobj!{T}(m::MosekLinearQuadraticModel, c::Array{T,1}) = MathProgBase.setobj!(m,convert(Array{Float64,1},c))
function MathProgBase.setobj!(m::MosekLinearQuadraticModel, c :: Array{Float64,1})
    n = min(length(c),m.numvar)
    Mosek.putclist(m.task,Int32[1:n;],c[1:n])
end


function MathProgBase.getconstrmatrix(m::MosekLinearQuadraticModel)
    numnz = sum(Int[ Mosek.getarownumnz(m.task,i) for i in 1:m.numcon ])
    asubi = Array{Int32}(numnz)
    asubj = Array{Int32}(numnz)
    aval  = Array{Float64}(numnz)

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

MathProgBase.addvar!(m::MosekLinearQuadraticModel, bl, bu, c) = MathProgBase.addvar!(m,convert(Float64,bl),convert(Float64,bu),convert(Float64,c))

function MathProgBase.addvar!(m::MosekLinearQuadraticModel,
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
    Mosek.putcj(m.task,m.numvar,c)
    push!(m.binvarflags,false)
end

MathProgBase.addvar!{T1,T2,T3,T4,T5}(m::MosekLinearQuadraticModel,
                                     subi::Array{T1,1},
                                     val ::Array{T2,1},
                                     bl  ::T3,
                                     bu  ::T4,
                                     c   ::T5) = MathProgBase.addvar!(m,convert(Array{Int32,1},subi),convert(Array{Float64,1},val),convert(Float64,bl),convert(Float64,bu),convert(Float64,c))
function MathProgBase.addvar!(m::MosekLinearQuadraticModel,
                              subi::Array{Int32,1},
                              val ::Array{Float64,1},
                              bl  ::Float64,
                              bu  ::Float64,
                              c   ::Float64)
    MathProgBase.addvar!(m,bl,bu,c)
    Mosek.putacol(m.task,m.numvar,subi,val)
end

MathProgBase.addconstr!{T1,T2,T3,T4}(m::MosekLinearQuadraticModel,
                                     subj::Array{T1,1},
                                     val ::Array{T2,1},
                                     bl  ::T3,
                                     bu  ::T4) = MathProgBase.addconstr!(m,convert(Array{Int32,1},subj),convert(Array{Float64,1},val),convert(Float64,bl),convert(Float64,bu))

function MathProgBase.addconstr!(m::MosekLinearQuadraticModel,
                                 subj::Array{Int32,1},
                                 val ::Array{Float64,1},
                                 bl  ::Float64,
                                 bu  ::Float64)
    m.numcon += 1
    push!(m.lincon,m.numcon)
    if bl > -Inf
        if bu < Inf
            if abs(bl-bu) < 1e-8
                push!(m.bkc,Mosek.MSK_BK_FX)
                push!(m.blc,bl)
                push!(m.buc,bl)
            else
                push!(m.bkc,Mosek.MSK_BK_RA)
                push!(m.blc,bl)
                push!(m.buc,bu)
            end
        else
            push!(m.bkc,Mosek.MSK_BK_LO)
            push!(m.blc,bl)
            push!(m.buc,Inf)
        end
    else
        if bu < Inf
            push!(m.bkc,Mosek.MSK_BK_UP)
            push!(m.blc,-Inf)
            push!(m.buc,bu)
        else
            push!(m.bkc,Mosek.MSK_BK_FR)
            push!(m.blc,-Inf)
            push!(m.buc,Inf)
        end
    end
    Mosek.appendcons(m.task,1);
    Mosek.putconbound(m.task,m.numcon,m.bkc[m.numcon],m.blc[m.numcon],m.buc[m.numcon])
    Mosek.putarow(m.task,m.numcon,subj,val)
end

MathProgBase.numlinconstr(m::MosekLinearQuadraticModel) = length(m.lincon)

MathProgBase.getobjval(m::MosekLinearQuadraticModel) = getobjval(m.task)

function MathProgBase.getsolution(m::MosekLinearQuadraticModel)
    sol = getsoldef(m.task)

    solsta = Mosek.getsolsta(m.task,sol)

    if solsta in [Mosek.MSK_SOL_STA_OPTIMAL,
                  Mosek.MSK_SOL_STA_PRIM_FEAS,
                  Mosek.MSK_SOL_STA_PRIM_AND_DUAL_FEAS,
                  Mosek.MSK_SOL_STA_NEAR_OPTIMAL,
                  Mosek.MSK_SOL_STA_NEAR_PRIM_FEAS,
                  Mosek.MSK_SOL_STA_NEAR_PRIM_AND_DUAL_FEAS,
                  Mosek.MSK_SOL_STA_INTEGER_OPTIMAL,
                  Mosek.MSK_SOL_STA_NEAR_INTEGER_OPTIMAL ]
        Mosek.getxx(m.task,sol)
  else
    throw(MosekMathProgModelError("No solution available"))
  end
end

function MathProgBase.getconstrsolution(m::MosekLinearQuadraticModel)
    sol = getsoldef(m.task)

    Mosek.getxc(m.task,sol)[m.lincon]
end

function MathProgBase.getreducedcosts(m::MosekLinearQuadraticModel)
    sol = getsoldef(m.task)
    solsta = Mosek.getsolsta(m.task,sol)

    if solsta in [Mosek.MSK_SOL_STA_OPTIMAL,
                  Mosek.MSK_SOL_STA_DUAL_FEAS,
                  Mosek.MSK_SOL_STA_PRIM_AND_DUAL_FEAS,
                  Mosek.MSK_SOL_STA_NEAR_OPTIMAL,
                  Mosek.MSK_SOL_STA_NEAR_DUAL_FEAS,
                  Mosek.MSK_SOL_STA_NEAR_PRIM_AND_DUAL_FEAS ]
        Mosek.getslx(m.task,sol) - Mosek.getsux(m.task,sol)
    else
        throw(MosekMathProgModelError("No solution available"))
    end
end

function MathProgBase.getconstrduals(m::MosekLinearQuadraticModel)
    sol = getsoldef(m.task)
    if sol == Mosek.MSK_SOL_ITG
        throw(Mosek.MosekMathProgModelError("Solution not available"))
    end
    solsta = Mosek.getsolsta(m.task,sol)

    if solsta in [Mosek.MSK_SOL_STA_OPTIMAL,
                  Mosek.MSK_SOL_STA_DUAL_FEAS,
                  Mosek.MSK_SOL_STA_PRIM_AND_DUAL_FEAS,
                  Mosek.MSK_SOL_STA_NEAR_OPTIMAL,
                  Mosek.MSK_SOL_STA_NEAR_DUAL_FEAS,
                  Mosek.MSK_SOL_STA_NEAR_PRIM_AND_DUAL_FEAS ]
        Mosek.gety(m.task,sol)[m.lincon]
    else
        throw(MosekMathProgModelError("Solution not available"))
    end
end

function MathProgBase.getbasis(m::MosekLinearQuadraticModel)
    if ! Mosek.solutiondef(m.task,Mosek.MSK_SOL_BAS)
        throw(Mosek.MosekMathProgModelError("Basis not available"))
    end
    sol = Mosek.MSK_SOL_BAS

    skx = Mosek.getskx(m.task,sol)
    skc = Mosek.getskc(m.task,sol)

    cbasis = [if     skx[i] == Mosek.MSK_SK_BAS :Basic
              elseif skx[i] == Mosek.MSK_SK_LO  :NonBasicAtLower
              elseif skx[i] == Mosek.MSK_SK_UP  :NonBasicAtUpper
              elseif skx[i] == Mosek.MSK_SK_FX  :NonBasicAtLower # or upper. Doesn't matter.
              else                        :SuperBasic
              end
              for i in 1:m.numvar ]
    rbasis = [if     skc[i] == Mosek.MSK_SK_BAS :Basic
              elseif skc[i] == Mosek.MSK_SK_LO  :NonBasicAtLower
              elseif skc[i] == Mosek.MSK_SK_UP  :NonBasicAtUpper
              elseif skc[i] == Mosek.MSK_SK_FX  :NonBasicAtLower # or upper. Doesn't matter.
              else                        :SuperBasic
              end
              for i in m.lincon]

    cbasis,rbasis
end

function MathProgBase.getunboundedray(m::MosekLinearQuadraticModel)
    soldef = getsoldef(m.task)

    solsta = Mosek.getsolsta(m.task,soldef)
    if solsta in [ Mosek.MSK_SOL_STA_DUAL_INFEAS_CER, Mosek.MSK_SOL_STA_NEAR_DUAL_INFEAS_CER ]
        Mosek.getxx(m.task,soldef)
    else
        throw(MosekMathProgModelError("No ray available"))
    end
end

function MathProgBase.getinfeasibilityray(m::MosekLinearQuadraticModel)
    sol = getsoldef(m.task)

    solsta = Mosek.getsolsta(m.task,sol)
    if solsta in [ Mosek.MSK_SOL_STA_PRIM_INFEAS_CER, Mosek.MSK_SOL_STA_NEAR_PRIM_INFEAS_CER ]
        Mosek.getsux(m.task,sol) - Mosek.getslx(m.task,sol)
    else
        throw(MosekMathProgModelError("No ray available"))
    end
end


MathProgBase.getrawsolver(m::MosekLinearQuadraticModel) = m.task

MathProgBase.getsimplexiter(m::MosekLinearQuadraticModel) = Mosek.getintinf(m.task,Mosek.MSK_IINFITEM_SIM_PRIMAL_ITER)+Mosek.getintinf(m.task,Mosek.MSK_IINFITEM_SIM_DUAL_ITER)+Mosek.getintinf(m.task,Mosek.MSK_IINFITEM_SIM_PRIMAL_DUAL_ITER)

MathProgBase.getbarrieriter(m::MosekLinearQuadraticModel) = Mosek.getintinf(m.task,Mosek.MSK_IINFITEM_SIM_INTPNT_ITER)

MathProgBase.setwarmstart!{T}(m::MosekLinearQuadraticModel, v::Array{T,1}) = MathProgBase.setwarmstart!(m,convert(Array{Float64,1},v))
function MathProgBase.setwarmstart!(m::MosekLinearQuadraticModel, v::Array{Float64,1})
    n = min(m.numvar,length(v))
    vals = Array{Float64}( n)
    vals[:] = v[1:n]

    nanidxs = find(isnan,vals)
    vals[nanidxs] = 0.0

    skx = Int32[ if isnan(vals[i]) Mosek.MSK_SK_UNK else Mosek.MSK_SK_BAS end for i in 1:n ]

    Mosek.putxxslice(m.task,Mosek.MSK_SOL_BAS,1,n+1,vals);
    Mosek.putskxslice(m.task,Mosek.MSK_SOL_BAS,1,n+1,skx);
end

function MathProgBase.optimize!(m::MosekLinearQuadraticModel)
    try
        m.lasttrm = Mosek.optimize(m.task)
        Mosek.solutionsummary(m.task,Mosek.MSK_STREAM_LOG)
    catch err
        if isa(err,MosekError)
            m.lasttrm = err.rcode
        else
            rethrow()
        end
    end
end

MathProgBase.status(m::MosekLinearQuadraticModel) = status(m.task,m.lasttrm)

MathProgBase.getobjbound(m::MosekLinearQuadraticModel) = Mosek.getdouinf(m.task,Mosek.MSK_DINF_MIO_OBJ_BOUND)

MathProgBase.getobjgap(m::MosekLinearQuadraticModel) = getobjgap(m.task)

MathProgBase.getsolvetime(m::MosekLinearQuadraticModel) = Mosek.getdouinf(m.task,Mosek.MSK_DINF_OPTIMIZER_TIME)

MathProgBase.getsense(m::MosekLinearQuadraticModel) = getsense(m.task)

MathProgBase.setsense!(m::MosekLinearQuadraticModel,sense) = setsense!(m.task,sense)

function MathProgBase.freemodel!(m::MosekLinearQuadraticModel)
    Mosek.deletetask(m.task)
    nothing
end

MathProgBase.numvar(m::MosekLinearQuadraticModel) = m.numvar
MathProgBase.numconstr(m::MosekLinearQuadraticModel) = length(m.lincon)

function MathProgBase.setvartype!(m::MosekLinearQuadraticModel,vtvec::Vector{Symbol})
    n = min(m.numvar,length(vtvec))
    if n > 0
        vts = Mosek.Variabletype[if     vt == :Cont Mosek.MSK_VAR_TYPE_CONT
                                 elseif vt == :Int  Mosek.MSK_VAR_TYPE_INT
                                 elseif vt == :Bin  Mosek.MSK_VAR_TYPE_INT
                                 else               Mosek.MSK_VAR_TYPE_CONT
                                 end
                    for vt in vtvec[1:n]]
        Mosek.putvartypelist(m.task,Int32[1:n;],vts)
        for i in find(vt -> vt == :Bin, vtvec[1:n])
            bl = max(m.blx[i],0.0)
            bu = min(m.bux[i],1.0)
            Mosek.putvarbound(m.task,i,Mosek.MSK_BK_RA,bl,bu)
        end

        # for all :Bin vars being changed to :Int or :Cont, restore original bounds
        for i in find(i -> (vtvec[i] == :Cont || vtvec[i] == :Int) && m.binvarflags[i], 1:n)
            Mosek.putvarbound(m.task,i,m.bkx[i],m.blx[i],m.bux[i])
        end

        for i in 1:n
            m.binvarflags[i] = vtvec[i] == :Bin
        end
    end
end

function MathProgBase.getvartype(m::MosekLinearQuadraticModel)
    mskvt = Mosek.getvartypelist(m.task,Int32[1:m.numvar;])
    [if mskvt[i] == Mosek.MSK_VAR_TYPE_INT
         if m.binvarflags[i]
             :Bin
         else
             :Int
         end
     else
         :Cont
     end
     for i in 1:m.numvar]
end


##############################################################
## Integer Programming
##############################################################

MathProgBase.getnodecount(m::MosekLinearQuadraticModel) = 0

##############################################################
## Quadratic
##############################################################

MathProgBase.numquadconstr(m::MosekLinearQuadraticModel) = length(m.quadcon)


MathProgBase.setquadobj!(m::MosekLinearQuadraticModel,subi,subj,valij) = MathProgBase.setquadobj!(m,convert(Array{Int32,1},subi),convert(Array{Int32,1},subj),convert(Array{Float64,1},valij))

# NOTE on data format: The matrix is specified by inputting only lower
# or upper triangular part. This means that whenever element (i,j) is
# added, (j,i) is implicitly added. Duplicates are added together
function MathProgBase.setquadobj!(m::MosekLinearQuadraticModel,
                                  subi  :: Array{Int32,1},
                                  subj  :: Array{Int32,1},
                                  valij :: Array{Float64,1})
    n = length(subi)
    let qsubi = subi[:],
        qsubj = subj[:]
        for i in 1:n
            if qsubi[i] < qsubj[i]
                tmp = qsubi[i]
                qsubi[i] = qsubj[i]
                qsubj[i] = tmp
            end
        end
        Mosek.putqobj(m.task,qsubi,qsubj,valij)
    end
end


function MathProgBase.addquadconstr!(m :: MosekLinearQuadraticModel,
                                     subj,
                                     valj,
                                     qsubi,
                                     qsubj,
                                     qvalij,
                                     sense  :: Char,
                                     bnd)
    MathProgBase.addquadconstr!(m,
                                convert(Array{Int32,1},subj),
                                convert(Array{Float64,1},valj),
                                convert(Array{Int32,1},qsubi),
                                convert(Array{Int32,1},qsubj),
                                convert(Array{Float64,1},qvalij),
                                sense,
                                convert(Float64,bnd))
end

function MathProgBase.addquadconstr!(m      :: MosekLinearQuadraticModel,
                                     subj   :: Array{Int32,1},
                                     valj   :: Array{Float64,1},
                                     qsubi  :: Array{Int32,1},
                                     qsubj  :: Array{Int32,1},
                                     qvalij :: Array{Float64,1},
                                     sense  :: Char,
                                     bnd    :: Float64)

    if     sense == '<'
        push!(m.bkc,Mosek.MSK_BK_UP)
    elseif sense == '>'
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
            elseif  qsubi[i] == qsubj[i]
                qval[i] *= 2
            end
        end
        Mosek.putqconk(m.task,m.numcon,qsubi,qsubj,qval)
    end

    Mosek.putarow(m.task,m.numcon,subj,valj)
end

function MathProgBase.getquadconstrsolution(m::MosekLinearQuadraticModel)
    sol = getsoldef(m.task)

    Mosek.getxc(sol)[m.quadcon]
end

function MathProgBase.getquadconstrduals(m::MosekLinearQuadraticModel)
    sol = getsoldef(m.task)
    solsta = Mosek.getsolsta(m.task,sol)

    if solsta in [Mosek.MSK_SOL_STA_OPTIMAL,
                  Mosek.MSK_SOL_STA_DUAL_FEAS,
                  Mosek.MSK_SOL_STA_PRIM_AND_DUAL_FEAS,
                  Mosek.MSK_SOL_STA_NEAR_OPTIMAL,
                  Mosek.MSK_SOL_STA_NEAR_DUAL_FEAS,
                  Mosek.MSK_SOL_STA_NEAR_PRIM_AND_DUAL_FEAS ]
        Mosek.gety(m.task)[m.quadcon]
    else
        throw(MosekMathProgModelError("Solution not available"))
    end
end

function MathProgBase.getquadinfeasibilityray(m::MosekLinearQuadraticModel)
    sol = getsoldef(m)
    solsta = getsolsta(m.task,sol)

    s = Mosek.getsux(m.task,sol) - Mosek.getslx(m.task,sol)
    if solsta in [ Mosek.MSK_SOL_STA_PRIM_INFEAS_CER, Mosek.MSK_SOL_STA_NEAR_PRIM_INFEAS_CER ]
        -s[m.quadcon]
    else
        throw(MosekMathProgModelError("No solution available"))
    end
end

function MathProgBase.getquadconstrRHS(m::MosekLinearQuadraticModel)
    m.blc[m.quadcon]
end

MathProgBase.setquadconstrRHS!(m::MosekLinearQuadraticModel, bnd) = MathProgBase.setquadconstrRHS!(m,convert(Array{Float64,1},bnd))
function MathProgBase.setquadconstrRHS!(m::MosekLinearQuadraticModel, bnd::Array{Float64,1})
    n = min(length(bnd),length(m.quadcon))
    m.blc[m.quadcon[1:n]] = bnd[1:n]

    Mosek.putconboundlist(m.task,m.quadcon[1:n],m.bkc,m.blc,m.buc)
end



##############################################################
## Nonlinear
#############################################################

MathProgBase.NonlinearModel(s::MosekSolver) = MosekNonlinearModel(MathProgBase.LinearQuadraticModel(s))

type CallbackData
    d::MathProgBase.AbstractNLPEvaluator
    numVar::Int
    numConstr::Int
    Ihess::Vector{Int32}
    Jhess::Vector{Int32}
    jac_colval::Vector{Int32} # Compressed sparse row form of Jacobian sparsity
    jac_rowstarts::Vector{Int32}
    jac_nzval::Vector{Float64} # storage for full jacobian (mosek may ask for a subset)
    jac_nnz_original::Int # nnz in NLP evaluator jacobian == length(Ijac)
    jac_idxmap::Vector{Int} # map from indices in NLP evaluator jacobian to mosek jacobian
    J_tmp::Vector{Float64} # storage for NLP evaluator jacobian
    g_tmp::Vector{Float64} # storage for constraint values
end

function msk_nl_getsp_wrapper_mpb(nlhandle::    Ptr{Void},
                                  numgrdobjnz:: Ptr{Int32}, # number of nonzeros in gradient of objective
                                  grdobjsub::   Ptr{Int32}, # subscripts of nonzeros in gradient of objective
                                  i_::          Int32,      # constraint index
                                  convali::     Ptr{Bool},  # 0/1 indicating whether constraint i is non-linear
                                  grdconinz::   Ptr{Int32}, # number of nonzeros in gradient of constraint i
                                  grdconisub::  Ptr{Int32}, # subscripts of nonzeros in gradient of constraint i
                                  yo::          Int32,      # 0/1 include objective in computation of hessian pattern
                                  numycnz::     Int32,      # number of constraints to include in computation of the hessian pattern
                                  ycsub::       Ptr{Int32}, # indexes of constraints to include in computation of the hessian pattern
                                  maxnumhesnz:: Int32,      # lengths of hessubi and hessubj
                                  numhesnz_::   Ptr{Int32}, # number of hessian nonzeros
                                  hessubi::     Ptr{Int32}, # column subscrips of hessian non-zeros
                                  hessubj::     Ptr{Int32}) # row subscripts of hessian non-zeros
    cb = unsafe_pointer_to_objref(nlhandle)::CallbackData

    i = i_+1

    if numgrdobjnz != C_NULL
        unsafe_store!(numgrdobjnz, convert(Int32,cb.numVar))
    end

    if grdobjsub != C_NULL
        grdobjsub_a = unsafe_wrap(Array{Int32,1},grdobjsub,(cb.numVar,))
        for i in 1:cb.numVar
            grdobjsub_a[i] = i-1
        end
    end

    if i <= cb.numConstr
        con_nnz = cb.jac_rowstarts[i+1] - cb.jac_rowstarts[i]
        if convali != C_NULL
            # can we only say that constraint is linear if we've added the linear part separately?
            if con_nnz > 0
                unsafe_store!(convali, convert(Int32,1))
            else
                unsafe_store!(convali, convert(Int32,0))
            end
        end

        if grdconinz != C_NULL
            unsafe_store!(grdconinz, convert(Int32, con_nnz))
        end

        if grdconisub != C_NULL
            if con_nnz > 0
                grdconisub_a = unsafe_wrap(Array{Int32,1},grdconisub,(con_nnz,))        
                grdconisub_a[1:con_nnz] = cb.jac_colval[cb.jac_rowstarts[i]:(cb.jac_rowstarts[i+1]-1)] - 1
            end
        end
    end

    hess_nnz = length(cb.Ihess)

    if numhesnz_ != C_NULL
        unsafe_store!(numhesnz_, convert(Int32, hess_nnz))
    end

    if hessubi != C_NULL && hessubj != C_NULL && maxnumhesnz >= hess_nnz
        hessubi_a = unsafe_wrap(Array{Int32,1},hessubi,(hess_nnz,))
        hessubj_a = unsafe_wrap(Array{Int32,1},hessubj,(hess_nnz,))

        for i in 1:hess_nnz
            hessubi_a[i] = cb.Ihess[i] - 1
            hessubj_a[i] = cb.Jhess[i] - 1
        end
    end

    return Int32(0)::Int32
end

function msk_nl_getva_wrapper_mpb(nlhandle    :: Ptr{Void},
                                  xx_         :: Ptr{Float64}, # input
                                  yo          :: Float64,
                                  yc_         :: Ptr{Float64}, # input, length = numcon
                                  objval      :: Ptr{Float64},
                                  numgrdobjnz :: Ptr{Int32},
                                  grdobjsub   :: Ptr{Int32},
                                  grdobjval   :: Ptr{Float64},
                                  numi_       :: Int32,
                                  subi_       :: Ptr{Int32},   # input
                                  conval      :: Ptr{Float64},
                                  grdconptrb_ :: Ptr{Int32},   # input
                                  grdconptre_ :: Ptr{Int32},   # input
                                  grdconsub_  :: Ptr{Int32},   # input
                                  grdconval_  :: Ptr{Float64},
                                  grdlag      :: Ptr{Float64},
                                  maxnumhesnz :: Int32,
                                  numhesnz    :: Ptr{Int32},
                                  hessubi     :: Ptr{Int32},
                                  hessubj     :: Ptr{Int32},
                                  hesval      :: Ptr{Float64})
    cb = unsafe_pointer_to_objref(nlhandle)::CallbackData

    numi = convert(Int,numi_)
    xx = unsafe_wrap(Array{Float64,1},xx_,(cb.numVar,))
    yc = unsafe_wrap(Array{Float64,1},yc_,(cb.numConstr,))
    subi = unsafe_wrap(Array{Int32,1},subi_,(numi,))

    if objval != C_NULL
        unsafe_store!(objval, MathProgBase.eval_f(cb.d, xx))
    end

    if numgrdobjnz != C_NULL
        unsafe_store!(numgrdobjnz, convert(Int32,cb.numVar))
    end

    if grdobjsub != C_NULL && grdobjval != C_NULL
        grdobjval_a = unsafe_wrap(Array{Float64,1},grdobjval,(cb.numVar,))
        grdobjsub_a = unsafe_wrap(Array{Int32,1},grdobjsub,(cb.numVar,))

        MathProgBase.eval_grad_f(cb.d, grdobjval_a, xx)
        
        for i in 1:cb.numVar      
            grdobjsub_a[i] = i-1
        end
    end

    if numi > 0 && conval != C_NULL
        conv   = unsafe_wrap(Array{Float64,1},conval,(numi,))
        MathProgBase.eval_g(cb.d, cb.g_tmp, xx)
        for i=1:numi
            conv[i] = cb.g_tmp[subi[i]+1]
        end
    end

    # do we need to compute Jacobian?
    if grdconval_ != C_NULL || grdlag != C_NULL
        MathProgBase.eval_jac_g(cb.d, cb.J_tmp, xx)
        fill!(cb.jac_nzval, 0.0)
        for i in 1:cb.jac_nnz_original
            cb.jac_nzval[cb.jac_idxmap[i]] += cb.J_tmp[i]
        end
    end

    if grdconval_ != C_NULL
        grdconptrb = unsafe_wrap(Array{Int32,1},grdconptrb_,(numi,))
        grdconptre = unsafe_wrap(Array{Int32,1},grdconptre_,(numi,))
        grdconsub = unsafe_wrap(Array{Int32,1},grdconsub_,(grdconptre[numi],))
        grdconval = unsafe_wrap(Array{Float64,1},grdconval_,(grdconptre[numi],))
        nz = 1
        for i=1:numi
            ptrb = grdconptrb[i]
            n    = grdconptre[i] - grdconptrb[i]
            con = subi[i]+1
            # make sure input matches the sparsity pattern we provided
            @assert n == (cb.jac_rowstarts[con+1]-cb.jac_rowstarts[con])
            for k in 1:n
                @assert (grdconsub[grdconptrb[i]+k]+1) ==  cb.jac_colval[cb.jac_rowstarts[con]+k-1]
            end
        end

        # extract subset of Jacobian that mosek asked for
        for i in 1:numi
            con = subi[i]+1
            grdconval[(grdconptrb[i]+1):grdconptre[i]] = view(cb.jac_nzval, cb.jac_rowstarts[con]:(cb.jac_rowstarts[con+1]-1))
        end
    end

    if grdlag != C_NULL
        # could use eval_jac_prod_t here, but just do a sparse matvec instead
        grdlag_a = unsafe_wrap(Array{Float64,1},grdlag,(cb.numVar,))
        MathProgBase.eval_grad_f(cb.d, grdlag_a, xx)
        scale!(grdlag_a, yo)
        Jmat = SparseMatrixCSC(cb.numVar,cb.numConstr,cb.jac_rowstarts,cb.jac_colval,cb.jac_nzval)
        A_mul_B!(-1.0, Jmat, yc, 1.0, grdlag_a)
    end

    nhesnz = length(cb.Ihess)

    if numhesnz != C_NULL
        unsafe_store!(numhesnz,convert(Int32,nhesnz))
    end

    if maxnumhesnz > 0 && hessubi != C_NULL && hessubj != C_NULL && hesval != C_NULL
        hessubi_a = unsafe_wrap(Array{Int32,1},hessubi,(nhesnz,))
        hessubj_a = unsafe_wrap(Array{Int32,1},hessubj,(nhesnz,))

        scale!(yc,-1)
        MathProgBase.eval_hesslag(cb.d,unsafe_wrap(Array{Float64,1},hesval,(nhesnz,)),xx,yo,yc)
        scale!(yc,-1)

        for i=1:length(hessubi_a)
            hessubi_a[i] = cb.Ihess[i]-1
            hessubj_a[i] = cb.Jhess[i]-1
        end
    end
    return convert(Int32,0)::Int32
end

function MathProgBase.loadproblem!(m::MosekNonlinearModel,
                                   numVar,
                                   numConstr,
                                   collb,
                                   colub,
                                   rowlb,
                                   rowub,
                                   sense::Symbol,
                                   d::MathProgBase.AbstractNLPEvaluator)
    MathProgBase.loadproblem!(m,
                              convert(Int32,numVar),
                              convert(Int32,numConstr),
                              convert(Array{Float64,1},collb),
                              convert(Array{Float64,1},colub),
                              convert(Array{Float64,1},rowlb),
                              convert(Array{Float64,1},rowub),
                              sense,d)
end

function MathProgBase.loadproblem!(m::MosekNonlinearModel,
                                   numVar::Int32,
                                   numConstr::Int32,
                                   collb::Array{Float64,1},
                                   colub::Array{Float64,1},
                                   rowlb::Array{Float64,1},
                                   rowub::Array{Float64,1},
                                   sense::Symbol,
                                   d::MathProgBase.AbstractNLPEvaluator)
    m = m.m

    if !(numVar == length(collb) == length(colub)) ||
       !(numConstr == length(rowlb) == length(rowub))
        throw(MosekMathProgModelError("Inconsistent data dimensions"))
    end

    Mosek.deletetask(m.task)
    m.task = Mosek.maketask(Mosek.msk_global_env)
    loadoptions!(m)

    Mosek.appendvars(m.task, numVar)
    Mosek.appendcons(m.task, numConstr)

    m.numvar = numVar
    m.numcon = numConstr

    (m.bkx,m.blx,m.bux) = makebounds(collb,colub)
    (m.bkc,m.blc,m.buc) = makebounds(rowlb,rowub)

    Mosek.putvarboundslice(m.task, 1, numVar+1, m.bkx, m.blx, m.bux)
    Mosek.putconboundslice(m.task, 1, numConstr+1, m.bkc, m.blc, m.buc)

    m.lincon = Int32[1:m.numcon;]
    m.quadcon = Int32[]

    # input bounds

    MathProgBase.setsense!(m, sense)

    MathProgBase.initialize(d, [:Grad, :Jac, :Hess])
    Ijac, Jjac = MathProgBase.jac_structure(d)
    Ihess, Jhess = MathProgBase.hesslag_structure(d)

    # Mosek needs jacobian sparsity in compressed sparse row format, which means that we need to 
    # reformat the indices while maintaining a map
    jac_nnz_original = length(Ijac)
    mergedindices = zeros(Int, jac_nnz_original)
    mergednnz = [0]
    mergedmap = zeros(Int, jac_nnz_original)

    # this is called for duplicate elements
    # assumption: at least one of the indices hasn't been combined before
    function combine(idx1,idx2)
        if mergedmap[idx1] == 0 && mergedmap[idx2] != 0
            mergednnz[1] += 1
            mergedmap[idx1] = idx2
            mergedindices[mergednnz[1]] = idx1
            return idx2
        else
            @assert mergedmap[idx2] == 0
            mergednnz[1] += 1
            mergedmap[idx2] = idx1
            mergedindices[mergednnz[1]] = idx2
            return idx1
        end
    end

    J = sparse(Jjac, Ijac, [i for i in 1:length(Ijac)], numVar, numConstr, combine)

    # map from original index into output index
    idxmap = zeros(Int, jac_nnz_original)
    for row in 1:numConstr
        for pos in J.colptr[row]:(J.colptr[row+1]-1)
            col = J.rowval[pos]
            origidx = J.nzval[pos] # this is the original index of this element
            idxmap[origidx] = pos
        end
    end
    for k in 1:mergednnz[1]
        origidx = mergedindices[k]
        mergedwith = mergedmap[origidx]
        @assert idxmap[origidx] == 0
        @assert idxmap[mergedwith] != 0
        idxmap[origidx] = idxmap[mergedwith]
    end
    for i in 1:jac_nnz_original
        @assert idxmap[i] != 0
    end

    cb = CallbackData(d, numVar, numConstr, Ihess, Jhess, J.rowval, J.colptr, zeros(J.colptr[numConstr+1]-1),
          jac_nnz_original, idxmap, zeros(jac_nnz_original), zeros(numConstr))

    nlgetsp = cfunction(msk_nl_getsp_wrapper_mpb,
                      Int32,
                      (Ptr{Void},Ptr{Int32},Ptr{Int32},Int32,Ptr{Bool},Ptr{Int32},Ptr{Int32},Int32,Int32,Ptr{Int32},Int32,Ptr{Int32},Ptr{Int32},Ptr{Int32}))
    nlgetva = cfunction(msk_nl_getva_wrapper_mpb,
                      Int32,
                      ( Ptr{Void}, # nlhandle
                        Ptr{Float64},Float64,Ptr{Float64}, # xx,yo,yc
                        Ptr{Float64},Ptr{Int32},Ptr{Int32},Ptr{Float64}, # objval,numgrdobjnz,grdobjsub,grdobjval
                        Int32,Ptr{Int32},Ptr{Float64}, # numi,subi,conval
                        Ptr{Int32},Ptr{Int32},Ptr{Int32},Ptr{Float64},Ptr{Float64},# grdconptrb,grdconptre,grdconsub,grdconval,grdlag
                        Int32,Ptr{Int32},Ptr{Int32},Ptr{Int32},Ptr{Float64}))# maxnumhesnz, numhesnz.hessubi,hessubj,hesval
    @Mosek.msk_ccall("putnlfunc",
             Int32, (Ptr{Void},Any,Ptr{Void},Ptr{Void}),
             m.task.task, cb, nlgetsp, nlgetva)
end










MathProgBase.getvarLB(m::MosekNonlinearModel)                            = MathProgBase.getvarLB(m.m)
MathProgBase.getvarUB(m::MosekNonlinearModel)                            = MathProgBase.getvarUB(m.m)
MathProgBase.getconstrLB(m::MosekNonlinearModel)                         = MathProgBase.getconstrLB(m.m)
MathProgBase.getconstrUB(m::MosekNonlinearModel)                         = MathProgBase.getconstrUB(m.m)
MathProgBase.setvarLB!(m::MosekNonlinearModel, bnd::Array{Float64,1})    = MathProgBase.setvarLB!(m.m,bnd)
MathProgBase.setvarUB!(m::MosekNonlinearModel, bnd::Array{Float64,1})    = MathProgBase.setvarUB!(m.m,bnd)
MathProgBase.setconstrLB!(m::MosekNonlinearModel, bnd::Array{Float64,1}) = MathProgBase.setconstrLB!(m.m,bnd)
MathProgBase.setconstrUB!(m::MosekNonlinearModel, bnd::Array{Float64,1}) = MathProgBase.setconstrUB!(m.m,bnd)
#MathProgBase.getobj(m::MosekNonlinearModel)                              = MathProgBase.getobj(m.m)
#MathProgBase.setobj!(m::MosekNonlinearModel, c :: Array{Float64,1})      = MathProgBase.setobj!(m,c)
#MathProgBase.getconstrmatrix(m::MosekNonlinearModel)                     = MathProgBase.getconstrmatrix(m.m)
#MathProgBase.addvar!(m::MosekNonlinearModel,bl::Float64,bu::Float64,c::Float64) = MathProgBase.addvar!(m.m)
#MathProgBase.addvar!(m::MosekNonlinearModel,subi::Array{Int32,1},val::Array{Float64,1},bl::Float64,bu::Float64,c::Float64) = MathProgBase.addvar!(m.m)
#MathProgBase.addconstr!(m::MosekNonlinearModel,subj::Array{Int32,1},val::Array{Float64,1},bl::Float64,bu::Float64) = MathProgBase.addconstr!(m.m)
#MathProgBase.numlinconstr(m::MosekNonlinearModel) = MathProgBase.numlinconstr(m.m)
MathProgBase.getobjval(m::MosekNonlinearModel)                           = MathProgBase.getobjval(m.m)
MathProgBase.getsolution(m::MosekNonlinearModel)                         = MathProgBase.getsolution(m.m)
MathProgBase.getconstrsolution(m::MosekNonlinearModel)                   = MathProgBase.getconstrsolution(m.m)
MathProgBase.getreducedcosts(m::MosekNonlinearModel)                     = MathProgBase.getreducedcosts(m.m)
MathProgBase.getconstrduals(m::MosekNonlinearModel)                      = MathProgBase.getconstrduals(m.m)
#MathProgBase.getbasis(m::MosekNonlinearModel)                            = MathProgBase.getbasis(m.m)
MathProgBase.getunboundedray(m::MosekNonlinearModel)                     = MathProgBase.getunboundedray(m.m)
MathProgBase.getinfeasibilityray(m::MosekNonlinearModel)                 = MathProgBase.getinfeasibilityray(m.m)
MathProgBase.getrawsolver(m::MosekNonlinearModel)                        = MathProgBase.getrawsolver(m.m)
#MathProgBase.getsimplexiter(m::MosekNonlinearModel)                      = MathProgBase.getsimplexiter(m.m)
MathProgBase.getbarrieriter(m::MosekNonlinearModel)                      = MathProgBase.getbarrieriter(m.m)
MathProgBase.setwarmstart!(m::MosekNonlinearModel, v::Array{Float64,1})  = MathProgBase.setwarmstart!(m.m, v)
MathProgBase.optimize!(m::MosekNonlinearModel)                           = MathProgBase.optimize!(m.m)
MathProgBase.status(m::MosekNonlinearModel)                              = MathProgBase.status(m.m)
#MathProgBase.getobjbound(m::MosekNonlinearModel) = MathProgBase.getobjbound(m.m)
#MathProgBase.getobjgap(m::MosekNonlinearModel) = MathProgBase.getobjgap(m.m)
MathProgBase.getsolvetime(m::MosekNonlinearModel)                        = MathProgBase.getsolvetime(m.m)
MathProgBase.getsense(m::MosekNonlinearModel)                            = MathProgBase.getsense(m.m)
MathProgBase.setsense!(m::MosekNonlinearModel,sense)                     = MathProgBase.setsense!(m.m)
MathProgBase.freemodel!(m::MosekNonlinearModel)                          = MathProgBase.freemodel!(m.m)
MathProgBase.numvar(m::MosekNonlinearModel)                              = MathProgBase.numvar(m.m)
MathProgBase.numconstr(m::MosekNonlinearModel)                           = MathProgBase.numconstr(m.m)
#MathProgBase.setvartype!(m::MosekNonlinearModel,vtvec::Vector{Symbol}) = MathProgBase.setvartype!(m,vtec)
#MathProgBase.getvartype(m::MosekNonlinearModel) = MathProgBase.getvartype(m.m)
