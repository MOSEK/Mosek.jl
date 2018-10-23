
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

mutable struct MosekNonlinearModel <: MathProgBase.AbstractLinearQuadraticModel
    m :: MosekLinearQuadraticModel
end

function MathProgBase.LinearQuadraticModel(s::Mosek.MosekSolver)
  r = MosekLinearQuadraticModel(Mosek.maketask(),
                                Bool[],
                                0,
                                0,
                                Int32[],
                                Float64[],
                                Float64[],

                                Int32[],
                                Float64[],
                                Float64[],

                                Int32[],
                                Int32[],

                                Mosek.MSK_RES_OK,

                                s.options)
    loadoptions!(r)
    r
end
##############################################################
## Linear
##############################################################

function MathProgBase.loadproblem!(m::MosekLinearQuadraticModel,
                                   A,
                                   collb :: Vector{T1},
                                   colub :: Vector{T2},
                                   obj   :: Vector{T3},
                                   rowlb :: Vector{T4},
                                   rowub :: Vector{T5},
                                   sense:: Symbol) where {T1,T2,T3,T4,T5}
    MathProgBase.loadproblem!(m,
                              convert(SparseMatrixCSC{Float64,Int},A),
                              convert(Vector{Float64},collb),
                              convert(Vector{Float64},colub),
                              convert(Vector{Float64},obj),
                              convert(Vector{Float64},rowlb),
                              convert(Vector{Float64},rowub),
                              sense)
end

function MathProgBase.loadproblem!(m::MosekLinearQuadraticModel,
                                   A::     SparseMatrixCSC{Float64,Int},
                                   collb:: Vector{Float64},
                                   colub:: Vector{Float64},
                                   obj::   Vector{Float64},
                                   rowlb:: Vector{Float64},
                                   rowub:: Vector{Float64},
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
    m.quadcon = Int32[]
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

        lincon = findall(i -> Mosek.getnumqconknz(m.task,i) == 0, 1:numcon)
        quadcon = findall(i -> Mosek.getnumqconknz(m.task,i) > 0, 1:numcon)

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

MathProgBase.setvarLB!(m::MosekLinearQuadraticModel, bnd::Vector{T}) where {T} = MathProgBase.setvarLB!(m,convert(Vector{Float64},bnd))
function MathProgBase.setvarLB!(m::MosekLinearQuadraticModel, bnd::Vector{Float64})
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
        idxs = convert(Vector{Int32},findall(v->v, m.binvarflags))
        bkx = Mosek.Boundkey[ Mosek.MSK_BK_RA for i in 1:length(idxs)]
        blx = Float64[ max(m.blx[i],0.0) for i in idxs ]
        bux = Float64[ min(m.bux[i],1.0) for i in idxs ]

        Mosek.putvarboundlist(m.task,idxs, bkx,blx,bux)
    end

    nothing
end

MathProgBase.setvarUB!(m::MosekLinearQuadraticModel, bnd::Vector{T}) where {T} = MathProgBase.setvarUB!(m,convert(Vector{Float64},bnd))
function MathProgBase.setvarUB!(m::MosekLinearQuadraticModel, bnd::Vector{Float64})
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
        idxs = convert(Vector{Int32},findall(v->v, m.binvarflags))
        bkx = Mosek.Boundkey[ Mosek.MSK_BK_RA for i in 1:length(idxs)]
        blx = Float64[ max(m.blx[i],0.0) for i in idxs ]
        bux = Float64[ min(m.bux[i],1.0) for i in idxs ]

        Mosek.putvarboundlist(m.task,idxs, bkx,blx,bux)
    end

    nothing
end


MathProgBase.setconstrLB!(m::MosekLinearQuadraticModel, bnd::Vector{T}) where {T} = MathProgBase.setconstrLB!(m,convert(Vector{Float64},bnd))
function MathProgBase.setconstrLB!(m::MosekLinearQuadraticModel, bnd::Vector{Float64})
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

MathProgBase.setconstrUB!(m::MosekLinearQuadraticModel, bnd::Vector{T}) where {T} = MathProgBase.setconstrUB!(m,convert(Vector{Float64},bnd))
function MathProgBase.setconstrUB!(m::MosekLinearQuadraticModel, bnd::Vector{Float64})
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

MathProgBase.setobj!(m::MosekLinearQuadraticModel, c::Vector{T}) where {T} = MathProgBase.setobj!(m,convert(Vector{Float64},c))
function MathProgBase.setobj!(m::MosekLinearQuadraticModel, c :: Vector{Float64})
    n = min(length(c),m.numvar)
    Mosek.putclist(m.task,Int32[1:n;],c[1:n])
end


function MathProgBase.getconstrmatrix(m::MosekLinearQuadraticModel)
    (ptrb,ptre,asubj,aval) = Mosek.getaslice(m.task, Mosek.MSK_ACC_CON,1,m.numcon+1)

    numnz = length(asubj)
    asubi = Vector{Int32}(undef,numnz)

    for i in 1:m.numcon
        asubi[ptrb[i]:ptre[i]-1] .= Int32(i)
    end

    sparse(asubi,asubj,aval,m.numcon,m.numvar)
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

MathProgBase.addvar!(m::MosekLinearQuadraticModel,
                     subi::Vector{T1},
                     val ::Vector{T2},
                     bl  ::T3,
                     bu  ::T4,
                     c   ::T5) where {T1,T2,T3,T4,T5} = MathProgBase.addvar!(m,convert(Vector{Int32},subi),convert(Vector{Float64},val),convert(Float64,bl),convert(Float64,bu),convert(Float64,c))
function MathProgBase.addvar!(m::MosekLinearQuadraticModel,
                              subi::Vector{Int32},
                              val ::Vector{Float64},
                              bl  ::Float64,
                              bu  ::Float64,
                              c   ::Float64)
    MathProgBase.addvar!(m,bl,bu,c)
    Mosek.putacol(m.task,m.numvar,subi,val)
end

MathProgBase.addconstr!(m::MosekLinearQuadraticModel,
                        subj::Vector{T1},
                        val ::Vector{T2},
                        bl  ::T3,
                        bu  ::T4) where {T1,T2,T3,T4} = MathProgBase.addconstr!(m,convert(Vector{Int32},subj),convert(Vector{Float64},val),convert(Float64,bl),convert(Float64,bu))

function MathProgBase.addconstr!(m::MosekLinearQuadraticModel,
                                 subj::Vector{Int32},
                                 val ::Vector{Float64},
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
                  Mosek.MSK_SOL_STA_INTEGER_OPTIMAL ]
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
                  Mosek.MSK_SOL_STA_PRIM_AND_DUAL_FEAS ]
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
                  Mosek.MSK_SOL_STA_PRIM_AND_DUAL_FEAS ]
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
    if solsta in [ Mosek.MSK_SOL_STA_DUAL_INFEAS_CER ]
        Mosek.getxx(m.task,soldef)
    else
        throw(MosekMathProgModelError("No ray available"))
    end
end

function MathProgBase.getinfeasibilityray(m::MosekLinearQuadraticModel)
    sol = getsoldef(m.task)

    solsta = Mosek.getsolsta(m.task,sol)
    if solsta in [ Mosek.MSK_SOL_STA_PRIM_INFEAS_CER ]
        Mosek.getslc(m.task,sol) - Mosek.getsuc(m.task,sol)
    else
        throw(MosekMathProgModelError("No ray available"))
    end
end


MathProgBase.getrawsolver(m::MosekLinearQuadraticModel) = m.task

MathProgBase.getsimplexiter(m::MosekLinearQuadraticModel) = Mosek.getintinf(m.task,Mosek.MSK_IINF_SIM_PRIMAL_ITER)+Mosek.getintinf(m.task,Mosek.MSK_IINF_SIM_DUAL_ITER)
MathProgBase.getbarrieriter(m::MosekLinearQuadraticModel) = Mosek.getintinf(m.task,Mosek.MSK_IINF_INTPNT_ITER)

MathProgBase.setwarmstart!(m::MosekLinearQuadraticModel, v::Vector{T}) where {T} = MathProgBase.setwarmstart!(m,convert(Vector{Float64},v))
function MathProgBase.setwarmstart!(m::MosekLinearQuadraticModel, v::Vector{Float64})
    n = min(m.numvar,length(v))
    vals = Vector{Float64}(undef,n)
    vals[:] = v[1:n]

    nanidxs = findall(isnan,vals)
    vals[nanidxs] = 0.0

    skx = Mosek.Stakey[ if isnan(vals[i]) Mosek.MSK_SK_UNK else Mosek.MSK_SK_BAS end for i in 1:n ]

    Mosek.putxxslice(m.task,Mosek.MSK_SOL_BAS,1,n+1,vals);
    Mosek.putskxslice(m.task,Mosek.MSK_SOL_BAS,1,n+1,skx);
end

function MathProgBase.optimize!(m::MosekLinearQuadraticModel)
    try
        m.lasttrm = Mosek.optimize(m.task)
        Mosek.solutionsummary(m.task,Mosek.MSK_STREAM_LOG)
    catch err
        if isa(err, Mosek.MosekError)
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
        for i in findall(vt -> vt == :Bin, vtvec[1:n])
            bl = max(m.blx[i],0.0)
            bu = min(m.bux[i],1.0)
            Mosek.putvarbound(m.task,i,Mosek.MSK_BK_RA,bl,bu)
        end

        # for all :Bin vars being changed to :Int or :Cont, restore original bounds
        for i in findall(i -> (vtvec[i] == :Cont || vtvec[i] == :Int) && m.binvarflags[i], 1:n)
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


MathProgBase.setquadobj!(m::MosekLinearQuadraticModel,subi,subj,valij) = MathProgBase.setquadobj!(m,convert(Vector{Int32},subi),convert(Vector{Int32},subj),convert(Vector{Float64},valij))

# NOTE on data format: The matrix is specified by inputting only lower
# or upper triangular part. This means that whenever element (i,j) is
# added, (j,i) is implicitly added. Duplicates are added together
function MathProgBase.setquadobj!(m::MosekLinearQuadraticModel,
                                  subi  :: Vector{Int32},
                                  subj  :: Vector{Int32},
                                  valij :: Vector{Float64})
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
                                convert(Vector{Int32},subj),
                                convert(Vector{Float64},valj),
                                convert(Vector{Int32},qsubi),
                                convert(Vector{Int32},qsubj),
                                convert(Vector{Float64},qvalij),
                                sense,
                                convert(Float64,bnd))
end

function MathProgBase.addquadconstr!(m      :: MosekLinearQuadraticModel,
                                     subj   :: Vector{Int32},
                                     valj   :: Vector{Float64},
                                     qsubi  :: Vector{Int32},
                                     qsubj  :: Vector{Int32},
                                     qvalij :: Vector{Float64},
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
                  Mosek.MSK_SOL_STA_PRIM_AND_DUAL_FEAS ]
        Mosek.gety(m.task)[m.quadcon]
    else
        throw(MosekMathProgModelError("Solution not available"))
    end
end

function MathProgBase.getquadinfeasibilityray(m::MosekLinearQuadraticModel)
    sol = getsoldef(m)
    solsta = getsolsta(m.task,sol)

    s = Mosek.getsux(m.task,sol) - Mosek.getslx(m.task,sol)
    if solsta in [ Mosek.MSK_SOL_STA_PRIM_INFEAS_CER ]
        -s[m.quadcon]
    else
        throw(MosekMathProgModelError("No solution available"))
    end
end

function MathProgBase.getquadconstrRHS(m::MosekLinearQuadraticModel)
    m.blc[m.quadcon]
end

MathProgBase.setquadconstrRHS!(m::MosekLinearQuadraticModel, bnd) = MathProgBase.setquadconstrRHS!(m,convert(Vector{Float64},bnd))
function MathProgBase.setquadconstrRHS!(m::MosekLinearQuadraticModel, bnd::Vector{Float64})
    n = min(length(bnd),length(m.quadcon))
    m.blc[m.quadcon[1:n]] = bnd[1:n]

    Mosek.putconboundlist(m.task,m.quadcon[1:n],m.bkc,m.blc,m.buc)
end



##############################################################
## Nonlinear hasbeen discontinued
#############################################################
