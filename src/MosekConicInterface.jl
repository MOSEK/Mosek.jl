const msk_accepted_cones = [:Free,
                            :Zero,
                            :NonNeg,
                            :NonPos,
                            :SOC,
                            :SOCRotated,
                            :SDP ]


MathProgBase.supportedcones(::Mosek.MosekSolver) = msk_accepted_cones



mutable struct MosekMathProgConicModel <: MathProgBase.AbstractConicModel
    task      :: Mosek.MSKtask

    # Length of the variable and constraint vector in the user model
    numvar    :: Int
    numcon    :: Int

    # varmap: Maps UserVarIndex -> MosekVarIndex. UserVarIndex is a
    # positive integer. MosekVarIndex is either positive or negative,
    # where positive integers map to linear variables, and negative
    # map to barvars. When varmap[UserVarIndex] < 0, it refers to a
    # Mosek barvar index, and barvarij refers to the linear element
    # index of the barvar.
    varmap   :: Vector{Int32}
    barvarij :: Vector{Int64}

    # varbk: The boundkeys of variables, index by UserVarIndex. These
    # are necessary when first setting a variable to Binary (thus
    # changing the bounds to MSK_BK_RA), then to Continuous.
    varbk      :: Vector{Mosek.Boundkey}
    # Flags indicating that these variables are binary (integer plus [0] bounded)
    binvarflag :: Vector{Bool}

    # Index of slack variables for conic constraints.
    #   conslack[UserConIndex] = 0 => No slack, constraint is linear
    #   conslack[UserConIndex] > 0 => Slack is linear variable (maps to MosekVarIndex)
    #   conslack[UserConIndex] < 0 => Slack is PSD variable (maps to MosekBarvarIndex).
    #                                 In this case barconij[UserConIndex] maps to the
    #                                 linear index of the element in barvar.
    conslack :: Vector{Int32}
    barconij :: Vector{Int64}
    conbk    :: Vector{Mosek.Boundkey}

    # last termination code, used for status(task)
    lasttrm :: Mosek.Rescode

    # Solver options
    options
end

function MathProgBase.ConicModel(s::Mosek.MosekSolver)
    m = MosekMathProgConicModel(Mosek.maketask(),
                                0,   # numvar
                                0,   # numcon

                                Vector{Int32}(0),  # varmap
                                Vector{Int64}(0),  # barvarij

                                Vector{Int32}(0),  # varbk

                                Vector{Bool}(0),   # binvarflag

                                Vector{Int32}(0),  # conslack
                                Vector{Int64}(0),  # barconij
                                Vector{Int32}(0),  # conbk
                                Mosek.MSK_RES_OK,
                                s.options)
    loadoptions!(m)

    m
end


function loadoptions!(m::MosekMathProgConicModel)
  loadoptions_internal!(m.task, m.options)
end


function MathProgBase.loadproblem!(m::MosekMathProgConicModel,
                                   c,
                                   A,
                                   b,
                                   constr_cones,
                                   var_cones)
    MathProgBase.loadproblem!(m,c,A,b,constr_cones,var_cones,:Min)
end

function MathProgBase.loadproblem!(m::MosekMathProgConicModel,
                                   c,
                                   A,
                                   b,
                                   constr_cones,
                                   var_cones,
                                   sense::Symbol)

    MathProgBase.loadproblem!(m,
                              convert(Array{Float64},c),
                              convert(SparseMatrixCSC{Float64,Int},A),
                              convert(Array{Float64,1},b),
                              convert(Array{Tuple{Symbol,Any},1},constr_cones),
                              convert(Array{Tuple{Symbol,Any},1},var_cones))
end

function MathProgBase.loadproblem!(m::MosekMathProgConicModel,
                                   c::Array{Float64},
                                   A::SparseMatrixCSC{Float64,Int},
                                   b::Array{Float64,1},
                                   constr_cones::Array{Tuple{Symbol,Any},1},
                                   var_cones   ::Array{Tuple{Symbol,Any},1},
                                   sense :: Symbol)
    # check data
    const N = length(c)
    const M = length(b)

    if M != A.m || N != A.n
        throw(MosekMathProgSolverInterface.MosekMathProgModelError("Invalid data dimensions"))
    end

    (numqcvar,numbarvar,numlinvarelm,numqcvarelm,numbarvarelm,totnumvar) = countcones(var_cones)
    (numqccon,numbarcon,numlinconelm,numqcconelm,numbarconelm,totnumcon) = countcones(constr_cones)

    # clear task data
    Mosek.resizetask(m.task,0,0,0,0,0)
    Mosek.putcfix(m.task,0.0)

    # allocate necessary variables and cons, reserve space for cones and barvars
    Mosek.appendvars(m.task,numlinvarelm+numqcvarelm+numqcconelm)
    Mosek.appendcons(m.task,totnumcon)
    Mosek.putmaxnumcone(m.task,numqcvar+numqccon)
    Mosek.putmaxnumbarvar(m.task,numbarvar+numbarcon)

    varmap     = Vector{Int32}(totnumvar) # nonnegative refer to linear vars, negative to barvars

    barvarij   = zeros(Int64,totnumvar)
    barvardim  = zeros(Int32,numbarvar+numbarcon)

    varbk      = Vector{Mosek.Boundkey}(totnumvar)

    linvarptr = 1
    barvarptr = 1

    let nvar = numlinvarelm+numqcvarelm
        varbkidx = 1
        bk = Vector{Mosek.Boundkey}(nvar)
        for (sym,idxs_) in var_cones
            idxs = coneidxstoarray(idxs_)

            n = length(idxs)
            if sym in [ :Free, :Zero, :NonPos, :NonNeg ]
                first = linvarptr
                last  = linvarptr+n-1
                linvarptr += n

                varmap[idxs] = first:last

                bk[first:last] =
                   if     sym == :Free   Mosek.MSK_BK_FR
                   elseif sym == :Zero   Mosek.MSK_BK_FX
                   elseif sym == :NonNeg Mosek.MSK_BK_LO
                   elseif sym == :NonPos Mosek.MSK_BK_UP
                   end

                varbk[idxs] = bk[first:last]
                for i in 1:length(idxs)
                    Mosek.putvarname(m.task,first+i-1,"x$(idxs[i])")
                end
                varbkidx += n
            elseif sym in [ :SOC, :SOCRotated ]
                first = linvarptr
                last  = linvarptr+n-1
                linvarptr += n

                varmap[idxs] = Int32[first:last;]

                bk[first:last] = Mosek.MSK_BK_FR
                if     sym == :SOC        Mosek.appendcone(m.task, Mosek.MSK_CT_QUAD,  0.0, [first:last;])
                elseif sym == :SOCRotated Mosek.appendcone(m.task, Mosek.MSK_CT_RQUAD, 0.0, [first:last;])
                end

                varbk[idxs] = Mosek.MSK_BK_FR
                varbkidx += n
            elseif sym == :SDP
                d = round(Int32,sqrt(.25+2*length(idxs))-0.5)
                trilsz = length(idxs)
                barvardim[barvarptr] = d
                Mosek.appendbarvars(m.task, Int32[d])
                varmap[idxs] = -barvarptr
                barvarij[idxs] = Int64[1:trilsz;]
                barvarptr += 1

                varbk[idxs] = Mosek.MSK_BK_FR
                varbkidx += trilsz
            end
        end
        bx = zeros(Float64,nvar)
        Mosek.putvarboundslice(m.task,Int32(1),Int32(nvar+1),bk,bx,bx)
    end

    conslack = zeros(Int32,M) # 0 means no slack, positive means linear var, negative means semidefinite slack
    barconij = zeros(Int64,M)
    conmap = Int32[1:M;]
    begin # add model constraints
        # Split A into linear and semidefinite columns
        nlinnz = 0
        for ci in 1:length(A.colptr)-1 # count linear nonzeros
            if (varmap[ci] >= 0) nlinnz += A.colptr[ci+1]-A.colptr[ci]  end
        end
        nbarnz = nnz(A)-nlinnz

        asubi = zeros(Int32,nlinnz)
        asubj = zeros(Int32,nlinnz)
        acof  = zeros(Float64,nlinnz)

        barasubi = zeros(Int32,nbarnz)
        barasubj = zeros(Int32,nbarnz)
        barvij   = zeros(Int64,nbarnz)
        baracof  = zeros(Float64,nbarnz)

        let ptr      = 1,
            barptr   = 1
            for ci in 1:A.n
                let n = A.colptr[ci+1]-A.colptr[ci]
                    if varmap[ci] > 0
                        asubi[ptr:ptr+n-1] = A.rowval[A.colptr[ci]:A.colptr[ci+1]-1]
                        asubj[ptr:ptr+n-1] = varmap[ci]
                        acof[ptr:ptr+n-1]  = A.nzval[A.colptr[ci]:A.colptr[ci+1]-1]

                        ptr += n
                    else
                        barasubi[barptr:barptr+n-1] = A.rowval[A.colptr[ci]:A.colptr[ci+1]-1]
                        barasubj[barptr:barptr+n-1] = -varmap[ci]
                        barvij[barptr:barptr+n-1]   = barvarij[ci]
                        baracof[barptr:barptr+n-1]  = A.nzval[A.colptr[ci]:A.colptr[ci+1]-1]

                        barptr += n
                    end
                end
            end
        end

        # Add linear part
        # NOTE: Since the conic API uses the form (b-Ax < K) we use -acof.

        let At = sparse(asubj,asubi,-acof,numlinvarelm+numqcvarelm,M)
            Mosek.putarowslice(m.task,1,M+1,At)
        end

        # Add sdp part
        if nbarnz > 0
            let perm = sortperm(barasubi),
                # NOTE on perm: by default sortperm is stable. Since we
                # barsubi/barsubj to be sorted by barsubj (by construction
                # from column-packed format), perm will be sorted
                # primarily by barasubi and secondarily by barasubj.
                nbarnz = length(barasubi)

                local k    = 1

                while k <= nbarnz
                    let i = barasubi[perm[k]],
                        j = barasubj[perm[k]]

                        let b = k
                            k += 1
                            while k <= nbarnz && barasubi[perm[k]] == i && barasubj[perm[k]] == j
                                k += 1
                            end

                            matidx =
                                let d = barvardim[j]
                                    ii,jj,vv = lintriltoijv(barvij[perm[b:k-1]],baracof[perm[b:k-1]],d)
                                    Mosek.appendsparsesymmat(m.task,barvardim[j], ii,jj,vv)
                                end
                            # NOTE: Since the conic API uses the form (b-Ax < K) we use the weight -1.0.
                            Mosek.putbaraij(m.task, i,j,Int64[matidx],Float64[-1.0])
                        end
                    end
                end
            end
        end

        # Add bounds and slacks
        conbk = Vector{Mosek.Boundkey}(M)
        let bk = conbk
            local conptr = 1

            for (sym,idxs_) in constr_cones
                idxs = coneidxstoarray(idxs_)
                local n = length(idxs)
                if sym in [ :Free, :Zero, :NonPos, :NonNeg ]
                    firstcon = conptr
                    lastcon  = conptr+n-1
                    conptr += n

                    conslack[idxs] = 0 # no slack

                    bk[idxs] =
                      if     sym == :Free   Mosek.MSK_BK_FR
                      elseif sym == :Zero   Mosek.MSK_BK_FX
                      elseif sym == :NonNeg Mosek.MSK_BK_LO
                      elseif sym == :NonPos Mosek.MSK_BK_UP
                      end
                elseif sym in [ :SOC, :SOCRotated ]
                    firstcon   = conptr
                    lastcon    = conptr+n-1
                    firstslack = linvarptr
                    lastslack  = linvarptr+n-1
                    conptr += n
                    linvarptr += n

                    conslack[idxs] = firstslack:lastslack # no slack
                    bk[idxs] = Mosek.MSK_BK_FX

                    # Append a variable vector s and make it conic
                    # Then add slacks to the rows: b-Ax - s = 0, s in C
                    local bx = zeros(Float64,n)
                    Mosek.putvarboundslice(m.task,Int32(firstslack),Int32(lastslack+1),Mosek.Boundkey[Mosek.MSK_BK_FR for i in 1:n],bx,bx)
                    Mosek.putaijlist(m.task,Int32[firstcon:lastcon;],Int32[firstslack:lastslack;],-ones(Float64,n))
                    if     sym == :SOC        Mosek.appendcone(m.task, Mosek.MSK_CT_QUAD,  0.0, Int32[firstslack:lastslack;])
                    elseif sym == :SOCRotated Mosek.appendcone(m.task, Mosek.MSK_CT_RQUAD, 0.0, Int32[firstslack:lastslack;])
                    end
                elseif sym == :SDP
                    firstcon   = conptr
                    lastcon    = conptr+n-1
                    barslackj  = barvarptr
                    d = floor(Int32,sqrt(.25+2*length(idxs))-0.5)

                    bk[firstcon:lastcon] = Mosek.MSK_BK_FX

                    barvardim[barvarptr] = d
                    Mosek.appendbarvars(m.task, Int32[d])

                    let i = firstcon
                        for vj in 1:d
                            for vi in vj:d
                                cof = (vj == vi) ? 1.0 : 1/sqrt(2)
                                const matidx = Mosek.appendsparsesymmat(m.task,d,Int32[vi],Int32[vj],Float64[cof])
                                Mosek.putbaraij(m.task,i,barslackj,Int64[matidx],Float64[-1.0])
                                barconij[i] = i-firstcon+1
                                i += 1
                            end
                        end
                    end
                    conslack[firstcon:lastcon] = -barvarptr

                    conptr += n
                    barvarptr += 1
                end
            end

            Mosek.putconboundslice(m.task,Int32(1),Int32(M+1),bk,-b,-b)
        end
    end


    # Input objective
    let lincidxs = find(j -> varmap[j] > 0 && abs(c[j]) > 1e-8,1:length(c)),
        numcnz   = length(lincidxs),
        barcidxs = find(j -> varmap[j] < 0 && abs(c[j]) > 1e-8,1:length(c)),
        numbarcnz = length(barcidxs)

        if numcnz > 0
            Mosek.putclist(m.task,varmap[lincidxs],c[lincidxs])
        end

        if numbarcnz > 0
            barvardim = Array{Int32}(numbarvar+numbarcon)
            n = numbarvar+numbarcon
            barptr = zeros(Int,n+1)
            for i in barcidxs
                barptr[1-varmap[i]] += 1
            end
            for i in 1:n
                barptr[i+1] += barptr[i]
            end

            barcsubi = Array{Int32}(numbarcnz)
            barcsubj = Array{Int32}(numbarcnz)
            barcval  = Array{Float64}(numbarcnz)
            for i in barcidxs
                j = -varmap[i]
                L = Mosek.getdimbarvarj(m.task,j)
                barvardim[j] = L
                ii,jj = lintriltoij(barvarij[i],L)
                barcsubi[barptr[j]+1] = ii
                barcsubj[barptr[j]+1] = jj
                if ii != jj
                    barcval[barptr[j]+1]  = c[i]/sqrt(2)
                else
                    barcval[barptr[j]+1]  = c[i]
                end
                barptr[j] += 1
            end

            for i in length(barptr)-1:-1:1
                barptr[i+1] = barptr[i]
            end
            barptr[1] = 0

            for j in 1:length(barptr)-1
                if barptr[j] < barptr[j+1]
                    pb = barptr[j]+1
                    pe = barptr[j+1]
                    d = barvardim[j]

                    const matidx = Mosek.appendsparsesymmat(m.task,
                                                            d,
                                                            barcsubi[pb:pe],
                                                            barcsubj[pb:pe],
                                                            barcval[pb:pe])
                    Mosek.putbarcj(m.task,j,Int64[matidx],Float64[1.0])
                end
            end
        end
        setsense!(m.task, sense)
    end


    m.varbk      = varbk
    m.numvar     = totnumvar # elements used in varmap
    m.varmap     = varmap
    m.barvarij   = barvarij
    m.binvarflag = fill(false,m.numvar)

    m.numcon     = totnumcon
    m.conslack   = conslack
    m.barconij   = barconij
    m.conbk      = conbk
end


function MathProgBase.setbvec!(m::MosekMathProgConicModel, b::Array{Float64,1})
    if length(b) != m.numcon
        throw(MosekMathProgSolverInterface.MosekMathProgModelError("Invalid b vector dimension"))
    end

    Mosek.putconboundslice(m.task,Int32(1),Int32(length(b)+1),m.conbk,-b,-b)
end

function MathProgBase.setbvec!(m::MosekMathProgConicModel, b)
    MathProgBase.setbvec!(m,collect(Float64,b))
end


function MathProgBase.writeproblem(m::MosekMathProgConicModel, filename::AbstractString)
    Mosek.writedata(m.task,filename)
end


function MathProgBase.getsolution(m::MosekMathProgConicModel)
    sol = getsoldef(m.task)
    xx = Mosek.getxx(m.task,sol)
    barx = [ Mosek.getbarxj(m.task,sol,j) for j in 1:Mosek.getnumbarvar(m.task) ]
    # rescale primal solution to svec form
    for j in 1:Mosek.getnumbarvar(m.task)
        L = Mosek.getdimbarvarj(m.task,j)
        r = 0
        for k in 1:L
            for i in k:L
                r += 1
                if i != k
                    barx[j][r] *= sqrt(2)
                end
            end
        end
    end

    Float64[ if (m.varmap[i] > 0) xx[m.varmap[i]] else barx[-m.varmap[i]][m.barvarij[i]] end
            for i in 1:m.numvar]
end

function MathProgBase.getvardual(m::MosekMathProgConicModel)
    sol = getsoldef(m.task)
    solsta = Mosek.getsolsta(m.task,sol)

    if sol == Mosek.MSK_SOL_BAS
        s = Mosek.getslx(m.task,sol) - Mosek.getsux(m.task,sol)

        Float64[s[m.varmap[i]] for i in 1:m.numvar]
    else
        s = Mosek.getslx(m.task,sol) - Mosek.getsux(m.task,sol) + Mosek.getsnx(m.task,sol)
        bars = [ Mosek.getbarsj(m.task,sol,j) for j in 1:Mosek.getnumbarvar(m.task) ]

        # rescale dual solution to svec form
        for j in 1:Mosek.getnumbarvar(m.task)
            L = Mosek.getdimbarvarj(m.task,j)
            r = 0
            for k in 1:L
                for i in k:L
                    r += 1
                    if i != k
                        bars[j][r] *= sqrt(2)
                    end
                end
            end
        end

        Float64[if (m.varmap[i] > 0) s[m.varmap[i]] else bars[-m.varmap[i]][m.barvarij[i]] end
                for i in 1:m.numvar]
    end

end


function getconstrsolution_internal(m::MosekMathProgConicModel)
    sol = getsoldef(m.task)

    xc = Mosek.getxc(m.task,sol)
    xx = Mosek.getxx(m.task,sol)
    barx = [ Mosek.getbarxj(j) for j in 1:Mosek.getnumbarvar(m.task) ]

    Float64[if     m.conslack[i] == 0 xc[i]
            elseif m.conslack[i] >  0 xx[m.conslack[i]]
            else                      barx[-m.conslack[i]][m.barconij[i]]
            end
            for i in 1:m.numcon]
end

function getvarduals_internal(m::MosekMathProgConicModel)
    sol = getsoldef(m.task)
    if sol == Mosek.MSK_SOL_ITG
        throw(Mosek.MosekMathProgModelError("No dual solution information available"))
    end

    s = Mosek.getslx(sol) - Mosek.getsux(sol) + Mosek.getsnx(sol)
    bars = [ Mosek.getbarsj(m.task,sol,j) for j in 1:Mosek.getnumbarvar(m.task) ]

    Float64[ if (m.varmap[i] > 0) s[m.varmap[i]] else bars[-m.varmap[i]][m.barvarij[i]] end
             for i in 1:m.numvar ]
end


function MathProgBase.getdual(m::MosekMathProgConicModel)
    sol = getsoldef(m.task)
    if sol == Mosek.MSK_SOL_ITG
        throw(Mosek.MosekMathProgModelError("No dual solution information available"))
    end

    y    = Mosek.gety(m.task,sol)
    snx  = (if sol == Mosek.MSK_SOL_ITR
              Mosek.getsnx(m.task,sol)
            else
              zeros(Float64,m.numvar)
            end)

    bars = [ Mosek.getbarsj(m.task,sol,j) for j in 1:Mosek.getnumbarvar(m.task) ]
    # rescale dual solution to svec form
    for j in 1:Mosek.getnumbarvar(m.task)
        L = Mosek.getdimbarvarj(m.task,j)
        r = 0
        for k in 1:L
            for i in k:L
                r += 1
                if i != k
                    bars[j][r] *= sqrt(2)
                end
            end
        end
    end

    Float64[if     m.conslack[i] == 0 y[i]
            elseif m.conslack[i] >  0 snx[m.conslack[i]]
            else                      bars[-m.conslack[i]][m.barconij[i]]
            end
            for i in 1:m.numcon ]
end





MathProgBase.getobjval(m::MosekMathProgConicModel) = getobjval(m.task)

function MathProgBase.optimize!(m::MosekMathProgConicModel)
    try
        m.lasttrm = Mosek.optimize(m.task)
        Mosek.solutionsummary(m.task,Mosek.MSK_STREAM_LOG)
    catch err
        m.lasttrm = err.rcode
    end
end

MathProgBase.status(m::MosekMathProgConicModel) =
begin
    status(m.task,m.lasttrm)
end

MathProgBase.setsense!(m::MosekMathProgConicModel, sense) = setsense!(m.task,sense)

MathProgBase.getobjbound(m::MosekMathProgConicModel) = Mosek.getdouinf(m.task,Mosek.MSK_DINF_MIO_OBJ_BOUND)

MathProgBase.getobjgap(m::MosekMathProgConicModel) = getobjgap(m::MosekMathProgConicModel)

MathProgBase.getsolvetime(m::MosekMathProgConicModel) = Mosek.getdouinf(m.task,Mosek.MSK_DINF_OPTIMIZER_TIME)

MathProgBase.getrawsolver(m::MosekMathProgConicModel) = m.task

MathProgBase.getsense(m::MosekMathProgConicModel) = getsense(m.task)
MathProgBase.numvar(m::MosekMathProgConicModel)    = m.numvar
MathProgBase.numconstr(m::MosekMathProgConicModel) = m.numcon

function MathProgBase.freemodel!(m::MosekMathProgConicModel)
    Mosek.deletetask(m.task)
    nothing
end


# NOTE: We simply disregard any integer SDP vars.
# NOTE: When setting :Bin, we *change* the domain of the variable to [0;1], irregardless what it was before.
function MathProgBase.setvartype!(m::MosekMathProgConicModel, intvarflag::Vector{Symbol})
    n = min(length(intvarflag),m.numvar)
    if n > 0
        all(x->in(x,[:Cont,:Int,:Bin]), intvarflag) || error("Invalid variable type present")

        idxs        = find(i -> m.varmap[i] > 0,1:n) # indexes into intvarflag for non-PSD vars

        newbk = Mosek.Boundkey[ if (intvarflag[i] == :Bin) Mosek.MSK_BK_RA else m.varbk[i] end for i in idxs ]
        newbl = Float64[0.0 for i in idxs ]
        newbu = Float64[if (intvarflag[i] == :Bin) 1.0 else 0.0 end for i in idxs ]
        newvt = Mosek.Variabletype[if (c == :Cont) Mosek.MSK_VAR_TYPE_CONT else Mosek.MSK_VAR_TYPE_INT end for c in intvarflag ]

        Mosek.putvartypelist(m.task,m.varmap[idxs],newvt)
        Mosek.putvarboundlist(m.task,m.varmap[idxs],newbk,newbl,newbu)

        m.binvarflag[idxs] = map(i -> intvarflag[i] == :Bin, idxs)
    end
end

function MathProgBase.getvartype(m::MosekMathProgConicModel)
    vartypes = [ if isbin :Bin else :Cont end for isbin in m.binvarflag[m.numvar]]
    idxs = find(i -> m.varmap[i] > 0 && vartype[i] == :Cont, 1:m.numvar)
    mskvartypes = Mosek.getvartypelist(m.task,m.varmap[idxs])
    intvaridxs = find(i -> mskvartypes[i] == Mosek.MSK_VAR_TYPE_INT, 1:length(mskvartypes))
    vartypes[intvaridxs] = :Int

    vartypes
end



# countcones :: Array{(Symbol,Tis),1} -> (Int,Int,Int,Int,Int,Int)
#
# Count number of elements in the cone product.
#
# Returns (numqcone,numsdpcone,numlin,numqconeelm,numsdpconeelm,vecsize)
# numqcone
#   Number of quadratic cones
# numsdpcone
#   Number of SDP cones
# numlin
#   Number of linear scalar elements
# numqconeelm
#   Total number of quadratic cone scalar elements
# numsdpconeelm
#   Total number of PSD cone scalar elements
# vecsize
#   Total number of scalar element (= numlin+numqconeelm+numsdpconeelm)
#
coneidxstoarray(idxs :: Int32) = Int[ convert(Int,idxs) ]
coneidxstoarray(idxs :: Int64) = Int[ convert(Int,idxs) ]
coneidxstoarray(idxs :: Array{Int,1}) = idxs
coneidxstoarray(idxs) = collect(Int,idxs)

function countcones{Tis}(cones :: Array{Tuple{Symbol,Tis},1})
    numlin        = 0 # linear and conic quadratic elements
    numsdpcone    = 0 # number of sdp cones
    numsdpconeelm = 0 # total number of elements in all sdp cones
    numqcone      = 0 # number of quadratic cones
    numqconeelm   = 0 # number of quadratic cone elements
    vecsize       = 0 # total number of elements (linear, conic and SDP)

    for (sym,idxs) in cones
        if ! (sym in msk_accepted_cones)
            throw(MosekMathProgModelError("Unsupported cone type"))
        end
        vecsize += length(idxs)

        if     sym == :SDP
            n = round(Int32,sqrt(.25+2*length(idxs))-0.5)
            if n*(n+1)/2 != size(idxs,1) # does not define the lower triangular part of a square matrix
                throw(MosekMathProgModelError("Invalid SDP cone definition"))
            end

            numsdpcone += 1
            numsdpconeelm += length(idxs)
        elseif sym in [ :SOC, :SOCRotated ]
            numqcone += 1
            numqconeelm += length(idxs)
        else
            numlin   += length(idxs)
        end
    end

    elmidxs = vcat([ coneidxstoarray(idxs) for (_,idxs) in cones ]...)
    sort!(elmidxs)

    # check for duplicated and missing elements
    for i in 2:vecsize
        if     elmidxs[i-1] == elmidxs[i]
            throw(MosekMathProgModelError("Invalid data: Intersecting cones"))
        elseif elmidxs[i-1] < elmidxs[i]-1
            throw(MosekMathProgModelError("Invalid data: Missing element in cone specification"))
        end
    end

    return (numqcone,numsdpcone,numlin,numqconeelm,numsdpconeelm,vecsize)
end


function arepeat{Tv}(a :: Array{Tv,1}, n :: Int)
    res = Array{Tv}(length(a)*n)
    m   = length(a)
    for i in 1:length(res):m
        res[i:i+m-1] = a
    end
    return res
end

function erepeat{Tv}(a :: Array{Tv,1}, n :: Int)
    res = Array{Tv}(length(a)*n)
    m   = length(a)
    for i in 0:length(a)-1
        res[i*n+1:(i+1)*n] = a[i]
    end
    return res
end


#internal
# Map linear index into column oriented lower triangular part of a
# square matrix to an (i,j) row,column index. It feels like there
# should be a closed term for computing i,j from L, but... :'(
function lintriltoij(L::Int64, n::Int32)
    let L = L-1
        local j = 0
        while L >= n-j
            L -= (n-j)
            j += 1
        end
        i = L+j
        (i+1,j+1)
    end
end

#internal
function lintriltoij(Ls::Array{Int64,1}, d::Int32)
    n = length(Ls)
    ii = Array{Int32}(length(Ls))
    jj = Array{Int32}(length(Ls))
    for (k,L) in enumerate(Ls)
        i,j = lintriltoij(L,d)
        ii[k] = i
        jj[k] = j
    end
    ii,jj
end

#internal
ijtolintril(i::Int32, j::Int32, d::Int32) =
    let i = int64(i-1),
        j = int64(j-1)
        if (i < j)
            (i*(2*d-i-1) >> 1)+j+1
        else
            (j*(2*d-j-1) >> 1)+i+1
        end
    end
#internal
ijtolintril(ii::Array{Int32,1}, jj::Array{Int32,1}, n::Int32) =
    map(i,j -> ijtolintril(i,j,n),ii,jj)

#internal
#  Parameters:
#
#  * d dimension of the matrix
#  * Ls List of linear indexes into matrix (inplicitly defines subi,subj)
#  * vs List of coefficient values
#
#  Toghether (Ls,vs) define subscripts and coefficients of the *full*
#  matrix. We convert this and return the lower triangular only on
#  (i,j,v)-form. Note that this means that all off-diagonal elements
#  in vs are multiplied by sqrt(2).
#
function lintriltoijv(Ls::Array{Int64,1}, vs::Array{Float64,1}, d::Int32)
    if length(Ls) == 0
        Array{Int32}(0),Array{Int32}(0),Array{Float64}(0)
    else
        perm = sortperm(Ls)
        # count unique
        nunique = 1
        for i in 2:length(perm)
            if Ls[perm[i-1]] < Ls[perm[i]]
                nunique += 1
            end
        end
        ii = zeros(Int32,  nunique)
        jj = zeros(Int32,  nunique)
        vv = zeros(Float64,nunique)

        let # NOTE: See https://github.com/JuliaLang/julia/issues/9134
            i,j = lintriltoij(Ls[perm[1]],d)
            ii[1] = i
            jj[1] = j
            if i == j
                vv[1] = vs[perm[1]]
            else
                vv[1] = vs[perm[1]]/sqrt(2)
            end
        end
        k = 1
        for i in 2:length(perm)
            if Ls[perm[i-1]] == Ls[perm[i]]
                if ii[k] == jj[k]
                    vv[k] += vs[perm[i]]
                else
                    vv[k] += vs[perm[i]]/sqrt(2)
                end
            else
                k += 1
                let # NOTE: See https://github.com/JuliaLang/julia/issues/9134
                    vi,vj = lintriltoij(Ls[perm[i]],d)
                    ii[k] = vi
                    jj[k] = vj
                    if vi == vj
                        vv[k] = vs[perm[i]]
                    else
                        vv[k] = vs[perm[i]]/sqrt(2)
                    end
                end
            end
        end
        ii,jj,vv
    end
end
