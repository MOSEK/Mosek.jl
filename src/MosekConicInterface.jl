
msk_accepted_cones = [:Free, 
                      :Zero,
                      :NonNeg,
                      :NonPos,
                      :SOC,
                      :SOCRotated,
                      :SDP ]



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
function countcones{Tis}(cones :: Array{@compat(Tuple{Symbol,Tis}),1})
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
            n = @compat(round(Int32,sqrt(.25+2*length(idxs))-0.5))
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
    

    elmidxs = vcat([ convert(Vector{Int},collect(idxs)) for (_,idxs) in cones ]...)
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
    res = Array(Tv,length(a)*n)
    m   = length(a)
    for i in 1:length(res):m
        res[i:i+m-1] = a
    end
    return res
end

function erepeat{Tv}(a :: Array{Tv,1}, n :: Int)
    res = Array(Tv,length(a)*n)
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
    ii = Array(Int32,length(Ls))
    jj = Array(Int32,length(Ls))
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
#  in vs are multiplied by 0.5.
# 
function lintriltoijv(Ls::Array{Int64,1}, vs::Array{Float64,1}, d::Int32)
    if length(Ls) == 0
        Array(Int32,0),Array(Int32,0),Array(Float64,0)
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
                vv[1] = 0.5 * vs[perm[1]]
            end
        end
        k = 1
        for i in 2:length(perm)
            if Ls[perm[i-1]] == Ls[perm[i]]
                if ii[k] == jj[k]
                    vv[k] += vs[perm[i]]
                else
                    vv[k] += 0.5*vs[perm[i]]
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
                        vv[k] = 0.5*vs[perm[i]]
                    end
                end
            end
        end
        ii,jj,vv
    end
end





function loadconicproblem!(m::MosekMathProgModel,c,A,b,constr_cones,var_cones) 
    loadconicproblem!(m,
                      convert(SparseMatrixCSC{Float64,Int},reshape(c,(length(c),1))),
                      convert(SparseMatrixCSC{Float64,Int},A),
                      convert(Array{Float64,1},b),
                      convert(Array{@compat(Tuple{Symbol,Any}),1},constr_cones),
                      convert(Array{@compat(Tuple{Symbol,Any}),1},var_cones))
end


function loadconicproblem!(m::MosekMathProgModel,
                           c::SparseMatrixCSC{Float64,Int},
                           A::SparseMatrixCSC{Float64,Int},
                           b::Array{Float64,1},
                           constr_cones::Array{@compat(Tuple{Symbol,Any}),1},
                           var_cones   ::Array{@compat(Tuple{Symbol,Any}),1})
    # check data
    const N = c.m
    const M = length(b)
    
    if M != A.m || N != A.n || c.n != 1
        throw(MosekMathProgModelError("Invalid data dimensions"))
    end

    (numqcvar,numbarvar,numlinvarelm,numqcvarelm,numbarvarelm,totnumvar) = countcones(var_cones)
    (numqccon,numbarcon,numlinconelm,numqcconelm,numbarconelm,totnumcon) = countcones(constr_cones)

    # clear task data
    putmaxnumvar(m.task,0)
    putmaxnumcon(m.task,0)
    putmaxnumbarvar(m.task,0)
    putcfix(m.task,0.0)

    # allocate necessary variables and cons, reserve space for cones and barvars
    appendvars(m.task,numlinvarelm+numqcvarelm + numqcconelm)
    appendcons(m.task,totnumcon)
    putmaxnumcone(m.task,numqcvar+numqccon)
    putmaxnumbarvar(m.task,numbarvar+numbarcon)

    m.probtype = 
      if     any(v -> first(v) in [:SDP],            var_cones)    MosekMathProgModel_SDP
      elseif any(v -> first(v) in [:SDP],            constr_cones) MosekMathProgModel_SDP
      elseif any(v -> first(v) in [:SOC,:SOCrotated],var_cones)    MosekMathProgModel_SOCP
      elseif any(v -> first(v) in [:SOC,:SOCrotated],constr_cones) MosekMathProgModel_SOCP
      else                                                         MosekMathProgModel_LINR
      end
    
    varmap     = Array(Int32,totnumvar) # nonnegative refer to linear vars, negative to barvars    

    barvarij   = zeros(Int64,totnumvar)
    barvardim  = zeros(Int32,numbarvar+numbarcon)

    linvarptr = 1
    barvarptr = 1
    
    let nvar = numlinvarelm+numqcvarelm
        bk = Array(Int32,nvar)        
        for (sym,idxs_) in var_cones
            idxs = convert(Vector{Int32},collect(idxs_))

            n = length(idxs)
            if sym in [ :Free, :Zero, :NonPos, :NonNeg ]             
                first = linvarptr
                last  = linvarptr+n-1
                linvarptr += n

                varmap[idxs] = first:last

                bk[first:last] =
                   if     sym == :Free   MSK_BK_FR
                   elseif sym == :Zero   MSK_BK_FX
                   elseif sym == :NonNeg MSK_BK_LO
                   elseif sym == :NonPos MSK_BK_UP
                   end            
            elseif sym in [ :SOC, :SOCRotated ]
                first = linvarptr
                last  = linvarptr+n-1
                linvarptr += n

                varmap[idxs] = Int32[first:last;]

                bk[first:last] = MSK_BK_FR
                if     sym == :SOC        appendcone(m.task, MSK_CT_QUAD,  0.0, idxs)
                elseif sym == :SOCRotated appendcone(m.task, MSK_CT_RQUAD, 0.0, idxs)
                end                
            elseif sym == :SDP
                d = @compat(round(Int32,sqrt(.25+2*length(idxs))-0.5))
                trilsz = length(idxs)
                barvardim[barvarptr] = d
                appendbarvars(m.task, Int32[d])
                varmap[idxs] = -barvarptr
                barvarij[idxs] = Int64[1:trilsz;]
                barvarptr += 1
            end
        end
        bx = zeros(Float64,nvar)
        putvarboundslice(m.task,1,nvar+1,bk,bx,bx)
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
            putarowslice(m.task,1,M+1,At)
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
                                    appendsparsesymmat(m.task,barvardim[j], ii,jj,vv)
                                end
                            # NOTE: Since the conic API uses the form (b-Ax < K) we use the weight -1.0.
                            putbaraij(m.task, i,j,Int64[matidx],Float64[-1.0])
                        end
                    end
                end
            end
        end

        # Add bounds and slacks
        let bk = Array(Int32,M)
            local conptr = 1

            for (sym,idxs_) in constr_cones
                idxs = convert(Vector{Int32},collect(idxs_))
                local n = length(idxs)
                if sym in [ :Free, :Zero, :NonPos, :NonNeg ]
                    firstcon = conptr
                    lastcon  = conptr+n-1
                    conptr += n

                    conslack[idxs] = 0 # no slack

                    bk[idxs] =
                      if     sym == :Free   MSK_BK_FR
                      elseif sym == :Zero   MSK_BK_FX
                      elseif sym == :NonNeg MSK_BK_LO
                      elseif sym == :NonPos MSK_BK_UP
                      end
                elseif sym in [ :SOC, :SOCRotated ]
                    firstcon   = conptr
                    lastcon    = conptr+n-1
                    firstslack = linvarptr
                    lastslack  = linvarptr+n-1
                    conptr += n
                    linvarptr += n

                    conslack[idxs] = firstslack:lastslack # no slack
                    bk[idxs] = MSK_BK_FX

                    # Append a variable vector s and make it conic
                    # Then add slacks to the rows: b-Ax - s = 0, s in C
                    local bx = zeros(Float64,n)
                    putvarboundslice(m.task,firstslack,lastslack+1,Int32[MSK_BK_FR for i in 1:n],bx,bx)
                    putaijlist(m.task,Int32[firstcon:lastcon;],Int32[firstslack:lastslack;],-ones(Float64,n))
                    if     sym == :SOC        appendcone(m.task, MSK_CT_QUAD,  0.0, Int32[firstslack:lastslack;])
                    elseif sym == :SOCRotated appendcone(m.task, MSK_CT_RQUAD, 0.0, Int32[firstslack:lastslack;])
                    end
                elseif sym == :SDP
                    firstcon   = conptr
                    lastcon    = conptr+n-1
                    barslackj  = barvarptr
                    d = int32(sqrt(.25+2*length(idxs))-0.5)

                    bk[firstcon:lastcon] = MSK_BK_FX

                    barvardim[barvarptr] = d
                    appendbarvars(m.task, Int32[d])

                    conptr += n
                    barvarptr += 1

                    let i = firstcon
                        for vj in 1:d
                            for vi in vj:d
                                cof = (if vj == vi 1.0 else 0.5 end )
                                const matidx = appendsparsesymmat(m.task,d,Int32[vi],Int32[vj],Float64[cof])
                                putbaraij(m.task,i,barslackj,Int64[matidx],Float64[-1.0])
                                barconij[i] = i-firstcon+1
                                i += 1
                            end
                        end
                    end
                    conmap[firstcon:lastcon] = -barvarptr
                end
            end

            putconboundslice(m.task,1,M+1,bk,-b,-b)
        end
    end


    # Input objective
    let numcnz    = count(j -> varmap[j] > 0,c.rowval),
        numbarcnz = length(c.rowval) - numcnz

        if numbarcnz == 0
            putclist(m.task,varmap[c.rowval],c.nzval)
        else
            lptr   = 1
            barptr = 1

            linidxs = find(j -> varmap[j] > 0, c.rowval)
            baridxs = find(j -> varmap[j] < 0, c.rowval)

            csub    = varmap[c.rowval[linidxs]]
            cval    = c.nzval[linidxs]
            barcsub = -varmap[c.rowval[baridxs]]
            barcij  = barvarij[c.rowval[baridxs]]
            barcval = c.nzval[baridxs]

            putclist(m.task,csub,cval)

            perm = sortperm(barcsub)
            k = 0
            while k < numbarcnz
                k += 1
                let b = k,
                    j = barcsub[perm[k]],
                    d = barvardim[j]
                    k += 1
                    while (k <= numbarcnz && barcsub[perm[k]] == j) k += 1 end

                    ii,jj,vv = lintriltoijv(barcij[perm[b:k-1]],barcval[perm[b:k-1]],d)
                    const matidx = appendsparsesymmat(m.task,d,ii,jj,vv)
                    putbarcj(m.task,j,Int64[matidx],Float64[1.0])
                end
            end
        end
        putobjsense(m.task, MSK_OBJECTIVE_SENSE_MINIMIZE)
    end

    m.probtype =
      if     numbarvar+numbarcon > 0 MosekMathProgModel_SDP
      elseif numqcvar+numqccon   > 0 MosekMathProgModel_SOCP
      else                           MosekMathProgModel_LINR
      end

    m.numvar     = totnumvar # elements used in varmap
    m.varmap     = varmap
    m.barvarij   = barvarij
    m.numbarvar  = numbarvar
    m.barvarmap  = Int32[] # barvars that appear in the linear variable are not accessable through the low level interface
    m.binvarflag = fill(false,m.numvar)

    m.numcon     = totnumcon
    m.conmap     = conmap
    m.conslack   = conslack
    m.barconij   = barconij
end


# this is the same as getconstrduals, except that it accepts cert. of primal infeasibility too
function getconicdual(m::MosekMathProgModel)
    soldef = getsoldef(m)
    if soldef < 0 throw(MosekMathProgModelError("No solution available")) end

    solsta = getsolsta(m.task,soldef)

    if solsta in [MSK_SOL_STA_OPTIMAL,
                  MSK_SOL_STA_DUAL_FEAS,
                  MSK_SOL_STA_PRIM_AND_DUAL_FEAS,
                  MSK_SOL_STA_NEAR_OPTIMAL,
                  MSK_SOL_STA_NEAR_DUAL_FEAS,
                  MSK_SOL_STA_NEAR_PRIM_AND_DUAL_FEAS,
                  # certificate, near certificate
                  MSK_SOL_STA_PRIM_INFEAS_CER,
                  MSK_SOL_STA_NEAR_PRIM_INFEAS_CER ]
        getcondual(m,soldef)
    else
        throw(MosekMathProgModelError("No solution available"))
    end
end
