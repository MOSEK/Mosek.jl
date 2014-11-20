function sparseToSparseTriple(mat::SparseMatrixCSC)
    if issym(mat)
        nnz = convert(Int64, (countnz(mat)+countnz(diag(mat))) / 2)
    elseif istriu(mat)
        nnz = nfilled(mat)
    else
        error("Matrix must be symmetric or upper triangular")
    end
    II = Array(Cint, nnz)
    JJ = Array(Cint, nnz)
    VV = Array(Cdouble, nnz)
    m, n = size(mat)
    k = 0
    colptr::Vector{Int64} = mat.colptr
    nzval::Vector{Float64} = mat.nzval

    for i = 1:n
        qi = convert(Cint, i)
        for j = colptr[i]:(colptr[i+1]-1)
            qj = convert(Cint, mat.rowval[j])
            if qi <= qj && nzval[j] != 0.0
                k += 1
                II[k] = qj
                JJ[k] = qi
                VV[k] = nzval[j]
            end
        end
    end
    return II,JJ,VV
end

function denseToSparseTriple(mat::Matrix)
    if issym(mat)
        nnz = convert(Int64, (countnz(mat)+countnz(diag(mat))) / 2)
        II = Array(Int32, nnz)
        JJ = Array(Int32, nnz)
        VV = Array(Float64, nnz)
        m, n = size(mat)
        cnt = 1
        for j in 1:m # get LOWER TRIANGULAR
            for i in j:n
                if mat[i,j] != 0.0
                    II[cnt] = i
                    JJ[cnt] = j
                    VV[cnt] = mat[i,j]
                    cnt += 1
                end
            end
        end
    elseif istriu(mat)
        nnz = nfilled(mat)
        II = Array(Int32, nnz)
        JJ = Array(Int32, nnz)
        VV = Array(Float64, nnz)
        m, n = size(mat)
        cnt = 1
        for i in 1:m # UPPER TRIANGULAR -> LOWER TRIANGULAR
            for j in i:n
                if mat[i,j] != 0.0
                    II[cnt] = j
                    JJ[cnt] = i
                    VV[cnt] = mat[i,j]
                    cnt += 1
                end
            end
        end
    # elseif istril(mat)
    #     nnz = nfilled(mat)
    #     II = Array(Int32, nnz)
    #     JJ = Array(Int32, nnz)
    #     VV = Array(Float64, nnz)
    #     m, n = size(mat)
    #     cnt = 1
    #     for j in 1:m # get LOWER TRIANGULAR
    #         for i in j:n
    #             if mat[i,j] != 0.0
    #                 II[cnt] = i
    #                 JJ[cnt] = j
    #                 VV[cnt] = mat[i,j]
    #                 cnt += 1
    #             end
    #         end
    #     end
    else
        error("Matrix must be symmetric or upper triangular")
    end
    return II,JJ,VV
end

function addsdpvar!(m::MosekMathProgModel, dim)   
    m.probtype = MosekMathProgModel_SDP
    appendbarvars(m.task, Cint[dim])
    const barvaridx = getnumbarvar(m.task)
    addUserBarvar(m,barvaridx)
    return int64(barvaridx)
end

function addsdpmatrix!(m::MosekMathProgModel, mat)
    if isa(mat, Matrix)
        II,JJ,VV = denseToSparseTriple(mat)
    elseif isa(mat, SparseMatrixCSC)
        II,JJ,VV = sparseToSparseTriple(mat)
    else
        II,JJ,VV = sparseToSparseTriple(sparse(mat))
    end
    idx = Mosek.appendsparsesymmat(m.task, size(mat,1), II, JJ, VV)
    return int64(idx)
end

# 
function addsdpconstr!(m::MosekMathProgModel, matvaridx, matcoefidx, scalidx, scalcoef, lb, ub)
    m.probtype = MosekMathProgModel_SDP
    appendcons(m.task,1)
    constridx = getnumcon(m.task)
    for i in 1:length(matvaridx)
        putbaraij(m.task, constridx, matvaridx[i], [matcoefidx[i]], [1])
    end

    if all(k -> m.varmap[k] > 0, scalidx)
        putarow(m.task,constridx,m.varmap[scalidx],scalcoef)
    else # SDP terms exist in linear part

        local idxs    = find(k -> m.varmap[k] > 0, scalidx)
        local baridxs = find(k -> m.varmap[k] < 0, scalidx)
        
        # add linear terms
        putarow(m.task,constridx,m.varmap[scalidx[idxs]],scalcoef[idxs])
        
        # add SDP terms
        let barcof   = scalcoef[idxs],
            barj     = m.varmap[scalidx[idxs]],
            barvarij = m.barvar[scalidx[idxs]]
            
            local perm = sortperm(barj)
            local np = length(perm)
            
            local p = 1
            while p <= np
                const b = p
                const j = barj[perm[p]]
                const dimbarvarj = getdimbarvarj(m.task,j)
                p += 1
                while (p <= np && barj[perm[p]] == j) p += 1 end
                
                const JJ = map(L -> int32(floor(n+0.5-sqrt((n+0.5)^2-2*L))), barvarij[perm[b:p-1]])
                const II = map(i -> int32(barvarij[perm[i+b-1]]-n*II[i]), 1:p-b-1)
                const VV = barcof[perm[b:p-1]]

                const matidx = appendsparsesymmat(m.task,dimbarvarj,II,JJ,VV)
                putbaraij(m.task, constridx, J, [matidx], [1.0])                
            end
            
        end
    end

    bk = getBoundsKey(lb, ub)
    putconbound(m.task,constridx,bk,lb,ub)

    userconidx = addUserCon(m,constridx)

    return int64(userconidx)
end

function setsdpobj!(m::MosekMathProgModel, matvaridx, matcoefidx)
    for (it,varidx) in enumerate(m.barvarmap[matvaridx])
        putbarcj(m.task, varidx, [matcoefidx[it]], [1])
    end
end

function getsdpsolution(m::MosekMathProgModel, idx)
    V = getbarxj(m.task, MSK_SOL_ITR, m.barvarmap[idx])
    n = convert(Int64, sqrt(8*length(V)+1)/2-1/2)
    cnt = 0
    A = Array(Float64,n,n)
    for j in 1:n
      cnt += 1
      A[j,j] = V[cnt]
      for i in (j+1):n
        cnt += 1
        A[i,j] = V[cnt]
        A[j,i] = V[cnt]
      end
    end
    return A
end

getsdpdual(m::MosekMathProgModel) = getconstrduals(m)

# hmm... missing getdualpsdvarsolution?
