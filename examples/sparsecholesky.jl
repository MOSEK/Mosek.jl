#
#   Copyright: Copyright (c) MOSEK ApS, Denmark. All rights reserved.
#
#   File:      sparsecholesky.jl
#
#   Purpose:   Demonstrate the sparse Cholesky factorization.
#
using Mosek

function printsparse(n          :: Int32,
                     perm       :: Vector{Int32},
                     diag       :: Vector{Float64},
                     lnzc       :: Vector{Int32},
                     lptrc      :: Vector{Int64},
                     lensubnval :: Int64,
                     lsubc      :: Vector{Int32},
                     lvalc      :: Vector{Float64})
    println("P = $(perm)")
    println("diag(D) = $(diag)")

    l = zeros(Int32,(n,n))
    for j in 1:n
        for i in lptrc[j]:lptrc[j]+lnzc[j]-1
            l[lsubc[i],j] = lvalc[i]
        end
    end
    println("L = $l")
end

function main()
    # Example from the manual
    # Observe that anzc, aptrc, asubc and avalc only specify the lower triangular part.
    n     = Int32(4)
    anzc  = Int32[4, 1, 1, 1]
    asubc = Int32[0, 1, 2, 3,
                     1,
                        2,
                           3]
    aptrc = Int64[0, 4, 5, 6]
    avalc = Float64[4.0, 1.0, 1.0, 1.0,
                         1.0,
                              1.0,
                                   1.0]
    b     = Float64[13.0, 3.0, 4.0, 5.0]

    (perm,
     diag,
     lnzc,
     lptrc,
     lensubnval,
     lsubc,
     lvalc) = computesparsecholesky(Int32(0),  #Mosek chooses number of threads
                                    Int32(1),  #Apply reordering heuristic
                                    1.0e-14,  #Singularity tolerance
                                    anzc, aptrc, asubc, avalc)

    printsparse(n, perm, diag, lnzc, lptrc, lensubnval, lsubc, lvalc)

    # Permuted b is stored as x.
    x = b[perm]

    # Compute  inv(L)*x.
    sparsetriangularsolvedense(MSK_TRANSPOSE_NO,  lnzc, lptrc, lsubc, lvalc, x)
    # Compute  inv(L^T)*x.
    sparsetriangularsolvedense(MSK_TRANSPOSE_YES, lnzc, lptrc, lsubc, lvalc, x)

    println("\nSolution A x = b, x = $([ x[j] for i in 1..n for j in 1..n if perm[j] == i ])")

    n     = Int32(3)
    anzc  = Int32[3, 2, 1]
    asub  = Int32[0, 1, 2, 1, 2, 2]
    aptr  = Int64[0, 3, 5, ]
    avalc = Float64[1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

    (perm,
     diag,
     lnzc,
     lptrc,
     lensubnval,
     lsubc,
     lvalc) = computesparsecholesky(0,        #Mosek chooses number of threads
                                    1,        #Apply reordering heuristic
                                    1.0e-14,  #Singularity tolerance
                                    anzc, aptrc, asubc, avalc)

    printsparse(n, perm, diag, lnzc, lptrc, lensubnval, lsubc, lvalc)
end


main()
