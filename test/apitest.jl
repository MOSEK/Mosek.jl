using Mosek, Test
using SparseArrays
using Printf

function test_lo1()
    printstream(msg::String) = nothing # print(msg)

    bkc = [MSK_BK_FX, MSK_BK_LO, MSK_BK_UP]
    blc = [30.0, 15.0, -Inf]
    buc = [30.0, +Inf, 25.0]
    bkx = [ MSK_BK_LO, MSK_BK_RA, MSK_BK_LO, MSK_BK_LO ]
    blx = [   0.0,  0.0,    0.0,    0.0]
    bux = [+Inf, 10.0, +Inf, +Inf]
    numvar = length(bkx)
    numcon = length(bkc)
    c = [ 3.0, 1.0, 5.0, 1.0 ]
    A = sparse([1, 2, 1, 2, 3, 1, 2, 2, 3],
               [1, 1, 2, 2, 2, 3, 3, 4, 4],
               [3.0, 2.0, 1.0, 1.0, 2.0, 2.0, 3.0, 1.0, 3.0 ],
               numcon,numvar)

    task = maketask()
    putstreamfunc(task,MSK_STREAM_LOG,printstream)
    putobjname(task,"lo1")
    appendcons(task,numcon)
    for i=1:numcon
        putconname(task,i,@sprintf("c%02d",i))
    end
    appendvars(task,numvar)
    for j=1:numvar
        putvarname(task,j,@sprintf("x%02d",j))
    end
    putclist(task,[1,2,3,4], c)

    putacolslice(task,
                 1, numvar+1,
                 A.colptr[1:numvar], A.colptr[2:numvar+1],
                 A.rowval,
                 A.nzval)
    putvarboundslice(task, 1, numvar+1, bkx,blx,bux)
    putconboundslice(task,1,numcon+1,bkc,blc,buc)
    putobjsense(task,MSK_OBJECTIVE_SENSE_MAXIMIZE)
    optimize(task,"mosek://solve.mosek.com:30080")
    solutionsummary(task,MSK_STREAM_MSG)

    solsta = getsolsta(task,MSK_SOL_BAS)
    prosta = getprosta(task,MSK_SOL_BAS)

    @test getsolsta(task,MSK_SOL_ITR) == MSK_SOL_STA_OPTIMAL
    @test solsta == MSK_SOL_STA_OPTIMAL
    @test prosta == MSK_PRO_STA_PRIM_AND_DUAL_FEAS

    if solsta in     [ MSK_SOL_STA_OPTIMAL ]
        xx = getxx(task,MSK_SOL_BAS)
        # check feasibility and optimality of solution
    end
end

function test_qo1()
    printstream(msg::String) = nothing # print(msg)
    bkc   = [ MSK_BK_LO ]
    blc   = [ 1.0 ]
    buc   = [ Inf ]
    bkx   = [ MSK_BK_LO, MSK_BK_LO, MSK_BK_LO ]
    blx   = [ 0.0,  0.0, 0.0 ]
    bux   = [ Inf,  Inf, Inf ]
    numvar = length(bkx)
    numcon = length(bkc)
    c     = [ 0.0, -1.0, 0.0 ]
    A     = sparse( [ 1, 1, 1 ],
                    [ 1, 2, 3 ],
                    [ 1.0, 1.0, 1.0 ],
                    numcon, numvar )
    task = maketask()
    putstreamfunc(task,MSK_STREAM_LOG,printstream)
    appendcons(task,numcon)
    appendvars(task,numvar)
    putclist(task,[1:numvar;],c)
    putvarboundslice(task,1,numvar+1,bkx,blx,bux)
    putacolslice(task,1,numvar+1,
                 A.colptr[1:numvar],A.colptr[2:numvar+1],
                 A.rowval,A.nzval)
    qsubi = [ 1,   2,    3,   3   ]
    qsubj = [ 1,   2,    1,   3   ]
    qval  = [ 2.0, 0.2, -1.0, 2.0 ]
    putqobj(task,qsubi,qsubj,qval)
    putobjsense(task,MSK_OBJECTIVE_SENSE_MINIMIZE)
    optimize(task,"mosek://solve.mosek.com:30080")
    solutionsummary(task,MSK_STREAM_MSG)

    prosta = getprosta(task,MSK_SOL_ITR)
    solsta = getsolsta(task,MSK_SOL_ITR)

    @test solsta == MSK_SOL_STA_OPTIMAL
    @test prosta == MSK_PRO_STA_PRIM_AND_DUAL_FEAS

    if solsta in     [ MSK_SOL_STA_OPTIMAL ]
        xx = getxx(task,MSK_SOL_ITR)
        # check feasibility and optimality of solution
    end
end

function test_qcqo1()
    printstream(msg::String) = nothing # print(msg)
    task = maketask()
    putstreamfunc(task,MSK_STREAM_LOG,m -> nothing)
    bkc   = [ MSK_BK_LO ]
    blc   = [ 1.0 ]
    buc   = [ Inf ]
    bkx   = [ MSK_BK_LO, MSK_BK_LO, MSK_BK_LO ]
    blx   = [ 0.0,  0.0, 0.0 ]
    bux   = [ Inf,  Inf, Inf ]
    c     = [ 0.0, -1.0, 0.0 ]
    asub  = [ 1 ,2, 3 ]
    aval  = [ 1.0, 1.0, 1.0 ]
    numvar = length(bkx)
    numcon = length(bkc)
    appendcons(task,numcon)
    appendvars(task,numvar)
    putcfix(task,0.0)
    putclist(task,[1:numvar;],c)
    putvarboundslice(task,1,numvar+1,bkx,blx,bux)
    putarow(task,1,asub,aval)
    putconbound(task,1,bkc[1],blc[1],buc[1])
    qsubi = [ 1,   2,    3,   3   ]
    qsubj = [ 1,   2,    1,   3   ]
    qval  = [ 2.0, 0.2, -1.0, 2.0 ]
    putqobj(task,qsubi,qsubj,qval)
    qsubi = [  1,    2,    3,   3   ]
    qsubj = [  1,    2,    3,   1   ]
    qval  = [ -2.0, -2.0, -0.2, 0.2 ]
    putqconk(task,1, qsubi,qsubj, qval)
    putobjsense(task,MSK_OBJECTIVE_SENSE_MINIMIZE)
    optimize(task,"mosek://solve.mosek.com:30080")
    solutionsummary(task,MSK_STREAM_MSG)

    prosta = getprosta(task,MSK_SOL_ITR)
    solsta = getsolsta(task,MSK_SOL_ITR)


    @test solsta in [MSK_SOL_STA_OPTIMAL ]
    @test prosta in [MSK_PRO_STA_PRIM_AND_DUAL_FEAS]

    if solsta in     [ MSK_SOL_STA_OPTIMAL ]
        xx = getxx(task,MSK_SOL_ITR)
        # check feasibility and optimality of solution
    end
end

function test_milo1()
    printstream(msg::String) = print(msg)  #print(msg)
    task = maketask()
    putstreamfunc(task,MSK_STREAM_LOG,printstream)
    bkc = [ MSK_BK_UP, MSK_BK_LO  ]
    blc = [      -Inf,      -4.0  ]
    buc = [     250.0,       Inf  ]
    bkx = [ MSK_BK_LO, MSK_BK_LO  ]
    blx = [       0.0,       0.0  ]
    bux = [       Inf,       Inf  ]
    c   = [       1.0,      0.64 ]
    A    = sparse( [ 1, 1, 2, 2],
                   [ 1, 2, 1, 2],
                   [ 50.0, 3.0, 31.0, -2.0])
    numvar = length(bkx)
    numcon = length(bkc)
    appendcons(task,numcon)
    appendvars(task,numvar)
    putclist(task,[1:numvar;],c)
    putvarboundslice(task,1,numvar+1,bkx,blx,bux)
    putacolslice(task,1,numvar+1, A.colptr[1:numvar],A.colptr[2:numvar+1],A.rowval,A.nzval)
    putconboundslice(task,1,numcon+1,bkc,blc,buc)
    putobjsense(task,MSK_OBJECTIVE_SENSE_MAXIMIZE)
    putvartypelist(task,[ 1, 2 ], [ MSK_VAR_TYPE_INT, MSK_VAR_TYPE_INT ])
    optimize(task,"mosek://solve.mosek.com:30080")
    solutionsummary(task,MSK_STREAM_MSG)
    prosta = getprosta(task,MSK_SOL_ITG)
    solsta = getsolsta(task,MSK_SOL_ITG)

    @test solsta in [MSK_SOL_STA_INTEGER_OPTIMAL]
    @test prosta in [MSK_PRO_STA_PRIM_FEAS]

    if solsta in     [ MSK_SOL_STA_INTEGER_OPTIMAL ]
        xx = getxx(task,MSK_SOL_ITG)
    end
end

function test_cqo()
        printstream(msg::String) = nothing # print(msg)
    callback(where,dinf,iinf,liinf) = 0
    bkc = [ MSK_BK_FX ]
    blc = [ 1.0 ]
    buc = [ 1.0 ]
    c   = [               0.0,              0.0,              0.0,
                          1.0,              1.0,              1.0 ]
    bkx = [ MSK_BK_LO,MSK_BK_LO,MSK_BK_LO,
            MSK_BK_FR,MSK_BK_FR,MSK_BK_FR ]
    blx = [               0.0,              0.0,              0.0,
                          -Inf,             -Inf,             -Inf ]
    bux = [               Inf,              Inf,              Inf,
                          Inf,              Inf,              Inf ]
    asub  = [ 1 ,2, 3 ]
    aval  = [ 1.0, 1.0, 1.0 ]
    numvar = length(bkx)
    numcon = length(bkc)
    task = maketask()
    putstreamfunc(task,MSK_STREAM_LOG,printstream)
    putcallbackfunc(task,callback)
    appendcons(task,numcon)
    appendvars(task,numvar)
    putclist(task,[1:6;],c)
    putvarboundslice(task,1,numvar+1,bkx,blx,bux)
    putarow(task,1,asub,aval)
    putconbound(task,1,bkc[1],blc[1],buc[1])
    appendcone(task,MSK_CT_QUAD, 0.0, [ 4, 1, 2 ])
    appendcone(task,MSK_CT_RQUAD, 0.0, [ 5, 6, 3 ])
    putobjsense(task,MSK_OBJECTIVE_SENSE_MINIMIZE)
    optimize(task,"mosek://solve.mosek.com:30080")
    solutionsummary(task,MSK_STREAM_MSG)
    prosta = getprosta(task,MSK_SOL_ITR)
    solsta = getsolsta(task,MSK_SOL_ITR)

    @test solsta in [MSK_SOL_STA_OPTIMAL]
    @test prosta in [MSK_PRO_STA_PRIM_AND_DUAL_FEAS]

    if solsta in     [ MSK_SOL_STA_OPTIMAL ]
        xx = getxx(task,MSK_SOL_ITR)
        # check feasibility and optimality of solution
    end
end

function test_sdo1()
    task = maketask()
    putstreamfunc(task,MSK_STREAM_LOG,m -> nothing)
    bkc = [MSK_BK_FX,
           MSK_BK_FX]
    blc = [1.0, 0.5]
    buc = [1.0, 0.5]
    A = sparse( [1,2,2],[1,2,3],[1.0, 1.0, 1.0])
    conesub = [1, 2, 3]
    barci = [1, 2, 2, 3, 3]
    barcj = [1, 1, 2, 2, 3]
    barcval = [2.0, 1.0, 2.0, 1.0, 2.0]
    barai   = Any[ [1, 2, 3],
                   [1, 2, 3, 2, 3, 3] ]
    baraj   = Any[ [1, 2, 3],
                   [1, 1, 1, 2, 2, 3] ]
    baraval = Any[ [1.0, 1.0, 1.0],
                   [1.0, 1.0, 1.0, 1.0, 1.0, 1.0] ]
    numvar = 3
    numcon = length(bkc)
    barvardim = [3]
    appendvars(task,numvar)
    appendcons(task,numcon)
    appendbarvars(task,barvardim)
    putcj(task, 1, 1.0)
    putvarboundslice(task,1,numvar+1,
                     [ MSK_BK_FR for i in 1:numvar ],
                     [ -Inf      for i in 1:numvar ],
                     [ +Inf      for i in 1:numvar ])
    putconboundslice(task,1,numcon+1, bkc,blc,buc)
    putacolslice(task,1,numvar+1,
                 A.colptr[1:numvar], A.colptr[2:numvar+1],
                 A.rowval,A.nzval)
    appendcone(task,MSK_CT_QUAD, 0.0, conesub)
    symc  = appendsparsesymmat(task,barvardim[1],
                               barci,
                               barcj,
                               barcval)

    syma0 = appendsparsesymmat(task,
                               barvardim[1],
                               barai[1],
                               baraj[1],+
                               baraval[1])
    syma1 = appendsparsesymmat(task,barvardim[1],
                               barai[2],
                               baraj[2],
                               baraval[2])
    putbarcj(task,1, [symc], [1.0])
    putbaraij(task,1, 1, [syma0], [1.0])
    putbaraij(task,2, 1, [syma1], [1.0])
    putobjsense(task,MSK_OBJECTIVE_SENSE_MINIMIZE)
    optimize(task,"mosek://solve.mosek.com:30080")
    solutionsummary(task,MSK_STREAM_MSG)
    prosta = getprosta(task,MSK_SOL_ITR)
    solsta = getsolsta(task,MSK_SOL_ITR)

    @test solsta in [MSK_SOL_STA_OPTIMAL]
    @test prosta in [MSK_PRO_STA_PRIM_AND_DUAL_FEAS]

    if solsta in     [ MSK_SOL_STA_OPTIMAL ]
        xx = getxx(task,MSK_SOL_ITR)
        # check feasibility and optimality of solution
    end
end

function test_removecones()
    task = maketask()
    appendvars(task, 7)
    appendcone(task, MSK_CT_QUAD, 0.0, [1, 2, 3])
    appendcone(task, MSK_CT_RQUAD, 0.0, [4, 5, 6, 7])
    @test getcone(task, 1) == (Mosek.MSK_CT_QUAD, 0.0, 3, Int32[1, 2, 3])
    @test getconeinfo(task, 1) == (Mosek.MSK_CT_QUAD, 0.0, 3)
    @test getcone(task, 2) == (Mosek.MSK_CT_RQUAD, 0.0, 4, Int32[4, 5, 6, 7])
    @test getconeinfo(task, 2) == (Mosek.MSK_CT_RQUAD, 0.0, 4)
    removecones(task, [1])
    info = getconeinfo(task, 1)
    # `info[1]` is Mosek.MSK_CT_QUAD
    #@test_broken info[1] == Mosek.MSK_CT_RQUAD
    @test info[2:end] == (0.0, 4)
    cone = getcone(task, 1)
    # `cone[1]` is Mosek.MSK_CT_QUAD
    #@test_broken cone[1] == Mosek.MSK_CT_RQUAD
    @test cone[2:end] == (0.0, 4, Int32[4, 5, 6, 7])
end

@testset "[apitest]" begin
    @testset "lo1" begin
        test_lo1()
    end

    @testset "qo1" begin
        test_qo1()
    end

    @testset "qcqo1" begin
        test_qcqo1()
    end

    @testset "milo1" begin
        test_milo1()
    end

    @testset "cqo" begin
        test_cqo()
    end

    @testset "sdo1" begin
        test_sdo1()
    end

    @testset "removecones" begin
        test_removecones()
    end

end
