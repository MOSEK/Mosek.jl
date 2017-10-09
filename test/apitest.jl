using Mosek, Base.Test

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
    optimize(task)
    solutionsummary(task,MSK_STREAM_MSG)

    solsta = getsolsta(task,MSK_SOL_BAS)
    prosta = getprosta(task,MSK_SOL_BAS)

    @test getsolsta(task,MSK_SOL_ITR) in [ MSK_SOL_STA_OPTIMAL, MSK_SOL_STA_NEAR_OPTIMAL ]
    @test solsta in [ MSK_SOL_STA_OPTIMAL, MSK_SOL_STA_NEAR_OPTIMAL ]
    @test prosta in [MSK_PRO_STA_PRIM_AND_DUAL_FEAS,MSK_PRO_STA_NEAR_PRIM_AND_DUAL_FEAS]

    if solsta in     [ MSK_SOL_STA_OPTIMAL, 
                       MSK_SOL_STA_NEAR_OPTIMAL ]
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
    optimize(task)
    solutionsummary(task,MSK_STREAM_MSG)

    prosta = getprosta(task,MSK_SOL_ITR)
    solsta = getsolsta(task,MSK_SOL_ITR)

    @test solsta in (MSK_SOL_STA_OPTIMAL, MSK_SOL_STA_NEAR_OPTIMAL)
    @test prosta in (MSK_PRO_STA_PRIM_AND_DUAL_FEAS,MSK_PRO_STA_NEAR_PRIM_AND_DUAL_FEAS)

    if solsta in     [ MSK_SOL_STA_OPTIMAL, 
                       MSK_SOL_STA_NEAR_OPTIMAL ]
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
    optimize(task)
    solutionsummary(task,MSK_STREAM_MSG)

    prosta = getprosta(task,MSK_SOL_ITR)
    solsta = getsolsta(task,MSK_SOL_ITR)


    @test solsta in (MSK_SOL_STA_OPTIMAL, MSK_SOL_STA_NEAR_OPTIMAL)
    @test prosta in (MSK_PRO_STA_PRIM_AND_DUAL_FEAS,MSK_PRO_STA_NEAR_PRIM_AND_DUAL_FEAS)

    if solsta in     [ MSK_SOL_STA_OPTIMAL, 
                       MSK_SOL_STA_NEAR_OPTIMAL ]
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
    optimize(task)
    solutionsummary(task,MSK_STREAM_MSG)
    prosta = getprosta(task,MSK_SOL_ITG)
    solsta = getsolsta(task,MSK_SOL_ITG)

    @test solsta in (MSK_SOL_STA_INTEGER_OPTIMAL, MSK_SOL_STA_NEAR_INTEGER_OPTIMAL)
    @test prosta in (MSK_PRO_STA_PRIM_FEAS,MSK_PRO_STA_NEAR_PRIM_FEAS)

    if solsta in     [ MSK_SOL_STA_INTEGER_OPTIMAL, 
                       MSK_SOL_STA_NEAR_INTEGER_OPTIMAL ]
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
    optimize(task)
    solutionsummary(task,MSK_STREAM_MSG)
    prosta = getprosta(task,MSK_SOL_ITR)
    solsta = getsolsta(task,MSK_SOL_ITR)

    @test solsta in (MSK_SOL_STA_OPTIMAL, MSK_SOL_STA_NEAR_OPTIMAL)
    @test prosta in (MSK_PRO_STA_PRIM_AND_DUAL_FEAS,MSK_PRO_STA_NEAR_PRIM_AND_DUAL_FEAS)

    if solsta in     [ MSK_SOL_STA_OPTIMAL, 
                       MSK_SOL_STA_NEAR_OPTIMAL ]
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
    optimize(task)
    solutionsummary(task,MSK_STREAM_MSG)
    prosta = getprosta(task,MSK_SOL_ITR)
    solsta = getsolsta(task,MSK_SOL_ITR)

    @test solsta in (MSK_SOL_STA_OPTIMAL, MSK_SOL_STA_NEAR_OPTIMAL)
    @test prosta in (MSK_PRO_STA_PRIM_AND_DUAL_FEAS,MSK_PRO_STA_NEAR_PRIM_AND_DUAL_FEAS)
    
    if solsta in     [ MSK_SOL_STA_OPTIMAL, 
                       MSK_SOL_STA_NEAR_OPTIMAL ]
        xx = getxx(task,MSK_SOL_ITR)
        # check feasibility and optimality of solution
    end
end

function test_nlo1()
    t = maketask()
    putstreamfunc(t,MSK_STREAM_LOG,m -> nothing)
    appendvars(t,3)
    appendcons(t,1)
    putvarbound(t,1,MSK_BK_LO,0.0, +Inf)  # x0
    putvarbound(t,2,MSK_BK_LO,0.0, +Inf)  # x1
    putvarbound(t,3,MSK_BK_LO,0.0, +Inf)  # x2
    putarow(t,1,[1,2,3],[1.0,1.0,1.0])   # x0 + x1 + x2
    putconbound(t,1,MSK_BK_FX,1.0, 1.0)  # = 1.0
    putclist(t,[1], [-1.0])

    function indexof(v,a)
        for i in 1:length(a)
            if a[i] == v return i
            end
        end
        return 0
    end
    evalobj(x :: Vector{Float64}) = log(x[2]+x[3])
    evalconi(x:: Vector{Float64},i:: Int32) = 0.0
    function grdobj(x :: Vector{Float64},sub:: Vector{Int32}, val:: Vector{Float64})
        for j in 1:length(sub)
            if sub[j] == 2 || sub[j] == 3
                val[j] = -1.0 / (x[2]+x[3])
            end
        end
    end
    function grdconi(x  :: Vector{Float64},
                     i  :: Int32,
                     sub:: Vector{Int32},
                     val:: Vector{Float64})
        none
    end
    function grdlag(x ::   Vector{Float64},
                    yo::   Float64,
                    yc::   Vector{Float64},
                    subi:: Vector{Int32},
                    val::  Vector{Float64})
        val[2] = - yo * 1.0 / (x[2]+x[3])
        val[3] = - yo * 1.0 / (x[2]+x[3])
    end 
    function heslag(x ::      Vector{Float64},
                    yo::      Float64,
                    yc::      Vector{Float64},
                    subi::    Vector{Int32},
                    hessubi:: Vector{Int32},
                    hessubj:: Vector{Int32},
                    hesval::  Vector{Float64})

        hessubi[1] = 2; hessubj[1] = 2; hesval[1] = (x[2]+x[3])^(-2)
        hessubi[2] = 3; hessubj[2] = 2; hesval[2] = (x[2]+x[3])^(-2)
        hessubi[3] = 3; hessubj[3] = 3; hesval[3] = (x[2]+x[3])^(-2)
    end

    putnlcallbacks(t,
                   [2,3], # subscripts of non-zeros in the gradient of the objective
                   Int[], # subscripts of non-zeros in the gradient of the constraints
                   [1,1], # rowptr for subscripts of non-zeros in the gradient of the constraints
                   [2,3,3], # hessubi
                   [2,2,3], # hessubj
                   evalobj,
                   evalconi,
                   grdlag,
                   grdobj,
                   grdconi,
                   heslag)

    optimize(t)

    solsta = getsolsta(t,MSK_SOL_ITR)
    prosta = getprosta(t,MSK_SOL_ITR)

    @test solsta in (MSK_SOL_STA_OPTIMAL, MSK_SOL_STA_NEAR_OPTIMAL)
    @test prosta in (MSK_PRO_STA_PRIM_AND_DUAL_FEAS,MSK_PRO_STA_NEAR_PRIM_AND_DUAL_FEAS)
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

    @testset "nlo1" begin
        test_nlo1()
    end
end

