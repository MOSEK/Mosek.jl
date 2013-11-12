using Mosek

function test_lo1()
    printstream(msg::String) = print(msg)

    bkc = [MSK_BK_FX MSK_BK_LO MSK_BK_UP]
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
    putacolslice(task,1,numvar+1,
                 A.colptr[1:numvar],A.colptr[2:numvar+1],
                 A.rowval,
                 A.nzval)
    putvarboundslice(task, 1, numvar+1, bkx,blx,bux)
    putconboundslice(task,1,numcon+1,bkc,blc,buc)
    putobjsense(task,MSK_OBJECTIVE_SENSE_MAXIMIZE)
    optimize(task)
    solutionsummary(task,MSK_STREAM_MSG)
    solsta = getsolsta(task,MSK_SOL_BAS)

    if solsta in     [ MSK_SOL_STA_OPTIMAL, 
                       MSK_SOL_STA_NEAR_OPTIMAL ]
      xx = getxx(task,MSK_SOL_BAS)
      # check feasibility and optimality of solution
      return true
    else
      return false
    end
end

function test_qo1()
    printstream(msg::String) = print(msg)
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
    putclist(task,[1:numvar],c)
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
    if solsta in     [ MSK_SOL_STA_OPTIMAL, 
                       MSK_SOL_STA_NEAR_OPTIMAL ]
      xx = getxx(task,MSK_SOL_ITR)
      # check feasibility and optimality of solution
      return true
    else
      return false
    end
end

function test_qcqo1()
    task = maketask()
    bkc   = [ MSK_BK_LO ]
    blc   = [ 1.0 ]
    buc   = [ Inf ]
    bkx   = [ MSK_BK_LO, MSK_BK_LO, MSK_BK_LO ]
    blx   = [ 0.0,  0.0, 0.0 ]
    bux   = [ Inf,  Inf, Inf ]
    c     = [ 0.0, -1.0, 0.0 ]
    asub  = [ 1 ,1, 1 ]
    aval  = [ 1.0, 1.0, 1.0 ]
    numvar = length(bkx)
    numcon = length(bkc)
    appendcons(task,numcon)
    appendvars(task,numvar)
    putcfix(task,0.0)
    putclist(task,[1:numvar],c)
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
    putqconk (task,1, qsubi,qsubj, qval) 
    putobjsense(task,MSK_OBJECTIVE_SENSE_MINIMIZE)
    optimize(task)
    solutionsummary(task,MSK_STREAM_MSG)
    prosta = getprosta(task,MSK_SOL_ITR)
    solsta = getsolsta(task,MSK_SOL_ITR)
    if solsta in     [ MSK_SOL_STA_OPTIMAL, 
                       MSK_SOL_STA_NEAR_OPTIMAL ]
      xx = getxx(task,MSK_SOL_ITR)
      # check feasibility and optimality of solution
      return true
    else
      return false
    end
end

function test_milo1()
    printstream(msg::String) = print(msg)
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
    putclist(task,[1:numvar],c)
    putvarboundslice(task,1,numvar+1,bkx,blx,bux)
    putacolslice(task,1,numvar+1, A.colptr[1:numvar],A.colptr[2:numvar+1],A.rowval,A.nzval)     
    putconboundslice(task,1,numcon+1,bkc,blc,buc)
    putobjsense(task,MSK_OBJECTIVE_SENSE_MAXIMIZE)
    putvartypelist(task,[ 1, 2 ], [ MSK_VAR_TYPE_INT, MSK_VAR_TYPE_INT ])
    optimize(task)
    solutionsummary(task,MSK_STREAM_MSG)
    prosta = getprosta(task,MSK_SOL_ITG)
    solsta = getsolsta(task,MSK_SOL_ITG)
    if solsta in     [ MSK_SOL_STA_INTEGER_OPTIMAL, 
                       MSK_SOL_STA_NEAR_INTEGER_OPTIMAL ]
      xx = getxx(task,MSK_SOL_ITG)
      # check feasibility and optimality of solution
      return true
    else
      return false
    end
end


function test_cqo1()
    printstream(msg::String) = print(msg)
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
    asub  = [ 1 ,1, 1 ]
    aval  = [ 1.0, 1.0, 1.0 ]
    numvar = length(bkx)
    numcon = length(bkc)
    task = maketask()
    putstreamfunc(task,MSK_STREAM_LOG,printstream)
    putcallbackfunc(task,callback)
    appendcons(task,numcon)
    appendvars(task,numvar)
    putclist(task,[1:6],c)
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
    
    if solsta in     [ MSK_SOL_STA_OPTIMAL, 
                       MSK_SOL_STA_NEAR_OPTIMAL ]
      xx = getxx(task,MSK_SOL_ITR)
      # check feasibility and optimality of solution
      return true
    else
      return false
    end
end
function test_sdo1()
    task = maketask()
    bkc = [MSK_BK_FX,
           MSK_BK_FX]
    blc = [1.0, 0.5]
    buc = [1.0, 0.5]
    A = sparse( [1,2,2],[1,2,3],[1.0, 1.0, 1.0])
    conesub = [1, 2, 3]
    barci = [1, 2, 2, 3, 3]
    barcj = [1, 1, 2, 2, 3]
    barcval = [2.0, 1.0, 2.0, 1.0, 2.0]
    barai   = { [1, 2, 3], 
                [1, 2, 3, 2, 3, 3] }
    baraj   = { [1, 2, 3]
                [1, 1, 1, 2, 2, 3] }
    baraval = { [1.0, 1.0, 1.0],
                [1.0, 1.0, 1.0, 1.0, 1.0, 1.0] }
    numvar = 3
    numcon = length(bkc)
    barvardim = [3]
    appendvars(task,numvar)
    appendcons(task,numcon)
    appendbarvars(task,barvardim)
    putcj(task, 1, 1.0)
    putvarboundslice(task,1,numvar+1,
                     [ MSK_BK_FR::Int32 for i in 1:numvar ],
                     [ -Inf             for i in 1:numvar ],
                     [ +Inf             for i in 1:numvar ])
    putconboundslice(task,1,numcon+1, bkc,blc,buc)
    putacolslice(task,1,numvar+1,
                 A.colptr[1:numvar], A.colptr[2:numvar+1],
                 A.rowval,A.nzval)
    appendcone(task,MSK_CT_QUAD, 0.0, conesub)
    symc  = appendsparsesymmat(task,barvardim[1], 
                               barci, 
                               barcj, 
                               barcval)
    syma0 = appendsparsesymmat(task,barvardim[1], 
                               barai[1], 
                               baraj[1], 
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
    
    if solsta in     [ MSK_SOL_STA_OPTIMAL, 
                       MSK_SOL_STA_NEAR_OPTIMAL ]
      xx = getxx(task,MSK_SOL_ITR)
      # check feasibility and optimality of solution
      return true
    else
      return false
    end
end

function test_nlo1()
    t = maketask()
    putstreamfunc(t,MSK_STREAM_LOG,msg -> print(msg))
    appendvars(t,3)
    appendcons(t,3)
    putvarbound(t,1,MSK_BK_RA,0.0, 10.0)  # x
    putvarbound(t,2,MSK_BK_RA,0.0, 10.0)  # y
    putvarbound(t,3,MSK_BK_LO,0.0, +Inf)  # t
    putarow(t,1,[1,2,3],[1.0, 2.0, -1.0]) # x + 2y - t
    putconbound(t,1,MSK_BK_UP,-Inf, 0.0)  # < 0.0
    putarow(t,2,[2],[-1.0])               # - y
    putconbound(t,2,MSK_BK_LO,0.0,+Inf)   # 0.0 <
    putarow(t,3,[1,2],[-10.0, -1.0])      # -10x - y
    putqconk(t,3,[1],[1],[2.0])           # x^2 
    putconbound(t,3,MSK_BK_UP,-Inf,-27.0) # < -27.0
    
    function indexof(v,a)
      for i in 1:length(a)
        if a[i] == v return i
        end
      end
      return 0
    end

    evalobj(x :: Array{Float64,1}) = x[3] ^ 2.3
    function evalconi(x:: Array{Float64,1},i:: Int32)
      if i == 2
        sqrt(10)*sqrt(x[1])
      else
        0.0
      end
    end

    function grdobj(x :: Array{Float64,1},sub:: Array{Int32,1}, val:: Array{Float64,1})
      for i=1:length(sub)
        if sub[i] == 3
          val[i] = 2.3 * x[3] ^ 1.3
        else
          val[i] = 0.0
        end
      end
    end
            
    function grdconi(x :: Array{Float64,1},
                     i:: Int32, 
                     sub:: Array{Int32,1}, 
                     val:: Array{Float64,1})
      if i == 2
        for k=1:length(sub)
          if sub[k] == 3
            val[k] = 0.5 * x[1] ^ (-0.5)
          else
            val[k] = 0.0
          end
        end
      else
        val[1:length(val)] = 0.0
      end
    end

    function grdlag(x ::   Array{Float64,1},
                    yo::   Float64,
                    yc::   Array{Float64,1},
                    subi:: Array{Int32,1},
                    val::  Array{Float64,1})
      val[3] = yo * 2.3 * x[3] ^ 1.3
      
      k = indexof(2,subi)
      val[1] += yc[k] * 0.5 * x[1] ^ (-0.5)
    end 
    
    function heslag(x ::      Array{Float64,1},
                    yo::      Float64,
                    yc::      Array{Float64,1},
                    subi::    Array{Int32,1},
                    hessubi:: Array{Int32,1},
                    hessubj:: Array{Int32,1},
                    hesval::  Array{Float64,1})
      # d^2/dt^2 (yo * t^2.3) = yo * 2.3 * 1.3 * t^0.3
      hessubi[1] = 3
      hessubj[1] = 3
      hesval[1]  = yo * 2.3 * 1.3 * x[3]^0.3
      
      # d^2/dx^2 (yc * 10^(1/2) * x^(1/2)) = yc * (-1/4) x^(-3/2)
      k = indexof(2,subi)      
      hessubi[2] = 1
      hessubj[2] = 1
      hesval[2]  = - yc[k] * 0.25 * x[1] ^ (-1.5)
    end

    putnlcallbacks(t,
                   [3], # subscripts of non-zeros in the gradient of the objective
                   [1], # subscripts of non-zeros in the gradient of the constraints
                   [1,1,2,2], # rowptr for subscripts of non-zeros in the gradient of the constraints
                   [1,3], # hessubi
                   [1,3], # hessubj
                   evalobj,
                   evalconi,
                   grdlag,
                   grdobj,
                   grdconi,
                   heslag)

    optimize(t)

    solutionsummary(t,MSK_STREAM_MSG)

    solsta = getsolsta(t,MSK_SOL_ITR)
    if solsta in     [ MSK_SOL_STA_OPTIMAL, 
                       MSK_SOL_STA_NEAR_OPTIMAL ]
      xx = getxx(task,MSK_SOL_ITR)
      # check feasibility and optimality of solution
      return true
    else
      return false
    end
end

assert(test_lo1() &&
       test_qo1() &&
       test_qcqo1() &&
       test_milo1() &&
       test_cqo1() && 
       test_sdo1())
# NOTE: nlo1 broken 

