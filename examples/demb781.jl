#
#  Copyright: Copyright (c) MOSEK ApS, Denmark. All rights reserved.
#
#  File:    demb781.jl
#
#  Purpose: Demonstrates how to solve a simple non-liner separable problem
#  using the SCopt interface for Python. Then problem is this:
#    Minimize   e^x2 + e^x3
#    Such that  e^x4 + e^x5                       <= 1
#                   x0 + x1 - x2                   = 0 
#                 - x0 - x1      - x3              = 0e+00           
#               0.5 x0                - x4         = 1.3862944
#                        x1                  - x5  = 0               
#              x0 ... x5 are unrestricted
#
##

using Mosek

numvar = 6 
numcon = 5 

bkc = [MSK_BK_UP,
       MSK_BK_FX,
       MSK_BK_FX,
       MSK_BK_FX,
       MSK_BK_FX]
blc = [ 0.0, 0.0, 0.0, 1.3862944, 0.0 ] 
buc = [ 1.0, 0.0, 0.0, 1.3862944, 0.0 ] 

bkx = fill(MSK_BK_FR, numvar)
blx = zeros(numvar) 
bux = zeros(numvar)

aptrb = [ 0, 0, 3, 6, 8 ] 
aptre = [ 0, 3, 6, 8, 10 ] 
asubi = [ 0, 1, 2, 3, 4 ] 
asubj = [ 0, 1, 2, 
          0, 1, 3, 
          0, 4, 
          1, 5 ]  
aval  = [  1.0,  1.0, -1.0, 
          -1.0, -1.0, -1.0, 
           0.5, -1.0, 
           1.0, -1.0 ]

printstream(msg::String) = print(msg)
task = maketask()
putstreamfunc(task,MSK_STREAM_LOG,printstream)

appendcons(task,numcon)
appendvars(task,numvar)
putobjsense(task,MSK_OBJECTIVE_SENSE_MINIMIZE)
putvarboundslice(task,1,numvar+1,bkx,blx,bux)
putconboundslice(task,1,numcon+1,bkc,blc,buc)

putarowlist(task,asubi+1,aptrb+1,aptre+1,asubj+1,aval)

opro  = [ MSK_OPR_EXP, MSK_OPR_EXP ] 
oprjo = Int32[ 2, 3 ] 
oprfo = [ 1.0, 1.0 ] 
oprgo = [ 1.0, 1.0 ] 
oprho = [ 0.0, 0.0 ] 


oprc  = [ MSK_OPR_EXP, MSK_OPR_EXP ] 
opric = Int32[ 0, 0 ] 
oprjc = Int32[ 4, 5 ] 
oprfc = [ 1.0, 1.0 ] 
oprgc = [ 1.0, 1.0 ] 
oprhc = [ 0.0, 0.0 ]


sch = scbegin(task, opro, oprjo, oprfo, oprgo, oprho, oprc,
              opric, oprjc, oprfc, oprgc, oprhc)
optimize(task)
scend(task, sch)
res = getsolutionslice(task,MSK_SOL_ITR,MSK_SOL_ITEM_XX,1,numvar)
println(res)