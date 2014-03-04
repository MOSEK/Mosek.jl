#   Purpose:   Demonstrate how to solve the following
#              optimization problem
#
#              minimize    f(x) + c^T x
#              subject to  l^c <= a x <= u^c
#                          l^x <=  x  <= u^x
#              where f is a nonlinear (convex) function.
#
#              The program solves the example:
#
#              minimize    - x_0 - ln(x_1 + x_2)
#              subject to  x_0 + x_1 + x_2 = 1
#                          x_0, x_1, x_2 >= 0.
#


using Mosek

t = maketask()
putstreamfunc(t,MSK_STREAM_LOG,msg -> print(msg))

appendvars(t,3)
appendcons(t,1)

putvarbound(t,1,MSK_BK_LO,0.0, +Inf)  # x0
putvarbound(t,2,MSK_BK_LO,0.0, +Inf)  # x1
putvarbound(t,3,MSK_BK_LO,0.0, +Inf)  # x2

putarow(t,1,[1,2,3],[1.0,1.0,1.0])   # x0 + x1 + x2
putconbound(t,1,MSK_BK_FX,1.0, 1.0)  # = 1.0

putclist(t,[1], [-1.0])

#define the evaluation functions
function indexof(v,a)
  for i in 1:length(a)
    if a[i] == v return i
    end
  end
  return 0
end

evalobj(x :: Array{Float64,1}) = log(x[2]+x[3])

evalconi(x:: Array{Float64,1},i:: Int32) = 0.0

function grdobj(x :: Array{Float64,1},sub:: Array{Int32,1}, val:: Array{Float64,1})  
  for j in 1:length(sub)    
    if sub[j] == 2 || sub[j] == 3
      val[j] = -1.0 / (x[2]+x[3])
    end
  end
end
        
function grdconi(x  :: Array{Float64,1},
                 i  :: Int32, 
                 sub:: Array{Int32,1}, 
                 val:: Array{Float64,1})
  none
end

function grdlag(x ::   Array{Float64,1},
                yo::   Float64,
                yc::   Array{Float64,1},
                subi:: Array{Int32,1},
                val::  Array{Float64,1})
  val[2] = - yo * 1.0 / (x[2]+x[3])
  val[3] = - yo * 1.0 / (x[2]+x[3])
end 


# NOTE: in this case we have no off-diagonal nonzeros in the hessian.
# If we do, then _either_ element (i,j) _or_ element (j,i) should be 
# specified.
function heslag(x ::      Array{Float64,1},
                yo::      Float64,
                yc::      Array{Float64,1},
                subi::    Array{Int32,1},
                hessubi:: Array{Int32,1},
                hessubj:: Array{Int32,1},
                hesval::  Array{Float64,1})

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

# Print a summary containing information
# about the solution for debugging purposes
solutionsummary(t,MSK_STREAM_MSG)

# Get status information about the solution
solsta = getsolsta(t,MSK_SOL_ITR)

if solsta in     [ MSK_SOL_STA_OPTIMAL, 
                   MSK_SOL_STA_NEAR_OPTIMAL ]
  xx = getxx(t,MSK_SOL_ITR)
  print("Optimal solution:\n")
  show(xx)
  println()
elseif solsta in [ MSK_SOL_STA_DUAL_INFEAS_CER,
                   MSK_SOL_STA_PRIM_INFEAS_CER,
                   MSK_SOL_STA_NEAR_DUAL_INFEAS_CER,
                   MSK_SOL_STA_NEAR_PRIM_INFEAS_CER ]
  print("Primal or dual infeasibility certificate found.\n")
elseif solsta == MSK_SOL_STA_UNKNOWN
  print("Unknown solution status\n")
else
  @printf("Other solution status (%d)\n",solsta)
end

