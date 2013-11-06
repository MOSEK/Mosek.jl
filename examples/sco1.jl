# Solve the problem
# 
#  minimize  t^2.3
#  such that t > x+2y
#            y < sqrt(x)*sqrt(10)
#            y > x^2-10x+27        
#            0 < x,y < 10
#            0 < t            
#  
# This translates to minimizing (x+2y)^2.3 in the region marked with "O"s below.
# The variables are bounded by a 10x10 rectangle starting at origo. 
#  
#  10 ++------------+#------------+-------------+------------#+---------*****
#     +             + #           +             +    x**0.5*sqrt(10)******* +
#     |               #                                 x**2-10*x+27 ###### |
#     |                #                             ******#                |
#   8 ++                #                      ******OOOOO#                ++
#     |                  #               ******OOOOOOOOOO#                  |
#     |                   #         ******OOOOOOOOOOOOOO#                   |
#     |                    #    *****OOOOOOOOOOOOOOOOOO#                    |
#   6 ++                   #****OOOOOOOOOOOOOOOOOOOOOOO#                   ++
#     |                 ****##OOOOOOOOOOOOOOOOOOOOOOO##                     |
#     |             ****     #OOOOOOOOOOOOOOOOOOOOOOO#                      |
#   4 ++         ***          ##OOOOOOOOOOOOOOOOOOO##                      ++
#     |       ****             ##OOOOOOOOOOOOOOOOO##                        |
#     |     ***                  ##OOOOOOOOOOOOO##                          |
#     |   ***                      ##OOOOOOOOO##                            |
#   2 ++ **                          #########                             ++
#     |**                                                                   |
#     |*                                                                    |
#     +*            +             +             +             +             +
#   0 *+------------+-------------+-------------+-------------+------------++
#     0             2             4             6             8             10
#
# NOTE: MOSEK can solve general convex problems (convex objective, convex constraints):
#   min f(x)
#   s.t.  g(x) < b
#         lbx < x < lby
# where f and g are convex functions. That is, each function must by convex on 
# the whole of its domain, and the domain is defined by the bound on the 
# variables it uses. This means that if
#   f(x) = - sqrt(x) + 0.5x - 4 
# then the variable bounds on x must be: 0 < x < b (for 0 < b < +inf). I particular, 
# adding a _constraint_ saying that 0 < x will _not_ change define the correct domain
# of x. Furethermore, MOSEK will not always be able to detect non-convexity, only if it 
# enters a region that produces a hessian that is not positive semi-definite.
# 
# NOTE: While it is possible to implement the whole f and g as non-linear functions,
# it much more efficient to split out the parts that MOSEK explicitly supports (the linear and
# the quadratic parts). The linear and quadratic ports can be inputted as for normal linear and
# quadratic problems.
#
# In the example below we split out the linear parts of f and g, so the problem we solve is:
#   min  f(t)
#   s.t.     x+2y -t < 0
#        0 < -y + g_1(x)
#            x^2 - 10x - y < -27
#   and
#        f(t) = t^2.3
#        g_1(x) = sqrt(10)*sqrt(x)

using Mosek

t = maketask()
putstreamfunc(t,MSK_STREAM_LOG,msg -> print(msg))

appendvars(t,3)
appendcons(t,3)

# define domains of x,y and t
putvarbound(t,1,MSK_BK_RA,0.0, 10.0)  # x
putvarbound(t,2,MSK_BK_RA,0.0, 10.0)  # y
putvarbound(t,3,MSK_BK_LO,0.0, +Inf)  # t

# x+2y-t < 0
putarow(t,1,[1,2,3],[1.0, 2.0, -1.0]) # x + 2y - t
putconbound(t,1,MSK_BK_UP,-Inf, 0.0)  # < 0.0

# 0 < sqrt(x)*sqrt(10) - y
putarow(t,2,[2],[-1.0])               # - y
putconbound(t,2,MSK_BK_LO,0.0,+Inf)   # 0.0 <

# x^2 -10x - y < -27
putarow(t,3,[1,2],[-10.0, -1.0])      # -10x - y
putqconk(t,3,[1],[1],[2.0])           # x^2 
putconbound(t,3,MSK_BK_UP,-Inf,-27.0) # < -27.0

#define the evaluation functions
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
  val[1] = yo * 2.3 * x[3] ^ 1.3
  
  k = indexof(2,subi)
  val[1] += yc[k] * 0.5 * x[1] ^ (-0.5)
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

# Print a summary containing information
# about the solution for debugging purposes
solutionsummary(t,MSK_STREAM_MSG)

# Get status information about the solution
solsta = getsolsta(t,MSK_SOL_ITR)

if solsta in     [ MSK_SOL_STA_OPTIMAL, 
                   MSK_SOL_STA_NEAR_OPTIMAL ]
  xx = getxx(t,MSK_SOL_ITR)
  print("Optimal solution:")
  print(xx)
elseif solsta in [ MSK_SOL_STA_DUAL_INFEAS_CER,
                   MSK_SOL_STA_PRIM_INFEAS_CER,
                   MSK_SOL_STA_NEAR_DUAL_INFEAS_CER,
                   MSK_SOL_STA_NEAR_PRIM_INFEAS_CER ]
  print("Primal or dual infeasibility certificate found.\n")
elseif solsta == MSK_SOL_STA_UNKNOWN
  print("Unknown solution status")
else
  @printf("Other solution status (%d)\n",solsta)
end
