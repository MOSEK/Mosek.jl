
#
# Copyright : Copyright (c) MOSEK ApS, Denmark. All rights reserved.
#
# File :      midual1.jl
#
#    Purpose:  Demonstrates how to compute dual values with
#              respect to a fixed integer solution for a MIO problem.
#
#              minimize    500 s1 + 300 s2 + 10 x1 + 14 x2
#              subject to  x1 + x2 >= 100
#                          0 <= x1 <= 70 s1
#                          0 <= x2 <= 80 s2
#                          s1, s2 - binary

using Mosek

midual1_ptf = """Task midual1
Objective
    Minimize + 10 x1 + 14 x2 + 500 s1 + 300 s2
Constraints
    'demand' [100] + x1 + x2
    'production1' [-inf;0] + x1 - 70 s1
    'production2' [-inf;0] + x2 - 80 s2
Variables
    x1 [0;+inf]
    x2 [0;+inf]
    s1 [0;1]
    s2 [0;1]
Integers
    s1 s2
"""

# The original mixed-integer problem
task = maketask()
readptfstring(task, midual1_ptf)

optimize(task)

if getprosta(task, MSK_SOL_ITG) != MSK_PRO_STA_PRIM_FEAS
  println("Unsuitable problem status, exiting")
  exit(-1)
end

xx = getxxslice(task, MSK_SOL_ITG, 1, 3)
println("x = $(xx[1]), $(xx[2])")

# Formulate the continuous fixed problem
fixTask = getfixedproblem(task)
optimize(fixTask)

if getprosta(fixTask, MSK_SOL_BAS) != MSK_PRO_STA_PRIM_AND_DUAL_FEAS
  println("Unsuitable problem status, exiting")
  exit(-1)
end

xfix = getxxslice(fixTask, MSK_SOL_BAS, 1, 3)
y = getyslice(fixTask, MSK_SOL_BAS, 1, 4)

println("xfix = $(xfix[1]), $(xfix[2])")
println("demand dual = $(y[1])")
println("production dual = $(y[2]), $(y[3])")
