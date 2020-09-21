# Example of the piecewise linear regration with Mosek, MosekTools, and JuMP
# See Mosek manual 9.2.3 Convex piecewise linear regression

using JuMP, Mosek, MosekTools
using RDatasets, Plots

# prepare data
data = dataset("datasets", "pressure");
Y = data.Temperature ./1000;    Y = Y[8:end]
X = data.Pressure ./1000;       X = X[8:end]
n = length(X)
k = 3
(M_lb, M_ub) = (-5, 5)

# build and solve the model
m = Model(Mosek.Optimizer)
@variable(m, t)
@variable(m, M_lb <= f[1:n] <= M_ub)
@variable(m, a[1:k])
@variable(m, b[1:k])
@variable(m, z[1:n,1:k], Bin)
@objective(m, Min,  t)
@constraint(m, [t;Y - f] in SecondOrderCone())
@constraint(m, cut_one[i=1:n,j=1:k], (a[j]*X[i]+b[j]) - f[i] <= M_ub * (1 - z[i,j]))
@constraint(m, cut_two[i=1:n,j=1:k], (a[j]*X[i]+b[j]) - f[i] >= M_lb * (1 - z[i,j]))
@constraint(m, sum_bin[i=1:n], sum(z[i,:]) == 1)
optimize!(m)
status=termination_status(model)

@show a_ = JuMP.value.(a)
@show b_ = JuMP.value.(b)

# make a plot
Plots.gr(ylims=(0,1))
Plots.scatter(X,Y,legend=false)
for j in 1:k
  f_z(x) = a_[j]*x + b_[j]
  Plots.plot!(X, f_z)
end
current()
