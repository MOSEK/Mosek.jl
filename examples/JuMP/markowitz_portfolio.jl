#
# In the following two examples V is the factored covariance matrix,
# i.e. Q = VV', μ is the vector of expected returns for a set of
# stocks, δ is the target expected return and S is a measure for the
# risk (specifically, the maximum variance).
#
# The interesting constraint here is the constraint t > x'Q x. We notice that
#   x' Q x = x'VVx = (x'V)^2 = || x' V ||^2
# And write
#   t > || x' V ||^2
# Or equivalently
#   (0.5, t, x' V) in RotatedSecondOrderCone
#
# Our initial total wealth is 1, thus sum(x)=1.
#

using JuMP
using MosekTools

# Find the position x with minimal risk such that we get at least thge
# expected return δ.
#
#
function portfolio_min_risk(solver, V :: Array{Float64,2}, μ :: Vector{Float64}, δ :: Float64)
    model = Model(with_optimizer(solver))
    n = size(μ,1)
    m = size(V,2)

    
    @variable(model, t >= 0)
    @variable(model, x[1:n] >= 0)

    @constraint(model, μ'x .>= δ)

    @constraint(model, sum(x) == 1.0)
    @constraint(model, [0.5; t; V' * x] in MOI.RotatedSecondOrderCone(m+2))

    @objective(model, Min, t)

    optimize!(model)

    value(t), [ value(item) for item in x ]
end

# Find the position x with maximum return such that the variance is at most σ.
function portfolio_max_return(solver, V :: Array{Float64,2}, μ :: Vector{Float64}, R :: Float64)
    model = Model(with_optimizer(solver))
    n = length(μ)
    m = size(V,2)
    
    @variable(model, t >= 0)
    @variable(model, x[1:n] >= 0)

    @constraint(model, sum(x) == 1.0)
    @constraint(model, [0.5; t; V' * x] in MOI.RotatedSecondOrderCone(m+2))
    @constraint(model, t <= R)

    @objective(model, Max, μ'x)

    optimize!(model)

    value(t), [ value(item) for item in x ]
end


# Find the position x that maximizes the Sharpe ratio S(x) =
# (μ'x-rf)/||Vx||. See MOSEK Cookbook 2018 section 3.3.4 for details.
#
# Here rf is the return of the risk-free asset.
#
# Note that in the implementation we minimize 1/S(x).
function portfolio_sharpe_ratio(solver, V :: Array{Float64,2}, μ :: Vector{Float64}, rf :: Float64)
    model = Model(with_optimizer(solver))
    n = length(μ)
    m = size(V,2)
    
    @variable(model, t >= 0)
    @variable(model, y[1:n] >= 0)
    @variable(model, z >= 0)

    @constraint(model, [t; V' * y] in MOI.SecondOrderCone(m+1))
    @constraint(model, sum(y) - z == 0.0)
    @constraint(model, (μ .- rf)'y == 1)

    @objective(model, Min, t)

    optimize!(model)

    zres = value(z)
    
    1.0/value(t), [ value(item)/zres for item in y ] 
end

#
# Find the position x with minimal risk, including a risk parity term,
# such that we get at least the expected return δ. See MOSEK Cookbook
# 2018 section 5.4.1 for details.
#
# Here c is the weight given to the parity term.
function portfolio_risk_parity(solver, V :: Array{Float64,2}, μ :: Vector{Float64}, δ :: Float64, c :: Float64)
    model = direct_model(solver())
    n = size(μ,1)
    m = size(V,2)
    
    @variable(model, t >= 0)
    @variable(model, x[1:n] >= 0)
    @variable(model, s[1:n] >= 0)

    @constraint(model, μ'x .>= δ)

    @constraint(model, sum(x) == 1.0)
    @constraint(model, [0.5; t; V' * x] in MOI.RotatedSecondOrderCone(m+2))
    for i in 1:n
        @constraint(model, [ s[i], 1.0, x[i] ] in MOI.ExponentialCone())
    end

    @objective(model, Min, t + c * sum(s))

    optimize!(model)

    value(t), [ value(item) for item in x ]
end








solver = Mosek.Optimizer
## Example data:

# Vector of expected returns
μ = [1.05, 1.3, 0.9, 1.0, 0.9, 1.5]

# lower bound on return of investment
δ = 1.1

# upper bound for standard deviation
σ = 0.018

# factor of covariance (we use the data-matrix)
V = [-0.0644  -0.0352  -0.0019   0.1015 ;
      0.1183  -0.0428  -0.0013  -0.0741 ;
     -0.0798  -0.1201   0.1040   0.0958 ;
     -0.0855  -0.1174   0.1042   0.0987 ;
     -0.0158  -0.0074   0.0110   0.0123 ;
      0.0636  -0.1341   0.0206   0.0497 ]


t1,x1 = portfolio_min_risk(solver,V,μ,δ)
t2,x2 = portfolio_max_return(solver,V,μ,σ)
t3,x3 = portfolio_sharpe_ratio(solver,V,μ,0.0)
t4,x4 = portfolio_risk_parity(solver,V,μ,δ,0.05)

println("Minimize risk with δ = $δ:")
println("  variance = $t1")
println("  x = $x1")
println("  expected return = $(μ'x1)")
println("")
println("Maximize expected return with variance <= $σ:")
println("  variance = $t2")
println("  x = $x2")
println("  expected return = $(μ'x2)")
println("")
println("Maximize Sharpe ratio:")
println("  Sharpe ratio = $t3")
println("  x = $x3")
println("  expected return = $(μ'x3)")
println("")
println("Minimize risk with parity:")
println("  risk = $t4")
println("  x = $x4")
println("  expected return = $(μ'x4)")
