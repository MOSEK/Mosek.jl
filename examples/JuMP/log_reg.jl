# Example of logistic regration with Mosek, MosekTools, and JuMP
# See Mosek manual 5.4.4 Logistic regression

using JuMP, Mosek, MosekTools
using RDatasets, Plots

# Iris dataset:  Our goal is to predict whether the iris species is versicolor
# using the sepal length and width and petal length and width.
iris = dataset("datasets", "iris");
Y = [species == "versicolor" ? 1.0 : -1.0 for species in iris.Species];
# Create feature matrix with one column for each feature
X = hcat(ones(size(iris, 1)), iris.SepalLength, iris.SepalWidth, iris.PetalLength, iris.PetalWidth);

# model size
n, p = size(X)
# regularizer
λ = 1e-2
# create model
model = Model(Mosek.Optimizer)
# variables
@variable(model, β[1:p])
@variable(model, β0)
# auxiliary variables for conic constraints
@variable(model, t[1:n])
@variable(model, r)
@variable(model, u̅[1:n])
@variable(model, u̲[1:n])
@variable(model, v̅[1:n])
@variable(model, v̲[1:n])

# minimize negative likelihood with a regularizer
@objective(model, Min, sum(t) + λ*r)
# conic constraints
for i in 1:n
    if Y[i] ==  1
        @constraint(model, u̅[i] + v̅[i] <= 1)
        @constraint(model, [-β0-X[i,:]'*β-t[i],1,u̅[i]] in MOI.ExponentialCone())
        @constraint(model, [-t[i],1,v̅[i]] in MOI.ExponentialCone())
    else
        @constraint(model, u̲[i] + v̲[i] <= 1)
        @constraint(model, [β0+X[i,:]'*β-t[i],1,u̲[i]] in MOI.ExponentialCone())
        @constraint(model, [-t[i],1,v̲[i]] in MOI.ExponentialCone())
    end
end
@constraint(model, [r; β] in SecondOrderCone())

# optimize
optimize!(model)
status=termination_status(model)
_β=JuMP.value.(β)
_β0=JuMP.value.(β0)
_t=sort!(JuMP.value.(t))

σ = 1 ./ (1 .+ exp.(-_β0.-X * _β))
if minimum(Y) == -1
  ŷ = (σ .> 0.5) .+ 2 .- 3 # bring the prediction in [-1, 1] range
else
  ŷ = (σ .> 0.5)
end
acc = count(Y .== ŷ) / n
print("accuracy: ", acc)

perm = sortperm(vec(X * _β));
plot(1:n, (Y[perm] .+ 1)/2, st=:scatter)
plot!(1:n, σ[perm])
