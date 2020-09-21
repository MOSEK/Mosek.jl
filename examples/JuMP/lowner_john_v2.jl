# Refined Löwner-John ellipsoid example with Mosek, MosekTools, and JuMP

# Computes the inner and outer Löwner-John ellipsoids around the feasible
# space of the JuMP model.

# Uses Polyhedra and CDDLib packages to express the polyhedral set of the JuMP
# model by inequalities and a convex hull of a finite set of points

using JuMP, Mosek, MosekTools
using LinearAlgebra, Polyhedra, CDDLib
using Plots

function lowner_john_inner_ellipsoid(A,b)
    """
    Löwner-John inner ellipsoid:
        the function takes H-representation (A,b) of a polytope and
        outputs the parameters of the maximum volume inscribed ellipsoid.

        polyhedral set S:
        S = { x ∈ Rn | ai' x ⩽ bi, i= 1,...,m }

        Proposition 3.7.1 [BenTalN01]
            ... Then the largest volume ellipsoid contained in S is
            E = { x = Z∗u + z∗ | uTu ≤1 }, where (Z∗,z∗) are from

            maximize    t                                           (1.0)
                s.t.    t ⩽ (DetZ)^1/n                              (1.1)
                        Z - PSD                                     (1.2)
                        || Z * ai ||_2 ≤ bi −ai' * z, i = 1,...,m   (1.3)
                        Z ∈ S^n, z ∈ Rn, t ∈ R                      (1.4)

        To model (1.1), refer to [Mosek, 6.2.3 Log-determinant]

            t ⩽ (DetZ)^1/n <==> [Z B; B' diag(B)] - PSD             (2.1)
                                B is lower triangular               (2.2)
                                t ⩽ (Π_i={1,...,n} Bii)^1/n         (2.3)

            where (2.3) maximizes the geometric mean of the Bii variables. It
            can be modelled using power cone [Mosek, 4.2.4 Geometric mean] and
            MOI.GeometricMeanCone available with the MathOptInterface.jl

    [BenTalN01]	    A. Ben-Tal and A. Nemirovski. Lectures on Modern Convex
                    Optimization: Analysis, Algorithms, and Engineering Applications.
                    MPS/SIAM Series on Optimization. SIAM, 2001.
    [Mosek]         https://docs.mosek.com/modeling-cookbook/index.html
    """
    # model size
    m,n = size(A)
    # model declaraton
    model = Model(optimizer_with_attributes(Mosek.Optimizer, "LOG" => 0))
    # variables
    @variable(model, t)
    @variable(model, Z[1:n,1:n], PSD)
    @variable(model, z[1:n])
    @variable(model, B[1:n,1:n])
    # maximize the volume
    @objective(model, Max, t)
    # SOC constraint (1.3)
    @constraint(model, SOC[i=1:m], [ b[i] - (A[i,:]' * z) ; (A[i,:]'*Z)'] in SecondOrderCone())
    # PSD constraint (2.1)
    @constraint(model, [Z LowerTriangular(B); LowerTriangular(B)' diagm(diag(B))] in PSDCone())
    # Power cone costraint (2.3)
    @constraint(model, [t; vec([B[i,i] for i in 1:n])] in MOI.GeometricMeanCone(n+1))
    # optimize
    optimize!(model)
    status=termination_status(model)
    CPU_time = MOI.get(model, MOI.SolveTime())
    # return solution
    return JuMP.value.(Z), JuMP.value.(z), CPU_time
end

function lownerjohn_outer(x)
    """
    Löwner-John outer ellipsoid:
        the function takes V-representation (x) of a polytope and
        outputs the parameters of the maximum volume inscribed ellipsoid.

        Polyhedral set S = Conv{x_1,...,x_n}

        Uses the same Mosek features as the inner ellipsoid function,
        refer to Proposition 3.7.2 [BenTalN01] for further details.
    """

    m,n = size(x)
    model = Model(optimizer_with_attributes(Mosek.Optimizer, "LOG" => 0))

    @variable(model,t >= 0)
    @variable(model,P[1:n,1:n],PSD)
    @variable(model,B[1:n,1:n])
    @variable(model,c[1:n])

    # (1, Px-c) in cone
    for i in 1:m
        @constraint(model,[1 ; (x[i,:]' * P)' - c] in MOI.SecondOrderCone(n+1))
    end

    @constraint(model, [P LowerTriangular(B); LowerTriangular(B)' diagm(diag(B))] in PSDCone())
    @constraint(model, [t; vec([B[i,i] for i in 1:n])] in MOI.GeometricMeanCone(n+1))

    @objective(model,Max,t)

    optimize!(model)
    status=termination_status(model)
    CPU_time = MOI.get(model, MOI.SolveTime())

    JuMP.value.(P), JuMP.value.(c), CPU_time
end

# create a JuMP model
function JuMP_model(n)
    model = Model(optimizer_with_attributes(Mosek.Optimizer, "LOG" => 0))
    @variable(model, x[1:n])
    @objective(model, Min, x'diagm(ones(n))'x)
    @constraint(model, con1[i=1:n], 0 <= x[i] <= 1)
    @constraint(model, ones(n)'x >= n/3)
    optimize!(model)
    return model
end

# model dimension
n = 2
# get JuMP model
model = JuMP_model(n)
@info("done exporting JuMP model")
# create a polyhedron from the feasible set of a JuMP model
poly = polyhedron(model, CDDLib.Library(:exact))
@info("done exporting the polyhedron")

# obtain H-representation of the polyhedron
hr = MixedMatHRep(hrep(poly))
A = hr.A
b = hr.b
@info("done obtaining H-representation")
# solve ellipsoid problem
Ci, di, CPU_time = lowner_john_inner_ellipsoid(A,b)
@info("done solving inner ellipsoid problem in $(CPU_time) seconds")

# obtain V-representation of the polyhedron
vr = vrep(poly)
vr = MixedMatVRep(vr)
p = vr.V
@info("done obtaining V-representation")
Po, co, CPU_time = lownerjohn_outer(p)
@info("done solving outter ellipsoid problem in $(CPU_time) seconds")
Poinv = Po^-1

# make a plot
if n == 2
    println("Plotting solution...")
    Plots.gr()
    Plots.plot([p[:,1] ; p[1,1]], [p[:,2] ; p[1,2]])
    Plots.plot!(t -> Ci[1,:]' * [ cos(t),sin(t) ] + di[1],
                t -> Ci[2,:]' * [ cos(t),sin(t) ] + di[2],
                0,2*pi)
    Plots.plot!(t -> Poinv[1,:]' * ( [cos(t),sin(t)] .+ co[1] ),
                t -> Poinv[2,:]' * ( [cos(t),sin(t)] .+ co[2] ),
                0,2*pi)
end
