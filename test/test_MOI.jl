using MathOptInterface, Base.Test,Mosek
# TODO using solvers?

const MOI = MathOptInterface

include(joinpath(Pkg.dir("MathOptInterface"),"test","contlinear.jl"))
@testset "Continuous linear problems" begin
    contlineartest(MosekSolver())
end

# include("contquadratic.jl")
# @testset "Continuous quadratic problems" begin
#     # contquadratictest(GurobiSolver())
# end

# include("contconic.jl")
# @testset "Continuous conic problems" begin
#     # contconictest(SCSSolver(verbose=0))
# end

# include("intlinear.jl")
# @testset "Mixed-integer linear problems" begin
#     # intlineartest(GLPKSolverMIP())
# end
