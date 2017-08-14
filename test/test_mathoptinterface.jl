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

include(joinpath(Pkg.dir("MathOptInterface"),"test","contconic.jl"))
@testset "Continuous conic problems" begin
    contconictest(MosekSolver())
end

include(joinpath(Pkg.dir("MathOptInterface"),"test","intlinear.jl"))
@testset "Mixed-integer linear problems" begin
    intlineartest(MosekSolver())
end
