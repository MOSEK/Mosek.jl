module TestMathOptInterface

using MathOptInterface, Base.Test, Mosek.MathOptInterfaceMosek
# TODO using solvers?

const MOI = MathOptInterface

include(joinpath(Pkg.dir("MathOptInterface"),"test","contlinear.jl"))
@testset "Continuous linear problems" begin
    contlineartest(MathOptMosekSolver(QUIET = true))
end

# include("contquadratic.jl")
# @testset "Continuous quadratic problems" begin
#     # contquadratictest(GurobiSolver())
# end

include(joinpath(Pkg.dir("MathOptInterface"),"test","contconic.jl"))
@testset "Continuous conic problems" begin
    contconictest(MathOptMosekSolver(QUIET = true))
end

include(joinpath(Pkg.dir("MathOptInterface"),"test","intlinear.jl"))
@testset "Mixed-integer linear problems" begin
    intlineartest(MathOptMosekSolver(QUIET = true))
end

end
