using Test, Mosek

@testset "Optimizer" begin
    err = ErrorException(
        "To use Mosek with JuMP (or MathOptInterface), you need to use the " *
        "package `MosekTools` (via `using MosekTools`). You may need to first" *
        " install it via `import Pkg; Pkg.add(\"MosekTools\")`.")
    @test_throws err Mosek.Optimizer()
end
