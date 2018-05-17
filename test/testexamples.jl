@testset "[distro examples]" begin
    @testset "lo1" begin
        include("../examples/lo1.jl")
    end
    @testset "cqo1" begin
        include("../examples/cqo1.jl")
    end
    @testset "qo1" begin
        include("../examples/qo1.jl")
    end
    @testset "qcqo1" begin
        include("../examples/qcqo1.jl")
    end
    @testset "milo1" begin
        include("../examples/milo1.jl")
    end
    @testset "sdo1" begin
        include("../examples/sdo1.jl")
    end
    @testset "eo1" begin
        include("../examples/eo1.jl")
    end
    @testset "pow1" begin
        include("../examples/pow1.jl")
    end
    @testset "solvebasis" begin
        include("../examples/solvebasis.jl")
    end
    @testset "portfolio" begin
        include("../examples/portfolio.jl")
    end
    @testset "production" begin
        include("../examples/production.jl")
    end
    @testset "callback" begin
        include("../examples/callback.jl")

    end
    @testset "feasrepairex1" begin
        push!(ARGS, "$(dirname(dirname(@__FILE__)))/examples/feasrepair.lp")
        include("../examples/feasrepairex1.jl")

    end
end
