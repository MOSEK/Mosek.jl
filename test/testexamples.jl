@testset "[distro examples]" begin
    # @testset "lo1" begin
    #     include("../examples/lo1.jl")
    # end
    # @testset "cqo1" begin
    #     include("../examples/cqo1.jl")
    # end
    # @testset "qo1" begin
    #     include("../examples/qo1.jl")
    # end
    # @testset "qcqo1" begin
    #     include("../examples/qcqo1.jl")
    # end
    # @testset "milo1" begin
    #     include("../examples/milo1.jl")
    # end
    # @testset "sdo1" begin
    #     include("../examples/sdo1.jl")
    # end
    # @testset "eo1" begin
    #     include("../examples/eo1.jl")
    # end
    # @testset "pow1" begin
    #     include("../examples/pow1.jl")
    # end
    # @testset "solvebasis" begin
    #     include("../examples/solvebasis.jl")
    # end
    # @testset "portfolio" begin
    #     include("../examples/portfolio.jl")
    # end
    # @testset "production" begin
    #     include("../examples/production.jl")
    # end
    # @testset "callback" begin
    #     include("../examples/callback.jl")

    # end

    @testset "[distro examples]" begin
        @testset "acc1" begin include("../examples/acc1.jl") end
        @testset "acc2" begin include("../examples/acc2.jl") end
        @testset "callback" begin include("../examples/callback.jl") end
        @testset "ceo1" begin include("../examples/ceo1.jl") end
        @testset "cqo1" begin include("../examples/cqo1.jl") end
        @testset "djc1" begin include("../examples/djc1.jl") end
        @testset "feasrepairex1" begin let ARGS = ["-"]; include("../examples/feasrepairex1.jl") end end
        @testset "gp1" begin include("../examples/gp1.jl") end
        @testset "helloworld" begin include("../examples/helloworld.jl") end
        @testset "lo1" begin include("../examples/lo1.jl") end
        @testset "lo2" begin include("../examples/lo2.jl") end
        @testset "logistics" begin include("../examples/logistics.jl") end
        @testset "mico1" begin include("../examples/mico1.jl") end
        @testset "milo1" begin include("../examples/milo1.jl") end
        @testset "mioinitsol" begin include("../examples/mioinitsol.jl") end
        #@testset "opt_server_async" begin include("../examples/opt_server_async.jl") end
        #@testset "opt_server_sync" begin include("../examples/opt_server_sync.jl") end
        @testset "parallel" begin include("../examples/parallel.jl") end
        @testset "parameters" begin include("../examples/parameters.jl") end
        @testset "portfolio_1_basic" begin include("../examples/portfolio_1_basic.jl") end
        @testset "portfolio_2_frontier" begin include("../examples/portfolio_2_frontier.jl") end
        @testset "portfolio_3_impact" begin include("../examples/portfolio_3_impact.jl") end
        @testset "portfolio_4_transcost" begin include("../examples/portfolio_4_transcost.jl") end
        @testset "portfolio_5_card" begin include("../examples/portfolio_5_card.jl") end
        @testset "portfolio_6_factor" begin include("../examples/portfolio_6_factor.jl") end
        @testset "portfolio_data" begin include("../examples/portfolio_data.jl") end
        @testset "pow1" begin include("../examples/pow1.jl") end
        @testset "production" begin include("../examples/production.jl") end
        @testset "qcqo1" begin include("../examples/qcqo1.jl") end
        @testset "qo1" begin include("../examples/qo1.jl") end
        @testset "reoptimization" begin include("../examples/reoptimization.jl") end
        @testset "sdo1" begin include("../examples/sdo1.jl") end
        @testset "sdo2" begin include("../examples/sdo2.jl") end
        #@testset "sdo_lmi" begin include("../examples/sdo_lmi.jl") end
        @testset "sensitivity" begin include("../examples/sensitivity.jl") end
        @testset "simple" begin include("../examples/simple.jl") end
        @testset "solutionquality" begin include("../examples/solutionquality.jl") end
        @testset "solvebasis" begin include("../examples/solvebasis.jl") end
        @testset "solvelinear" begin include("../examples/solvelinear.jl") end
        #@testset "concurrent1" begin include("../examples/concurrent1.jl") end
    end
end
