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
        @testset "acc1" begin include("acc1.jl") end
        @testset "acc2" begin include("acc2.jl") end
        #@testset "callback" begin include("callback.jl") end
        @testset "ceo1" begin include("ceo1.jl") end
        @testset "concurrent1" begin include("concurrent1.jl") end
        @testset "cqo1" begin include("cqo1.jl") end
        @testset "djc1" begin include("djc1.jl") end
        @testset "feasrepairex1" begin include("feasrepairex1.jl") end
        @testset "gp1" begin include("gp1.jl") end
        @testset "helloworld" begin include("helloworld.jl") end
        @testset "lo1" begin include("lo1.jl") end
        @testset "lo2" begin include("lo2.jl") end
        @testset "logistics" begin include("logistics.jl") end
        @testset "mico1" begin include("mico1.jl") end
        @testset "milo1" begin include("milo1.jl") end
        @testset "mioinitsol" begin include("mioinitsol.jl") end
        #@testset "opt_server_async" begin include("opt_server_async.jl") end
        #@testset "opt_server_sync" begin include("opt_server_sync.jl") end
        @testset "parallel" begin include("parallel.jl") end
        @testset "parameters" begin include("parameters.jl") end
        @testset "portfolio_1_basic" begin include("portfolio_1_basic.jl") end
        @testset "portfolio_2_frontier" begin include("portfolio_2_frontier.jl") end
        @testset "portfolio_3_impact" begin include("portfolio_3_impact.jl") end
        @testset "portfolio_4_transcost" begin include("portfolio_4_transcost.jl") end
        @testset "portfolio_5_card" begin include("portfolio_5_card.jl") end
        @testset "portfolio_6_factor" begin include("portfolio_6_factor.jl") end
        @testset "portfolio_data" begin include("portfolio_data.jl") end
        @testset "pow1" begin include("pow1.jl") end
        @testset "production" begin include("production.jl") end
        @testset "qcqo1" begin include("qcqo1.jl") end
        @testset "qo1" begin include("qo1.jl") end
        @testset "reoptimization" begin include("reoptimization.jl") end
        @testset "sdo1" begin include("sdo1.jl") end
        @testset "sdo2" begin include("sdo2.jl") end
        @testset "sdo_lmi" begin include("sdo_lmi.jl") end
        @testset "sensitivity" begin include("sensitivity.jl") end
        @testset "simple" begin include("simple.jl") end
        @testset "solutionquality" begin include("solutionquality.jl") end
        @testset "solvebasis" begin include("solvebasis.jl") end
        @testset "solvelinear" begin include("solvelinear.jl") end
    end

end
