module MosekExtraMathProgTests

using MathProgBase
using MathProgBase.SolverInterface
using Mosek
using Base.Test


solver = MosekSolver(QUIET=true)


@testset "[mathprogextra]" begin

    # @testset "addelmtest" begin
    #     objtol = 1e-7
    #     primaltol = 1e-6
    #     m = MathProgBase.model(Mosek.MosekSolver())

    #     loadproblem!(m,
    #                  sparse([1,1,1],[1,2,3],[1.0,1.0,1.0]), # A
    #                  Float64[-1.0,-1.0,-1.0],Float64[1.0,1.0,1.0], # lbx,ubx
    #                  zeros(Float64,3), # c
    #                  Float64[1.0], Float64[1.0], # lbc, ubc
    #                  :Max)

    #     j0  = MathProgBase.addvar!(m, 0.0, Inf, 0.0)
    #     i0  = MathProgBase.addconstr!(m, [4], [1.0], 1.0, 1.0)
    #     k0 = MathProgBase.addquadconstr!(m, [], [], [1,2,3], [1,2,3], [-1.0,1.0,1.0], '<', 0.0)

    #     j1  = MathProgBase.addvar!(m, 0.0, Inf, 0.0)
    #     i1  = MathProgBase.addconstr!(m, [5], [1.0], 1.0, 1.0)
    #     k1 = MathProgBase.addquadconstr!(m, [], [], [1,3], [2,3], [-1.0,1.0], '<', 0.0)

    #     println("j0 = ",j0,", i1 = ",i0,", k0 = ",k0)
    #     println("j1 = ",j1,", i1 = ",i1,", k1 = ",k1)

    #     Mosek.writedata(m.task,"mathprogtestextra.opf")
    #     println("End")

    #     #MathProgBase.optimize!(m)
    #     #stat=MathProgBase.status(m) # = :Unbounded
    # end

    @testset "dualsigntest" begin
        duals=true
        # Problem 4 - lo1 from MOSEK docs
        # Property: All duals are non-zero
        # min      [ 3 1 5 1 ] * [ x y z w ]
        #  st [35]-[ 3 1 2   ]               ZERO
        #     [15]-[ 2 1 3 1 ] * [ x y z w ] NONPOS
        #     [25]-[   2   3 ]               NONNEG
        #     [10]-[   1     ]               NONPOS
        #     x,y,z,w > 0
        # Test: Signs of duals for ZERO, NONNEG, NONPOS
        println("Problem 4")
        duals = true
        if duals
            m = MathProgBase.ConicModel(solver)
            MathProgBase.loadproblem!(m,
                                      [ -3.0, -1.0, -5.0, -1.0 ],
                                      [  3.0   1.0   2.0   0.0 ;
                                       2.0   1.0   3.0   1.0 ;
                                       0.0   2.0   0.0   3.0 ;
                                       0.0   1.0   0.0   0.0 ],
                                      [ 30.0, 15.0, 25.0, 10.0 ],
            Any[(:Zero,1),(:NonPos,2),(:NonNeg,3),(:NonNeg,4)],
            [(:NonNeg,1:4)])
            MathProgBase.optimize!(m)
            @test MathProgBase.status(m) == :Optimal
            d = MathProgBase.getdual(m)
            # d[1] is free

            @test d[2] <= 1e-6
            @test d[3] >= -1e-6
            @test d[4] >= -1e-6
        end
    end

    # @testset "cqo1test1" begin
    #     s=MosekSolver()
    #     duals=true
    #     # Problem cqo1 from MOSEK docs
    #     # min  x4 + 0.5 * x5 + x6
    #     # s.t. x1 + x2 + 2 x3= 1
    #     #      x4^2  > x1^2 + x2^2
    #     #      x5*x6 > x3^2
    #     #      x1,x2,x3,x4,x5,x6 > 0

    #     # Input as conic, add conic constraints as quadratic constraints
    #     let m = MathProgBase.Conic(s)
    #         MathProgBase.loadproblem!(m,
    #                                   #  x1  x2  x3  x4  x5  x6
    #                                   [ 0.0 0.0 0.0 1.0 0.5 1.0 ], # c
    #                                   [ 1.0 1.0 2.0 0.0 0.0 0.0 ], # A
    #                                   [ 1.0 ], # b
    #                                   [ (:Zero,1) ],   # constr domain
    #                                   [ (:NonNeg,1:6) ]) # var    domain
    #         MathProgBase.addquadconstr!(m, [],[], Int32[ 4,1,2 ], Int32[4,1,2], Float64[-1.0, 1.0, 1.0], '<', 0.0 )
    #         MathProgBase.addquadconstr!(m, [],[], Int32[ 5,3   ], Int32[6,3  ], Float64[-1.0, 1.0     ], '<', 0.0 )

    #         MathProgBase.optimize!(m)

    #         @test MathProgBase.status(m) == :Optimal

    #         pobj = MathProgBase.getobjval(m)

    #         @test isapprox(pobj, 7.07106782e-01, atol=1e-6)

    #         xx = MathProgBase.getsolution(m)
    #         xc = MathProgBase.getconstrsolution(m)

    #         @test abs(pobj-dot(xx,[ 0.0,0.0,0.0,1.0,0.5,1.0 ])) < 1e-6
    #     end
    # end
    # @testset "cqo1test2" begin
    #     s=MosekSolver()
    #     duals=true

    #     # Input as linear, add conic constraints as quadratic constraints
    #     let m = MathProgBase.LinearQuadraticModel(s)
    #         MathProgBase.loadproblem!(m,
    #                                   [ 1.0 1.0 2.0 0.0 0.0 0.0 ], # A
    #                                   [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ], # blx
    #                                   [ Inf, Inf, Inf, Inf, Inf, Inf ], # bux
    #                                   [ 0.0, 0.0, 0.0, 1.0, 0.5, 1.0 ], # c
    #                                   [ 1.0 ], # blc
    #                                   [ 1.0 ], # buc
    #                                   :Min )
    #         MathProgBase.addquadconstr!(m, [],[], Int32[ 4,1,2 ], Int32[4,1,2], Float64[-1.0, 1.0, 1.0], '<', 0.0 )
    #         MathProgBase.addquadconstr!(m, [],[], Int32[ 5,3   ], Int32[6,3  ], Float64[-1.0, 1.0     ], '<', 0.0 )

    #         MathProgBase.optimize!(m)

    #         @test MathProgBase.status(m) == :Optimal

    #         pobj = MathProgBase.getobjval(m)

    #         @test isapprox(pobj, 7.07106782e-01, atol=1e-6)

    #         xx = MathProgBase.getsolution(m)
    #         xc = MathProgBase.getconstrsolution(m)
    #         @test abs(pobj-dot(xx,[ 0.0,0.0,0.0,1.0,0.5,1.0 ])) < 1e-6

    #     end
    # end

    # @testset "cqo1test4" begin
    #     s=MosekSolver()
    #     duals=true
    #     # Input as linear, add conic constraints as quadratic constraints
    #     let m = MathProgBase.model(s)
    #         MathProgBase.loadproblem!(m,
    #                                   [ 1.0 1.0 2.0 0.0 0.0 0.0 ], # A
    #                                   [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ], # blx
    #                                   [ Inf, Inf, Inf, Inf, Inf, Inf ], # bux
    #                                   [ 0.0, 0.0, 0.0, 1.0, 0.5, 1.0 ], # c
    #                                   [ 1.0 ], # blc
    #                                   [ 1.0 ], # buc
    #         :Min )
    #         MathProgBase.addquadconstr!(m, [],[], Int32[ 2,1,4 ], Int32[2,1,4], Float64[ 1.0, 1.0,-1.0], "<", 0.0 )
    #         MathProgBase.addquadconstr!(m, [],[], Int32[ 3,5   ], Int32[3,6  ], Float64[ 1.0,-1.0     ], "<", 0.0 )

    #         MathProgBase.optimize!(m)
            
    #         @test MathProgBase.status(m) == :Optimal

    #         pobj = MathProgBase.getobjval(m)
            
    #         @test isapprox(pobj, 7.07106782e-01, atol=1e-6)
            
    #         xx = MathProgBase.getsolution(m)
    #         xc = MathProgBase.getconstrsolution(m)

    #         @test abs(pobj-dot(xx,[ 0.0,0.0,0.0,1.0,0.5,1.0 ])) < 1e-6
    #     end
    # end
    
    
    @testset "qcqo1test" begin
        duals=true
        # Problem qcqo1 from MOSEK docs
        # min     - x2      + x1^2 -     x1*x3+ 0.1 x2^2 +     x3^2

        # s.t. x1 + x2 + x3 - x1^2 + 0.2 x1*x3 -    x2^2 - 0.1 x3^2 >= 1.0
        #      x1,x2,x3 >= 0
        
        m = MathProgBase.LinearQuadraticModel(solver)
        MathProgBase.loadproblem!(m,
                                  zeros(Float64,(0,3)), # A
                                  [ 0.0, 0.0, 0.0 ], # blx
                                  [ Inf, Inf, Inf ], # bux
                                  [ 0.0, -1.0, 0.0 ], # c
                                  Float64[], # blc
                                  Float64[], # buc
                                  :Min)
        MathProgBase.setquadobj!(m,[1,1,2,3],[1,3,2,3], [2.0,-1.0,0.2,2.0]) # specify lower triangular
        MathProgBase.addquadconstr!(m,  # specify quadratid _terms_ (not triangular)
                                    [1,2,3], [1.0,1.0,1.0], # linear part
                                    [1,1,2,3], [1,3,2,3], [-1.0,0.2,-1.0,-0.1],
                                    '>',
                                    1.0)
        MathProgBase.optimize!(m)
        
        @test MathProgBase.status(m) == :Optimal

        pobj = MathProgBase.getobjval(m)
        
        @test isapprox(pobj, -4.91766743e-01, atol=1e-6)
        
        xx = MathProgBase.getsolution(m)
        
        @test isapprox(pobj, -xx[2]+xx[1]^2-xx[1]*xx[3]+0.1*xx[2]^2+xx[3]^2, atol=1e-6)

        #xc = MathProgBase.getquadconstrsolution(m)
        #@test_approx_eq_eps xc[1] (xx[1]+xx[2]+xx[3] - xx[1]^2 + 0.2*xx[1]*xx[3] - xx[2]^2 - 0.1 * xx[3]^2) 1e-6
    end

    @testset "modifymodeltest" begin
        mmin = MathProgBase.LinearQuadraticModel(solver)
        mmax = MathProgBase.LinearQuadraticModel(solver)

        MathProgBase.loadproblem!(mmin,
                     [  1.0   0.0   0.0   0.0 ;
                        0.0   1.0   0.0   0.0 ;
                        0.0   0.0   1.0   0.0 ;
                        0.0   0.0   0.0   1.0 ],
                     [ -Inf, -Inf, -Inf, -Inf ], # blx
                     [  Inf,  Inf,  Inf,  Inf ], # bux
                     [  1.0,  1.0,  1.0,  1.0 ], # c
                     [  0.0,  0.0,  0.0,  0.0 ], # blc
                     [  0.0,  0.0,  0.0,  0.0 ], # buc
                     :Min)

        MathProgBase.loadproblem!(mmax,
                     [  1.0   0.0   0.0   0.0 ;
                        0.0   1.0   0.0   0.0 ;
                        0.0   0.0   1.0   0.0 ;
                        0.0   0.0   0.0   1.0 ],
                     [ -Inf, -Inf, -Inf, -Inf ], # blx
                     [  Inf,  Inf,  Inf,  Inf ], # bux
                     [  1.0,  1.0,  1.0,  1.0 ], # c
                     [  0.0,  0.0,  0.0,  0.0 ], # blc
                     [  0.0,  0.0,  0.0,  0.0 ], # buc
                     :Max)

        writeproblem(mmax,"testmax_orig.opf")
        # test: modify buc
        setconstrUB!(mmax,[ 1.0, 2.0, 3.0, 4.0 ]);
        MathProgBase.optimize!(mmax)
        @test MathProgBase.status(mmax) == :Optimal
        xx = MathProgBase.getsolution(mmax)

        @test isapprox(xx[1], 1.0, atol=1e-8)
        @test isapprox(xx[2], 2.0, atol=1e-8)
        @test isapprox(xx[3], 3.0, atol=1e-8)
        @test isapprox(xx[4], 4.0, atol=1e-8)

        writeproblem(mmax,"testmax.opf")
        #writeproblem(mmax,"testmin.task")

        # test: modify blc
        setconstrLB!(mmin,[ -1.0, -2.0, -3.0, -4.0 ]);
        MathProgBase.optimize!(mmin)
        @test MathProgBase.status(mmin) == :Optimal
        xx = MathProgBase.getsolution(mmin)

        @test isapprox(xx[1], -1.0, atol=1e-8)
        @test isapprox(xx[2], -2.0, atol=1e-8)
        @test isapprox(xx[3], -3.0, atol=1e-8)
        @test isapprox(xx[4], -4.0, atol=1e-8)

        # test: modify bux
        setvarUB!(mmax,[ 0.5, 1.5, 2.5, 3.5 ]);
        MathProgBase.optimize!(mmax)
        @test MathProgBase.status(mmax) == :Optimal
        xx = MathProgBase.getsolution(mmax)

        @test isapprox(xx[1], 0.5, atol=1e-8)
        @test isapprox(xx[2], 1.5, atol=1e-8)
        @test isapprox(xx[3], 2.5, atol=1e-8)
        @test isapprox(xx[4], 3.5, atol=1e-8)

        # test: modify blx
        setvarLB!(mmin,[ -0.5, -1.5, -2.5, -3.5 ]);
        MathProgBase.optimize!(mmin)
        @test MathProgBase.status(mmin) == :Optimal
        xx = MathProgBase.getsolution(mmin)

        @test isapprox(xx[1], -0.5, atol=1e-8)
        @test isapprox(xx[2], -1.5, atol=1e-8)
        @test isapprox(xx[3], -2.5, atol=1e-8)
        @test isapprox(xx[4], -3.5, atol=1e-8)
    end
end

end # module
