using Mosek, Mosek.Ext, Test
import Pkg

@testset "[show]" begin
    for fname in [ "25fv47.task", "lj-inner.task"]
        maketask(filename=joinpath(Pkg.dir("Mosek"),"test",fname)) do t
            open(joinpath(Pkg.dir("Mosek"),"test",join("show-",fname)),"w") do f
                println(f,t)

                n = getnumvar(t)
                m = getnumcon(t)
                numcone = getnumcone(t)
                numbarvar = getnumbarvar(t)

                if n > 0
                    for j in 1:n
                        println(f,t[Var(j)])
                    end
                end

                if m > 0
                    for j in 1:m
                        println(f,t[Con(j)])
                    end
                end

                if numcone > 0
                    for j in 1:numcone
                        println(f,t[Cone(j)])
                    end
                end

                if numbarvar > 0
                    for j in 1:getnumbarvar(t)
                        println(f,t[Barvar(j)])
                    end
                end
                
                for which in Soltype
                    if solutiondef(t,which)
                        sol = t[Sol(which)]

                        println(f,sol)

                        if n > 0
                            for j in 1:n
                                println(f,sol[Var(j)])
                            end
                        end

                        if m > 0
                            for j in 1:m
                                println(f,sol[Con(j)])
                            end
                        end

                        if numbarvar > 0
                            for j in 1:getnumbarvar(t)
                                println(f,sol[Barvar(j)])
                            end
                        end
                    end
                end    
            end
        end
    end
end
