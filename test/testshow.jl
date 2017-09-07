using Mosek, Mosek.Ext

maketask(filename=joinpath(Pkg.dir("Mosek"),"test","25fv47.opf")) do t
    println(t)

    n = getnumvar(t)
    m = getnumcon(t)
    numcone = getnumcone(t)
    numbarvar = getnumbarvar(t)

    if n > 0
        for j in 1:n
            println(t[Var(j)])
        end
    end

    if m > 0
        for j in 1:m
            println(t[Con(j)])
        end
    end

    if numcone > 0
        for j in 1:numcone
            println(t[Cone(j)])
        end
    end

    if numbarvar > 0
        for j in 1:getnumbarvar(t)
            println(t[Barvar(j)])
        end
    end
end
