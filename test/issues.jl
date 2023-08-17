using Mosek
using Test


issue188() =
    maketask() do task
        appendvars(task,6)
        appendcone(task,MSK_CT_QUAD,0.0,[1,2,3])
        appendcone(task,MSK_CT_RQUAD,0.0,[4,5,6])

        removecones(task,[1])
        (ct,cpar,nummem) = getconeinfo(task,1)
        
        @test ct == MSK_CT_RQUAD
    end

gitlab_2145() = 
    maketask() do task
        appendvars(task,1)
        domidx = appendrdomain(task,3)
        appendafes(task,1)
        appendacc(task,domidx,Int64[1,1,1],nothing)
    end

@testset "[Github issues]" begin
    issue188()
end

@testset "[Gitlab issues]" begin
    gitlab_2145()
end
