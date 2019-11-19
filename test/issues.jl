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



@testset "[Github issues]" begin
    issue188()
end