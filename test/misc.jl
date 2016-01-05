
module TestMiscSDP

import Mosek
using JuMP,FactCheck

facts("[sdp] Robust uncertainty example") do
    solver = Mosek.MosekSolver(QUIET=true)

    context("With solver $(typeof(solver))") do
        include(joinpath(Pkg.dir("JuMP"), "test","data","robust_uncertainty.jl"))
        R = 1
        d = 3
        𝛿 = 0.05
        ɛ = 0.05
        N = ceil((2+2log(2/𝛿))^2) + 1

        Γ1(𝛿,N) = (R/sqrt(N))*(2+sqrt(2*log(1/𝛿)))
        Γ2(𝛿,N) = (2R^2/sqrt(N))*(2+sqrt(2*log(2/𝛿)))

        #for d in [3,5,8]; context("d = $d") do
        for d in [3]; context("d = $d") do
            μhat = μhats[d]
            M = Ms[d]
            Σhat = 1/(d-1)*(M-ones(d)*μhat')'*(M-ones(d)*μhat')

            m = Model(solver=solver)

            @defVar(m, Σ[1:d,1:d], SDP)
            @defVar(m, u[1:d])
            @defVar(m, μ[1:d])

            @defVar(m, t1 >= 0)
            @defVar(m, L1[1:d])
            @addConstraint(m, L1 .== (μ-μhat))
            @addConstraint(m, sum{L1[i]^2, i=1:d} <= t1^2)
            @addConstraint(m, t1 <= Γ1(𝛿/2,N))

            @defVar(m, t2 >= 0)
            @defVar(m, L2[1:d,1:d])
            @addConstraint(m, L2 .== (Σ-Σhat))
            @addConstraint(m, sum{L2[i,j]^2, i=1:d, j=1:d} <= t2^2)
            @addConstraint(m, t2 <= Γ2(𝛿/2,N))

            A = [(1-ɛ)/ɛ (u-μ)';
                 (u-μ)     Σ   ]
            @addSDPConstraint(m, A >= 0)

            c = cs[d]
            @setObjective(m, Max, dot(c,u))

            stat = solve(m)

            object = getObjectiveValue(m)
            exact = dot(μhat,c) + Γ1(𝛿/2,N)*norm(c) + sqrt((1-ɛ)/ɛ)*sqrt(dot(c,(Σhat+Γ2(𝛿/2,N)*eye(d,d))*c))
            @fact stat --> :Optimal
            @fact abs(object - exact) --> roughly(0, 1e-5)


            resΣ  = getValue(Σ)
            resu  = getValue(u)
            resμ  = getValue(μ)
            rest1 = getValue(t1)
            resL1 = getValue(L1)
            rest2 = getValue(t2)
            resL2 = getValue(L2)
            
            println("object = $object, exact = $exact")
            println("object val : dot(c,u) : $(dot(c,resu))")

            println("L1 == μ-μhat : $(resL1 - (resμ-μhat))")
            println("sum{L1[i]^2, i=1:d} <= t1^2: $((sum(resL1.^2)) - (rest1.^2))")
            println("t1 <= Γ1(𝛿/2,N) : $(rest1 - (Γ1(𝛿/2,N)))")
            println("L2 .== (Σ-Σhat) : $(resL2 - (resΣ-Σhat))")
            println("sum{L2[i,j]^2, i=1:d, j=1:d} <= t2^2 : $(sum(resL2.^2)-rest2^2)")
            #@addConstraint(m, sum{L2[i,j]^2, i=1:d, j=1:d} <= t2^2)
            println("t2 <= Γ2(𝛿/2,N) : $(rest2 - Γ2(𝛿/2,N))")
            #@addConstraint(m, t2 <= Γ2(𝛿/2,N))
            println("---------- PSD con")
            resA = [(1-ɛ)/ɛ     (resu-resμ)';
                    (resu-resμ) resΣ   ]

            println("Eigvals A: $(eigvals(resA))")
            println("---------- vars")
            println("Σ PSD. Eigvals : $(eigvals(resΣ))")
            println("u free")
            println("μ free")
            println("t1 >= 0: $(rest1)")
            println("L1 free")
            println("t2 >= 0: $(rest2)")
            println("L2 free")

            intmodel = getInternalModel(m)
            Mosek.writedata(intmodel.task,"problem-$d.task")
        end
        end
    end

end
end

module TestX
using Mosek, MathProgBase

I = [1,1,2,2,3,3,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,32,32,33,33,34,34,35,36,37,38,39,40]
J = [14,10,15,11,16,12,13,18,1,21,2,24,3,19,2,22,4,25,5,20,3,23,5,26,6,17,13,17,13,14,15,16,17,18,19,20,21,22,23,24,25,26,7,10,8,11,9,12,1,2,3,4,5,6]
V = [1.0,-1.0,1.0,-1.0,1.0,-1.0,1.0,1.0,-1.0,1.0,-1.0,1.0,-1.0,1.0,-1.0,1.0,-1.0,1.0,-1.0,1.0,-1.0,1.0,-1.0,1.0,-1.0,1.0,1.0,1.0,-1.0,1.0,1.0,1.0,-1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,-1.0,1.0,-1.0,1.0,-1.0,1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0]
b = [-0.6088331477422722,-0.02852197818671054,-0.9016074145524158,0.4999165215049682,-0.025485683857737418,-0.07276320066249353,0.056565030662414396,-0.07276320066249353,-0.3319954143963322,0.3286736785015988,0.056565030662414396,0.3286736785015988,-0.43682537007384287,1.0516057442061533,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,18.999999999999996,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
var_cones = Any[(:SDP,1:6),(:NonNeg,[13,17]),(:Free,[7,8,9,10,11,12,14,15,16,18,19,20,21,22,23,24,25,26])]
con_cones = Any[(:NonNeg,[4,14]),(:NonPos,[15,16]),(:Zero,[1,2,3,5,6,7,8,9,10,11,12,13]),(:SOC,17:20),(:SOC,21:30),(:SDP,31:40)]
f = [-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,0.2269228958615089,-0.26647992165696055,-0.35104470584603914,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0]
numconstr = 40
numvar = 26
A = sparse(I,J,V,numconstr,numvar)

m = MathProgBase.ConicModel(MosekSolver(LOG=0))

MathProgBase.loadproblem!(m, f, A, b, con_cones, var_cones)
writedata(m.task,"problem-back1.task")
println("varbk before : $(m.varbk) $(m.varmap)")
MathProgBase.setvartype!(m, fill(:Cont,numvar))
println("varbk after  : $(m.varbk) $(m.varmap)")
writedata(m.task,"problem-back2.task")




MathProgBase.optimize!(m)


@assert MathProgBase.status(m) == :Optimal

obj = MathProgBase.getobjval(m)
trueobj = -2.7231533727479005
@show trueobj
@show obj
@show abs(obj-trueobj)

writedata(m.task,"problem-back.task")

end
