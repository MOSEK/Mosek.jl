using JuMP
using Base.Test

include(joinpath(Pkg.dir("JuMP"),"test","solvers.jl"))

include(joinpath(Pkg.dir("JuMP"),"test","model.jl"));        length(   lp_solvers) == 0 && warn("Model tests not run!")
include(joinpath(Pkg.dir("JuMP"),"test","probmod.jl"));      length(   lp_solvers) == 0 && warn("Prob. mod. tests not run!")
include(joinpath(Pkg.dir("JuMP"),"test","callback.jl"));     length( lazy_solvers) == 0 && warn("Callback tests not run!")
include(joinpath(Pkg.dir("JuMP"),"test","qcqpmodel.jl"));    length( quad_solvers) == 0 && warn("Quadratic tests not run!")
#include(joinpath(Pkg.dir("JuMP"),"test","nonlinear.jl"));    length(  nlp_solvers) == 0 && warn("Nonlinear tests not run!")
#                                                    length(minlp_solvers) == 0 && warn("Mixed-integer Nonlinear tests not run!")
include(joinpath(Pkg.dir("JuMP"),"test","sdp.jl"));          length(  sdp_solvers) == 0 && warn("Semidefinite tests not run!")
include(joinpath(Pkg.dir("JuMP"),"test","socduals.jl"));     length(conic_solvers_with_duals) == 0 && warn("Conic solvers with duals tests not run!")
