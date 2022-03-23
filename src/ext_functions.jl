module Ext
using ..Mosek
using Printf

struct VarByIndex index :: Int32 end
struct VarByName name :: String end
struct BarvarByIndex index :: Int64 end
struct BarvarByName name :: String end
struct ConByIndex index :: Int32 end
struct ConByName  name  :: String end
struct ConeByIndex index :: Int32 end
struct ConeByName  name  :: String end

struct Sol which :: Mosek.Soltype end

struct Solution
    t::Mosek.Task
    which :: Mosek.Soltype
end

struct VariableSolution
    t::Mosek.Task
    which :: Mosek.Soltype
    index :: Int32
end

struct ConstraintSolution
    t::Mosek.Task
    which :: Mosek.Soltype
    index :: Int32
end

struct BarvarSolution
    t::Mosek.Task
    which :: Mosek.Soltype
    index :: Int64
end

struct ObjSolution
    t::Mosek.Task
    which :: Mosek.Soltype
end

struct Variable
    t::Mosek.Task
    index :: Int32
end

struct SemidefiniteVariable
    t::Mosek.Task
    index :: Int32
end

struct Constraint
    t::Mosek.Task
    index :: Int32
end


struct ConeConstraint
    t::Mosek.Task
    index :: Int32
end

struct Obj
end

struct Objective
    t::Mosek.Task
end

struct Symmat
    index :: Int
end

Base.getindex(t::Mosek.Task, sol :: Sol) = Solution(t,sol.which)

Base.getindex(t::Mosek.Task, index :: VarByIndex) = Variable(t,index.index)
function Base.getindex(t::Mosek.Task, name :: VarByName)
    ok,index = getvarnameindex(t,name.name)
    if ok > 0
        getindex(t,VarByIndex(index))
    else
        throw(KeyError("'$name' not found"))
    end
end

Base.getindex(t::Mosek.Task, index :: BarvarByIndex) = SemidefiniteVariable(t,index.index)
function Base.getindex(t::Mosek.Task, name :: BarvarByName)
    ok,index = getbarvarnameindex(t,name.name)
    if ok > 0
        getindex(t,BarvarByIndex(index))
    else
        throw(KeyError("'$name' not found"))
    end
end

Base.getindex(t::Mosek.Task, index :: ConByIndex) = Constraint(t,index.index)
function Base.getindex(t::Mosek.Task, name :: ConByName)
    ok,index = getconnameindex(t,name.name)
    if ok > 0
        getindex(t,ConByIndex(index))
    else
        throw(KeyError("'$name' not found"))
    end
end


Base.getindex(t::Mosek.Task, index :: ConeByIndex) = ConeConstraint(t,index.index)
function Base.getindex(t::Mosek.Task, name :: ConeByName)
    ok,index = getconenameindex(t,name.name)
    if ok > 0
        getindex(t,ConeByIndex(index))
    else
        throw(KeyError("'$name' not found"))
    end
end

Base.getindex(t::Mosek.Task, index :: Obj) = Objective(t)
function Base.getindex(t::Mosek.Task, index :: Symmat)
    dim,nz,tp = getsymmatinfo(t,index.index)
    subi,subj,valij = getsparsesymmat(t,index.index)
    for i in 1:nz
        if subi[i] != subj[i]
            push!(subi,subj[i])
            push!(subj,subi[i])
            push!(valij,valij[i])
        end
    end
    sparse(subi,subj,valij,dim,dim)
end

Base.getindex(s::Solution, index :: VarByIndex) = VariableSolution(s.t,s.which,index.index)
Base.getindex(s::Solution, index :: ConByIndex) = ConstraintSolution(s.t,s.which,index.index)
Base.getindex(s::Solution, index :: BarvarByIndex) = BarvarSolution(s.t,s.which,index.index)
Base.getindex(s::Solution, index :: Obj) = ObjSolution(s.t,s.which)

function Base.show(f::IO, var :: Variable)
    bk,bl,bu = getvarbound(var.t,var.index)
    name = getvarname(var.t,var.index)

    if bk == MSK_BK_FX
        print(f,"Variable('$name' = $bl)")
    else
        blstr =
            if     bk == MSK_BK_LO || bk == MSK_BK_RA "[$bl"
            else   "]-Inf"
            end
        bustr =
            if     bk == MSK_BK_UP || bk == MSK_BK_RA "$bu]"
            else   "+Inf["
            end

        print(f,"Variable('$name' ∈ $blstr;$bustr)")
    end
end

function Base.show(f::IO, var :: SemidefiniteVariable)
    if var.index > 0 && var.index <= getnumbarvar(var.t)
        dim,ns,tp = getsymmatinfo(var.t,var.index)
        name = mkbarvarname(var.t,var.index)

        print(f,"SemidefiniteVariable('$name' ∈ 𝒞_S($dim))")
    else
        error("Invalid reference")
    end
end

function Base.show(f::IO, cone :: ConeConstraint)
    name = getconename(cone.t,cone.index)
    ct,cp,nummem,submem = getcone(cone.t,cone.index)

    dom =
        if     ct == MSK_CT_ZERO "𝒞_0"
        elseif ct == MSK_CT_QUAD "𝒞_q"
        elseif ct == MSK_CT_RQUAD "𝒞_qr"
        elseif ct == MSK_CT_PPOW "𝒞_pow{$cp}"
        elseif ct == MSK_CT_PEXP "𝒞_exp"
        else "?"
        end

    varnames = join(map(
        j -> let varname = getvarname(cone.t,j)
               if length(varname) == 0
                 "#x$j"
               else
                 Mosek.escapename(varname)
               end
             end,
        submem),",")

    print(f,"ConeByIndex('$name': ($varnames) ∈ $dom($nummem))")
end

mkvarname(t,j) =
    let n = getvarname(t,j)
        if length(n) == 0
            "#x$j"
        else
            Mosek.escapename(n)
        end
    end

mkconname(t,j) =
    let n = getconname(t,j)
        if length(n) == 0
            "#c$j"
        else
            Mosek.escapename(n)
        end
    end

mkbarvarname(t,j) =
    let n = getbarvarname(t,j)
        if length(n) == 0
            "#X̄$j"
        else
            Mosek.escapename(n)
        end
    end

Base.show(f::IO,con :: Objective) = Base.show(f,con,20)
#Base.showall(f::IO,con :: Objective) = Base.show(f,con,0)

function Base.show(f::IO, con :: Objective, limit :: Int)
    t = con.t
    name = getobjname(t)

    print(f,"Objective('$name': ")

    termslimit = if limit > 0 limit else typemax(Int) end
    termsomitted = 0
    nterms = 0
    for j in 1:getnumvar(t)
        cj = getcj(t,j)
        if cj < 0 || cj > 0
            if nterms < termslimit
                varname = Mosek.escapename(getvarname(t,j))
                if length(varname) == 0
                    varname = "#x$j"
                end
                print(f,"$(Mosek.fmtcof(cj)) $varname ")
                nterms += 1
            else
                termsomitted += 1
            end
        end
    end

    if getnumbarvar(t) > 0
        _,barcidx = getbarcsparsity(t)

        for idx in barcidx
            if nterms < termslimit
                (barcj,num,sub,w) = getbarcidx(t,idx)

                barvarname = mkbarvarname(t,barcj)

                if num == 1
                    print(f,"$(Mosek.fmtcof(w[1])) #MX$(sub[1]) ⋅ $barvarname) ")
                elseif num > 1
                    print(f,"(")
                    printf(f,"$(Mosek.fmtcof(w[1])) #MX$(sub[1])")
                    print(f,") ⋅ $barvarname ")
                end
                nterms += 1
            else
                termsomitted += 1
            end
        end
    end

    numqobjnz = getnumqobjnz(t)
    if numqobjnz > 0
        (nnz,qcsubi,qcsubj,qcval) = getqobj(t)
        for k in 1:nnz
            if nterms < termslimit
                if qcsubi[k] == qcsubj[k]
                    ni = mkvarname(t,qcsubi[k])
                    print(f,"$(Mosek.fmtcof(qcval[k])) $ni ^ 2")
                else
                    ni = mkvarname(t,qcsubi[k])
                    nj = mkvarname(t,qcsubj[k])
                    print(f,"$(Mosek.fmtcof(qcval[k])) $ni ⋅ $nj")
                end
                nterms += 1
            else
                termsomitted += 1
            end
        end
    end

    if termsomitted > 0
        print(f," ...($termsomitted terms omitted)")
    end

    print(f,")")
end


function Base.show(f::IO, con :: Constraint)
    t = con.t
    i = con.index
    bk,bl,bu = getconbound(t,i)
    name = getconname(t,i)

    print(f,"Constraint('$name': ")
    if bk == MSK_BK_FX
        print(f,"$bl = ")
    else
        blstr =
            if     bk == MSK_BK_LO || bk == MSK_BK_RA "[$bl"
            else   "]-Inf"
            end
        bustr =
            if     bk == MSK_BK_UP || bk == MSK_BK_RA "$bu]"
            else   "+Inf["
            end

        print(f,"$blstr;$bustr ∋ ")
    end

    (nzi,subi,vali) = getarow(t,Int32(i))

    for k in 1:nzi
        varname = Mosek.escapename(getvarname(t,subi[k]))
        if length(varname) == 0
            varname = "#x$(subi[k])"
        end
        print(f,"$(Mosek.fmtcof(vali[k])) $varname ")
    end

    numqconnz = getnumqconknz(t,i)
    if numqconnz > 0
        (nnz,qcsubi,qcsubj,qcval) = getqconk(t,i)
        for k in 1:nnz
            if qcsubi[k] == qcsubj[k]
                ni = mkvarname(t,qcsubi[k])
                print(f,"$(Mosek.fmtcof(qcval[k])) $ni ^ 2")
            else
                ni = mkvarname(t,qcsubi[k])
                nj = mkvarname(t,qcsubj[k])
                print(f,"$(Mosek.fmtcof(qcval[k])) $ni ⋅ $nj")
            end
        end
    end



    _,baraidx = getbarasparsity(t)
    if length(baraidx) > 0
        # count
        numbarnzi = count(idx -> getbaraidxij(t,idx)[1] == i, baraidx)
        idxs = Vector{Int}(undef,numbarnzi)
        k = 0
        for idx in baraidx
            (barai,baraj) = getbaraidxij(t,idx)
            if barai == i
                k += 1
                idxs[k] = idx
            end
        end

        for idx in idxs
            (barai,baraj,num,sub,w) = getbaraidx(t,idx)
            barvarname = mkbarvarname(t,baraj)

            if num == 1
                print(f,"$(Mosek.fmtcof(w[1])) #MX$(sub[1]) ⋅ $barvarname ")
            elseif num > 1
                print(f,"(")
                printf(f,"$(Mosek.fmtcof(w[1])) #MX$(sub[1])")
                print(f,") ⋅ $barvarname ")
            end
        end
    end

    print(f,")")
end




#Base.showall(f::IO, sol :: Solution) = Base.show(f,sol,0)
Base.show(f::IO, sol :: Solution) = Base.show(f,sol,20)

solstainfo(solsta) =
    if     solsta == MSK_SOL_STA_DUAL_FEAS                "DualFeasible",false,true
    elseif solsta == MSK_SOL_STA_DUAL_ILLPOSED_CER        "DualIllposedCertificate",true,false
    elseif solsta == MSK_SOL_STA_DUAL_INFEAS_CER          "DualInfeasibilityCertificate",true,false
    elseif solsta == MSK_SOL_STA_INTEGER_OPTIMAL          "IntegerOptimal",true,false
    elseif solsta == MSK_SOL_STA_OPTIMAL                  "Optimal",true,true
    elseif solsta == MSK_SOL_STA_PRIM_AND_DUAL_FEAS       "PrimalAndDualFeasible",true,true
    elseif solsta == MSK_SOL_STA_PRIM_FEAS                "PrimalFeasible",true,false
    elseif solsta == MSK_SOL_STA_PRIM_ILLPOSED_CER        "PrimalIllposedCertificate",false,true
    elseif solsta == MSK_SOL_STA_PRIM_INFEAS_CER          "PrimakInfeasibleCertificate",false,true
    else "Unknown",true,true
    end

function Base.show(f::IO, sol :: Solution, limit :: Int)
    t = sol.t
    if 0 != solutiondef(t,sol.which)
        solname =
            if     sol.which == MSK_SOL_ITG "Integer Solution"
            elseif sol.which == MSK_SOL_BAS "Basis Solution"
            else                            "Interior Solution"
            end
        numvar = getnumvar(t)
        numcon = getnumcon(t)
        numbarvar = getnumbarvar(t)
        prosta,solsta,skc,skx,skn,xc,xx,y,slc,suc,slx,sux,snx = getsolution(t,sol.which)

        solstaname,pdef,ddef = solstainfo(solsta)

        println(f,"$solname, status = $solstaname")

        limitnumvar = if limit == 0 numvar else min(numvar,limit) end
        limitnumcon = if limit == 0 numcon else min(numcon,limit) end

        if pdef && ddef
            pobj = getprimalobj(t,sol.which)
            dobj = getdualobj(t,sol.which)
            println(f,"    Objective: $pobj | $dobj")
        elseif pdef
            pobj = getprimalobj(t,sol.which)
            println(f,"    Objective: $pobj | -")
        elseif ddef
            dobj = getdualobj(t,sol.which)
            println(f,"    Objective: - | $dobj")
        else
            println(f,"    Objective: - | -")
        end

        if numvar > 0
            println(f,"    Variable solution")

            if pdef && ddef
                @printf(f,"        %-20s  %13s  %13s  %13s  %13s\n","name","level","dual lower","dual upper","dual conic")
                for j in 1:limitnumvar
                    name = mkvarname(t,j)
                    @printf(f,"        %-20s: %13.4e  %13.4e  %13.4e  %13.4e\n",name,xx[j],slx[j],sux[j],snx[j])
                end
            elseif pdef
                @printf(f,"        %-20s  %13s  -  -  -\n","name","level")

                for j in 1:limitnumvar
                    name = mkvarname(t,j)
                    @printf(f,"        %-20s: %13.4e  -  -  -\n",name,xx[j])
                end
            elseif ddef
                @printf(f,"        %-20s    %13s  %13s  %13s\n","name","dual lower","dual upper","dual conic")
                for j in 1:limitnumvar
                    name = mkvarname(t,j)
                    @printf(f,"        %-20s: - %13.4e  %13.4e  %13.4e\n",name,slx[j],sux[j],snx[j])
                end
            end
            if limitnumvar < numvar
                println(f,"        ... ($(numvar-limitnumvar) variables omitted)")
            end
        end

        if numbarvar > 0
            println(f,"    PSD Variable solution")
            for k in 1:numbarvar
                dim = getdimbarvarj(t,k)
                markatrow = dim >> 1
                barx = getbarxj(t,sol.which,k)
                bars = getbarsj(t,sol.which,k)
                name = mkbarvarname(t,k)

                println(f,"        $name: Symmetric $dim × $dim")
                if pdef
                    px = 1
                    for i in 1:dim
                        if i == markatrow print(f,"            X̄ = |")
                        else              print(f,"                |")
                        end
                        if i > 1
                            for k in 1:i-1 print(f,"           ") end
                        end
                        for j in i:dim
                            @printf(f," %10.2e",barx[px])
                            px += 1
                        end
                        println(f," |")
                    end
                end
                if ddef
                    if pdef println(f) end

                    ps = 1

                    for i in 1:dim
                        if i == markatrow print(f,"            S̄ = |")
                        else              print(f,"                |")
                        end

                        if i > 1
                            for k in 1:i-1 print(f,"           ") end
                        end
                        for j in i:dim
                            @printf(f," %10.2e",bars[ps])
                            ps += 1
                        end

                        println(f," |")
                    end
                end
            end
        end

        if numcon > 0
            println(f,"    Constraint solution")
            if pdef && ddef
                @printf(f,"        %-20s  %13s  %13s  %13s  %13s\n","name","level","dual lower","dual upper","y")
                for j in 1:limitnumcon
                    name = mkconname(t,j)
                    @printf(f,"        %-20s: %13.4e  %13.4e  %13.4e  %13.4e\n",name,xc[j],slc[j],suc[j],y[j]) # 412
                end
            elseif pdef
                @printf(f,"        %-20s  %13s  -  -  -\n","name","level")
                for j in 1:limitnumcon
                    name = mkconname(t,j)
                    @printf(f,"        %-20s: %13.4e  -  -  -\n",name,xc[j])
                end
            elseif ddef
                @printf(f,"        %-20s    %13s  %13s  %13s\n","name","dual lower","dual upper","y")
                for j in 1:limitnumcon
                    name = mkconname(t,j)
                    @printf(f,"        %-20s: - %13.4e  %13.4e  %13e\n",name,slc[j],suc[j],y[j])
                end
            end
            if limitnumcon < numcon
                println(f,"        ... ($(numcon-limitnumcon) constraints omitted)")
            end
        end

    else
        error("Solution not defined")
    end
end

function Base.values(sol::VariableSolution)
    t = sol.t
    xx =
        try
            getxxslice(t,sol.which,sol.index,sol.index+1)[1]
        catch
            NaN
        end
    slx =
        try
            getslxslice(t,sol.which,sol.index,sol.index+1)[1]
        catch
            NaN
        end
    sux =
        try
            getsuxslice(t,sol.which,sol.index,sol.index+1)[1]
        catch
            NaN
        end
    snx =
        try
            getsnxslice(t,sol.which,sol.index,sol.index+1)[1]
        catch
            NaN
        end

    xx,slx,sux,snx
end

function Base.show(f::IO,sol::VariableSolution)
    name = mkvarname(sol.t,sol.index)
    xx,slx,sux,snx = Base.values(sol)
    print(f,"$name: $xx, dual lower: $slx, dual upper: $sux, dual conic: $snx")
end


function Base.values(sol::ObjSolution)
    pdef,pobj =
        try
            true,getprimalobj(sol.t,sol.which)
        catch
            false,NaN
        end

    ddef,dobj =
        try
            true,getdualobj(sol.t,sol.which)
        catch
            false,NaN
        end
    pobj,dobj
end

function Base.show(f::IO,sol::ObjSolution)
    pdef,pobj =
        try
            true,getprimalobj(sol.t,sol.which)
        catch
            false,NaN
        end

    ddef,dobj =
        try
            true,getdualobj(sol.t,sol.which)
        catch
            false,NaN
        end
    if pdef && ddef
        print(f,"Objective($pobj,$dobj)")
    elseif pdef
        print(f,"Objective($pobj,-)")
    elseif ddef
        print(f,"Objective(-,$dobj)")
    else
        print(f,"Objective(-,-)")
    end
end


function Base.values(sol::BarvarSolution)
    t = sol.t
    dim = getdimbarvarj(t,sol.index)
    barx =
        try
            v = getbarxj(t,sol.which,sol.index)
            res = Array{Float64,2}(undef,(dim,dim))
            k = 0
            for i in 1:dim
                k += 1
                res[i,i] = v[k]
                for j in i+1:dim
                    k += 1
                    res[i,j] = v[k]
                    res[j,i] = v[k]
                end
            end
            res
        catch
            fill(NaN,(dim,dim))
        end
    bars =
        try
            v = getbarsj(t,sol.which,sol.index)
            res = Array{Float64,2}(undef,(dim,dim))
            k = 0
            for i in 1:dim
                k += 1
                res[i,i] = v[k]
                for j in i+1:dim
                    k += 1
                    res[i,j] = v[k]
                    res[j,i] = v[k]
                end
            end
            res
        catch
            fill(NaN,(dim,dim))
        end
    barx,bars
end

function Base.show(f::IO,sol::BarvarSolution)
    t = sol.t
    name = mkvarname(t,sol.index)
    dim = getdimbarvarj(t,sol.index)
    markatrow = dim >> 1
    barx =
        try
            getbarxj(t,sol.which,sol.index)
        catch
            fill(NaN,(dim*(dim+1)) >> 1)
        end
    bars =
        try
            getbarsj(t,sol.which,sol.index)
        catch
            fill(NaN,(dim*(dim+1)) >> 1)
        end


    println(f,"$name: Symmetric $dim × $dim")
    px = 1
    ps = 1
    for i in 1:dim
        if i == markatrow print(f,"  X̄ = |")
        else              print(f,"      |")
        end
        if i > 1
            for k in 1:i-1 print(f,"           ") end
        end
        for j in i:dim
            @printf(f," %10.2e",barx[px])
            px += 1
        end

        println(f," |")
    end
    println(f)

    for i in 1:dim
        if i == markatrow print(f,"  S̄ = |")
        else              print(f,"      |")
        end
        if i > 1
            for k in 1:i-1 print(f,"           ") end
        end
        for j in i:dim
            @printf(f," %10.2e",bars[ps])
            ps += 1
        end

        println(f," |")
    end

end


function Base.values(sol::ConstraintSolution)
    t = sol.t
    xc =
        try
            getxcslice(t,sol.which,sol.index,sol.index+1)[1]
        catch
            NaN
        end
    slc =
        try
            getslcslice(t,sol.which,sol.index,sol.index+1)[1]
        catch
            NaN
        end
    suc =
        try
            getsucslice(t,sol.which,sol.index,sol.index+1)[1]
        catch
            NaN
        end
    y   =
        try
            getyslice(t,sol.which,sol.index,sol.index+1)[1]
        catch
            NaN
        end
    xc,slc,suc,y
end

function Base.show(f::IO,sol::ConstraintSolution)
    name = mkconname(sol.t,sol.index)
    xc,slc,suc,y = Base.values(sol)
    print(f,"$name: $xc, dual lower: $slc, dual upper: $suc, dual conic: $y")
end






function findvars(t::Mosek.Task, name::String)
    res = VarByIndex[]
    for j in 1:getnumvar(t)
        if getvarname(t,j) == name
            push!(res,VarByIndex(j))
        end
    end
    res
end

function findvars(t::Mosek.Task, regex::Regex)
    res = VarByIndex[]
    for j in 1:getnumvar(t)
        if ismatch(regex,getvarname(t,j))
            push!(res,VarByIndex(j))
        end
    end
    res
end



function findcons(t::Mosek.Task, name::String)
    res = ConByIndex[]
    for j in 1:getnumcon(t)
        if getconname(t,j) == name
            push!(res,ConByIndex(j))
        end
    end
    res
end

function findcons(t::Mosek.Task, regex::Regex)
    res = ConByIndex[]
    for j in 1:getnumcon(t)
        if ismatch(regex,getconname(t,j))
            push!(res,ConByIndex(j))
        end
    end
    res
end

function findcones(t::Mosek.Task, name::String)
    res = ConeByIndex[]
    for j in 1:getnumcone(t)
        if getconename(t,j) == name
            push!(res,ConeByIndex(j))
        end
    end
    res
end

function findcones(t::Mosek.Task, regex::Regex)
    res = ConeByIndex[]
    for j in 1:getnumcone(t)
        if ismatch(regex,getconename(t,j))
            push!(res,ConeByIndex(j))
        end
    end
    res
end

Var(index::I) where {I <: Integer} = VarByIndex(convert(Int32,index))
Var(name::String) = VarByName(name)

Barvar(index::I) where {I <: Integer} = BarvarByIndex(convert(Int32,index))
Barvar(name::String) = BarvarByName(name)

Con(index::I) where {I <: Integer} = ConByIndex(convert(Int32,index))
Con(name::String) = ConByName(name)

Cone(index::I) where {I <: Integer} = ConeByIndex(convert(Int32,index))
Cone(name::String) = ConeByName(name)

###########

function Base.getindex(t::Mosek.Task, ri :: UnitRange{T}, rj :: Colon) where { T <: Integer }
    first = Int32(ri.start)
    last  = Int32(ri.stop+1)
    m = last-first
    n = getnumvar(t)

    subi,subj,valij = getarowslicetrip(t,first,last)
    sparse(subi-(first-1),subj,valij,m,n)
end

function Base.getindex(t::Mosek.Task, ri :: T, rj :: Colon) where { T <: Integer }
    i  = Int32(ri)
    n = getnumvar(t)

    subj,valij = getarow(t,i)
    sparsevec(subj,valij,n)
end

function Base.getindex(t::Mosek.Task, ri :: Colon, rj :: UnitRange{T}) where { T <: Integer }
    first = Int32(ri.start)
    last  = Int32(ri.stop+1)
    m = getnumcon(t)
    n = last-first

    subi,subj,valij = getacolslicetrip(t,first,last)
    sparse(subi,subj,valij,m,n)
end

function Base.getindex(t::Mosek.Task, ri :: Colon, rj :: T) where { T <: Integer }
    j = Int32(rj)
    m = getnumcon(t)

    subi,valij = getacol(t,j)
    sparsevector(subi,valij,m)
end

function Base.getindex(t::Mosek.Task, ri :: Colon, rj :: Colon)
    subi,subj,valij = getarowslicetrip(t,1,getnumcon(t)+1)
    sparse(subi,subj,valij,getnumcon(t),getnumvar(t))
end

Base.getindex(t::Mosek.Task, ri :: UnitRange{T}, rj :: UnitRange{T}) where { T <: Integer } = getindex(t,ri,Colon())[:,rj]

###########


export
    Var,
    Barvar,
    Con,
    Cone,
    Obj,
    Symmat,
    Barvar,
    Sol,
    findvars,
    findcons,
    findcones
end
