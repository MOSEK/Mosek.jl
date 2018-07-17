using Mosek




escapename(n::String) =
    replace(n,r"(^#)","_#")

fmtcof(c::Float64) =
    if     c < 0.0
        if     c < -1.0 || c > -1.0
             "- $(-c)"
        else
            "-"
        end
    elseif c > 0.0
        if c < 1.0 || c > 1.0
            "+ $c"
        else
            "+"
        end
    else
        "0"
    end
        
        

write(t::Mosek.Task, filename::String) =
    open(filename,"w") do f
        write(t,f)
    end

function write(t::Mosek.Task, f::IO)
    numvar    = getnumvar(t)
    numbarvar = getnumbarvar(t)
    numcon    = getnumcon(t) 
    numcone   = getnumcone(t)
    numanz    = getnumanz(t)

    varnames    = String[ escapename(getvarname(t,i)) for i in 1:numvar] 
    barvarnames = String[ escapename(getbarvarname(t,i)) for i in 1:numbarvar] 
    connames    = String[ escapename(getconname(t,i)) for i in 1:numcon]
    conenames   = String[ escapename(getconename(t,i)) for i in 1:numcone]

    for i in 1:numvar
        if length(varnames[i]) == 0
            varnames[i] = "#x$i"
        end
    end

    for i in 1:numcon
        if length(connames[i]) == 0
            connames[i] = "#c$i"
        end
    end

    for i in 1:numcone
        if length(conenames[i]) == 0
            conenames[i] = "#k$i"
        end
    end

    for i in 1:numbarvar
        if length(barvarnames[i]) == 0
            barvarnames[i] = "#X̂$i"
        end
    end

    maxvarnamelen = max(map(n -> length(n),varnames)...)
    maxconnamelen = max(map(n -> length(n),connames)...)
    
    taskname = gettaskname(t)
    
    println(f,"Task $(escape_string(taskname))")

    _,barcidx = getbarcsparsity(t)

    _,baraidx = getbarasparsity(t)
    baraptrb = zeros(Int,numcon+1)
    for idx in baraidx
        (barai,baraj) = getbaraidxij(t,idx)
        baraptrb[barai+1] += 1
    end
    for i in 1:numcon
        baraptrb[i+1] += baraptrb[i]
    end
    baraptrb += 1
    
    baraptre = copy(baraptrb)
    baraidx_idx = Vector{Int}(length(baraidx))
    for idx in baraidx
        (barai,baraj) = getbaraidxij(t,idx)
        baraidx_idx[baraptre[barai]] = idx
        baraptre[barai] += 1
    end
        
    # Objective
    objsense = getobjsense(t)
    objname  = getobjname(t)

    if objsense == MSK_OBJECTIVE_SENSE_MAXIMIZE
        println(f,"Maximize $(escape_string(objname))")
    else
        println(f,"Minimze $(escape_string(objname))")
    end
    c = getc(t)
    print(f,"    ")
    for j in 1:numvar
        if c[j] < 0.0 || c[j] > 0.0
            print(f,"$(fmtcof(c[j])) $(varnames[j]) ")
        end
    end

    # Bar entries go here
    
    cfix = getcfix(t)
    if cfix < 0 || cfix > 0
        print(f,fmtcof(cfix))
    end
    println()

    # Constraints
    if numcon+numcone > 0
        println(f,"Subject to")

        for i in 1:numcon
            (nzi,subi,vali) = getarow(t,Int32(i))
            (bk,bl,bu) = getconbound(t,Int32(i))

            
            print(f,"    $(connames[i]): ")
            if bk == MSK_BK_LO || bk == MSK_BK_RA
                print(f,"$bl < ")
            end
            for k in 1:nzi
                print(f,"$(fmtcof(vali[k])) $(varnames[subi[k]]) ")
            end

            for k in baraptrb[i]:baraptrb[i+1]-1
                idx = baraidx_idx[k]
                (barai,baraj,num,sub,w) = getbaraidx(t,idx)
                if num == 1
                    print(f,"$(fmtcof(w[1])) #MX$(sub[1]) ⋅ $(barvarnames[baraj]) ")
                elseif num > 1
                    print(f,"(")
                    printf(f,"$(fmtcof(w[1])) #MX$(sub[1])")                    
                    print(f,") ⋅ $(barvarnames[baraj]) ")
                end
            end

            if bk == MSK_BK_UP || bk == MSK_BK_RA
                println(f,"< $bu")
            elseif bk == MSK_BK_FX
                println(f,"= $bu")
            end
            
        end
    
        for k in 1:numcone
            (ct,conepar,nummem,submem) = getcone(t,k)
            dom =
                if ct == MSK_CT_QUAD "Q"
                elseif ct == MSK_CT_RQUAD "Q_r"
                else "?"
                end

            print(f,"    $(conenames[k]): ($(varnames[submem[1]])")
            for i in 2:nummem
                print(f,",$(varnames[submem[i]])")
            end
            println(f,") ∈ $dom(nummem)")
        end
    end

    # Variable bounds
    if numvar > 0
        println(f,"Variables")
        for j in 1:numvar
            bk,bl,bu = getvarbound(t,j)

            if bk == MSK_BK_FX
                println(f,"    $(varnames[j]) = $bl")
            else
                blstr =
                    if     bk == MSK_BK_LO || bk == MSK_BK_RA "[$bl"
                    else   "]-Inf"
                    end
                bustr = 
                    if     bk == MSK_BK_UP || bk == MSK_BK_RA "$bu]"
                    else   "+Inf["
                    end

                println(f,"    $(varnames[j]) ∈ $blstr;$bustr")
            end
        end

        for j in 1:numbarvar
            dim = getdimbarvarj(t,j)
            println(f,"    $(barvarnames[j]) ∈ S($dim)")
        end

    end

    # Matsto
    numsymmat = getnumsymmat(t)
    if numsymmat > 0
        println(f,"Symmetric matrixes")
        for i in 1:numsymmat
            dim,nz,tp = getsymmatinfo(t,i)
            (msubi,msubj,mvalij) = getsparsesymmat(t,i)
            print(f,"    #MX$i ($dim): [ " )
            if nz > 0                
                print(f,"($(msubi[1]),$(msubj[1]),$(mvalij[1]))")
                for i in 2:nz
                    print(f,",($(msubi[i]),$(msubj[i]),$(mvalij[i]))")
                end
            end
            println(f," ]")
        end
    end
end


maketask() do t
    readdata(t,ARGS[1])
    write(t,STDOUT)
end
