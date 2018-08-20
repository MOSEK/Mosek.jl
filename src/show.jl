
escapename(n::String) =
    replace(n,r"(^#)" => "_#")

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

#Base.showall(f::IO, t::Task) = showlimited(f,t,0)
Base.show(f::IO, t::Task) = showlimited(f,t,20)

Base.show(io::IOContext, t::Task) =
    if get(io, :short, false)
        showlimited(io,t,20)
    else
        showlimited(io,t,0)
    end

function showlimited(f::IO, t::Mosek.Task, limit :: Int)
    numvar    = getnumvar(t)
    numbarvar = getnumbarvar(t)
    numcon    = getnumcon(t) 
    numcone   = getnumcone(t)
    numanz    = getnumanz(t)
    numqobjnz = getnumqobjnz(t)
    numqconnz = if numcon > 0 sum(i -> getnumqconknz(t,i),1:numcon) else 0 end

    limitnumcon = if limit > 0 min(limit,numcon) else numcon end
    limitnumvar = if limit > 0 min(limit,numvar) else numvar end
    limitnumcone = if limit > 0 min(limit,numcone) else numcone end
    limitnumbarvar = if limit > 0 min(limit,numbarvar) else numbarvar end

    termslimit = numvar+numbarvar+numqobjnz+numqconnz
    if limit > 0 termslimit = min(limit,termslimit) end

    if numvar == 0 && numcon == 0 && numbarvar == 0
        println(f,"Empty Task")
    else
        varnames = String[ escapename(getvarname(t,i)) for i in 1:numvar] 
        barvarnames = String[ escapename(getbarvarname(t,i)) for i in 1:numbarvar] 
        connames = String[ escapename(getconname(t,i)) for i in 1:numcon]
        conenames = String[ escapename(getconename(t,i)) for i in 1:numcone]

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
                barvarnames[i] = "#X̄$i"
            end
        end

        maxvarnamelen = if numvar > 0 max(map(n -> length(n),varnames)...) else 0 end
        maxconnamelen = if numcon > 0 max(map(n -> length(n),connames)...) else 0 end

        taskname = gettaskname(t)

        println(f,"Task '$(escape_string(taskname))'")

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
        baraptrb .+= 1

        baraptre = copy(baraptrb)
        baraidx_idx = Vector{Int}(undef,length(baraidx))
        for idx in baraidx
            (barai,baraj) = getbaraidxij(t,idx)
            baraidx_idx[baraptre[barai]] = idx
            baraptre[barai] += 1
        end

        # Objective
        objsense = getobjsense(t)
        objname  = escapename(getobjname(t))
        if objname == ""
            objname = "#obj"
        end

        if objsense == MSK_OBJECTIVE_SENSE_MAXIMIZE
            println(f,"    Maximize")
        else
            println(f,"    Minimize")
        end
        c = getc(t)
        print(f,"        $objname: ")
        terms = 0
        omittedterms = 0
        for j in 1:numvar
            if c[j] < 0.0 || c[j] > 0.0
               if terms < termslimit
                   print(f,"$(fmtcof(c[j])) $(varnames[j]) ")
                   terms += 1
               else
                   omittedterms += 1
               end
            end
        end


        for idx in barcidx
            if terms < termslimit
                (barcj,num,sub,w) = getbarcidx(t,idx)
                if num == 1
                    print(f,"$(fmtcof(w[1])) #MX$(sub[1]) ⋅ $(barvarnames[barcj]) ")
                elseif num > 1
                    print(f,"(")
                    printf(f,"$(fmtcof(w[1])) #MX$(sub[1])")
                    print(f,") ⋅ $(barvarnames[barcj]) ")
                end
                terms += 1
            else
                omittedterms += 1
            end
        end

        numqobjnz = getnumqobjnz(t)
        if numqobjnz > 0
            (nnz,qcsubi,qcsubj,qcval) = getqobj(t)
            for k in 1:nnz
                if terms < limitterms
                    if qcsubi[k] == qcsubj[k]
                        print(f,"$(Mosek.fmtcof(qcval[k])) $(varnames[qcsubi[k]]) ^ 2")
                    else
                        ni = varnames[qcsubi[k]]
                        nj = varnames[qcsubj[k]]
                        print(f,"$(Mosek.fmtcof(qcval[k])) $ni ⋅ $nj")
                    end
                    terms += 1
                else
                    omittedterms += 1
                end
            end
        end

        if omittedterms > 0
            print(f," ... ($omittedterms terms omitted)")
        end

        cfix = getcfix(t)
        if cfix < 0 || cfix > 0
            print(f,fmtcof(cfix))
        end
        println(f)

        # Constraints
        if numcon+numcone > 0
            println(f,"    Subject to")

            for i in 1:limitnumcon
                (nzi,subi,vali) = getarow(t,Int32(i))
                (bk,bl,bu) = getconbound(t,Int32(i))

                
                print(f,"        $(connames[i]): ")
                if bk == MSK_BK_LO || bk == MSK_BK_RA
                    print(f,"$bl < ")
                end
                terms = 0
                omittedterms = 0
                for k in 1:nzi
                    if terms < termslimit
                        print(f,"$(fmtcof(vali[k])) $(varnames[subi[k]]) ")
                        terms += 1
                    else
                        omittedterms += 1
                    end
                end

                for k in baraptrb[i]:baraptrb[i+1]-1
                    if terms < termslimit
                        idx = baraidx_idx[k]
                        (barai,baraj,num,sub,w) = getbaraidx(t,idx)
                        if num == 1
                            print(f,"$(fmtcof(w[1])) #MX$(sub[1]) ⋅ $(barvarnames[baraj]) ")
                        elseif num > 1
                            print(f,"(")
                            printf(f,"$(fmtcof(w[1])) #MX$(sub[1])")                    
                            print(f,") ⋅ $(barvarnames[baraj]) ")
                        end
                        terms += 1
                    else
                        omittedterms += 1
                    end
                end

                numqconnz = getnumqconknz(t,i)
                if numqconnz > 0
                    (nnz,qcsubi,qcsubj,qcval) = getqconk(t,i)
                    for k in 1:nnz
                        if terms < limitterms
                            if qcsubi[k] == qcsubj[k]
                                print(f,"$(Mosek.fmtcof(qcval[k])) $(varnames[qcsubi[k]]) ^ 2")
                            else
                                ni = varnames[qcsubi[k]]
                                nj = varnames[qcsubj[k]]
                                print(f,"$(Mosek.fmtcof(qcval[k])) $ni ⋅ $nj")
                            end
                            terms += 1
                        else
                            omittedterms += 1
                        end
                    end
                end

                if omittedterms > 0
                    print(f," ... ($omittedterms terms omitted)")
                end

                if bk == MSK_BK_UP || bk == MSK_BK_RA
                    println(f,"< $bu")
                elseif bk == MSK_BK_FX
                    println(f,"= $bu")
                else
                    println(f,"")
                end
            end
            if limitnumcon < numcon
                println(f,"        ... ($(numcon-limitnumcon) constraints omitted)")
            end

            if numcone > 0
                for k in 1:limitnumcone
                    (ct,conepar,nummem,submem) = getcone(t,k)
                    dom =
                        if     ct == MSK_CT_ZERO "𝒞_0"
                        elseif ct == MSK_CT_QUAD "𝒞_q"
                        elseif ct == MSK_CT_RQUAD "𝒞_qr"
                        elseif ct == MSK_CT_PPOW "𝒞_pow{$conepar}"
                        elseif ct == MSK_CT_PEXP "𝒞_exp"
                        else "?{$ct}"
                        end

                    
                    print(f,"        $(conenames[k]): (")

                    if nummem > 0
                        print(f,"$(varnames[submem[1]])")
                        for i in 2:nummem
                            print(f,",$(varnames[submem[i]])")
                        end
                    end
                    println(f,") ∈ $dom($nummem)")
                end

                if limitnumcone < numcone
                    println(f,"        ... ($(numcone-limitnumcone) cones omitted)")
                end
            end
        end

        # Variable bounds
        if numvar+numbarvar > 0
            println(f,"    Variables")
            if numvar > 0
                for j in 1:limitnumvar
                    bk,bl,bu = getvarbound(t,j)

                    if bk == MSK_BK_FX
                        println(f,"        $(varnames[j]) = $bl")
                    else
                        blstr =
                            if     bk == MSK_BK_LO || bk == MSK_BK_RA "[$bl"
                            else   "]-Inf"
                            end
                        bustr =
                            if     bk == MSK_BK_UP || bk == MSK_BK_RA "$bu]"
                            else   "+Inf["
                            end

                        println(f,"        $(varnames[j]) ∈ $blstr;$bustr")
                    end
                end
            end
            if limitnumvar < numvar
                println(f,"        ... ($(numvar-limitnumvar) variable bounds omitted)")
            end


            for j in 1:limitnumbarvar
                dim = getdimbarvarj(t,j)
                println(f,"        $(barvarnames[j]) ∈ S($dim)")
            end
            if limitnumbarvar < numbarvar
                println(f,"        ... ($(numbarvar-limitnumbarvar) semidefinite variables omitted)")
            end
        end

        # Matsto
        numsymmat = getnumsymmat(t)
        limitnumsymmat = if limit > 0 min(limit,numsymmat) else numsymmat end
        if numsymmat > 0
            println(f,"    Symmetric matrixes")
            for i in 1:limitnumsymmat
                dim,nz,tp = getsymmatinfo(t,i)
                (msubi,msubj,mvalij) = getsparsesymmat(t,i)
                print(f,"        #MX$i ($dim): [ " )
                if nz > 0
                    limitnz = if limit > 0 min(limit,nz) else nz end
                    print(f,"($(msubi[1]),$(msubj[1]),$(mvalij[1]))")
                    for i in 2:limitnz
                        print(f,",($(msubi[i]),$(msubj[i]),$(mvalij[i]))")
                    end
                    if limitnz < nz
                        print(f,"... ($(nz-limitnz) nonzeros omitted)")
                    end
                end
                println(f," ]")
            end
            if limitnumsymmat < numsymmat
                println(f,"        ... ($(numsymmat-limitnumsymmat) matrixes omitted)")
            end
            
        end
    end
end
