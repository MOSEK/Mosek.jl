# Contents of this file is generated. Do not edit by hand
# Target: Mosek 10.1.11
macro MSK_analyzeproblem(task,whichstream)
  quote
     local res = disable_sigint(()->ccall((:MSK_analyzeproblem,libmosek),Int32,(Ptr{Nothing},Int32,),$(esc(task)),$(esc(whichstream))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_analyzenames(task,whichstream,nametype)
  quote
     local res = disable_sigint(()->ccall((:MSK_analyzenames,libmosek),Int32,(Ptr{Nothing},Int32,Int32,),$(esc(task)),$(esc(whichstream)),$(esc(nametype))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_analyzesolution(task,whichstream,whichsol)
  quote
     local res = disable_sigint(()->ccall((:MSK_analyzesolution,libmosek),Int32,(Ptr{Nothing},Int32,Int32,),$(esc(task)),$(esc(whichstream)),$(esc(whichsol))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_initbasissolve(task,basis)
  quote
     local res = disable_sigint(()->ccall((:MSK_initbasissolve,libmosek),Int32,(Ptr{Nothing},Ptr{Int32},),$(esc(task)),$(esc(basis))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_solvewithbasis(task,transp,numnz,sub,val,numnzout)
  quote
     local res = disable_sigint(()->ccall((:MSK_solvewithbasis,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Ptr{Int32},Ptr{Float64},Ref{Int32},),$(esc(task)),$(esc(transp)),$(esc(numnz)),$(esc(sub)),$(esc(val)),$(esc(numnzout))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_basiscond(task,nrmbasis,nrminvbasis)
  quote
     local res = disable_sigint(()->ccall((:MSK_basiscond,libmosek),Int32,(Ptr{Nothing},Ref{Float64},Ref{Float64},),$(esc(task)),$(esc(nrmbasis)),$(esc(nrminvbasis))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_appendcons(task,num)
  quote
     local res = disable_sigint(()->ccall((:MSK_appendcons,libmosek),Int32,(Ptr{Nothing},Int32,),$(esc(task)),$(esc(num))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_appendvars(task,num)
  quote
     local res = disable_sigint(()->ccall((:MSK_appendvars,libmosek),Int32,(Ptr{Nothing},Int32,),$(esc(task)),$(esc(num))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_removecons(task,num,subset)
  quote
     local res = disable_sigint(()->ccall((:MSK_removecons,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{Int32},),$(esc(task)),$(esc(num)),$(esc(subset))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_removevars(task,num,subset)
  quote
     local res = disable_sigint(()->ccall((:MSK_removevars,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{Int32},),$(esc(task)),$(esc(num)),$(esc(subset))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_removebarvars(task,num,subset)
  quote
     local res = disable_sigint(()->ccall((:MSK_removebarvars,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{Int32},),$(esc(task)),$(esc(num)),$(esc(subset))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_removecones(task,num,subset)
  quote
     local res = disable_sigint(()->ccall((:MSK_removecones,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{Int32},),$(esc(task)),$(esc(num)),$(esc(subset))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_appendbarvars(task,num,dim)
  quote
     local res = disable_sigint(()->ccall((:MSK_appendbarvars,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{Int32},),$(esc(task)),$(esc(num)),$(esc(dim))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_appendcone(task,ct,conepar,nummem,submem)
  quote
     local res = disable_sigint(()->ccall((:MSK_appendcone,libmosek),Int32,(Ptr{Nothing},Int32,Float64,Int32,Ptr{Int32},),$(esc(task)),$(esc(ct)),$(esc(conepar)),$(esc(nummem)),$(esc(submem))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_appendconeseq(task,ct,conepar,nummem,j)
  quote
     local res = disable_sigint(()->ccall((:MSK_appendconeseq,libmosek),Int32,(Ptr{Nothing},Int32,Float64,Int32,Int32,),$(esc(task)),$(esc(ct)),$(esc(conepar)),$(esc(nummem)),$(esc(j))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_appendconesseq(task,num,ct,conepar,nummem,j)
  quote
     local res = disable_sigint(()->ccall((:MSK_appendconesseq,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{Int32},Ptr{Float64},Ptr{Int32},Int32,),$(esc(task)),$(esc(num)),$(esc(ct)),$(esc(conepar)),$(esc(nummem)),$(esc(j))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_bktostr(task,bk,str)
  quote
     local res = disable_sigint(()->ccall((:MSK_bktostr,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{UInt8},),$(esc(task)),$(esc(bk)),$(esc(str))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_chgconbound(task,i,lower,finite,value)
  quote
     local res = disable_sigint(()->ccall((:MSK_chgconbound,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Int32,Float64,),$(esc(task)),$(esc(i)),$(esc(lower)),$(esc(finite)),$(esc(value))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_chgvarbound(task,j,lower,finite,value)
  quote
     local res = disable_sigint(()->ccall((:MSK_chgvarbound,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Int32,Float64,),$(esc(task)),$(esc(j)),$(esc(lower)),$(esc(finite)),$(esc(value))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_conetypetostr(task,ct,str)
  quote
     local res = disable_sigint(()->ccall((:MSK_conetypetostr,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{UInt8},),$(esc(task)),$(esc(ct)),$(esc(str))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_deletetask(task)
  quote
     local res = disable_sigint(()->ccall((:MSK_deletetask,libmosek),Int32,(Ref{Ptr{Nothing}},),$(esc(task))))
     if res != 0
       throw(MosekError(res,""))
     end
     nothing
  end
end
macro MSK_freetask(task,buffer)
  quote
     disable_sigint(()->ccall((:MSK_freetask,libmosek),Cvoid,(Ptr{Nothing},Ptr{Cvoid},),$(esc(task)),$(esc(buffer))))
  end
end
macro MSK_freedbgtask(task,buffer,file,line)
  quote
     disable_sigint(()->ccall((:MSK_freedbgtask,libmosek),Cvoid,(Ptr{Nothing},Ptr{Cvoid},Ptr{UInt8},UInt32,),$(esc(task)),$(esc(buffer)),$(esc(file)),$(esc(line))))
  end
end
macro MSK_getaij(task,i,j,aij)
  quote
     local res = disable_sigint(()->ccall((:MSK_getaij,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Ref{Float64},),$(esc(task)),$(esc(i)),$(esc(j)),$(esc(aij))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getapiecenumnz(task,firsti,lasti,firstj,lastj,numnz)
  quote
     local res = disable_sigint(()->ccall((:MSK_getapiecenumnz,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Int32,Int32,Ref{Int32},),$(esc(task)),$(esc(firsti)),$(esc(lasti)),$(esc(firstj)),$(esc(lastj)),$(esc(numnz))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getacolnumnz(task,i,nzj)
  quote
     local res = disable_sigint(()->ccall((:MSK_getacolnumnz,libmosek),Int32,(Ptr{Nothing},Int32,Ref{Int32},),$(esc(task)),$(esc(i)),$(esc(nzj))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getacol(task,j,nzj,subj,valj)
  quote
     local res = disable_sigint(()->ccall((:MSK_getacol,libmosek),Int32,(Ptr{Nothing},Int32,Ref{Int32},Ptr{Int32},Ptr{Float64},),$(esc(task)),$(esc(j)),$(esc(nzj)),$(esc(subj)),$(esc(valj))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getacolslice(task,first,last,maxnumnz,ptrb,ptre,sub,val)
  quote
     local res = disable_sigint(()->ccall((:MSK_getacolslice,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Int32,Ptr{Int32},Ptr{Int32},Ptr{Int32},Ptr{Float64},),$(esc(task)),$(esc(first)),$(esc(last)),$(esc(maxnumnz)),$(esc(ptrb)),$(esc(ptre)),$(esc(sub)),$(esc(val))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getacolslice64(task,first,last,maxnumnz,ptrb,ptre,sub,val)
  quote
     local res = disable_sigint(()->ccall((:MSK_getacolslice64,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Int64,Ptr{Int64},Ptr{Int64},Ptr{Int32},Ptr{Float64},),$(esc(task)),$(esc(first)),$(esc(last)),$(esc(maxnumnz)),$(esc(ptrb)),$(esc(ptre)),$(esc(sub)),$(esc(val))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getarownumnz(task,i,nzi)
  quote
     local res = disable_sigint(()->ccall((:MSK_getarownumnz,libmosek),Int32,(Ptr{Nothing},Int32,Ref{Int32},),$(esc(task)),$(esc(i)),$(esc(nzi))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getarow(task,i,nzi,subi,vali)
  quote
     local res = disable_sigint(()->ccall((:MSK_getarow,libmosek),Int32,(Ptr{Nothing},Int32,Ref{Int32},Ptr{Int32},Ptr{Float64},),$(esc(task)),$(esc(i)),$(esc(nzi)),$(esc(subi)),$(esc(vali))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getacolslicenumnz(task,first,last,numnz)
  quote
     local res = disable_sigint(()->ccall((:MSK_getacolslicenumnz,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Ref{Int32},),$(esc(task)),$(esc(first)),$(esc(last)),$(esc(numnz))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getacolslicenumnz64(task,first,last,numnz)
  quote
     local res = disable_sigint(()->ccall((:MSK_getacolslicenumnz64,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Ref{Int64},),$(esc(task)),$(esc(first)),$(esc(last)),$(esc(numnz))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getarowslicenumnz(task,first,last,numnz)
  quote
     local res = disable_sigint(()->ccall((:MSK_getarowslicenumnz,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Ref{Int32},),$(esc(task)),$(esc(first)),$(esc(last)),$(esc(numnz))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getarowslicenumnz64(task,first,last,numnz)
  quote
     local res = disable_sigint(()->ccall((:MSK_getarowslicenumnz64,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Ref{Int64},),$(esc(task)),$(esc(first)),$(esc(last)),$(esc(numnz))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getarowslice(task,first,last,maxnumnz,ptrb,ptre,sub,val)
  quote
     local res = disable_sigint(()->ccall((:MSK_getarowslice,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Int32,Ptr{Int32},Ptr{Int32},Ptr{Int32},Ptr{Float64},),$(esc(task)),$(esc(first)),$(esc(last)),$(esc(maxnumnz)),$(esc(ptrb)),$(esc(ptre)),$(esc(sub)),$(esc(val))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getarowslice64(task,first,last,maxnumnz,ptrb,ptre,sub,val)
  quote
     local res = disable_sigint(()->ccall((:MSK_getarowslice64,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Int64,Ptr{Int64},Ptr{Int64},Ptr{Int32},Ptr{Float64},),$(esc(task)),$(esc(first)),$(esc(last)),$(esc(maxnumnz)),$(esc(ptrb)),$(esc(ptre)),$(esc(sub)),$(esc(val))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getatrip(task,maxnumnz,subi,subj,val)
  quote
     local res = disable_sigint(()->ccall((:MSK_getatrip,libmosek),Int32,(Ptr{Nothing},Int64,Ptr{Int32},Ptr{Int32},Ptr{Float64},),$(esc(task)),$(esc(maxnumnz)),$(esc(subi)),$(esc(subj)),$(esc(val))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getarowslicetrip(task,first,last,maxnumnz,subi,subj,val)
  quote
     local res = disable_sigint(()->ccall((:MSK_getarowslicetrip,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Int64,Ptr{Int32},Ptr{Int32},Ptr{Float64},),$(esc(task)),$(esc(first)),$(esc(last)),$(esc(maxnumnz)),$(esc(subi)),$(esc(subj)),$(esc(val))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getacolslicetrip(task,first,last,maxnumnz,subi,subj,val)
  quote
     local res = disable_sigint(()->ccall((:MSK_getacolslicetrip,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Int64,Ptr{Int32},Ptr{Int32},Ptr{Float64},),$(esc(task)),$(esc(first)),$(esc(last)),$(esc(maxnumnz)),$(esc(subi)),$(esc(subj)),$(esc(val))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getconbound(task,i,bk,bl,bu)
  quote
     local res = disable_sigint(()->ccall((:MSK_getconbound,libmosek),Int32,(Ptr{Nothing},Int32,Ref{Int32},Ref{Float64},Ref{Float64},),$(esc(task)),$(esc(i)),$(esc(bk)),$(esc(bl)),$(esc(bu))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getvarbound(task,i,bk,bl,bu)
  quote
     local res = disable_sigint(()->ccall((:MSK_getvarbound,libmosek),Int32,(Ptr{Nothing},Int32,Ref{Int32},Ref{Float64},Ref{Float64},),$(esc(task)),$(esc(i)),$(esc(bk)),$(esc(bl)),$(esc(bu))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getconboundslice(task,first,last,bk,bl,bu)
  quote
     local res = disable_sigint(()->ccall((:MSK_getconboundslice,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Ptr{Int32},Ptr{Float64},Ptr{Float64},),$(esc(task)),$(esc(first)),$(esc(last)),$(esc(bk)),$(esc(bl)),$(esc(bu))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getvarboundslice(task,first,last,bk,bl,bu)
  quote
     local res = disable_sigint(()->ccall((:MSK_getvarboundslice,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Ptr{Int32},Ptr{Float64},Ptr{Float64},),$(esc(task)),$(esc(first)),$(esc(last)),$(esc(bk)),$(esc(bl)),$(esc(bu))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getcj(task,j,cj)
  quote
     local res = disable_sigint(()->ccall((:MSK_getcj,libmosek),Int32,(Ptr{Nothing},Int32,Ref{Float64},),$(esc(task)),$(esc(j)),$(esc(cj))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getc(task,c)
  quote
     local res = disable_sigint(()->ccall((:MSK_getc,libmosek),Int32,(Ptr{Nothing},Ptr{Float64},),$(esc(task)),$(esc(c))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getcallbackfunc(task,func,handle)
  quote
     local res = disable_sigint(()->ccall((:MSK_getcallbackfunc,libmosek),Int32,(Ptr{Nothing},Ref{Ptr{Cvoid}},Ref{Any},),$(esc(task)),$(esc(func)),$(esc(handle))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getcfix(task,cfix)
  quote
     local res = disable_sigint(()->ccall((:MSK_getcfix,libmosek),Int32,(Ptr{Nothing},Ref{Float64},),$(esc(task)),$(esc(cfix))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getcone(task,k,ct,conepar,nummem,submem)
  quote
     local res = disable_sigint(()->ccall((:MSK_getcone,libmosek),Int32,(Ptr{Nothing},Int32,Ref{Int32},Ref{Float64},Ref{Int32},Ptr{Int32},),$(esc(task)),$(esc(k)),$(esc(ct)),$(esc(conepar)),$(esc(nummem)),$(esc(submem))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getconeinfo(task,k,ct,conepar,nummem)
  quote
     local res = disable_sigint(()->ccall((:MSK_getconeinfo,libmosek),Int32,(Ptr{Nothing},Int32,Ref{Int32},Ref{Float64},Ref{Int32},),$(esc(task)),$(esc(k)),$(esc(ct)),$(esc(conepar)),$(esc(nummem))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getclist(task,num,subj,c)
  quote
     local res = disable_sigint(()->ccall((:MSK_getclist,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{Int32},Ptr{Float64},),$(esc(task)),$(esc(num)),$(esc(subj)),$(esc(c))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getcslice(task,first,last,c)
  quote
     local res = disable_sigint(()->ccall((:MSK_getcslice,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Ptr{Float64},),$(esc(task)),$(esc(first)),$(esc(last)),$(esc(c))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getdouinf(task,whichdinf,dvalue)
  quote
     local res = disable_sigint(()->ccall((:MSK_getdouinf,libmosek),Int32,(Ptr{Nothing},Int32,Ref{Float64},),$(esc(task)),$(esc(whichdinf)),$(esc(dvalue))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getdouparam(task,param,parvalue)
  quote
     local res = disable_sigint(()->ccall((:MSK_getdouparam,libmosek),Int32,(Ptr{Nothing},Int32,Ref{Float64},),$(esc(task)),$(esc(param)),$(esc(parvalue))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getdualobj(task,whichsol,dualobj)
  quote
     local res = disable_sigint(()->ccall((:MSK_getdualobj,libmosek),Int32,(Ptr{Nothing},Int32,Ref{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(dualobj))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getenv(task,env)
  quote
     local res = disable_sigint(()->ccall((:MSK_getenv,libmosek),Int32,(Ptr{Nothing},Ref{Ptr{Nothing}},),$(esc(task)),$(esc(env))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getinfindex(task,inftype,infname,infindex)
  quote
     local res = disable_sigint(()->ccall((:MSK_getinfindex,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{UInt8},Ref{Int32},),$(esc(task)),$(esc(inftype)),$(esc(infname)),$(esc(infindex))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getinfmax(task,inftype,infmax)
  quote
     local res = disable_sigint(()->ccall((:MSK_getinfmax,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{Int32},),$(esc(task)),$(esc(inftype)),$(esc(infmax))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getinfname(task,inftype,whichinf,infname)
  quote
     local res = disable_sigint(()->ccall((:MSK_getinfname,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Ptr{UInt8},),$(esc(task)),$(esc(inftype)),$(esc(whichinf)),$(esc(infname))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getintinf(task,whichiinf,ivalue)
  quote
     local res = disable_sigint(()->ccall((:MSK_getintinf,libmosek),Int32,(Ptr{Nothing},Int32,Ref{Int32},),$(esc(task)),$(esc(whichiinf)),$(esc(ivalue))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getlintinf(task,whichliinf,ivalue)
  quote
     local res = disable_sigint(()->ccall((:MSK_getlintinf,libmosek),Int32,(Ptr{Nothing},Int32,Ref{Int64},),$(esc(task)),$(esc(whichliinf)),$(esc(ivalue))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getintparam(task,param,parvalue)
  quote
     local res = disable_sigint(()->ccall((:MSK_getintparam,libmosek),Int32,(Ptr{Nothing},Int32,Ref{Int32},),$(esc(task)),$(esc(param)),$(esc(parvalue))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getmaxnamelen(task,maxlen)
  quote
     local res = disable_sigint(()->ccall((:MSK_getmaxnamelen,libmosek),Int32,(Ptr{Nothing},Ref{Int32},),$(esc(task)),$(esc(maxlen))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getmaxnumanz(task,maxnumanz)
  quote
     local res = disable_sigint(()->ccall((:MSK_getmaxnumanz,libmosek),Int32,(Ptr{Nothing},Ref{Int32},),$(esc(task)),$(esc(maxnumanz))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getmaxnumanz64(task,maxnumanz)
  quote
     local res = disable_sigint(()->ccall((:MSK_getmaxnumanz64,libmosek),Int32,(Ptr{Nothing},Ref{Int64},),$(esc(task)),$(esc(maxnumanz))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getmaxnumcon(task,maxnumcon)
  quote
     local res = disable_sigint(()->ccall((:MSK_getmaxnumcon,libmosek),Int32,(Ptr{Nothing},Ref{Int32},),$(esc(task)),$(esc(maxnumcon))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getmaxnumvar(task,maxnumvar)
  quote
     local res = disable_sigint(()->ccall((:MSK_getmaxnumvar,libmosek),Int32,(Ptr{Nothing},Ref{Int32},),$(esc(task)),$(esc(maxnumvar))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getnadouinf(task,infitemname,dvalue)
  quote
     local res = disable_sigint(()->ccall((:MSK_getnadouinf,libmosek),Int32,(Ptr{Nothing},Ptr{UInt8},Ref{Float64},),$(esc(task)),$(esc(infitemname)),$(esc(dvalue))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getnadouparam(task,paramname,parvalue)
  quote
     local res = disable_sigint(()->ccall((:MSK_getnadouparam,libmosek),Int32,(Ptr{Nothing},Ptr{UInt8},Ref{Float64},),$(esc(task)),$(esc(paramname)),$(esc(parvalue))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getnaintinf(task,infitemname,ivalue)
  quote
     local res = disable_sigint(()->ccall((:MSK_getnaintinf,libmosek),Int32,(Ptr{Nothing},Ptr{UInt8},Ref{Int32},),$(esc(task)),$(esc(infitemname)),$(esc(ivalue))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getnaintparam(task,paramname,parvalue)
  quote
     local res = disable_sigint(()->ccall((:MSK_getnaintparam,libmosek),Int32,(Ptr{Nothing},Ptr{UInt8},Ref{Int32},),$(esc(task)),$(esc(paramname)),$(esc(parvalue))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getbarvarnamelen(task,i,len)
  quote
     local res = disable_sigint(()->ccall((:MSK_getbarvarnamelen,libmosek),Int32,(Ptr{Nothing},Int32,Ref{Int32},),$(esc(task)),$(esc(i)),$(esc(len))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getbarvarname(task,i,sizename,name)
  quote
     local res = disable_sigint(()->ccall((:MSK_getbarvarname,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Ptr{UInt8},),$(esc(task)),$(esc(i)),$(esc(sizename)),$(esc(name))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getbarvarnameindex(task,somename,asgn,index)
  quote
     local res = disable_sigint(()->ccall((:MSK_getbarvarnameindex,libmosek),Int32,(Ptr{Nothing},Ptr{UInt8},Ref{Int32},Ref{Int32},),$(esc(task)),$(esc(somename)),$(esc(asgn)),$(esc(index))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putconname(task,i,name)
  quote
     local res = disable_sigint(()->ccall((:MSK_putconname,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{UInt8},),$(esc(task)),$(esc(i)),$(esc(name))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putvarname(task,j,name)
  quote
     local res = disable_sigint(()->ccall((:MSK_putvarname,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{UInt8},),$(esc(task)),$(esc(j)),$(esc(name))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putconename(task,j,name)
  quote
     local res = disable_sigint(()->ccall((:MSK_putconename,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{UInt8},),$(esc(task)),$(esc(j)),$(esc(name))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putbarvarname(task,j,name)
  quote
     local res = disable_sigint(()->ccall((:MSK_putbarvarname,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{UInt8},),$(esc(task)),$(esc(j)),$(esc(name))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putdomainname(task,domidx,name)
  quote
     local res = disable_sigint(()->ccall((:MSK_putdomainname,libmosek),Int32,(Ptr{Nothing},Int64,Ptr{UInt8},),$(esc(task)),$(esc(domidx)),$(esc(name))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putdjcname(task,djcidx,name)
  quote
     local res = disable_sigint(()->ccall((:MSK_putdjcname,libmosek),Int32,(Ptr{Nothing},Int64,Ptr{UInt8},),$(esc(task)),$(esc(djcidx)),$(esc(name))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putaccname(task,accidx,name)
  quote
     local res = disable_sigint(()->ccall((:MSK_putaccname,libmosek),Int32,(Ptr{Nothing},Int64,Ptr{UInt8},),$(esc(task)),$(esc(accidx)),$(esc(name))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getvarnamelen(task,i,len)
  quote
     local res = disable_sigint(()->ccall((:MSK_getvarnamelen,libmosek),Int32,(Ptr{Nothing},Int32,Ref{Int32},),$(esc(task)),$(esc(i)),$(esc(len))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getvarname(task,j,sizename,name)
  quote
     local res = disable_sigint(()->ccall((:MSK_getvarname,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Ptr{UInt8},),$(esc(task)),$(esc(j)),$(esc(sizename)),$(esc(name))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getconnamelen(task,i,len)
  quote
     local res = disable_sigint(()->ccall((:MSK_getconnamelen,libmosek),Int32,(Ptr{Nothing},Int32,Ref{Int32},),$(esc(task)),$(esc(i)),$(esc(len))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getconname(task,i,sizename,name)
  quote
     local res = disable_sigint(()->ccall((:MSK_getconname,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Ptr{UInt8},),$(esc(task)),$(esc(i)),$(esc(sizename)),$(esc(name))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getconnameindex(task,somename,asgn,index)
  quote
     local res = disable_sigint(()->ccall((:MSK_getconnameindex,libmosek),Int32,(Ptr{Nothing},Ptr{UInt8},Ref{Int32},Ref{Int32},),$(esc(task)),$(esc(somename)),$(esc(asgn)),$(esc(index))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getvarnameindex(task,somename,asgn,index)
  quote
     local res = disable_sigint(()->ccall((:MSK_getvarnameindex,libmosek),Int32,(Ptr{Nothing},Ptr{UInt8},Ref{Int32},Ref{Int32},),$(esc(task)),$(esc(somename)),$(esc(asgn)),$(esc(index))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getconenamelen(task,i,len)
  quote
     local res = disable_sigint(()->ccall((:MSK_getconenamelen,libmosek),Int32,(Ptr{Nothing},Int32,Ref{Int32},),$(esc(task)),$(esc(i)),$(esc(len))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getconename(task,i,sizename,name)
  quote
     local res = disable_sigint(()->ccall((:MSK_getconename,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Ptr{UInt8},),$(esc(task)),$(esc(i)),$(esc(sizename)),$(esc(name))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getconenameindex(task,somename,asgn,index)
  quote
     local res = disable_sigint(()->ccall((:MSK_getconenameindex,libmosek),Int32,(Ptr{Nothing},Ptr{UInt8},Ref{Int32},Ref{Int32},),$(esc(task)),$(esc(somename)),$(esc(asgn)),$(esc(index))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getdomainnamelen(task,domidx,len)
  quote
     local res = disable_sigint(()->ccall((:MSK_getdomainnamelen,libmosek),Int32,(Ptr{Nothing},Int64,Ref{Int32},),$(esc(task)),$(esc(domidx)),$(esc(len))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getdomainname(task,domidx,sizename,name)
  quote
     local res = disable_sigint(()->ccall((:MSK_getdomainname,libmosek),Int32,(Ptr{Nothing},Int64,Int32,Ptr{UInt8},),$(esc(task)),$(esc(domidx)),$(esc(sizename)),$(esc(name))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getdjcnamelen(task,djcidx,len)
  quote
     local res = disable_sigint(()->ccall((:MSK_getdjcnamelen,libmosek),Int32,(Ptr{Nothing},Int64,Ref{Int32},),$(esc(task)),$(esc(djcidx)),$(esc(len))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getdjcname(task,djcidx,sizename,name)
  quote
     local res = disable_sigint(()->ccall((:MSK_getdjcname,libmosek),Int32,(Ptr{Nothing},Int64,Int32,Ptr{UInt8},),$(esc(task)),$(esc(djcidx)),$(esc(sizename)),$(esc(name))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getaccnamelen(task,accidx,len)
  quote
     local res = disable_sigint(()->ccall((:MSK_getaccnamelen,libmosek),Int32,(Ptr{Nothing},Int64,Ref{Int32},),$(esc(task)),$(esc(accidx)),$(esc(len))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getaccname(task,accidx,sizename,name)
  quote
     local res = disable_sigint(()->ccall((:MSK_getaccname,libmosek),Int32,(Ptr{Nothing},Int64,Int32,Ptr{UInt8},),$(esc(task)),$(esc(accidx)),$(esc(sizename)),$(esc(name))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getnastrparam(task,paramname,sizeparamname,len,parvalue)
  quote
     local res = disable_sigint(()->ccall((:MSK_getnastrparam,libmosek),Int32,(Ptr{Nothing},Ptr{UInt8},Int32,Ref{Int32},Ptr{UInt8},),$(esc(task)),$(esc(paramname)),$(esc(sizeparamname)),$(esc(len)),$(esc(parvalue))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getnumanz(task,numanz)
  quote
     local res = disable_sigint(()->ccall((:MSK_getnumanz,libmosek),Int32,(Ptr{Nothing},Ref{Int32},),$(esc(task)),$(esc(numanz))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getnumanz64(task,numanz)
  quote
     local res = disable_sigint(()->ccall((:MSK_getnumanz64,libmosek),Int32,(Ptr{Nothing},Ref{Int64},),$(esc(task)),$(esc(numanz))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getnumcon(task,numcon)
  quote
     local res = disable_sigint(()->ccall((:MSK_getnumcon,libmosek),Int32,(Ptr{Nothing},Ref{Int32},),$(esc(task)),$(esc(numcon))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getnumcone(task,numcone)
  quote
     local res = disable_sigint(()->ccall((:MSK_getnumcone,libmosek),Int32,(Ptr{Nothing},Ref{Int32},),$(esc(task)),$(esc(numcone))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getnumconemem(task,k,nummem)
  quote
     local res = disable_sigint(()->ccall((:MSK_getnumconemem,libmosek),Int32,(Ptr{Nothing},Int32,Ref{Int32},),$(esc(task)),$(esc(k)),$(esc(nummem))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getnumintvar(task,numintvar)
  quote
     local res = disable_sigint(()->ccall((:MSK_getnumintvar,libmosek),Int32,(Ptr{Nothing},Ref{Int32},),$(esc(task)),$(esc(numintvar))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getnumparam(task,partype,numparam)
  quote
     local res = disable_sigint(()->ccall((:MSK_getnumparam,libmosek),Int32,(Ptr{Nothing},Int32,Ref{Int32},),$(esc(task)),$(esc(partype)),$(esc(numparam))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getnumqconknz(task,k,numqcnz)
  quote
     local res = disable_sigint(()->ccall((:MSK_getnumqconknz,libmosek),Int32,(Ptr{Nothing},Int32,Ref{Int32},),$(esc(task)),$(esc(k)),$(esc(numqcnz))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getnumqconknz64(task,k,numqcnz)
  quote
     local res = disable_sigint(()->ccall((:MSK_getnumqconknz64,libmosek),Int32,(Ptr{Nothing},Int32,Ref{Int64},),$(esc(task)),$(esc(k)),$(esc(numqcnz))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getnumqobjnz(task,numqonz)
  quote
     local res = disable_sigint(()->ccall((:MSK_getnumqobjnz,libmosek),Int32,(Ptr{Nothing},Ref{Int32},),$(esc(task)),$(esc(numqonz))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getnumqobjnz64(task,numqonz)
  quote
     local res = disable_sigint(()->ccall((:MSK_getnumqobjnz64,libmosek),Int32,(Ptr{Nothing},Ref{Int64},),$(esc(task)),$(esc(numqonz))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getnumvar(task,numvar)
  quote
     local res = disable_sigint(()->ccall((:MSK_getnumvar,libmosek),Int32,(Ptr{Nothing},Ref{Int32},),$(esc(task)),$(esc(numvar))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getnumbarvar(task,numbarvar)
  quote
     local res = disable_sigint(()->ccall((:MSK_getnumbarvar,libmosek),Int32,(Ptr{Nothing},Ref{Int32},),$(esc(task)),$(esc(numbarvar))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getmaxnumbarvar(task,maxnumbarvar)
  quote
     local res = disable_sigint(()->ccall((:MSK_getmaxnumbarvar,libmosek),Int32,(Ptr{Nothing},Ref{Int32},),$(esc(task)),$(esc(maxnumbarvar))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getdimbarvarj(task,j,dimbarvarj)
  quote
     local res = disable_sigint(()->ccall((:MSK_getdimbarvarj,libmosek),Int32,(Ptr{Nothing},Int32,Ref{Int32},),$(esc(task)),$(esc(j)),$(esc(dimbarvarj))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getlenbarvarj(task,j,lenbarvarj)
  quote
     local res = disable_sigint(()->ccall((:MSK_getlenbarvarj,libmosek),Int32,(Ptr{Nothing},Int32,Ref{Int64},),$(esc(task)),$(esc(j)),$(esc(lenbarvarj))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getobjname(task,sizeobjname,objname)
  quote
     local res = disable_sigint(()->ccall((:MSK_getobjname,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{UInt8},),$(esc(task)),$(esc(sizeobjname)),$(esc(objname))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getobjnamelen(task,len)
  quote
     local res = disable_sigint(()->ccall((:MSK_getobjnamelen,libmosek),Int32,(Ptr{Nothing},Ref{Int32},),$(esc(task)),$(esc(len))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getparamname(task,partype,param,parname)
  quote
     local res = disable_sigint(()->ccall((:MSK_getparamname,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Ptr{UInt8},),$(esc(task)),$(esc(partype)),$(esc(param)),$(esc(parname))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getparammax(task,partype,parammax)
  quote
     local res = disable_sigint(()->ccall((:MSK_getparammax,libmosek),Int32,(Ptr{Nothing},Int32,Ref{Int32},),$(esc(task)),$(esc(partype)),$(esc(parammax))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getprimalobj(task,whichsol,primalobj)
  quote
     local res = disable_sigint(()->ccall((:MSK_getprimalobj,libmosek),Int32,(Ptr{Nothing},Int32,Ref{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(primalobj))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getprobtype(task,probtype)
  quote
     local res = disable_sigint(()->ccall((:MSK_getprobtype,libmosek),Int32,(Ptr{Nothing},Ref{Int32},),$(esc(task)),$(esc(probtype))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getqconk64(task,k,maxnumqcnz,numqcnz,qcsubi,qcsubj,qcval)
  quote
     local res = disable_sigint(()->ccall((:MSK_getqconk64,libmosek),Int32,(Ptr{Nothing},Int32,Int64,Ref{Int64},Ptr{Int32},Ptr{Int32},Ptr{Float64},),$(esc(task)),$(esc(k)),$(esc(maxnumqcnz)),$(esc(numqcnz)),$(esc(qcsubi)),$(esc(qcsubj)),$(esc(qcval))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getqconk(task,k,maxnumqcnz,numqcnz,qcsubi,qcsubj,qcval)
  quote
     local res = disable_sigint(()->ccall((:MSK_getqconk,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Ref{Int32},Ptr{Int32},Ptr{Int32},Ptr{Float64},),$(esc(task)),$(esc(k)),$(esc(maxnumqcnz)),$(esc(numqcnz)),$(esc(qcsubi)),$(esc(qcsubj)),$(esc(qcval))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getqobj(task,maxnumqonz,numqonz,qosubi,qosubj,qoval)
  quote
     local res = disable_sigint(()->ccall((:MSK_getqobj,libmosek),Int32,(Ptr{Nothing},Int32,Ref{Int32},Ptr{Int32},Ptr{Int32},Ptr{Float64},),$(esc(task)),$(esc(maxnumqonz)),$(esc(numqonz)),$(esc(qosubi)),$(esc(qosubj)),$(esc(qoval))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getqobj64(task,maxnumqonz,numqonz,qosubi,qosubj,qoval)
  quote
     local res = disable_sigint(()->ccall((:MSK_getqobj64,libmosek),Int32,(Ptr{Nothing},Int64,Ref{Int64},Ptr{Int32},Ptr{Int32},Ptr{Float64},),$(esc(task)),$(esc(maxnumqonz)),$(esc(numqonz)),$(esc(qosubi)),$(esc(qosubj)),$(esc(qoval))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getqobjij(task,i,j,qoij)
  quote
     local res = disable_sigint(()->ccall((:MSK_getqobjij,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Ref{Float64},),$(esc(task)),$(esc(i)),$(esc(j)),$(esc(qoij))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getsolution(task,whichsol,problemsta,solutionsta,skc,skx,skn,xc,xx,y,slc,suc,slx,sux,snx)
  quote
     local res = disable_sigint(()->ccall((:MSK_getsolution,libmosek),Int32,(Ptr{Nothing},Int32,Ref{Int32},Ref{Int32},Ptr{Int32},Ptr{Int32},Ptr{Int32},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(problemsta)),$(esc(solutionsta)),$(esc(skc)),$(esc(skx)),$(esc(skn)),$(esc(xc)),$(esc(xx)),$(esc(y)),$(esc(slc)),$(esc(suc)),$(esc(slx)),$(esc(sux)),$(esc(snx))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getsolutionnew(task,whichsol,problemsta,solutionsta,skc,skx,skn,xc,xx,y,slc,suc,slx,sux,snx,doty)
  quote
     local res = disable_sigint(()->ccall((:MSK_getsolutionnew,libmosek),Int32,(Ptr{Nothing},Int32,Ref{Int32},Ref{Int32},Ptr{Int32},Ptr{Int32},Ptr{Int32},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(problemsta)),$(esc(solutionsta)),$(esc(skc)),$(esc(skx)),$(esc(skn)),$(esc(xc)),$(esc(xx)),$(esc(y)),$(esc(slc)),$(esc(suc)),$(esc(slx)),$(esc(sux)),$(esc(snx)),$(esc(doty))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getsolsta(task,whichsol,solutionsta)
  quote
     local res = disable_sigint(()->ccall((:MSK_getsolsta,libmosek),Int32,(Ptr{Nothing},Int32,Ref{Int32},),$(esc(task)),$(esc(whichsol)),$(esc(solutionsta))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getprosta(task,whichsol,problemsta)
  quote
     local res = disable_sigint(()->ccall((:MSK_getprosta,libmosek),Int32,(Ptr{Nothing},Int32,Ref{Int32},),$(esc(task)),$(esc(whichsol)),$(esc(problemsta))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getskc(task,whichsol,skc)
  quote
     local res = disable_sigint(()->ccall((:MSK_getskc,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{Int32},),$(esc(task)),$(esc(whichsol)),$(esc(skc))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getskx(task,whichsol,skx)
  quote
     local res = disable_sigint(()->ccall((:MSK_getskx,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{Int32},),$(esc(task)),$(esc(whichsol)),$(esc(skx))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getskn(task,whichsol,skn)
  quote
     local res = disable_sigint(()->ccall((:MSK_getskn,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{Int32},),$(esc(task)),$(esc(whichsol)),$(esc(skn))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getxc(task,whichsol,xc)
  quote
     local res = disable_sigint(()->ccall((:MSK_getxc,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(xc))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getxx(task,whichsol,xx)
  quote
     local res = disable_sigint(()->ccall((:MSK_getxx,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(xx))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_gety(task,whichsol,y)
  quote
     local res = disable_sigint(()->ccall((:MSK_gety,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(y))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getslc(task,whichsol,slc)
  quote
     local res = disable_sigint(()->ccall((:MSK_getslc,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(slc))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getaccdoty(task,whichsol,accidx,doty)
  quote
     local res = disable_sigint(()->ccall((:MSK_getaccdoty,libmosek),Int32,(Ptr{Nothing},Int32,Int64,Ptr{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(accidx)),$(esc(doty))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getaccdotys(task,whichsol,doty)
  quote
     local res = disable_sigint(()->ccall((:MSK_getaccdotys,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(doty))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_evaluateacc(task,whichsol,accidx,activity)
  quote
     local res = disable_sigint(()->ccall((:MSK_evaluateacc,libmosek),Int32,(Ptr{Nothing},Int32,Int64,Ptr{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(accidx)),$(esc(activity))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_evaluateaccs(task,whichsol,activity)
  quote
     local res = disable_sigint(()->ccall((:MSK_evaluateaccs,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(activity))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getsuc(task,whichsol,suc)
  quote
     local res = disable_sigint(()->ccall((:MSK_getsuc,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(suc))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getslx(task,whichsol,slx)
  quote
     local res = disable_sigint(()->ccall((:MSK_getslx,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(slx))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getsux(task,whichsol,sux)
  quote
     local res = disable_sigint(()->ccall((:MSK_getsux,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(sux))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getsnx(task,whichsol,snx)
  quote
     local res = disable_sigint(()->ccall((:MSK_getsnx,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(snx))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getskcslice(task,whichsol,first,last,skc)
  quote
     local res = disable_sigint(()->ccall((:MSK_getskcslice,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Int32,Ptr{Int32},),$(esc(task)),$(esc(whichsol)),$(esc(first)),$(esc(last)),$(esc(skc))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getskxslice(task,whichsol,first,last,skx)
  quote
     local res = disable_sigint(()->ccall((:MSK_getskxslice,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Int32,Ptr{Int32},),$(esc(task)),$(esc(whichsol)),$(esc(first)),$(esc(last)),$(esc(skx))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getxcslice(task,whichsol,first,last,xc)
  quote
     local res = disable_sigint(()->ccall((:MSK_getxcslice,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Int32,Ptr{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(first)),$(esc(last)),$(esc(xc))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getxxslice(task,whichsol,first,last,xx)
  quote
     local res = disable_sigint(()->ccall((:MSK_getxxslice,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Int32,Ptr{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(first)),$(esc(last)),$(esc(xx))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getyslice(task,whichsol,first,last,y)
  quote
     local res = disable_sigint(()->ccall((:MSK_getyslice,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Int32,Ptr{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(first)),$(esc(last)),$(esc(y))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getslcslice(task,whichsol,first,last,slc)
  quote
     local res = disable_sigint(()->ccall((:MSK_getslcslice,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Int32,Ptr{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(first)),$(esc(last)),$(esc(slc))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getsucslice(task,whichsol,first,last,suc)
  quote
     local res = disable_sigint(()->ccall((:MSK_getsucslice,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Int32,Ptr{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(first)),$(esc(last)),$(esc(suc))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getslxslice(task,whichsol,first,last,slx)
  quote
     local res = disable_sigint(()->ccall((:MSK_getslxslice,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Int32,Ptr{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(first)),$(esc(last)),$(esc(slx))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getsuxslice(task,whichsol,first,last,sux)
  quote
     local res = disable_sigint(()->ccall((:MSK_getsuxslice,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Int32,Ptr{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(first)),$(esc(last)),$(esc(sux))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getsnxslice(task,whichsol,first,last,snx)
  quote
     local res = disable_sigint(()->ccall((:MSK_getsnxslice,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Int32,Ptr{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(first)),$(esc(last)),$(esc(snx))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getbarxj(task,whichsol,j,barxj)
  quote
     local res = disable_sigint(()->ccall((:MSK_getbarxj,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Ptr{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(j)),$(esc(barxj))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getbarxslice(task,whichsol,first,last,slicesize,barxslice)
  quote
     local res = disable_sigint(()->ccall((:MSK_getbarxslice,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Int32,Int64,Ptr{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(first)),$(esc(last)),$(esc(slicesize)),$(esc(barxslice))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getbarsj(task,whichsol,j,barsj)
  quote
     local res = disable_sigint(()->ccall((:MSK_getbarsj,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Ptr{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(j)),$(esc(barsj))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getbarsslice(task,whichsol,first,last,slicesize,barsslice)
  quote
     local res = disable_sigint(()->ccall((:MSK_getbarsslice,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Int32,Int64,Ptr{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(first)),$(esc(last)),$(esc(slicesize)),$(esc(barsslice))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putskc(task,whichsol,skc)
  quote
     local res = disable_sigint(()->ccall((:MSK_putskc,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{Int32},),$(esc(task)),$(esc(whichsol)),$(esc(skc))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putskx(task,whichsol,skx)
  quote
     local res = disable_sigint(()->ccall((:MSK_putskx,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{Int32},),$(esc(task)),$(esc(whichsol)),$(esc(skx))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putxc(task,whichsol,xc)
  quote
     local res = disable_sigint(()->ccall((:MSK_putxc,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(xc))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putxx(task,whichsol,xx)
  quote
     local res = disable_sigint(()->ccall((:MSK_putxx,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(xx))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_puty(task,whichsol,y)
  quote
     local res = disable_sigint(()->ccall((:MSK_puty,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(y))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putslc(task,whichsol,slc)
  quote
     local res = disable_sigint(()->ccall((:MSK_putslc,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(slc))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putsuc(task,whichsol,suc)
  quote
     local res = disable_sigint(()->ccall((:MSK_putsuc,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(suc))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putslx(task,whichsol,slx)
  quote
     local res = disable_sigint(()->ccall((:MSK_putslx,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(slx))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putsux(task,whichsol,sux)
  quote
     local res = disable_sigint(()->ccall((:MSK_putsux,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(sux))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putsnx(task,whichsol,sux)
  quote
     local res = disable_sigint(()->ccall((:MSK_putsnx,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(sux))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putaccdoty(task,whichsol,accidx,doty)
  quote
     local res = disable_sigint(()->ccall((:MSK_putaccdoty,libmosek),Int32,(Ptr{Nothing},Int32,Int64,Ptr{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(accidx)),$(esc(doty))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putskcslice(task,whichsol,first,last,skc)
  quote
     local res = disable_sigint(()->ccall((:MSK_putskcslice,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Int32,Ptr{Int32},),$(esc(task)),$(esc(whichsol)),$(esc(first)),$(esc(last)),$(esc(skc))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putskxslice(task,whichsol,first,last,skx)
  quote
     local res = disable_sigint(()->ccall((:MSK_putskxslice,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Int32,Ptr{Int32},),$(esc(task)),$(esc(whichsol)),$(esc(first)),$(esc(last)),$(esc(skx))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putxcslice(task,whichsol,first,last,xc)
  quote
     local res = disable_sigint(()->ccall((:MSK_putxcslice,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Int32,Ptr{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(first)),$(esc(last)),$(esc(xc))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putxxslice(task,whichsol,first,last,xx)
  quote
     local res = disable_sigint(()->ccall((:MSK_putxxslice,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Int32,Ptr{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(first)),$(esc(last)),$(esc(xx))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putyslice(task,whichsol,first,last,y)
  quote
     local res = disable_sigint(()->ccall((:MSK_putyslice,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Int32,Ptr{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(first)),$(esc(last)),$(esc(y))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putslcslice(task,whichsol,first,last,slc)
  quote
     local res = disable_sigint(()->ccall((:MSK_putslcslice,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Int32,Ptr{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(first)),$(esc(last)),$(esc(slc))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putsucslice(task,whichsol,first,last,suc)
  quote
     local res = disable_sigint(()->ccall((:MSK_putsucslice,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Int32,Ptr{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(first)),$(esc(last)),$(esc(suc))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putslxslice(task,whichsol,first,last,slx)
  quote
     local res = disable_sigint(()->ccall((:MSK_putslxslice,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Int32,Ptr{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(first)),$(esc(last)),$(esc(slx))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putsuxslice(task,whichsol,first,last,sux)
  quote
     local res = disable_sigint(()->ccall((:MSK_putsuxslice,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Int32,Ptr{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(first)),$(esc(last)),$(esc(sux))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putsnxslice(task,whichsol,first,last,snx)
  quote
     local res = disable_sigint(()->ccall((:MSK_putsnxslice,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Int32,Ptr{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(first)),$(esc(last)),$(esc(snx))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putbarxj(task,whichsol,j,barxj)
  quote
     local res = disable_sigint(()->ccall((:MSK_putbarxj,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Ptr{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(j)),$(esc(barxj))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putbarsj(task,whichsol,j,barsj)
  quote
     local res = disable_sigint(()->ccall((:MSK_putbarsj,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Ptr{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(j)),$(esc(barsj))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getpviolcon(task,whichsol,num,sub,viol)
  quote
     local res = disable_sigint(()->ccall((:MSK_getpviolcon,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Ptr{Int32},Ptr{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(num)),$(esc(sub)),$(esc(viol))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getpviolvar(task,whichsol,num,sub,viol)
  quote
     local res = disable_sigint(()->ccall((:MSK_getpviolvar,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Ptr{Int32},Ptr{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(num)),$(esc(sub)),$(esc(viol))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getpviolbarvar(task,whichsol,num,sub,viol)
  quote
     local res = disable_sigint(()->ccall((:MSK_getpviolbarvar,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Ptr{Int32},Ptr{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(num)),$(esc(sub)),$(esc(viol))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getpviolcones(task,whichsol,num,sub,viol)
  quote
     local res = disable_sigint(()->ccall((:MSK_getpviolcones,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Ptr{Int32},Ptr{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(num)),$(esc(sub)),$(esc(viol))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getpviolacc(task,whichsol,numaccidx,accidxlist,viol)
  quote
     local res = disable_sigint(()->ccall((:MSK_getpviolacc,libmosek),Int32,(Ptr{Nothing},Int32,Int64,Ptr{Int64},Ptr{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(numaccidx)),$(esc(accidxlist)),$(esc(viol))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getpvioldjc(task,whichsol,numdjcidx,djcidxlist,viol)
  quote
     local res = disable_sigint(()->ccall((:MSK_getpvioldjc,libmosek),Int32,(Ptr{Nothing},Int32,Int64,Ptr{Int64},Ptr{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(numdjcidx)),$(esc(djcidxlist)),$(esc(viol))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getdviolcon(task,whichsol,num,sub,viol)
  quote
     local res = disable_sigint(()->ccall((:MSK_getdviolcon,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Ptr{Int32},Ptr{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(num)),$(esc(sub)),$(esc(viol))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getdviolvar(task,whichsol,num,sub,viol)
  quote
     local res = disable_sigint(()->ccall((:MSK_getdviolvar,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Ptr{Int32},Ptr{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(num)),$(esc(sub)),$(esc(viol))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getdviolbarvar(task,whichsol,num,sub,viol)
  quote
     local res = disable_sigint(()->ccall((:MSK_getdviolbarvar,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Ptr{Int32},Ptr{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(num)),$(esc(sub)),$(esc(viol))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getdviolcones(task,whichsol,num,sub,viol)
  quote
     local res = disable_sigint(()->ccall((:MSK_getdviolcones,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Ptr{Int32},Ptr{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(num)),$(esc(sub)),$(esc(viol))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getdviolacc(task,whichsol,numaccidx,accidxlist,viol)
  quote
     local res = disable_sigint(()->ccall((:MSK_getdviolacc,libmosek),Int32,(Ptr{Nothing},Int32,Int64,Ptr{Int64},Ptr{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(numaccidx)),$(esc(accidxlist)),$(esc(viol))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getsolutioninfo(task,whichsol,pobj,pviolcon,pviolvar,pviolbarvar,pviolcone,pviolitg,dobj,dviolcon,dviolvar,dviolbarvar,dviolcone)
  quote
     local res = disable_sigint(()->ccall((:MSK_getsolutioninfo,libmosek),Int32,(Ptr{Nothing},Int32,Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(pobj)),$(esc(pviolcon)),$(esc(pviolvar)),$(esc(pviolbarvar)),$(esc(pviolcone)),$(esc(pviolitg)),$(esc(dobj)),$(esc(dviolcon)),$(esc(dviolvar)),$(esc(dviolbarvar)),$(esc(dviolcone))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getsolutioninfonew(task,whichsol,pobj,pviolcon,pviolvar,pviolbarvar,pviolcone,pviolacc,pvioldjc,pviolitg,dobj,dviolcon,dviolvar,dviolbarvar,dviolcone,dviolacc)
  quote
     local res = disable_sigint(()->ccall((:MSK_getsolutioninfonew,libmosek),Int32,(Ptr{Nothing},Int32,Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(pobj)),$(esc(pviolcon)),$(esc(pviolvar)),$(esc(pviolbarvar)),$(esc(pviolcone)),$(esc(pviolacc)),$(esc(pvioldjc)),$(esc(pviolitg)),$(esc(dobj)),$(esc(dviolcon)),$(esc(dviolvar)),$(esc(dviolbarvar)),$(esc(dviolcone)),$(esc(dviolacc))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getdualsolutionnorms(task,whichsol,nrmy,nrmslc,nrmsuc,nrmslx,nrmsux,nrmsnx,nrmbars)
  quote
     local res = disable_sigint(()->ccall((:MSK_getdualsolutionnorms,libmosek),Int32,(Ptr{Nothing},Int32,Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(nrmy)),$(esc(nrmslc)),$(esc(nrmsuc)),$(esc(nrmslx)),$(esc(nrmsux)),$(esc(nrmsnx)),$(esc(nrmbars))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getprimalsolutionnorms(task,whichsol,nrmxc,nrmxx,nrmbarx)
  quote
     local res = disable_sigint(()->ccall((:MSK_getprimalsolutionnorms,libmosek),Int32,(Ptr{Nothing},Int32,Ref{Float64},Ref{Float64},Ref{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(nrmxc)),$(esc(nrmxx)),$(esc(nrmbarx))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getsolutionslice(task,whichsol,solitem,first,last,values)
  quote
     local res = disable_sigint(()->ccall((:MSK_getsolutionslice,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Int32,Int32,Ptr{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(solitem)),$(esc(first)),$(esc(last)),$(esc(values))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getreducedcosts(task,whichsol,first,last,redcosts)
  quote
     local res = disable_sigint(()->ccall((:MSK_getreducedcosts,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Int32,Ptr{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(first)),$(esc(last)),$(esc(redcosts))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getstrparam(task,param,maxlen,len,parvalue)
  quote
     local res = disable_sigint(()->ccall((:MSK_getstrparam,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Ref{Int32},Ptr{UInt8},),$(esc(task)),$(esc(param)),$(esc(maxlen)),$(esc(len)),$(esc(parvalue))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getstrparamlen(task,param,len)
  quote
     local res = disable_sigint(()->ccall((:MSK_getstrparamlen,libmosek),Int32,(Ptr{Nothing},Int32,Ref{Int32},),$(esc(task)),$(esc(param)),$(esc(len))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getsymbcon(task,i,sizevalue,name,value)
  quote
     local res = disable_sigint(()->ccall((:MSK_getsymbcon,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Ptr{UInt8},Ref{Int32},),$(esc(task)),$(esc(i)),$(esc(sizevalue)),$(esc(name)),$(esc(value))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_gettasknamelen(task,len)
  quote
     local res = disable_sigint(()->ccall((:MSK_gettasknamelen,libmosek),Int32,(Ptr{Nothing},Ref{Int32},),$(esc(task)),$(esc(len))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_gettaskname(task,sizetaskname,taskname)
  quote
     local res = disable_sigint(()->ccall((:MSK_gettaskname,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{UInt8},),$(esc(task)),$(esc(sizetaskname)),$(esc(taskname))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getvartype(task,j,vartype)
  quote
     local res = disable_sigint(()->ccall((:MSK_getvartype,libmosek),Int32,(Ptr{Nothing},Int32,Ref{Int32},),$(esc(task)),$(esc(j)),$(esc(vartype))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getvartypelist(task,num,subj,vartype)
  quote
     local res = disable_sigint(()->ccall((:MSK_getvartypelist,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{Int32},Ptr{Int32},),$(esc(task)),$(esc(num)),$(esc(subj)),$(esc(vartype))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_inputdata(task,maxnumcon,maxnumvar,numcon,numvar,c,cfix,aptrb,aptre,asub,aval,bkc,blc,buc,bkx,blx,bux)
  quote
     local res = disable_sigint(()->ccall((:MSK_inputdata,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Int32,Int32,Ptr{Float64},Float64,Ptr{Int32},Ptr{Int32},Ptr{Int32},Ptr{Float64},Ptr{Int32},Ptr{Float64},Ptr{Float64},Ptr{Int32},Ptr{Float64},Ptr{Float64},),$(esc(task)),$(esc(maxnumcon)),$(esc(maxnumvar)),$(esc(numcon)),$(esc(numvar)),$(esc(c)),$(esc(cfix)),$(esc(aptrb)),$(esc(aptre)),$(esc(asub)),$(esc(aval)),$(esc(bkc)),$(esc(blc)),$(esc(buc)),$(esc(bkx)),$(esc(blx)),$(esc(bux))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_inputdata64(task,maxnumcon,maxnumvar,numcon,numvar,c,cfix,aptrb,aptre,asub,aval,bkc,blc,buc,bkx,blx,bux)
  quote
     local res = disable_sigint(()->ccall((:MSK_inputdata64,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Int32,Int32,Ptr{Float64},Float64,Ptr{Int64},Ptr{Int64},Ptr{Int32},Ptr{Float64},Ptr{Int32},Ptr{Float64},Ptr{Float64},Ptr{Int32},Ptr{Float64},Ptr{Float64},),$(esc(task)),$(esc(maxnumcon)),$(esc(maxnumvar)),$(esc(numcon)),$(esc(numvar)),$(esc(c)),$(esc(cfix)),$(esc(aptrb)),$(esc(aptre)),$(esc(asub)),$(esc(aval)),$(esc(bkc)),$(esc(blc)),$(esc(buc)),$(esc(bkx)),$(esc(blx)),$(esc(bux))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_isdouparname(task,parname,param)
  quote
     local res = disable_sigint(()->ccall((:MSK_isdouparname,libmosek),Int32,(Ptr{Nothing},Ptr{UInt8},Ref{Int32},),$(esc(task)),$(esc(parname)),$(esc(param))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_isintparname(task,parname,param)
  quote
     local res = disable_sigint(()->ccall((:MSK_isintparname,libmosek),Int32,(Ptr{Nothing},Ptr{UInt8},Ref{Int32},),$(esc(task)),$(esc(parname)),$(esc(param))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_isstrparname(task,parname,param)
  quote
     local res = disable_sigint(()->ccall((:MSK_isstrparname,libmosek),Int32,(Ptr{Nothing},Ptr{UInt8},Ref{Int32},),$(esc(task)),$(esc(parname)),$(esc(param))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_linkfiletotaskstream(task,whichstream,filename,append)
  quote
     local res = disable_sigint(()->ccall((:MSK_linkfiletotaskstream,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{UInt8},Int32,),$(esc(task)),$(esc(whichstream)),$(esc(filename)),$(esc(append))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_linkfunctotaskstream(task,whichstream,handle,func)
  quote
     local res = disable_sigint(()->ccall((:MSK_linkfunctotaskstream,libmosek),Int32,(Ptr{Nothing},Int32,Any,Ptr{Cvoid},),$(esc(task)),$(esc(whichstream)),$(esc(handle)),$(esc(func))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_unlinkfuncfromtaskstream(task,whichstream)
  quote
     local res = disable_sigint(()->ccall((:MSK_unlinkfuncfromtaskstream,libmosek),Int32,(Ptr{Nothing},Int32,),$(esc(task)),$(esc(whichstream))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_clonetask(task,clonedtask)
  quote
     local res = disable_sigint(()->ccall((:MSK_clonetask,libmosek),Int32,(Ptr{Nothing},Ref{Ptr{Nothing}},),$(esc(task)),$(esc(clonedtask))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_primalrepair(task,wlc,wuc,wlx,wux)
  quote
     local res = disable_sigint(()->ccall((:MSK_primalrepair,libmosek),Int32,(Ptr{Nothing},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},),$(esc(task)),$(esc(wlc)),$(esc(wuc)),$(esc(wlx)),$(esc(wux))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_infeasibilityreport(task,whichstream,whichsol)
  quote
     local res = disable_sigint(()->ccall((:MSK_infeasibilityreport,libmosek),Int32,(Ptr{Nothing},Int32,Int32,),$(esc(task)),$(esc(whichstream)),$(esc(whichsol))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_toconic(task)
  quote
     local res = disable_sigint(()->ccall((:MSK_toconic,libmosek),Int32,(Ptr{Nothing},),$(esc(task))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_optimize(task)
  quote
     local res = disable_sigint(()->ccall((:MSK_optimize,libmosek),Int32,(Ptr{Nothing},),$(esc(task))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_optimizetrm(task,trmcode)
  quote
     local res = disable_sigint(()->ccall((:MSK_optimizetrm,libmosek),Int32,(Ptr{Nothing},Ref{Int32},),$(esc(task)),$(esc(trmcode))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_printparam(task)
  quote
     local res = disable_sigint(()->ccall((:MSK_printparam,libmosek),Int32,(Ptr{Nothing},),$(esc(task))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_probtypetostr(task,probtype,str)
  quote
     local res = disable_sigint(()->ccall((:MSK_probtypetostr,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{UInt8},),$(esc(task)),$(esc(probtype)),$(esc(str))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_prostatostr(task,problemsta,str)
  quote
     local res = disable_sigint(()->ccall((:MSK_prostatostr,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{UInt8},),$(esc(task)),$(esc(problemsta)),$(esc(str))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putresponsefunc(task,responsefunc,handle)
  quote
     local res = disable_sigint(()->ccall((:MSK_putresponsefunc,libmosek),Int32,(Ptr{Nothing},Ptr{Cvoid},Any,),$(esc(task)),$(esc(responsefunc)),$(esc(handle))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_commitchanges(task)
  quote
     local res = disable_sigint(()->ccall((:MSK_commitchanges,libmosek),Int32,(Ptr{Nothing},),$(esc(task))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getatruncatetol(task,tolzero)
  quote
     local res = disable_sigint(()->ccall((:MSK_getatruncatetol,libmosek),Int32,(Ptr{Nothing},Ptr{Float64},),$(esc(task)),$(esc(tolzero))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putatruncatetol(task,tolzero)
  quote
     local res = disable_sigint(()->ccall((:MSK_putatruncatetol,libmosek),Int32,(Ptr{Nothing},Float64,),$(esc(task)),$(esc(tolzero))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putaij(task,i,j,aij)
  quote
     local res = disable_sigint(()->ccall((:MSK_putaij,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Float64,),$(esc(task)),$(esc(i)),$(esc(j)),$(esc(aij))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putaijlist(task,num,subi,subj,valij)
  quote
     local res = disable_sigint(()->ccall((:MSK_putaijlist,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{Int32},Ptr{Int32},Ptr{Float64},),$(esc(task)),$(esc(num)),$(esc(subi)),$(esc(subj)),$(esc(valij))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putaijlist64(task,num,subi,subj,valij)
  quote
     local res = disable_sigint(()->ccall((:MSK_putaijlist64,libmosek),Int32,(Ptr{Nothing},Int64,Ptr{Int32},Ptr{Int32},Ptr{Float64},),$(esc(task)),$(esc(num)),$(esc(subi)),$(esc(subj)),$(esc(valij))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putacol(task,j,nzj,subj,valj)
  quote
     local res = disable_sigint(()->ccall((:MSK_putacol,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Ptr{Int32},Ptr{Float64},),$(esc(task)),$(esc(j)),$(esc(nzj)),$(esc(subj)),$(esc(valj))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putarow(task,i,nzi,subi,vali)
  quote
     local res = disable_sigint(()->ccall((:MSK_putarow,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Ptr{Int32},Ptr{Float64},),$(esc(task)),$(esc(i)),$(esc(nzi)),$(esc(subi)),$(esc(vali))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putarowslice(task,first,last,ptrb,ptre,asub,aval)
  quote
     local res = disable_sigint(()->ccall((:MSK_putarowslice,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Ptr{Int32},Ptr{Int32},Ptr{Int32},Ptr{Float64},),$(esc(task)),$(esc(first)),$(esc(last)),$(esc(ptrb)),$(esc(ptre)),$(esc(asub)),$(esc(aval))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putarowslice64(task,first,last,ptrb,ptre,asub,aval)
  quote
     local res = disable_sigint(()->ccall((:MSK_putarowslice64,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Ptr{Int64},Ptr{Int64},Ptr{Int32},Ptr{Float64},),$(esc(task)),$(esc(first)),$(esc(last)),$(esc(ptrb)),$(esc(ptre)),$(esc(asub)),$(esc(aval))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putarowlist(task,num,sub,ptrb,ptre,asub,aval)
  quote
     local res = disable_sigint(()->ccall((:MSK_putarowlist,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{Int32},Ptr{Int32},Ptr{Int32},Ptr{Int32},Ptr{Float64},),$(esc(task)),$(esc(num)),$(esc(sub)),$(esc(ptrb)),$(esc(ptre)),$(esc(asub)),$(esc(aval))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putarowlist64(task,num,sub,ptrb,ptre,asub,aval)
  quote
     local res = disable_sigint(()->ccall((:MSK_putarowlist64,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{Int32},Ptr{Int64},Ptr{Int64},Ptr{Int32},Ptr{Float64},),$(esc(task)),$(esc(num)),$(esc(sub)),$(esc(ptrb)),$(esc(ptre)),$(esc(asub)),$(esc(aval))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putacolslice(task,first,last,ptrb,ptre,asub,aval)
  quote
     local res = disable_sigint(()->ccall((:MSK_putacolslice,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Ptr{Int32},Ptr{Int32},Ptr{Int32},Ptr{Float64},),$(esc(task)),$(esc(first)),$(esc(last)),$(esc(ptrb)),$(esc(ptre)),$(esc(asub)),$(esc(aval))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putacolslice64(task,first,last,ptrb,ptre,asub,aval)
  quote
     local res = disable_sigint(()->ccall((:MSK_putacolslice64,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Ptr{Int64},Ptr{Int64},Ptr{Int32},Ptr{Float64},),$(esc(task)),$(esc(first)),$(esc(last)),$(esc(ptrb)),$(esc(ptre)),$(esc(asub)),$(esc(aval))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putacollist(task,num,sub,ptrb,ptre,asub,aval)
  quote
     local res = disable_sigint(()->ccall((:MSK_putacollist,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{Int32},Ptr{Int32},Ptr{Int32},Ptr{Int32},Ptr{Float64},),$(esc(task)),$(esc(num)),$(esc(sub)),$(esc(ptrb)),$(esc(ptre)),$(esc(asub)),$(esc(aval))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putacollist64(task,num,sub,ptrb,ptre,asub,aval)
  quote
     local res = disable_sigint(()->ccall((:MSK_putacollist64,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{Int32},Ptr{Int64},Ptr{Int64},Ptr{Int32},Ptr{Float64},),$(esc(task)),$(esc(num)),$(esc(sub)),$(esc(ptrb)),$(esc(ptre)),$(esc(asub)),$(esc(aval))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putbaraij(task,i,j,num,sub,weights)
  quote
     local res = disable_sigint(()->ccall((:MSK_putbaraij,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Int64,Ptr{Int64},Ptr{Float64},),$(esc(task)),$(esc(i)),$(esc(j)),$(esc(num)),$(esc(sub)),$(esc(weights))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putbaraijlist(task,num,subi,subj,alphaptrb,alphaptre,matidx,weights)
  quote
     local res = disable_sigint(()->ccall((:MSK_putbaraijlist,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{Int32},Ptr{Int32},Ptr{Int64},Ptr{Int64},Ptr{Int64},Ptr{Float64},),$(esc(task)),$(esc(num)),$(esc(subi)),$(esc(subj)),$(esc(alphaptrb)),$(esc(alphaptre)),$(esc(matidx)),$(esc(weights))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putbararowlist(task,num,subi,ptrb,ptre,subj,nummat,matidx,weights)
  quote
     local res = disable_sigint(()->ccall((:MSK_putbararowlist,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{Int32},Ptr{Int64},Ptr{Int64},Ptr{Int32},Ptr{Int64},Ptr{Int64},Ptr{Float64},),$(esc(task)),$(esc(num)),$(esc(subi)),$(esc(ptrb)),$(esc(ptre)),$(esc(subj)),$(esc(nummat)),$(esc(matidx)),$(esc(weights))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getnumbarcnz(task,nz)
  quote
     local res = disable_sigint(()->ccall((:MSK_getnumbarcnz,libmosek),Int32,(Ptr{Nothing},Ref{Int64},),$(esc(task)),$(esc(nz))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getnumbaranz(task,nz)
  quote
     local res = disable_sigint(()->ccall((:MSK_getnumbaranz,libmosek),Int32,(Ptr{Nothing},Ref{Int64},),$(esc(task)),$(esc(nz))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getbarcsparsity(task,maxnumnz,numnz,idxj)
  quote
     local res = disable_sigint(()->ccall((:MSK_getbarcsparsity,libmosek),Int32,(Ptr{Nothing},Int64,Ref{Int64},Ptr{Int64},),$(esc(task)),$(esc(maxnumnz)),$(esc(numnz)),$(esc(idxj))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getbarasparsity(task,maxnumnz,numnz,idxij)
  quote
     local res = disable_sigint(()->ccall((:MSK_getbarasparsity,libmosek),Int32,(Ptr{Nothing},Int64,Ref{Int64},Ptr{Int64},),$(esc(task)),$(esc(maxnumnz)),$(esc(numnz)),$(esc(idxij))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getbarcidxinfo(task,idx,num)
  quote
     local res = disable_sigint(()->ccall((:MSK_getbarcidxinfo,libmosek),Int32,(Ptr{Nothing},Int64,Ref{Int64},),$(esc(task)),$(esc(idx)),$(esc(num))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getbarcidxj(task,idx,j)
  quote
     local res = disable_sigint(()->ccall((:MSK_getbarcidxj,libmosek),Int32,(Ptr{Nothing},Int64,Ref{Int32},),$(esc(task)),$(esc(idx)),$(esc(j))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getbarcidx(task,idx,maxnum,j,num,sub,weights)
  quote
     local res = disable_sigint(()->ccall((:MSK_getbarcidx,libmosek),Int32,(Ptr{Nothing},Int64,Int64,Ref{Int32},Ref{Int64},Ptr{Int64},Ptr{Float64},),$(esc(task)),$(esc(idx)),$(esc(maxnum)),$(esc(j)),$(esc(num)),$(esc(sub)),$(esc(weights))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getbaraidxinfo(task,idx,num)
  quote
     local res = disable_sigint(()->ccall((:MSK_getbaraidxinfo,libmosek),Int32,(Ptr{Nothing},Int64,Ref{Int64},),$(esc(task)),$(esc(idx)),$(esc(num))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getbaraidxij(task,idx,i,j)
  quote
     local res = disable_sigint(()->ccall((:MSK_getbaraidxij,libmosek),Int32,(Ptr{Nothing},Int64,Ref{Int32},Ref{Int32},),$(esc(task)),$(esc(idx)),$(esc(i)),$(esc(j))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getbaraidx(task,idx,maxnum,i,j,num,sub,weights)
  quote
     local res = disable_sigint(()->ccall((:MSK_getbaraidx,libmosek),Int32,(Ptr{Nothing},Int64,Int64,Ref{Int32},Ref{Int32},Ref{Int64},Ptr{Int64},Ptr{Float64},),$(esc(task)),$(esc(idx)),$(esc(maxnum)),$(esc(i)),$(esc(j)),$(esc(num)),$(esc(sub)),$(esc(weights))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getnumbarcblocktriplets(task,num)
  quote
     local res = disable_sigint(()->ccall((:MSK_getnumbarcblocktriplets,libmosek),Int32,(Ptr{Nothing},Ref{Int64},),$(esc(task)),$(esc(num))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putbarcblocktriplet(task,num,subj,subk,subl,valjkl)
  quote
     local res = disable_sigint(()->ccall((:MSK_putbarcblocktriplet,libmosek),Int32,(Ptr{Nothing},Int64,Ptr{Int32},Ptr{Int32},Ptr{Int32},Ptr{Float64},),$(esc(task)),$(esc(num)),$(esc(subj)),$(esc(subk)),$(esc(subl)),$(esc(valjkl))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getbarcblocktriplet(task,maxnum,num,subj,subk,subl,valjkl)
  quote
     local res = disable_sigint(()->ccall((:MSK_getbarcblocktriplet,libmosek),Int32,(Ptr{Nothing},Int64,Ref{Int64},Ptr{Int32},Ptr{Int32},Ptr{Int32},Ptr{Float64},),$(esc(task)),$(esc(maxnum)),$(esc(num)),$(esc(subj)),$(esc(subk)),$(esc(subl)),$(esc(valjkl))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putbarablocktriplet(task,num,subi,subj,subk,subl,valijkl)
  quote
     local res = disable_sigint(()->ccall((:MSK_putbarablocktriplet,libmosek),Int32,(Ptr{Nothing},Int64,Ptr{Int32},Ptr{Int32},Ptr{Int32},Ptr{Int32},Ptr{Float64},),$(esc(task)),$(esc(num)),$(esc(subi)),$(esc(subj)),$(esc(subk)),$(esc(subl)),$(esc(valijkl))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getnumbarablocktriplets(task,num)
  quote
     local res = disable_sigint(()->ccall((:MSK_getnumbarablocktriplets,libmosek),Int32,(Ptr{Nothing},Ref{Int64},),$(esc(task)),$(esc(num))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getbarablocktriplet(task,maxnum,num,subi,subj,subk,subl,valijkl)
  quote
     local res = disable_sigint(()->ccall((:MSK_getbarablocktriplet,libmosek),Int32,(Ptr{Nothing},Int64,Ref{Int64},Ptr{Int32},Ptr{Int32},Ptr{Int32},Ptr{Int32},Ptr{Float64},),$(esc(task)),$(esc(maxnum)),$(esc(num)),$(esc(subi)),$(esc(subj)),$(esc(subk)),$(esc(subl)),$(esc(valijkl))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putmaxnumafe(task,maxnumafe)
  quote
     local res = disable_sigint(()->ccall((:MSK_putmaxnumafe,libmosek),Int32,(Ptr{Nothing},Int64,),$(esc(task)),$(esc(maxnumafe))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getnumafe(task,numafe)
  quote
     local res = disable_sigint(()->ccall((:MSK_getnumafe,libmosek),Int32,(Ptr{Nothing},Ref{Int64},),$(esc(task)),$(esc(numafe))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_appendafes(task,num)
  quote
     local res = disable_sigint(()->ccall((:MSK_appendafes,libmosek),Int32,(Ptr{Nothing},Int64,),$(esc(task)),$(esc(num))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putafefentry(task,afeidx,varidx,value)
  quote
     local res = disable_sigint(()->ccall((:MSK_putafefentry,libmosek),Int32,(Ptr{Nothing},Int64,Int32,Float64,),$(esc(task)),$(esc(afeidx)),$(esc(varidx)),$(esc(value))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putafefentrylist(task,numentr,afeidx,varidx,val)
  quote
     local res = disable_sigint(()->ccall((:MSK_putafefentrylist,libmosek),Int32,(Ptr{Nothing},Int64,Ptr{Int64},Ptr{Int32},Ptr{Float64},),$(esc(task)),$(esc(numentr)),$(esc(afeidx)),$(esc(varidx)),$(esc(val))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_emptyafefrow(task,afeidx)
  quote
     local res = disable_sigint(()->ccall((:MSK_emptyafefrow,libmosek),Int32,(Ptr{Nothing},Int64,),$(esc(task)),$(esc(afeidx))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_emptyafefcol(task,varidx)
  quote
     local res = disable_sigint(()->ccall((:MSK_emptyafefcol,libmosek),Int32,(Ptr{Nothing},Int32,),$(esc(task)),$(esc(varidx))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_emptyafefrowlist(task,numafeidx,afeidx)
  quote
     local res = disable_sigint(()->ccall((:MSK_emptyafefrowlist,libmosek),Int32,(Ptr{Nothing},Int64,Ptr{Int64},),$(esc(task)),$(esc(numafeidx)),$(esc(afeidx))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_emptyafefcollist(task,numvaridx,varidx)
  quote
     local res = disable_sigint(()->ccall((:MSK_emptyafefcollist,libmosek),Int32,(Ptr{Nothing},Int64,Ptr{Int32},),$(esc(task)),$(esc(numvaridx)),$(esc(varidx))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putafefrow(task,afeidx,numnz,varidx,val)
  quote
     local res = disable_sigint(()->ccall((:MSK_putafefrow,libmosek),Int32,(Ptr{Nothing},Int64,Int32,Ptr{Int32},Ptr{Float64},),$(esc(task)),$(esc(afeidx)),$(esc(numnz)),$(esc(varidx)),$(esc(val))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putafefrowlist(task,numafeidx,afeidx,numnzrow,ptrrow,lenidxval,varidx,val)
  quote
     local res = disable_sigint(()->ccall((:MSK_putafefrowlist,libmosek),Int32,(Ptr{Nothing},Int64,Ptr{Int64},Ptr{Int32},Ptr{Int64},Int64,Ptr{Int32},Ptr{Float64},),$(esc(task)),$(esc(numafeidx)),$(esc(afeidx)),$(esc(numnzrow)),$(esc(ptrrow)),$(esc(lenidxval)),$(esc(varidx)),$(esc(val))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putafefcol(task,varidx,numnz,afeidx,val)
  quote
     local res = disable_sigint(()->ccall((:MSK_putafefcol,libmosek),Int32,(Ptr{Nothing},Int32,Int64,Ptr{Int64},Ptr{Float64},),$(esc(task)),$(esc(varidx)),$(esc(numnz)),$(esc(afeidx)),$(esc(val))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getafefrownumnz(task,afeidx,numnz)
  quote
     local res = disable_sigint(()->ccall((:MSK_getafefrownumnz,libmosek),Int32,(Ptr{Nothing},Int64,Ref{Int32},),$(esc(task)),$(esc(afeidx)),$(esc(numnz))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getafefnumnz(task,numnz)
  quote
     local res = disable_sigint(()->ccall((:MSK_getafefnumnz,libmosek),Int32,(Ptr{Nothing},Ref{Int64},),$(esc(task)),$(esc(numnz))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getafefrow(task,afeidx,numnz,varidx,val)
  quote
     local res = disable_sigint(()->ccall((:MSK_getafefrow,libmosek),Int32,(Ptr{Nothing},Int64,Ref{Int32},Ptr{Int32},Ptr{Float64},),$(esc(task)),$(esc(afeidx)),$(esc(numnz)),$(esc(varidx)),$(esc(val))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getafeftrip(task,afeidx,varidx,val)
  quote
     local res = disable_sigint(()->ccall((:MSK_getafeftrip,libmosek),Int32,(Ptr{Nothing},Ptr{Int64},Ptr{Int32},Ptr{Float64},),$(esc(task)),$(esc(afeidx)),$(esc(varidx)),$(esc(val))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putafebarfentry(task,afeidx,barvaridx,numterm,termidx,termweight)
  quote
     local res = disable_sigint(()->ccall((:MSK_putafebarfentry,libmosek),Int32,(Ptr{Nothing},Int64,Int32,Int64,Ptr{Int64},Ptr{Float64},),$(esc(task)),$(esc(afeidx)),$(esc(barvaridx)),$(esc(numterm)),$(esc(termidx)),$(esc(termweight))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putafebarfentrylist(task,numafeidx,afeidx,barvaridx,numterm,ptrterm,lenterm,termidx,termweight)
  quote
     local res = disable_sigint(()->ccall((:MSK_putafebarfentrylist,libmosek),Int32,(Ptr{Nothing},Int64,Ptr{Int64},Ptr{Int32},Ptr{Int64},Ptr{Int64},Int64,Ptr{Int64},Ptr{Float64},),$(esc(task)),$(esc(numafeidx)),$(esc(afeidx)),$(esc(barvaridx)),$(esc(numterm)),$(esc(ptrterm)),$(esc(lenterm)),$(esc(termidx)),$(esc(termweight))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putafebarfrow(task,afeidx,numentr,barvaridx,numterm,ptrterm,lenterm,termidx,termweight)
  quote
     local res = disable_sigint(()->ccall((:MSK_putafebarfrow,libmosek),Int32,(Ptr{Nothing},Int64,Int32,Ptr{Int32},Ptr{Int64},Ptr{Int64},Int64,Ptr{Int64},Ptr{Float64},),$(esc(task)),$(esc(afeidx)),$(esc(numentr)),$(esc(barvaridx)),$(esc(numterm)),$(esc(ptrterm)),$(esc(lenterm)),$(esc(termidx)),$(esc(termweight))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_emptyafebarfrow(task,afeidx)
  quote
     local res = disable_sigint(()->ccall((:MSK_emptyafebarfrow,libmosek),Int32,(Ptr{Nothing},Int64,),$(esc(task)),$(esc(afeidx))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_emptyafebarfrowlist(task,numafeidx,afeidxlist)
  quote
     local res = disable_sigint(()->ccall((:MSK_emptyafebarfrowlist,libmosek),Int32,(Ptr{Nothing},Int64,Ptr{Int64},),$(esc(task)),$(esc(numafeidx)),$(esc(afeidxlist))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putafebarfblocktriplet(task,numtrip,afeidx,barvaridx,subk,subl,valkl)
  quote
     local res = disable_sigint(()->ccall((:MSK_putafebarfblocktriplet,libmosek),Int32,(Ptr{Nothing},Int64,Ptr{Int64},Ptr{Int32},Ptr{Int32},Ptr{Int32},Ptr{Float64},),$(esc(task)),$(esc(numtrip)),$(esc(afeidx)),$(esc(barvaridx)),$(esc(subk)),$(esc(subl)),$(esc(valkl))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getafebarfnumblocktriplets(task,numtrip)
  quote
     local res = disable_sigint(()->ccall((:MSK_getafebarfnumblocktriplets,libmosek),Int32,(Ptr{Nothing},Ref{Int64},),$(esc(task)),$(esc(numtrip))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getafebarfblocktriplet(task,maxnumtrip,numtrip,afeidx,barvaridx,subk,subl,valkl)
  quote
     local res = disable_sigint(()->ccall((:MSK_getafebarfblocktriplet,libmosek),Int32,(Ptr{Nothing},Int64,Ref{Int64},Ptr{Int64},Ptr{Int32},Ptr{Int32},Ptr{Int32},Ptr{Float64},),$(esc(task)),$(esc(maxnumtrip)),$(esc(numtrip)),$(esc(afeidx)),$(esc(barvaridx)),$(esc(subk)),$(esc(subl)),$(esc(valkl))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getafebarfnumrowentries(task,afeidx,numentr)
  quote
     local res = disable_sigint(()->ccall((:MSK_getafebarfnumrowentries,libmosek),Int32,(Ptr{Nothing},Int64,Ref{Int32},),$(esc(task)),$(esc(afeidx)),$(esc(numentr))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getafebarfrowinfo(task,afeidx,numentr,numterm)
  quote
     local res = disable_sigint(()->ccall((:MSK_getafebarfrowinfo,libmosek),Int32,(Ptr{Nothing},Int64,Ref{Int32},Ref{Int64},),$(esc(task)),$(esc(afeidx)),$(esc(numentr)),$(esc(numterm))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getafebarfrow(task,afeidx,barvaridx,ptrterm,numterm,termidx,termweight)
  quote
     local res = disable_sigint(()->ccall((:MSK_getafebarfrow,libmosek),Int32,(Ptr{Nothing},Int64,Ptr{Int32},Ptr{Int64},Ptr{Int64},Ptr{Int64},Ptr{Float64},),$(esc(task)),$(esc(afeidx)),$(esc(barvaridx)),$(esc(ptrterm)),$(esc(numterm)),$(esc(termidx)),$(esc(termweight))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putafeg(task,afeidx,g)
  quote
     local res = disable_sigint(()->ccall((:MSK_putafeg,libmosek),Int32,(Ptr{Nothing},Int64,Float64,),$(esc(task)),$(esc(afeidx)),$(esc(g))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putafeglist(task,numafeidx,afeidx,g)
  quote
     local res = disable_sigint(()->ccall((:MSK_putafeglist,libmosek),Int32,(Ptr{Nothing},Int64,Ptr{Int64},Ptr{Float64},),$(esc(task)),$(esc(numafeidx)),$(esc(afeidx)),$(esc(g))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getafeg(task,afeidx,g)
  quote
     local res = disable_sigint(()->ccall((:MSK_getafeg,libmosek),Int32,(Ptr{Nothing},Int64,Ref{Float64},),$(esc(task)),$(esc(afeidx)),$(esc(g))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getafegslice(task,first,last,g)
  quote
     local res = disable_sigint(()->ccall((:MSK_getafegslice,libmosek),Int32,(Ptr{Nothing},Int64,Int64,Ptr{Float64},),$(esc(task)),$(esc(first)),$(esc(last)),$(esc(g))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putafegslice(task,first,last,slice)
  quote
     local res = disable_sigint(()->ccall((:MSK_putafegslice,libmosek),Int32,(Ptr{Nothing},Int64,Int64,Ptr{Float64},),$(esc(task)),$(esc(first)),$(esc(last)),$(esc(slice))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putmaxnumdjc(task,maxnumdjc)
  quote
     local res = disable_sigint(()->ccall((:MSK_putmaxnumdjc,libmosek),Int32,(Ptr{Nothing},Int64,),$(esc(task)),$(esc(maxnumdjc))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getnumdjc(task,num)
  quote
     local res = disable_sigint(()->ccall((:MSK_getnumdjc,libmosek),Int32,(Ptr{Nothing},Ref{Int64},),$(esc(task)),$(esc(num))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getdjcnumdomain(task,djcidx,numdomain)
  quote
     local res = disable_sigint(()->ccall((:MSK_getdjcnumdomain,libmosek),Int32,(Ptr{Nothing},Int64,Ref{Int64},),$(esc(task)),$(esc(djcidx)),$(esc(numdomain))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getdjcnumdomaintot(task,numdomaintot)
  quote
     local res = disable_sigint(()->ccall((:MSK_getdjcnumdomaintot,libmosek),Int32,(Ptr{Nothing},Ref{Int64},),$(esc(task)),$(esc(numdomaintot))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getdjcnumafe(task,djcidx,numafe)
  quote
     local res = disable_sigint(()->ccall((:MSK_getdjcnumafe,libmosek),Int32,(Ptr{Nothing},Int64,Ref{Int64},),$(esc(task)),$(esc(djcidx)),$(esc(numafe))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getdjcnumafetot(task,numafetot)
  quote
     local res = disable_sigint(()->ccall((:MSK_getdjcnumafetot,libmosek),Int32,(Ptr{Nothing},Ref{Int64},),$(esc(task)),$(esc(numafetot))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getdjcnumterm(task,djcidx,numterm)
  quote
     local res = disable_sigint(()->ccall((:MSK_getdjcnumterm,libmosek),Int32,(Ptr{Nothing},Int64,Ref{Int64},),$(esc(task)),$(esc(djcidx)),$(esc(numterm))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getdjcnumtermtot(task,numtermtot)
  quote
     local res = disable_sigint(()->ccall((:MSK_getdjcnumtermtot,libmosek),Int32,(Ptr{Nothing},Ref{Int64},),$(esc(task)),$(esc(numtermtot))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putmaxnumacc(task,maxnumacc)
  quote
     local res = disable_sigint(()->ccall((:MSK_putmaxnumacc,libmosek),Int32,(Ptr{Nothing},Int64,),$(esc(task)),$(esc(maxnumacc))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getnumacc(task,num)
  quote
     local res = disable_sigint(()->ccall((:MSK_getnumacc,libmosek),Int32,(Ptr{Nothing},Ref{Int64},),$(esc(task)),$(esc(num))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_appendacc(task,domidx,numafeidx,afeidxlist,b)
  quote
     local res = disable_sigint(()->ccall((:MSK_appendacc,libmosek),Int32,(Ptr{Nothing},Int64,Int64,Ptr{Int64},Ptr{Float64},),$(esc(task)),$(esc(domidx)),$(esc(numafeidx)),$(esc(afeidxlist)),$(esc(b))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_appendaccs(task,numaccs,domidxs,numafeidx,afeidxlist,b)
  quote
     local res = disable_sigint(()->ccall((:MSK_appendaccs,libmosek),Int32,(Ptr{Nothing},Int64,Ptr{Int64},Int64,Ptr{Int64},Ptr{Float64},),$(esc(task)),$(esc(numaccs)),$(esc(domidxs)),$(esc(numafeidx)),$(esc(afeidxlist)),$(esc(b))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_appendaccseq(task,domidx,numafeidx,afeidxfirst,b)
  quote
     local res = disable_sigint(()->ccall((:MSK_appendaccseq,libmosek),Int32,(Ptr{Nothing},Int64,Int64,Int64,Ptr{Float64},),$(esc(task)),$(esc(domidx)),$(esc(numafeidx)),$(esc(afeidxfirst)),$(esc(b))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_appendaccsseq(task,numaccs,domidxs,numafeidx,afeidxfirst,b)
  quote
     local res = disable_sigint(()->ccall((:MSK_appendaccsseq,libmosek),Int32,(Ptr{Nothing},Int64,Ptr{Int64},Int64,Int64,Ptr{Float64},),$(esc(task)),$(esc(numaccs)),$(esc(domidxs)),$(esc(numafeidx)),$(esc(afeidxfirst)),$(esc(b))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putacc(task,accidx,domidx,numafeidx,afeidxlist,b)
  quote
     local res = disable_sigint(()->ccall((:MSK_putacc,libmosek),Int32,(Ptr{Nothing},Int64,Int64,Int64,Ptr{Int64},Ptr{Float64},),$(esc(task)),$(esc(accidx)),$(esc(domidx)),$(esc(numafeidx)),$(esc(afeidxlist)),$(esc(b))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putacclist(task,numaccs,accidxs,domidxs,numafeidx,afeidxlist,b)
  quote
     local res = disable_sigint(()->ccall((:MSK_putacclist,libmosek),Int32,(Ptr{Nothing},Int64,Ptr{Int64},Ptr{Int64},Int64,Ptr{Int64},Ptr{Float64},),$(esc(task)),$(esc(numaccs)),$(esc(accidxs)),$(esc(domidxs)),$(esc(numafeidx)),$(esc(afeidxlist)),$(esc(b))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putaccb(task,accidx,lengthb,b)
  quote
     local res = disable_sigint(()->ccall((:MSK_putaccb,libmosek),Int32,(Ptr{Nothing},Int64,Int64,Ptr{Float64},),$(esc(task)),$(esc(accidx)),$(esc(lengthb)),$(esc(b))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putaccbj(task,accidx,j,bj)
  quote
     local res = disable_sigint(()->ccall((:MSK_putaccbj,libmosek),Int32,(Ptr{Nothing},Int64,Int64,Float64,),$(esc(task)),$(esc(accidx)),$(esc(j)),$(esc(bj))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getaccdomain(task,accidx,domidx)
  quote
     local res = disable_sigint(()->ccall((:MSK_getaccdomain,libmosek),Int32,(Ptr{Nothing},Int64,Ref{Int64},),$(esc(task)),$(esc(accidx)),$(esc(domidx))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getaccn(task,accidx,n)
  quote
     local res = disable_sigint(()->ccall((:MSK_getaccn,libmosek),Int32,(Ptr{Nothing},Int64,Ref{Int64},),$(esc(task)),$(esc(accidx)),$(esc(n))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getaccntot(task,n)
  quote
     local res = disable_sigint(()->ccall((:MSK_getaccntot,libmosek),Int32,(Ptr{Nothing},Ref{Int64},),$(esc(task)),$(esc(n))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getaccafeidxlist(task,accidx,afeidxlist)
  quote
     local res = disable_sigint(()->ccall((:MSK_getaccafeidxlist,libmosek),Int32,(Ptr{Nothing},Int64,Ptr{Int64},),$(esc(task)),$(esc(accidx)),$(esc(afeidxlist))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getaccb(task,accidx,b)
  quote
     local res = disable_sigint(()->ccall((:MSK_getaccb,libmosek),Int32,(Ptr{Nothing},Int64,Ptr{Float64},),$(esc(task)),$(esc(accidx)),$(esc(b))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getaccs(task,domidxlist,afeidxlist,b)
  quote
     local res = disable_sigint(()->ccall((:MSK_getaccs,libmosek),Int32,(Ptr{Nothing},Ptr{Int64},Ptr{Int64},Ptr{Float64},),$(esc(task)),$(esc(domidxlist)),$(esc(afeidxlist)),$(esc(b))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getaccfnumnz(task,accfnnz)
  quote
     local res = disable_sigint(()->ccall((:MSK_getaccfnumnz,libmosek),Int32,(Ptr{Nothing},Ref{Int64},),$(esc(task)),$(esc(accfnnz))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getaccftrip(task,frow,fcol,fval)
  quote
     local res = disable_sigint(()->ccall((:MSK_getaccftrip,libmosek),Int32,(Ptr{Nothing},Ptr{Int64},Ptr{Int32},Ptr{Float64},),$(esc(task)),$(esc(frow)),$(esc(fcol)),$(esc(fval))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getaccgvector(task,g)
  quote
     local res = disable_sigint(()->ccall((:MSK_getaccgvector,libmosek),Int32,(Ptr{Nothing},Ptr{Float64},),$(esc(task)),$(esc(g))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getaccbarfnumblocktriplets(task,numtrip)
  quote
     local res = disable_sigint(()->ccall((:MSK_getaccbarfnumblocktriplets,libmosek),Int32,(Ptr{Nothing},Ref{Int64},),$(esc(task)),$(esc(numtrip))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getaccbarfblocktriplet(task,maxnumtrip,numtrip,acc_afe,bar_var,blk_row,blk_col,blk_val)
  quote
     local res = disable_sigint(()->ccall((:MSK_getaccbarfblocktriplet,libmosek),Int32,(Ptr{Nothing},Int64,Ref{Int64},Ptr{Int64},Ptr{Int32},Ptr{Int32},Ptr{Int32},Ptr{Float64},),$(esc(task)),$(esc(maxnumtrip)),$(esc(numtrip)),$(esc(acc_afe)),$(esc(bar_var)),$(esc(blk_row)),$(esc(blk_col)),$(esc(blk_val))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_appenddjcs(task,num)
  quote
     local res = disable_sigint(()->ccall((:MSK_appenddjcs,libmosek),Int32,(Ptr{Nothing},Int64,),$(esc(task)),$(esc(num))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putdjc(task,djcidx,numdomidx,domidxlist,numafeidx,afeidxlist,b,numterms,termsizelist)
  quote
     local res = disable_sigint(()->ccall((:MSK_putdjc,libmosek),Int32,(Ptr{Nothing},Int64,Int64,Ptr{Int64},Int64,Ptr{Int64},Ptr{Float64},Int64,Ptr{Int64},),$(esc(task)),$(esc(djcidx)),$(esc(numdomidx)),$(esc(domidxlist)),$(esc(numafeidx)),$(esc(afeidxlist)),$(esc(b)),$(esc(numterms)),$(esc(termsizelist))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putdjcslice(task,idxfirst,idxlast,numdomidx,domidxlist,numafeidx,afeidxlist,b,numterms,termsizelist,termsindjc)
  quote
     local res = disable_sigint(()->ccall((:MSK_putdjcslice,libmosek),Int32,(Ptr{Nothing},Int64,Int64,Int64,Ptr{Int64},Int64,Ptr{Int64},Ptr{Float64},Int64,Ptr{Int64},Ptr{Int64},),$(esc(task)),$(esc(idxfirst)),$(esc(idxlast)),$(esc(numdomidx)),$(esc(domidxlist)),$(esc(numafeidx)),$(esc(afeidxlist)),$(esc(b)),$(esc(numterms)),$(esc(termsizelist)),$(esc(termsindjc))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getdjcdomainidxlist(task,djcidx,domidxlist)
  quote
     local res = disable_sigint(()->ccall((:MSK_getdjcdomainidxlist,libmosek),Int32,(Ptr{Nothing},Int64,Ptr{Int64},),$(esc(task)),$(esc(djcidx)),$(esc(domidxlist))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getdjcafeidxlist(task,djcidx,afeidxlist)
  quote
     local res = disable_sigint(()->ccall((:MSK_getdjcafeidxlist,libmosek),Int32,(Ptr{Nothing},Int64,Ptr{Int64},),$(esc(task)),$(esc(djcidx)),$(esc(afeidxlist))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getdjcb(task,djcidx,b)
  quote
     local res = disable_sigint(()->ccall((:MSK_getdjcb,libmosek),Int32,(Ptr{Nothing},Int64,Ptr{Float64},),$(esc(task)),$(esc(djcidx)),$(esc(b))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getdjctermsizelist(task,djcidx,termsizelist)
  quote
     local res = disable_sigint(()->ccall((:MSK_getdjctermsizelist,libmosek),Int32,(Ptr{Nothing},Int64,Ptr{Int64},),$(esc(task)),$(esc(djcidx)),$(esc(termsizelist))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getdjcs(task,domidxlist,afeidxlist,b,termsizelist,numterms)
  quote
     local res = disable_sigint(()->ccall((:MSK_getdjcs,libmosek),Int32,(Ptr{Nothing},Ptr{Int64},Ptr{Int64},Ptr{Float64},Ptr{Int64},Ptr{Int64},),$(esc(task)),$(esc(domidxlist)),$(esc(afeidxlist)),$(esc(b)),$(esc(termsizelist)),$(esc(numterms))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putconbound(task,i,bkc,blc,buc)
  quote
     local res = disable_sigint(()->ccall((:MSK_putconbound,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Float64,Float64,),$(esc(task)),$(esc(i)),$(esc(bkc)),$(esc(blc)),$(esc(buc))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putconboundlist(task,num,sub,bkc,blc,buc)
  quote
     local res = disable_sigint(()->ccall((:MSK_putconboundlist,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{Int32},Ptr{Int32},Ptr{Float64},Ptr{Float64},),$(esc(task)),$(esc(num)),$(esc(sub)),$(esc(bkc)),$(esc(blc)),$(esc(buc))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putconboundlistconst(task,num,sub,bkc,blc,buc)
  quote
     local res = disable_sigint(()->ccall((:MSK_putconboundlistconst,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{Int32},Int32,Float64,Float64,),$(esc(task)),$(esc(num)),$(esc(sub)),$(esc(bkc)),$(esc(blc)),$(esc(buc))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putconboundslice(task,first,last,bkc,blc,buc)
  quote
     local res = disable_sigint(()->ccall((:MSK_putconboundslice,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Ptr{Int32},Ptr{Float64},Ptr{Float64},),$(esc(task)),$(esc(first)),$(esc(last)),$(esc(bkc)),$(esc(blc)),$(esc(buc))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putconboundsliceconst(task,first,last,bkc,blc,buc)
  quote
     local res = disable_sigint(()->ccall((:MSK_putconboundsliceconst,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Int32,Float64,Float64,),$(esc(task)),$(esc(first)),$(esc(last)),$(esc(bkc)),$(esc(blc)),$(esc(buc))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putvarbound(task,j,bkx,blx,bux)
  quote
     local res = disable_sigint(()->ccall((:MSK_putvarbound,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Float64,Float64,),$(esc(task)),$(esc(j)),$(esc(bkx)),$(esc(blx)),$(esc(bux))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putvarboundlist(task,num,sub,bkx,blx,bux)
  quote
     local res = disable_sigint(()->ccall((:MSK_putvarboundlist,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{Int32},Ptr{Int32},Ptr{Float64},Ptr{Float64},),$(esc(task)),$(esc(num)),$(esc(sub)),$(esc(bkx)),$(esc(blx)),$(esc(bux))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putvarboundlistconst(task,num,sub,bkx,blx,bux)
  quote
     local res = disable_sigint(()->ccall((:MSK_putvarboundlistconst,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{Int32},Int32,Float64,Float64,),$(esc(task)),$(esc(num)),$(esc(sub)),$(esc(bkx)),$(esc(blx)),$(esc(bux))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putvarboundslice(task,first,last,bkx,blx,bux)
  quote
     local res = disable_sigint(()->ccall((:MSK_putvarboundslice,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Ptr{Int32},Ptr{Float64},Ptr{Float64},),$(esc(task)),$(esc(first)),$(esc(last)),$(esc(bkx)),$(esc(blx)),$(esc(bux))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putvarboundsliceconst(task,first,last,bkx,blx,bux)
  quote
     local res = disable_sigint(()->ccall((:MSK_putvarboundsliceconst,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Int32,Float64,Float64,),$(esc(task)),$(esc(first)),$(esc(last)),$(esc(bkx)),$(esc(blx)),$(esc(bux))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putcallbackfunc(task,func,handle)
  quote
     local res = disable_sigint(()->ccall((:MSK_putcallbackfunc,libmosek),Int32,(Ptr{Nothing},Ptr{Cvoid},Any,),$(esc(task)),$(esc(func)),$(esc(handle))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putcfix(task,cfix)
  quote
     local res = disable_sigint(()->ccall((:MSK_putcfix,libmosek),Int32,(Ptr{Nothing},Float64,),$(esc(task)),$(esc(cfix))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putcj(task,j,cj)
  quote
     local res = disable_sigint(()->ccall((:MSK_putcj,libmosek),Int32,(Ptr{Nothing},Int32,Float64,),$(esc(task)),$(esc(j)),$(esc(cj))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putobjsense(task,sense)
  quote
     local res = disable_sigint(()->ccall((:MSK_putobjsense,libmosek),Int32,(Ptr{Nothing},Int32,),$(esc(task)),$(esc(sense))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getobjsense(task,sense)
  quote
     local res = disable_sigint(()->ccall((:MSK_getobjsense,libmosek),Int32,(Ptr{Nothing},Ref{Int32},),$(esc(task)),$(esc(sense))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putclist(task,num,subj,val)
  quote
     local res = disable_sigint(()->ccall((:MSK_putclist,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{Int32},Ptr{Float64},),$(esc(task)),$(esc(num)),$(esc(subj)),$(esc(val))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putcslice(task,first,last,slice)
  quote
     local res = disable_sigint(()->ccall((:MSK_putcslice,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Ptr{Float64},),$(esc(task)),$(esc(first)),$(esc(last)),$(esc(slice))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putbarcj(task,j,num,sub,weights)
  quote
     local res = disable_sigint(()->ccall((:MSK_putbarcj,libmosek),Int32,(Ptr{Nothing},Int32,Int64,Ptr{Int64},Ptr{Float64},),$(esc(task)),$(esc(j)),$(esc(num)),$(esc(sub)),$(esc(weights))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putcone(task,k,ct,conepar,nummem,submem)
  quote
     local res = disable_sigint(()->ccall((:MSK_putcone,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Float64,Int32,Ptr{Int32},),$(esc(task)),$(esc(k)),$(esc(ct)),$(esc(conepar)),$(esc(nummem)),$(esc(submem))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putmaxnumdomain(task,maxnumdomain)
  quote
     local res = disable_sigint(()->ccall((:MSK_putmaxnumdomain,libmosek),Int32,(Ptr{Nothing},Int64,),$(esc(task)),$(esc(maxnumdomain))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getnumdomain(task,numdomain)
  quote
     local res = disable_sigint(()->ccall((:MSK_getnumdomain,libmosek),Int32,(Ptr{Nothing},Ref{Int64},),$(esc(task)),$(esc(numdomain))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_appendrplusdomain(task,n,domidx)
  quote
     local res = disable_sigint(()->ccall((:MSK_appendrplusdomain,libmosek),Int32,(Ptr{Nothing},Int64,Ref{Int64},),$(esc(task)),$(esc(n)),$(esc(domidx))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_appendrminusdomain(task,n,domidx)
  quote
     local res = disable_sigint(()->ccall((:MSK_appendrminusdomain,libmosek),Int32,(Ptr{Nothing},Int64,Ref{Int64},),$(esc(task)),$(esc(n)),$(esc(domidx))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_appendrdomain(task,n,domidx)
  quote
     local res = disable_sigint(()->ccall((:MSK_appendrdomain,libmosek),Int32,(Ptr{Nothing},Int64,Ref{Int64},),$(esc(task)),$(esc(n)),$(esc(domidx))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_appendrzerodomain(task,n,domidx)
  quote
     local res = disable_sigint(()->ccall((:MSK_appendrzerodomain,libmosek),Int32,(Ptr{Nothing},Int64,Ref{Int64},),$(esc(task)),$(esc(n)),$(esc(domidx))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_appendquadraticconedomain(task,n,domidx)
  quote
     local res = disable_sigint(()->ccall((:MSK_appendquadraticconedomain,libmosek),Int32,(Ptr{Nothing},Int64,Ref{Int64},),$(esc(task)),$(esc(n)),$(esc(domidx))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_appendrquadraticconedomain(task,n,domidx)
  quote
     local res = disable_sigint(()->ccall((:MSK_appendrquadraticconedomain,libmosek),Int32,(Ptr{Nothing},Int64,Ref{Int64},),$(esc(task)),$(esc(n)),$(esc(domidx))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_appendprimalexpconedomain(task,domidx)
  quote
     local res = disable_sigint(()->ccall((:MSK_appendprimalexpconedomain,libmosek),Int32,(Ptr{Nothing},Ref{Int64},),$(esc(task)),$(esc(domidx))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_appenddualexpconedomain(task,domidx)
  quote
     local res = disable_sigint(()->ccall((:MSK_appenddualexpconedomain,libmosek),Int32,(Ptr{Nothing},Ref{Int64},),$(esc(task)),$(esc(domidx))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_appendprimalgeomeanconedomain(task,n,domidx)
  quote
     local res = disable_sigint(()->ccall((:MSK_appendprimalgeomeanconedomain,libmosek),Int32,(Ptr{Nothing},Int64,Ref{Int64},),$(esc(task)),$(esc(n)),$(esc(domidx))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_appenddualgeomeanconedomain(task,n,domidx)
  quote
     local res = disable_sigint(()->ccall((:MSK_appenddualgeomeanconedomain,libmosek),Int32,(Ptr{Nothing},Int64,Ref{Int64},),$(esc(task)),$(esc(n)),$(esc(domidx))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_appendprimalpowerconedomain(task,n,nleft,alpha,domidx)
  quote
     local res = disable_sigint(()->ccall((:MSK_appendprimalpowerconedomain,libmosek),Int32,(Ptr{Nothing},Int64,Int64,Ptr{Float64},Ref{Int64},),$(esc(task)),$(esc(n)),$(esc(nleft)),$(esc(alpha)),$(esc(domidx))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_appenddualpowerconedomain(task,n,nleft,alpha,domidx)
  quote
     local res = disable_sigint(()->ccall((:MSK_appenddualpowerconedomain,libmosek),Int32,(Ptr{Nothing},Int64,Int64,Ptr{Float64},Ref{Int64},),$(esc(task)),$(esc(n)),$(esc(nleft)),$(esc(alpha)),$(esc(domidx))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_appendsvecpsdconedomain(task,n,domidx)
  quote
     local res = disable_sigint(()->ccall((:MSK_appendsvecpsdconedomain,libmosek),Int32,(Ptr{Nothing},Int64,Ref{Int64},),$(esc(task)),$(esc(n)),$(esc(domidx))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getdomaintype(task,domidx,domtype)
  quote
     local res = disable_sigint(()->ccall((:MSK_getdomaintype,libmosek),Int32,(Ptr{Nothing},Int64,Ref{Int32},),$(esc(task)),$(esc(domidx)),$(esc(domtype))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getdomainn(task,domidx,n)
  quote
     local res = disable_sigint(()->ccall((:MSK_getdomainn,libmosek),Int32,(Ptr{Nothing},Int64,Ref{Int64},),$(esc(task)),$(esc(domidx)),$(esc(n))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getpowerdomaininfo(task,domidx,n,nleft)
  quote
     local res = disable_sigint(()->ccall((:MSK_getpowerdomaininfo,libmosek),Int32,(Ptr{Nothing},Int64,Ref{Int64},Ref{Int64},),$(esc(task)),$(esc(domidx)),$(esc(n)),$(esc(nleft))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getpowerdomainalpha(task,domidx,alpha)
  quote
     local res = disable_sigint(()->ccall((:MSK_getpowerdomainalpha,libmosek),Int32,(Ptr{Nothing},Int64,Ptr{Float64},),$(esc(task)),$(esc(domidx)),$(esc(alpha))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_appendsparsesymmat(task,dim,nz,subi,subj,valij,idx)
  quote
     local res = disable_sigint(()->ccall((:MSK_appendsparsesymmat,libmosek),Int32,(Ptr{Nothing},Int32,Int64,Ptr{Int32},Ptr{Int32},Ptr{Float64},Ref{Int64},),$(esc(task)),$(esc(dim)),$(esc(nz)),$(esc(subi)),$(esc(subj)),$(esc(valij)),$(esc(idx))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_appendsparsesymmatlist(task,num,dims,nz,subi,subj,valij,idx)
  quote
     local res = disable_sigint(()->ccall((:MSK_appendsparsesymmatlist,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{Int32},Ptr{Int64},Ptr{Int32},Ptr{Int32},Ptr{Float64},Ptr{Int64},),$(esc(task)),$(esc(num)),$(esc(dims)),$(esc(nz)),$(esc(subi)),$(esc(subj)),$(esc(valij)),$(esc(idx))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getsymmatinfo(task,idx,dim,nz,mattype)
  quote
     local res = disable_sigint(()->ccall((:MSK_getsymmatinfo,libmosek),Int32,(Ptr{Nothing},Int64,Ref{Int32},Ref{Int64},Ref{Int32},),$(esc(task)),$(esc(idx)),$(esc(dim)),$(esc(nz)),$(esc(mattype))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getnumsymmat(task,num)
  quote
     local res = disable_sigint(()->ccall((:MSK_getnumsymmat,libmosek),Int32,(Ptr{Nothing},Ref{Int64},),$(esc(task)),$(esc(num))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getsparsesymmat(task,idx,maxlen,subi,subj,valij)
  quote
     local res = disable_sigint(()->ccall((:MSK_getsparsesymmat,libmosek),Int32,(Ptr{Nothing},Int64,Int64,Ptr{Int32},Ptr{Int32},Ptr{Float64},),$(esc(task)),$(esc(idx)),$(esc(maxlen)),$(esc(subi)),$(esc(subj)),$(esc(valij))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putdouparam(task,param,parvalue)
  quote
     local res = disable_sigint(()->ccall((:MSK_putdouparam,libmosek),Int32,(Ptr{Nothing},Int32,Float64,),$(esc(task)),$(esc(param)),$(esc(parvalue))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putintparam(task,param,parvalue)
  quote
     local res = disable_sigint(()->ccall((:MSK_putintparam,libmosek),Int32,(Ptr{Nothing},Int32,Int32,),$(esc(task)),$(esc(param)),$(esc(parvalue))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putmaxnumcon(task,maxnumcon)
  quote
     local res = disable_sigint(()->ccall((:MSK_putmaxnumcon,libmosek),Int32,(Ptr{Nothing},Int32,),$(esc(task)),$(esc(maxnumcon))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putmaxnumcone(task,maxnumcone)
  quote
     local res = disable_sigint(()->ccall((:MSK_putmaxnumcone,libmosek),Int32,(Ptr{Nothing},Int32,),$(esc(task)),$(esc(maxnumcone))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getmaxnumcone(task,maxnumcone)
  quote
     local res = disable_sigint(()->ccall((:MSK_getmaxnumcone,libmosek),Int32,(Ptr{Nothing},Ref{Int32},),$(esc(task)),$(esc(maxnumcone))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putmaxnumvar(task,maxnumvar)
  quote
     local res = disable_sigint(()->ccall((:MSK_putmaxnumvar,libmosek),Int32,(Ptr{Nothing},Int32,),$(esc(task)),$(esc(maxnumvar))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putmaxnumbarvar(task,maxnumbarvar)
  quote
     local res = disable_sigint(()->ccall((:MSK_putmaxnumbarvar,libmosek),Int32,(Ptr{Nothing},Int32,),$(esc(task)),$(esc(maxnumbarvar))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putmaxnumanz(task,maxnumanz)
  quote
     local res = disable_sigint(()->ccall((:MSK_putmaxnumanz,libmosek),Int32,(Ptr{Nothing},Int64,),$(esc(task)),$(esc(maxnumanz))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putmaxnumqnz(task,maxnumqnz)
  quote
     local res = disable_sigint(()->ccall((:MSK_putmaxnumqnz,libmosek),Int32,(Ptr{Nothing},Int64,),$(esc(task)),$(esc(maxnumqnz))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getmaxnumqnz(task,maxnumqnz)
  quote
     local res = disable_sigint(()->ccall((:MSK_getmaxnumqnz,libmosek),Int32,(Ptr{Nothing},Ref{Int32},),$(esc(task)),$(esc(maxnumqnz))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getmaxnumqnz64(task,maxnumqnz)
  quote
     local res = disable_sigint(()->ccall((:MSK_getmaxnumqnz64,libmosek),Int32,(Ptr{Nothing},Ref{Int64},),$(esc(task)),$(esc(maxnumqnz))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putnadouparam(task,paramname,parvalue)
  quote
     local res = disable_sigint(()->ccall((:MSK_putnadouparam,libmosek),Int32,(Ptr{Nothing},Ptr{UInt8},Float64,),$(esc(task)),$(esc(paramname)),$(esc(parvalue))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putnaintparam(task,paramname,parvalue)
  quote
     local res = disable_sigint(()->ccall((:MSK_putnaintparam,libmosek),Int32,(Ptr{Nothing},Ptr{UInt8},Int32,),$(esc(task)),$(esc(paramname)),$(esc(parvalue))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putnastrparam(task,paramname,parvalue)
  quote
     local res = disable_sigint(()->ccall((:MSK_putnastrparam,libmosek),Int32,(Ptr{Nothing},Ptr{UInt8},Ptr{UInt8},),$(esc(task)),$(esc(paramname)),$(esc(parvalue))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putobjname(task,objname)
  quote
     local res = disable_sigint(()->ccall((:MSK_putobjname,libmosek),Int32,(Ptr{Nothing},Ptr{UInt8},),$(esc(task)),$(esc(objname))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putparam(task,parname,parvalue)
  quote
     local res = disable_sigint(()->ccall((:MSK_putparam,libmosek),Int32,(Ptr{Nothing},Ptr{UInt8},Ptr{UInt8},),$(esc(task)),$(esc(parname)),$(esc(parvalue))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putqcon(task,numqcnz,qcsubk,qcsubi,qcsubj,qcval)
  quote
     local res = disable_sigint(()->ccall((:MSK_putqcon,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{Int32},Ptr{Int32},Ptr{Int32},Ptr{Float64},),$(esc(task)),$(esc(numqcnz)),$(esc(qcsubk)),$(esc(qcsubi)),$(esc(qcsubj)),$(esc(qcval))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putqconk(task,k,numqcnz,qcsubi,qcsubj,qcval)
  quote
     local res = disable_sigint(()->ccall((:MSK_putqconk,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Ptr{Int32},Ptr{Int32},Ptr{Float64},),$(esc(task)),$(esc(k)),$(esc(numqcnz)),$(esc(qcsubi)),$(esc(qcsubj)),$(esc(qcval))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putqobj(task,numqonz,qosubi,qosubj,qoval)
  quote
     local res = disable_sigint(()->ccall((:MSK_putqobj,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{Int32},Ptr{Int32},Ptr{Float64},),$(esc(task)),$(esc(numqonz)),$(esc(qosubi)),$(esc(qosubj)),$(esc(qoval))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putqobjij(task,i,j,qoij)
  quote
     local res = disable_sigint(()->ccall((:MSK_putqobjij,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Float64,),$(esc(task)),$(esc(i)),$(esc(j)),$(esc(qoij))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putsolution(task,whichsol,skc,skx,skn,xc,xx,y,slc,suc,slx,sux,snx)
  quote
     local res = disable_sigint(()->ccall((:MSK_putsolution,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{Int32},Ptr{Int32},Ptr{Int32},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(skc)),$(esc(skx)),$(esc(skn)),$(esc(xc)),$(esc(xx)),$(esc(y)),$(esc(slc)),$(esc(suc)),$(esc(slx)),$(esc(sux)),$(esc(snx))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putsolutionnew(task,whichsol,skc,skx,skn,xc,xx,y,slc,suc,slx,sux,snx,doty)
  quote
     local res = disable_sigint(()->ccall((:MSK_putsolutionnew,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{Int32},Ptr{Int32},Ptr{Int32},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},),$(esc(task)),$(esc(whichsol)),$(esc(skc)),$(esc(skx)),$(esc(skn)),$(esc(xc)),$(esc(xx)),$(esc(y)),$(esc(slc)),$(esc(suc)),$(esc(slx)),$(esc(sux)),$(esc(snx)),$(esc(doty))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putconsolutioni(task,i,whichsol,sk,x,sl,su)
  quote
     local res = disable_sigint(()->ccall((:MSK_putconsolutioni,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Int32,Float64,Float64,Float64,),$(esc(task)),$(esc(i)),$(esc(whichsol)),$(esc(sk)),$(esc(x)),$(esc(sl)),$(esc(su))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putvarsolutionj(task,j,whichsol,sk,x,sl,su,sn)
  quote
     local res = disable_sigint(()->ccall((:MSK_putvarsolutionj,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Int32,Float64,Float64,Float64,Float64,),$(esc(task)),$(esc(j)),$(esc(whichsol)),$(esc(sk)),$(esc(x)),$(esc(sl)),$(esc(su)),$(esc(sn))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putsolutionyi(task,i,whichsol,y)
  quote
     local res = disable_sigint(()->ccall((:MSK_putsolutionyi,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Float64,),$(esc(task)),$(esc(i)),$(esc(whichsol)),$(esc(y))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putstrparam(task,param,parvalue)
  quote
     local res = disable_sigint(()->ccall((:MSK_putstrparam,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{UInt8},),$(esc(task)),$(esc(param)),$(esc(parvalue))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_puttaskname(task,taskname)
  quote
     local res = disable_sigint(()->ccall((:MSK_puttaskname,libmosek),Int32,(Ptr{Nothing},Ptr{UInt8},),$(esc(task)),$(esc(taskname))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putvartype(task,j,vartype)
  quote
     local res = disable_sigint(()->ccall((:MSK_putvartype,libmosek),Int32,(Ptr{Nothing},Int32,Int32,),$(esc(task)),$(esc(j)),$(esc(vartype))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putvartypelist(task,num,subj,vartype)
  quote
     local res = disable_sigint(()->ccall((:MSK_putvartypelist,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{Int32},Ptr{Int32},),$(esc(task)),$(esc(num)),$(esc(subj)),$(esc(vartype))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_readdata(task,filename)
  quote
     local res = disable_sigint(()->ccall((:MSK_readdata,libmosek),Int32,(Ptr{Nothing},Ptr{UInt8},),$(esc(task)),$(esc(filename))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_readdatacb(task,hread,h,format,compress,path)
  quote
     local res = disable_sigint(()->ccall((:MSK_readdatacb,libmosek),Int32,(Ptr{Nothing},Ptr{Cvoid},Any,Int32,Int32,Ptr{UInt8},),$(esc(task)),$(esc(hread)),$(esc(h)),$(esc(format)),$(esc(compress)),$(esc(path))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_writedatahandle(task,func,handle,format,compress)
  quote
     local res = disable_sigint(()->ccall((:MSK_writedatahandle,libmosek),Int32,(Ptr{Nothing},Ptr{Cvoid},Any,Int32,Int32,),$(esc(task)),$(esc(func)),$(esc(handle)),$(esc(format)),$(esc(compress))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_readdataformat(task,filename,format,compress)
  quote
     local res = disable_sigint(()->ccall((:MSK_readdataformat,libmosek),Int32,(Ptr{Nothing},Ptr{UInt8},Int32,Int32,),$(esc(task)),$(esc(filename)),$(esc(format)),$(esc(compress))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_readdataautoformat(task,filename)
  quote
     local res = disable_sigint(()->ccall((:MSK_readdataautoformat,libmosek),Int32,(Ptr{Nothing},Ptr{UInt8},),$(esc(task)),$(esc(filename))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_readparamfile(task,filename)
  quote
     local res = disable_sigint(()->ccall((:MSK_readparamfile,libmosek),Int32,(Ptr{Nothing},Ptr{UInt8},),$(esc(task)),$(esc(filename))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_readsolution(task,whichsol,filename)
  quote
     local res = disable_sigint(()->ccall((:MSK_readsolution,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{UInt8},),$(esc(task)),$(esc(whichsol)),$(esc(filename))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_readjsonsol(task,filename)
  quote
     local res = disable_sigint(()->ccall((:MSK_readjsonsol,libmosek),Int32,(Ptr{Nothing},Ptr{UInt8},),$(esc(task)),$(esc(filename))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_readsummary(task,whichstream)
  quote
     local res = disable_sigint(()->ccall((:MSK_readsummary,libmosek),Int32,(Ptr{Nothing},Int32,),$(esc(task)),$(esc(whichstream))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_resizetask(task,maxnumcon,maxnumvar,maxnumcone,maxnumanz,maxnumqnz)
  quote
     local res = disable_sigint(()->ccall((:MSK_resizetask,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Int32,Int64,Int64,),$(esc(task)),$(esc(maxnumcon)),$(esc(maxnumvar)),$(esc(maxnumcone)),$(esc(maxnumanz)),$(esc(maxnumqnz))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_checkmemtask(task,file,line)
  quote
     local res = disable_sigint(()->ccall((:MSK_checkmemtask,libmosek),Int32,(Ptr{Nothing},Ptr{UInt8},Int32,),$(esc(task)),$(esc(file)),$(esc(line))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getmemusagetask(task,meminuse,maxmemuse)
  quote
     local res = disable_sigint(()->ccall((:MSK_getmemusagetask,libmosek),Int32,(Ptr{Nothing},Ref{Int64},Ref{Int64},),$(esc(task)),$(esc(meminuse)),$(esc(maxmemuse))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_setdefaults(task)
  quote
     local res = disable_sigint(()->ccall((:MSK_setdefaults,libmosek),Int32,(Ptr{Nothing},),$(esc(task))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_sktostr(task,sk,str)
  quote
     local res = disable_sigint(()->ccall((:MSK_sktostr,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{UInt8},),$(esc(task)),$(esc(sk)),$(esc(str))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_solstatostr(task,solutionsta,str)
  quote
     local res = disable_sigint(()->ccall((:MSK_solstatostr,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{UInt8},),$(esc(task)),$(esc(solutionsta)),$(esc(str))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_solutiondef(task,whichsol,isdef)
  quote
     local res = disable_sigint(()->ccall((:MSK_solutiondef,libmosek),Int32,(Ptr{Nothing},Int32,Ref{Int32},),$(esc(task)),$(esc(whichsol)),$(esc(isdef))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_deletesolution(task,whichsol)
  quote
     local res = disable_sigint(()->ccall((:MSK_deletesolution,libmosek),Int32,(Ptr{Nothing},Int32,),$(esc(task)),$(esc(whichsol))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_onesolutionsummary(task,whichstream,whichsol)
  quote
     local res = disable_sigint(()->ccall((:MSK_onesolutionsummary,libmosek),Int32,(Ptr{Nothing},Int32,Int32,),$(esc(task)),$(esc(whichstream)),$(esc(whichsol))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_solutionsummary(task,whichstream)
  quote
     local res = disable_sigint(()->ccall((:MSK_solutionsummary,libmosek),Int32,(Ptr{Nothing},Int32,),$(esc(task)),$(esc(whichstream))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_updatesolutioninfo(task,whichsol)
  quote
     local res = disable_sigint(()->ccall((:MSK_updatesolutioninfo,libmosek),Int32,(Ptr{Nothing},Int32,),$(esc(task)),$(esc(whichsol))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_optimizersummary(task,whichstream)
  quote
     local res = disable_sigint(()->ccall((:MSK_optimizersummary,libmosek),Int32,(Ptr{Nothing},Int32,),$(esc(task)),$(esc(whichstream))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_strtoconetype(task,str,conetype)
  quote
     local res = disable_sigint(()->ccall((:MSK_strtoconetype,libmosek),Int32,(Ptr{Nothing},Ptr{UInt8},Ref{Int32},),$(esc(task)),$(esc(str)),$(esc(conetype))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_strtosk(task,str,sk)
  quote
     local res = disable_sigint(()->ccall((:MSK_strtosk,libmosek),Int32,(Ptr{Nothing},Ptr{UInt8},Ref{Int32},),$(esc(task)),$(esc(str)),$(esc(sk))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_whichparam(task,parname,partype,param)
  quote
     local res = disable_sigint(()->ccall((:MSK_whichparam,libmosek),Int32,(Ptr{Nothing},Ptr{UInt8},Ref{Int32},Ref{Int32},),$(esc(task)),$(esc(parname)),$(esc(partype)),$(esc(param))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_writedata(task,filename)
  quote
     local res = disable_sigint(()->ccall((:MSK_writedata,libmosek),Int32,(Ptr{Nothing},Ptr{UInt8},),$(esc(task)),$(esc(filename))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_writetask(task,filename)
  quote
     local res = disable_sigint(()->ccall((:MSK_writetask,libmosek),Int32,(Ptr{Nothing},Ptr{UInt8},),$(esc(task)),$(esc(filename))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_writebsolution(task,filename,compress)
  quote
     local res = disable_sigint(()->ccall((:MSK_writebsolution,libmosek),Int32,(Ptr{Nothing},Ptr{UInt8},Int32,),$(esc(task)),$(esc(filename)),$(esc(compress))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_writebsolutionhandle(task,func,handle,compress)
  quote
     local res = disable_sigint(()->ccall((:MSK_writebsolutionhandle,libmosek),Int32,(Ptr{Nothing},Ptr{Cvoid},Any,Int32,),$(esc(task)),$(esc(func)),$(esc(handle)),$(esc(compress))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_readbsolution(task,filename,compress)
  quote
     local res = disable_sigint(()->ccall((:MSK_readbsolution,libmosek),Int32,(Ptr{Nothing},Ptr{UInt8},Int32,),$(esc(task)),$(esc(filename)),$(esc(compress))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_writesolutionfile(task,filename)
  quote
     local res = disable_sigint(()->ccall((:MSK_writesolutionfile,libmosek),Int32,(Ptr{Nothing},Ptr{UInt8},),$(esc(task)),$(esc(filename))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_readsolutionfile(task,filename)
  quote
     local res = disable_sigint(()->ccall((:MSK_readsolutionfile,libmosek),Int32,(Ptr{Nothing},Ptr{UInt8},),$(esc(task)),$(esc(filename))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_readtask(task,filename)
  quote
     local res = disable_sigint(()->ccall((:MSK_readtask,libmosek),Int32,(Ptr{Nothing},Ptr{UInt8},),$(esc(task)),$(esc(filename))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_readopfstring(task,data)
  quote
     local res = disable_sigint(()->ccall((:MSK_readopfstring,libmosek),Int32,(Ptr{Nothing},Ptr{UInt8},),$(esc(task)),$(esc(data))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_readlpstring(task,data)
  quote
     local res = disable_sigint(()->ccall((:MSK_readlpstring,libmosek),Int32,(Ptr{Nothing},Ptr{UInt8},),$(esc(task)),$(esc(data))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_readjsonstring(task,data)
  quote
     local res = disable_sigint(()->ccall((:MSK_readjsonstring,libmosek),Int32,(Ptr{Nothing},Ptr{UInt8},),$(esc(task)),$(esc(data))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_readptfstring(task,data)
  quote
     local res = disable_sigint(()->ccall((:MSK_readptfstring,libmosek),Int32,(Ptr{Nothing},Ptr{UInt8},),$(esc(task)),$(esc(data))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_writeparamfile(task,filename)
  quote
     local res = disable_sigint(()->ccall((:MSK_writeparamfile,libmosek),Int32,(Ptr{Nothing},Ptr{UInt8},),$(esc(task)),$(esc(filename))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getinfeasiblesubproblem(task,whichsol,inftask)
  quote
     local res = disable_sigint(()->ccall((:MSK_getinfeasiblesubproblem,libmosek),Int32,(Ptr{Nothing},Int32,Ref{Ptr{Nothing}},),$(esc(task)),$(esc(whichsol)),$(esc(inftask))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_writesolution(task,whichsol,filename)
  quote
     local res = disable_sigint(()->ccall((:MSK_writesolution,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{UInt8},),$(esc(task)),$(esc(whichsol)),$(esc(filename))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_writejsonsol(task,filename)
  quote
     local res = disable_sigint(()->ccall((:MSK_writejsonsol,libmosek),Int32,(Ptr{Nothing},Ptr{UInt8},),$(esc(task)),$(esc(filename))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_primalsensitivity(task,numi,subi,marki,numj,subj,markj,leftpricei,rightpricei,leftrangei,rightrangei,leftpricej,rightpricej,leftrangej,rightrangej)
  quote
     local res = disable_sigint(()->ccall((:MSK_primalsensitivity,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{Int32},Ptr{Int32},Int32,Ptr{Int32},Ptr{Int32},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},),$(esc(task)),$(esc(numi)),$(esc(subi)),$(esc(marki)),$(esc(numj)),$(esc(subj)),$(esc(markj)),$(esc(leftpricei)),$(esc(rightpricei)),$(esc(leftrangei)),$(esc(rightrangei)),$(esc(leftpricej)),$(esc(rightpricej)),$(esc(leftrangej)),$(esc(rightrangej))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_sensitivityreport(task,whichstream)
  quote
     local res = disable_sigint(()->ccall((:MSK_sensitivityreport,libmosek),Int32,(Ptr{Nothing},Int32,),$(esc(task)),$(esc(whichstream))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_dualsensitivity(task,numj,subj,leftpricej,rightpricej,leftrangej,rightrangej)
  quote
     local res = disable_sigint(()->ccall((:MSK_dualsensitivity,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{Int32},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},),$(esc(task)),$(esc(numj)),$(esc(subj)),$(esc(leftpricej)),$(esc(rightpricej)),$(esc(leftrangej)),$(esc(rightrangej))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_getlasterror64(task,lastrescode,sizelastmsg,lastmsglen,lastmsg)
  quote
     local res = disable_sigint(()->ccall((:MSK_getlasterror64,libmosek),Int32,(Ptr{Nothing},Ref{Int32},Int64,Ref{Int64},Ptr{UInt8},),$(esc(task)),$(esc(lastrescode)),$(esc(sizelastmsg)),$(esc(lastmsglen)),$(esc(lastmsg))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_optimizermt(task,address,accesstoken,trmcode)
  quote
     local res = disable_sigint(()->ccall((:MSK_optimizermt,libmosek),Int32,(Ptr{Nothing},Ptr{UInt8},Ptr{UInt8},Ref{Int32},),$(esc(task)),$(esc(address)),$(esc(accesstoken)),$(esc(trmcode))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_asyncoptimize(task,address,accesstoken,token)
  quote
     local res = disable_sigint(()->ccall((:MSK_asyncoptimize,libmosek),Int32,(Ptr{Nothing},Ptr{UInt8},Ptr{UInt8},Ptr{UInt8},),$(esc(task)),$(esc(address)),$(esc(accesstoken)),$(esc(token))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_asyncstop(task,address,accesstoken,token)
  quote
     local res = disable_sigint(()->ccall((:MSK_asyncstop,libmosek),Int32,(Ptr{Nothing},Ptr{UInt8},Ptr{UInt8},Ptr{UInt8},),$(esc(task)),$(esc(address)),$(esc(accesstoken)),$(esc(token))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_asyncpoll(task,address,accesstoken,token,respavailable,resp,trm)
  quote
     local res = disable_sigint(()->ccall((:MSK_asyncpoll,libmosek),Int32,(Ptr{Nothing},Ptr{UInt8},Ptr{UInt8},Ptr{UInt8},Ref{Int32},Ref{Int32},Ref{Int32},),$(esc(task)),$(esc(address)),$(esc(accesstoken)),$(esc(token)),$(esc(respavailable)),$(esc(resp)),$(esc(trm))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_asyncgetresult(task,address,accesstoken,token,respavailable,resp,trm)
  quote
     local res = disable_sigint(()->ccall((:MSK_asyncgetresult,libmosek),Int32,(Ptr{Nothing},Ptr{UInt8},Ptr{UInt8},Ptr{UInt8},Ref{Int32},Ref{Int32},Ref{Int32},),$(esc(task)),$(esc(address)),$(esc(accesstoken)),$(esc(token)),$(esc(respavailable)),$(esc(resp)),$(esc(trm))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_putoptserverhost(task,host)
  quote
     local res = disable_sigint(()->ccall((:MSK_putoptserverhost,libmosek),Int32,(Ptr{Nothing},Ptr{UInt8},),$(esc(task)),$(esc(host))))
     if res != 0
       throw(MosekError(res,getlasterrormsg($(esc(task)))))
     end
     nothing
  end
end
macro MSK_optimizebatch(env,israce,maxtime,numthreads,numtask,task,trmcode,rcode)
  quote
     local res = disable_sigint(()->ccall((:MSK_optimizebatch,libmosek),Int32,(Ptr{Nothing},Int32,Float64,Int32,Int64,Ptr{Ptr{Nothing}},Ptr{Int32},Ptr{Int32},),$(esc(env)),$(esc(israce)),$(esc(maxtime)),$(esc(numthreads)),$(esc(numtask)),$(esc(task)),$(esc(trmcode)),$(esc(rcode))))
     if res != 0
       throw(MosekError(res,""))
     end
     nothing
  end
end
macro MSK_callbackcodetostr(code,callbackcodestr)
  quote
     local res = disable_sigint(()->ccall((:MSK_callbackcodetostr,libmosek),Int32,(Int32,Ptr{UInt8},),$(esc(code)),$(esc(callbackcodestr))))
     if res != 0
       throw(MosekError(res,""))
     end
     nothing
  end
end
macro MSK_checkoutlicense(env,feature)
  quote
     local res = disable_sigint(()->ccall((:MSK_checkoutlicense,libmosek),Int32,(Ptr{Nothing},Int32,),$(esc(env)),$(esc(feature))))
     if res != 0
       throw(MosekError(res,""))
     end
     nothing
  end
end
macro MSK_checkinlicense(env,feature)
  quote
     local res = disable_sigint(()->ccall((:MSK_checkinlicense,libmosek),Int32,(Ptr{Nothing},Int32,),$(esc(env)),$(esc(feature))))
     if res != 0
       throw(MosekError(res,""))
     end
     nothing
  end
end
macro MSK_checkinall(env)
  quote
     local res = disable_sigint(()->ccall((:MSK_checkinall,libmosek),Int32,(Ptr{Nothing},),$(esc(env))))
     if res != 0
       throw(MosekError(res,""))
     end
     nothing
  end
end
macro MSK_expirylicenses(env,expiry)
  quote
     local res = disable_sigint(()->ccall((:MSK_expirylicenses,libmosek),Int32,(Ptr{Nothing},Ref{Int64},),$(esc(env)),$(esc(expiry))))
     if res != 0
       throw(MosekError(res,""))
     end
     nothing
  end
end
macro MSK_resetexpirylicenses(env)
  quote
     local res = disable_sigint(()->ccall((:MSK_resetexpirylicenses,libmosek),Int32,(Ptr{Nothing},),$(esc(env))))
     if res != 0
       throw(MosekError(res,""))
     end
     nothing
  end
end
macro MSK_getbuildinfo(buildstate,builddate)
  quote
     local res = disable_sigint(()->ccall((:MSK_getbuildinfo,libmosek),Int32,(Ptr{UInt8},Ptr{UInt8},),$(esc(buildstate)),$(esc(builddate))))
     if res != 0
       throw(MosekError(res,""))
     end
     nothing
  end
end
macro MSK_getresponseclass(r,rc)
  quote
     local res = disable_sigint(()->ccall((:MSK_getresponseclass,libmosek),Int32,(Int32,Ref{Int32},),$(esc(r)),$(esc(rc))))
     if res != 0
       throw(MosekError(res,""))
     end
     nothing
  end
end
macro MSK_deleteenv(env)
  quote
     local res = disable_sigint(()->ccall((:MSK_deleteenv,libmosek),Int32,(Ref{Ptr{Nothing}},),$(esc(env))))
     if res != 0
       throw(MosekError(res,""))
     end
     nothing
  end
end
macro MSK_echointro(env,longver)
  quote
     local res = disable_sigint(()->ccall((:MSK_echointro,libmosek),Int32,(Ptr{Nothing},Int32,),$(esc(env)),$(esc(longver))))
     if res != 0
       throw(MosekError(res,""))
     end
     nothing
  end
end
macro MSK_freeenv(env,buffer)
  quote
     disable_sigint(()->ccall((:MSK_freeenv,libmosek),Cvoid,(Ptr{Nothing},Ptr{Cvoid},),$(esc(env)),$(esc(buffer))))
  end
end
macro MSK_freedbgenv(env,buffer,file,line)
  quote
     disable_sigint(()->ccall((:MSK_freedbgenv,libmosek),Cvoid,(Ptr{Nothing},Ptr{Cvoid},Ptr{UInt8},UInt32,),$(esc(env)),$(esc(buffer)),$(esc(file)),$(esc(line))))
  end
end
macro MSK_getcodedesc(code,symname,str)
  quote
     local res = disable_sigint(()->ccall((:MSK_getcodedesc,libmosek),Int32,(Int32,Ptr{UInt8},Ptr{UInt8},),$(esc(code)),$(esc(symname)),$(esc(str))))
     if res != 0
       throw(MosekError(res,""))
     end
     nothing
  end
end
macro MSK_getsymbcondim(env,num,maxlen)
  quote
     local res = disable_sigint(()->ccall((:MSK_getsymbcondim,libmosek),Int32,(Ptr{Nothing},Ref{Int32},Ref{CSize},),$(esc(env)),$(esc(num)),$(esc(maxlen))))
     if res != 0
       throw(MosekError(res,""))
     end
     nothing
  end
end
macro MSK_rescodetostr(res,str)
  quote
     local res = disable_sigint(()->ccall((:MSK_rescodetostr,libmosek),Int32,(Int32,Ptr{UInt8},),$(esc(res)),$(esc(str))))
     if res != 0
       throw(MosekError(res,""))
     end
     nothing
  end
end
macro MSK_iinfitemtostr(item,str)
  quote
     local res = disable_sigint(()->ccall((:MSK_iinfitemtostr,libmosek),Int32,(Int32,Ptr{UInt8},),$(esc(item)),$(esc(str))))
     if res != 0
       throw(MosekError(res,""))
     end
     nothing
  end
end
macro MSK_dinfitemtostr(item,str)
  quote
     local res = disable_sigint(()->ccall((:MSK_dinfitemtostr,libmosek),Int32,(Int32,Ptr{UInt8},),$(esc(item)),$(esc(str))))
     if res != 0
       throw(MosekError(res,""))
     end
     nothing
  end
end
macro MSK_liinfitemtostr(item,str)
  quote
     local res = disable_sigint(()->ccall((:MSK_liinfitemtostr,libmosek),Int32,(Int32,Ptr{UInt8},),$(esc(item)),$(esc(str))))
     if res != 0
       throw(MosekError(res,""))
     end
     nothing
  end
end
macro MSK_getversion(major,minor,revision)
  quote
     local res = disable_sigint(()->ccall((:MSK_getversion,libmosek),Int32,(Ref{Int32},Ref{Int32},Ref{Int32},),$(esc(major)),$(esc(minor)),$(esc(revision))))
     if res != 0
       throw(MosekError(res,""))
     end
     nothing
  end
end
macro MSK_checkversion(env,major,minor,revision)
  quote
     local res = disable_sigint(()->ccall((:MSK_checkversion,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Int32,),$(esc(env)),$(esc(major)),$(esc(minor)),$(esc(revision))))
     if res != 0
       throw(MosekError(res,""))
     end
     nothing
  end
end
macro MSK_linkfiletoenvstream(env,whichstream,filename,append)
  quote
     local res = disable_sigint(()->ccall((:MSK_linkfiletoenvstream,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{UInt8},Int32,),$(esc(env)),$(esc(whichstream)),$(esc(filename)),$(esc(append))))
     if res != 0
       throw(MosekError(res,""))
     end
     nothing
  end
end
macro MSK_linkfunctoenvstream(env,whichstream,handle,func)
  quote
     local res = disable_sigint(()->ccall((:MSK_linkfunctoenvstream,libmosek),Int32,(Ptr{Nothing},Int32,Any,Ptr{Cvoid},),$(esc(env)),$(esc(whichstream)),$(esc(handle)),$(esc(func))))
     if res != 0
       throw(MosekError(res,""))
     end
     nothing
  end
end
macro MSK_unlinkfuncfromenvstream(env,whichstream)
  quote
     local res = disable_sigint(()->ccall((:MSK_unlinkfuncfromenvstream,libmosek),Int32,(Ptr{Nothing},Int32,),$(esc(env)),$(esc(whichstream))))
     if res != 0
       throw(MosekError(res,""))
     end
     nothing
  end
end
macro MSK_makeenv(env,dbgfile)
  quote
     local res = disable_sigint(()->ccall((:MSK_makeenv,libmosek),Int32,(Ref{Ptr{Nothing}},Ptr{UInt8},),$(esc(env)),$(esc(dbgfile))))
     if res != 0
       throw(MosekError(res,""))
     end
     nothing
  end
end
macro MSK_putlicensedebug(env,licdebug)
  quote
     local res = disable_sigint(()->ccall((:MSK_putlicensedebug,libmosek),Int32,(Ptr{Nothing},Int32,),$(esc(env)),$(esc(licdebug))))
     if res != 0
       throw(MosekError(res,""))
     end
     nothing
  end
end
macro MSK_putlicensecode(env,code)
  quote
     local res = disable_sigint(()->ccall((:MSK_putlicensecode,libmosek),Int32,(Ptr{Nothing},Ptr{Int32},),$(esc(env)),$(esc(code))))
     if res != 0
       throw(MosekError(res,""))
     end
     nothing
  end
end
macro MSK_putlicensewait(env,licwait)
  quote
     local res = disable_sigint(()->ccall((:MSK_putlicensewait,libmosek),Int32,(Ptr{Nothing},Int32,),$(esc(env)),$(esc(licwait))))
     if res != 0
       throw(MosekError(res,""))
     end
     nothing
  end
end
macro MSK_putlicensepath(env,licensepath)
  quote
     local res = disable_sigint(()->ccall((:MSK_putlicensepath,libmosek),Int32,(Ptr{Nothing},Ptr{UInt8},),$(esc(env)),$(esc(licensepath))))
     if res != 0
       throw(MosekError(res,""))
     end
     nothing
  end
end
macro MSK_maketask(env,maxnumcon,maxnumvar,task)
  quote
     local res = disable_sigint(()->ccall((:MSK_maketask,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Ref{Ptr{Nothing}},),$(esc(env)),$(esc(maxnumcon)),$(esc(maxnumvar)),$(esc(task))))
     if res != 0
       throw(MosekError(res,""))
     end
     nothing
  end
end
macro MSK_makeemptytask(env,task)
  quote
     local res = disable_sigint(()->ccall((:MSK_makeemptytask,libmosek),Int32,(Ptr{Nothing},Ref{Ptr{Nothing}},),$(esc(env)),$(esc(task))))
     if res != 0
       throw(MosekError(res,""))
     end
     nothing
  end
end
macro MSK_putexitfunc(env,exitfunc,handle)
  quote
     local res = disable_sigint(()->ccall((:MSK_putexitfunc,libmosek),Int32,(Ptr{Nothing},Ptr{Cvoid},Any,),$(esc(env)),$(esc(exitfunc)),$(esc(handle))))
     if res != 0
       throw(MosekError(res,""))
     end
     nothing
  end
end
macro MSK_checkmemenv(env,file,line)
  quote
     local res = disable_sigint(()->ccall((:MSK_checkmemenv,libmosek),Int32,(Ptr{Nothing},Ptr{UInt8},Int32,),$(esc(env)),$(esc(file)),$(esc(line))))
     if res != 0
       throw(MosekError(res,""))
     end
     nothing
  end
end
macro MSK_axpy(env,n,alpha,x,y)
  quote
     local res = disable_sigint(()->ccall((:MSK_axpy,libmosek),Int32,(Ptr{Nothing},Int32,Float64,Ptr{Float64},Ptr{Float64},),$(esc(env)),$(esc(n)),$(esc(alpha)),$(esc(x)),$(esc(y))))
     if res != 0
       throw(MosekError(res,""))
     end
     nothing
  end
end
macro MSK_dot(env,n,x,y,xty)
  quote
     local res = disable_sigint(()->ccall((:MSK_dot,libmosek),Int32,(Ptr{Nothing},Int32,Ptr{Float64},Ptr{Float64},Ref{Float64},),$(esc(env)),$(esc(n)),$(esc(x)),$(esc(y)),$(esc(xty))))
     if res != 0
       throw(MosekError(res,""))
     end
     nothing
  end
end
macro MSK_gemv(env,transa,m,n,alpha,a,x,beta,y)
  quote
     local res = disable_sigint(()->ccall((:MSK_gemv,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Int32,Float64,Ptr{Float64},Ptr{Float64},Float64,Ptr{Float64},),$(esc(env)),$(esc(transa)),$(esc(m)),$(esc(n)),$(esc(alpha)),$(esc(a)),$(esc(x)),$(esc(beta)),$(esc(y))))
     if res != 0
       throw(MosekError(res,""))
     end
     nothing
  end
end
macro MSK_gemm(env,transa,transb,m,n,k,alpha,a,b,beta,c)
  quote
     local res = disable_sigint(()->ccall((:MSK_gemm,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Int32,Int32,Int32,Float64,Ptr{Float64},Ptr{Float64},Float64,Ptr{Float64},),$(esc(env)),$(esc(transa)),$(esc(transb)),$(esc(m)),$(esc(n)),$(esc(k)),$(esc(alpha)),$(esc(a)),$(esc(b)),$(esc(beta)),$(esc(c))))
     if res != 0
       throw(MosekError(res,""))
     end
     nothing
  end
end
macro MSK_syrk(env,uplo,trans,n,k,alpha,a,beta,c)
  quote
     local res = disable_sigint(()->ccall((:MSK_syrk,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Int32,Int32,Float64,Ptr{Float64},Float64,Ptr{Float64},),$(esc(env)),$(esc(uplo)),$(esc(trans)),$(esc(n)),$(esc(k)),$(esc(alpha)),$(esc(a)),$(esc(beta)),$(esc(c))))
     if res != 0
       throw(MosekError(res,""))
     end
     nothing
  end
end
macro MSK_computesparsecholesky(env,numthreads,ordermethod,tolsingular,n,anzc,aptrc,asubc,avalc,perm,diag,lnzc,lptrc,lensubnval,lsubc,lvalc)
  quote
     local res = disable_sigint(()->ccall((:MSK_computesparsecholesky,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Float64,Int32,Ptr{Int32},Ptr{Int64},Ptr{Int32},Ptr{Float64},Ref{Ptr{Int32}},Ref{Ptr{Float64}},Ref{Ptr{Int32}},Ref{Ptr{Int64}},Ref{Int64},Ref{Ptr{Int32}},Ref{Ptr{Float64}},),$(esc(env)),$(esc(numthreads)),$(esc(ordermethod)),$(esc(tolsingular)),$(esc(n)),$(esc(anzc)),$(esc(aptrc)),$(esc(asubc)),$(esc(avalc)),$(esc(perm)),$(esc(diag)),$(esc(lnzc)),$(esc(lptrc)),$(esc(lensubnval)),$(esc(lsubc)),$(esc(lvalc))))
     if res != 0
       throw(MosekError(res,""))
     end
     nothing
  end
end
macro MSK_sparsetriangularsolvedense(env,transposed,n,lnzc,lptrc,lensubnval,lsubc,lvalc,b)
  quote
     local res = disable_sigint(()->ccall((:MSK_sparsetriangularsolvedense,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Ptr{Int32},Ptr{Int64},Int64,Ptr{Int32},Ptr{Float64},Ptr{Float64},),$(esc(env)),$(esc(transposed)),$(esc(n)),$(esc(lnzc)),$(esc(lptrc)),$(esc(lensubnval)),$(esc(lsubc)),$(esc(lvalc)),$(esc(b))))
     if res != 0
       throw(MosekError(res,""))
     end
     nothing
  end
end
macro MSK_potrf(env,uplo,n,a)
  quote
     local res = disable_sigint(()->ccall((:MSK_potrf,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Ptr{Float64},),$(esc(env)),$(esc(uplo)),$(esc(n)),$(esc(a))))
     if res != 0
       throw(MosekError(res,""))
     end
     nothing
  end
end
macro MSK_syeig(env,uplo,n,a,w)
  quote
     local res = disable_sigint(()->ccall((:MSK_syeig,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Ptr{Float64},Ptr{Float64},),$(esc(env)),$(esc(uplo)),$(esc(n)),$(esc(a)),$(esc(w))))
     if res != 0
       throw(MosekError(res,""))
     end
     nothing
  end
end
macro MSK_syevd(env,uplo,n,a,w)
  quote
     local res = disable_sigint(()->ccall((:MSK_syevd,libmosek),Int32,(Ptr{Nothing},Int32,Int32,Ptr{Float64},Ptr{Float64},),$(esc(env)),$(esc(uplo)),$(esc(n)),$(esc(a)),$(esc(w))))
     if res != 0
       throw(MosekError(res,""))
     end
     nothing
  end
end
macro MSK_licensecleanup()
  quote
     local res = disable_sigint(()->ccall((:MSK_licensecleanup,libmosek),Int32,()))
     if res != 0
       throw(MosekError(res,""))
     end
     nothing
  end
end
