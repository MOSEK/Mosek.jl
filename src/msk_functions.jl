# Contents of this file is generated. Do not edit by hand!
# MOSEK 7.0.0.92

export
  analyzenames,
  analyzeproblem,
  analyzesolution,
  appendbarvars,
  appendcone,
  appendconeseq,
  appendconesseq,
  appendcons,
  appendsparsesymmat,
  appendstat,
  appendvars,
  basiscond,
  bktostr,
  callbackcodetostr,
  checkconvexity,
  checkmem,
  chgbound,
  chgconbound,
  chgvarbound,
  commitchanges,
  conetypetostr,
  deletesolution,
  dualsensitivity,
  getacol,
  getacolnumnz,
  getacolslicetrip,
  getaij,
  getapiecenumnz,
  getarow,
  getarownumnz,
  getarowslicetrip,
  getaslice,
  getaslicenumnz,
  getbarablocktriplet,
  getbaraidx,
  getbaraidxij,
  getbaraidxinfo,
  getbarasparsity,
  getbarcblocktriplet,
  getbarcidx,
  getbarcidxinfo,
  getbarcidxj,
  getbarcsparsity,
  getbarsj,
  getbarvarname,
  getbarvarnameindex,
  getbarvarnamelen,
  getbarxj,
  getbound,
  getboundslice,
  getc,
  getcfix,
  getcj,
  getconbound,
  getconboundslice,
  getcone,
  getconeinfo,
  getconename,
  getconenameindex,
  getconenamelen,
  getconname,
  getconnameindex,
  getconnamelen,
  getcslice,
  getdbi,
  getdcni,
  getdeqi,
  getdimbarvarj,
  getdouinf,
  getdouparam,
  getdualobj,
  getdviolbarvar,
  getdviolcon,
  getdviolcones,
  getdviolvar,
  getinfeasiblesubproblem,
  getinfname,
  getinti,
  getintinf,
  getintparam,
  getlenbarvarj,
  getlintinf,
  getmaxnumanz,
  getmaxnumbarvar,
  getmaxnumcon,
  getmaxnumcone,
  getmaxnumqnz,
  getmaxnumvar,
  getmemusage,
  getnadouinf,
  getnadouparam,
  getnaintinf,
  getnaintparam,
  getnastrparam,
  getnumanz,
  getnumanz64,
  getnumbarablocktriplets,
  getnumbaranz,
  getnumbarcblocktriplets,
  getnumbarcnz,
  getnumbarvar,
  getnumcon,
  getnumcone,
  getnumconemem,
  getnumintvar,
  getnumparam,
  getnumqconknz,
  getnumqconknz64,
  getnumqobjnz,
  getnumsymmat,
  getnumvar,
  getobjname,
  getobjnamelen,
  getobjsense,
  getparamname,
  getpbi,
  getpcni,
  getpeqi,
  getprimalobj,
  getprobtype,
  getprosta,
  getpviolbarvar,
  getpviolcon,
  getpviolcones,
  getpviolvar,
  getqconk,
  getqobj,
  getqobj64,
  getqobjij,
  getreducedcosts,
  getskc,
  getskcslice,
  getskx,
  getskxslice,
  getslc,
  getslcslice,
  getslx,
  getslxslice,
  getsnx,
  getsnxslice,
  getsolsta,
  getsolution,
  getsolutioni,
  getsolutioninf,
  getsolutioninfo,
  getsolutionslice,
  getsparsesymmat,
  getstrparam,
  getstrparamlen,
  getsuc,
  getsucslice,
  getsux,
  getsuxslice,
  getsymmatinfo,
  gettaskname,
  gettasknamelen,
  getvarbound,
  getvarboundslice,
  getvarbranchdir,
  getvarbranchpri,
  getvarname,
  getvarnameindex,
  getvarnamelen,
  getvartype,
  getvartypelist,
  getxc,
  getxcslice,
  getxx,
  getxxslice,
  gety,
  getyslice,
  initbasissolve,
  inputdata,
  isdouparname,
  isintparname,
  isstrparname,
  linkfiletostream,
  onesolutionsummary,
  optimizeconcurrent,
  optimizersummary,
  optimize,
  primalrepair,
  primalsensitivity,
  printdata,
  printparam,
  putacol,
  putacollist,
  putacolslice,
  putaij,
  putaijlist,
  putaijlist64,
  putarow,
  putarowlist,
  putarowslice,
  putbarablocktriplet,
  putbaraij,
  putbarcblocktriplet,
  putbarcj,
  putbarsj,
  putbarvarname,
  putbarxj,
  putbound,
  putboundlist,
  putboundslice,
  putcfix,
  putcj,
  putclist,
  putconbound,
  putconboundlist,
  putconboundslice,
  putcone,
  putconename,
  putconname,
  putcslice,
  putdouparam,
  putintparam,
  putmaxnumanz,
  putmaxnumbarvar,
  putmaxnumcon,
  putmaxnumcone,
  putmaxnumqnz,
  putmaxnumvar,
  putnadouparam,
  putnaintparam,
  putnastrparam,
  putobjname,
  putobjsense,
  putparam,
  putqcon,
  putqconk,
  putqobj,
  putqobjij,
  putskc,
  putskcslice,
  putskx,
  putskxslice,
  putslc,
  putslcslice,
  putslx,
  putslxslice,
  putsnx,
  putsnxslice,
  putsolution,
  putsolutioni,
  putsolutionyi,
  putstrparam,
  putsuc,
  putsucslice,
  putsux,
  putsuxslice,
  puttaskname,
  putvarbound,
  putvarboundlist,
  putvarboundslice,
  putvarbranchorder,
  putvarname,
  putvartype,
  putvartypelist,
  putxc,
  putxcslice,
  putxx,
  putxxslice,
  puty,
  putyslice,
  readbranchpriorities,
  readdata,
  readdataformat,
  readparamfile,
  readsolution,
  readsummary,
  readtask,
  relaxprimal,
  removebarvars,
  removecones,
  removecons,
  removevars,
  resizetask,
  sensitivityreport,
  setdefaults,
  solutiondef,
  solutionsummary,
  solvewithbasis,
  startstat,
  stopstat,
  strtoconetype,
  strtosk,
  updatesolutioninfo,
  writebranchpriorities,
  writedata,
  writeparamfile,
  writesolution,
  writetask,
  checkinlicense,
  checkoutlicense,
  echointro,
  getcodedesc,
  getversion,
  licensecleanup,
  linkfiletostream,
  putdllpath,
  putkeepdlls,
  putlicensecode,
  putlicensedebug,
  putlicensepath,
  putlicensewait

function analyzenames(task_:: MSKtask,whichstream_:: Int32,nametype_:: Int32)
  res = @msk_ccall( "analyzenames",Int32,(Ptr{Void},Int32,Int32,),task_.task,whichstream_,nametype_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

function analyzeproblem(task_:: MSKtask,whichstream_:: Int32)
  res = @msk_ccall( "analyzeproblem",Int32,(Ptr{Void},Int32,),task_.task,whichstream_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

function analyzesolution(task_:: MSKtask,whichstream_:: Int32,whichsol_:: Int32)
  res = @msk_ccall( "analyzesolution",Int32,(Ptr{Void},Int32,Int32,),task_.task,whichstream_,whichsol_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

appendbarvars(task:: MSKtask,dim:: Array) = appendbarvars(convert(MSKtask,task),convert(Array{Int32},dim))
function appendbarvars(task_:: MSKtask,dim_:: Array{Int32})
  num_ = minimum([ length(dim_) ])
  __tmp_var_0 = if (typeof(dim_) != Array{Int32}) convert(Array{Int32},dim_) else dim_ end
  res = @msk_ccall( "appendbarvars",Int32,(Ptr{Void},Int32,Ptr{Int32},),task_.task,num_,__tmp_var_0)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

appendcone(task:: MSKtask,conetype:: Int32,conepar:: Float64,submem:: Array) = appendcone(convert(MSKtask,task),convert(Int32,conetype),convert(Float64,conepar),convert(Array{Int32},submem))
function appendcone(task_:: MSKtask,conetype_:: Int32,conepar_:: Float64,submem_:: Array{Int32})
  nummem_ = minimum([ length(submem_) ])
  __tmp_var_0 = if (typeof(submem_) != Array{Int32}) convert(Array{Int32},submem_) else submem_ end
  res = @msk_ccall( "appendcone",Int32,(Ptr{Void},Int32,Float64,Int32,Ptr{Int32},),task_.task,conetype_,conepar_,nummem_,__tmp_var_0-1)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

appendconeseq(task:: MSKtask,conetype:: Int32,conepar:: Float64,nummem,j) = appendconeseq(convert(MSKtask,task),convert(Int32,conetype),convert(Float64,conepar),convert(Int32,nummem),convert(Int32,j))
function appendconeseq(task_:: MSKtask,conetype_:: Int32,conepar_:: Float64,nummem_:: Int32,j_:: Int32)
  res = @msk_ccall( "appendconeseq",Int32,(Ptr{Void},Int32,Float64,Int32,Int32,),task_.task,conetype_,conepar_,nummem_,j_-1)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

appendconesseq(task:: MSKtask,conetype:: Array{Int32},conepar:: Array{Float64},nummem:: Array,j) = appendconesseq(convert(MSKtask,task),convert(Array{Int32},conetype),convert(Array{Float64},conepar),convert(Array{Int32},nummem),convert(Int32,j))
function appendconesseq(task_:: MSKtask,conetype_:: Array{Int32},conepar_:: Array{Float64},nummem_:: Array{Int32},j_:: Int32)
  num_ = minimum([ length(conetype_),length(conepar_),length(nummem_) ])
  __tmp_var_0 = if (typeof(nummem_) != Array{Int32}) convert(Array{Int32},nummem_) else nummem_ end
  res = @msk_ccall( "appendconesseq",Int32,(Ptr{Void},Int32,Ptr{Int32},Ptr{Float64},Ptr{Int32},Int32,),task_.task,num_,conetype_,conepar_,__tmp_var_0,j_-1)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

appendcons(task:: MSKtask,num) = appendcons(convert(MSKtask,task),convert(Int32,num))
function appendcons(task_:: MSKtask,num_:: Int32)
  res = @msk_ccall( "appendcons",Int32,(Ptr{Void},Int32,),task_.task,num_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

appendsparsesymmat(task:: MSKtask,dim,subi:: Array,subj:: Array,valij:: Array{Float64}) = appendsparsesymmat(convert(MSKtask,task),convert(Int32,dim),convert(Array{Int32},subi),convert(Array{Int32},subj),convert(Array{Float64},valij))
function appendsparsesymmat(task_:: MSKtask,dim_:: Int32,subi_:: Array{Int32},subj_:: Array{Int32},valij_:: Array{Float64})
  nz_ = minimum([ length(subi_),length(subj_),length(valij_) ])
  __tmp_var_0 = if (typeof(subi_) != Array{Int32}) convert(Array{Int32},subi_) else subi_ end
  __tmp_var_1 = if (typeof(subj_) != Array{Int32}) convert(Array{Int32},subj_) else subj_ end
  idx_ = Array(Int64,(1,))
  res = @msk_ccall( "appendsparsesymmat",Int32,(Ptr{Void},Int32,Int64,Ptr{Int32},Ptr{Int32},Ptr{Float64},Ptr{Int64},),task_.task,dim_,nz_,__tmp_var_0-1,__tmp_var_1-1,valij_,idx_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int64,idx_[1]+1))
end

function appendstat(task_:: MSKtask)
  res = @msk_ccall( "appendstat",Int32,(Ptr{Void},),task_.task)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

appendvars(task:: MSKtask,num) = appendvars(convert(MSKtask,task),convert(Int32,num))
function appendvars(task_:: MSKtask,num_:: Int32)
  res = @msk_ccall( "appendvars",Int32,(Ptr{Void},Int32,),task_.task,num_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

function basiscond(task_:: MSKtask)
  nrmbasis_ = Array(Float64,(1,))
  nrminvbasis_ = Array(Float64,(1,))
  res = @msk_ccall( "basiscond",Int32,(Ptr{Void},Ptr{Float64},Ptr{Float64},),task_.task,nrmbasis_,nrminvbasis_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Float64,nrmbasis_[1]),convert(Float64,nrminvbasis_[1]))
end

function bktostr(task_:: MSKtask,bk_:: Int32)
  str_ = zeros(Uint8,MSK_MAX_STR_LEN)
  res = @msk_ccall( "bktostr",Int32,(Ptr{Void},Int32,Ptr{Uint8},),task_.task,bk_,str_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (bytestring(str_))
end

function callbackcodetostr(code_:: Int32)
  callbackcodestr_ = zeros(Uint8,MSK_MAX_STR_LEN)
  res = @msk_ccall( "callbackcodetostr",Int32,(Int32,Ptr{Uint8},),code_,callbackcodestr_)
  if res != 0
    throw (MosekError(res,""))
  end
  (bytestring(callbackcodestr_))
end

function checkconvexity(task_:: MSKtask)
  res = @msk_ccall( "checkconvexity",Int32,(Ptr{Void},),task_.task)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

checkmem(task:: MSKtask,file:: String,line) = checkmem(convert(MSKtask,task),convert(String,file),convert(Int32,line))
function checkmem(task_:: MSKtask,file_:: String,line_:: Int32)
  res = @msk_ccall( "checkmemtask",Int32,(Ptr{Void},Ptr{Uint8},Int32,),task_.task,bytestring(file_),line_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

chgbound(task:: MSKtask,accmode:: Int32,i,lower,finite,value:: Float64) = chgbound(convert(MSKtask,task),convert(Int32,accmode),convert(Int32,i),convert(Int32,lower),convert(Int32,finite),convert(Float64,value))
function chgbound(task_:: MSKtask,accmode_:: Int32,i_:: Int32,lower_:: Int32,finite_:: Int32,value_:: Float64)
  res = @msk_ccall( "chgbound",Int32,(Ptr{Void},Int32,Int32,Int32,Int32,Float64,),task_.task,accmode_,i_-1,lower_,finite_,value_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

chgconbound(task:: MSKtask,i,lower,finite,value:: Float64) = chgconbound(convert(MSKtask,task),convert(Int32,i),convert(Int32,lower),convert(Int32,finite),convert(Float64,value))
function chgconbound(task_:: MSKtask,i_:: Int32,lower_:: Int32,finite_:: Int32,value_:: Float64)
  res = @msk_ccall( "chgconbound",Int32,(Ptr{Void},Int32,Int32,Int32,Float64,),task_.task,i_-1,lower_,finite_,value_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

chgvarbound(task:: MSKtask,j,lower,finite,value:: Float64) = chgvarbound(convert(MSKtask,task),convert(Int32,j),convert(Int32,lower),convert(Int32,finite),convert(Float64,value))
function chgvarbound(task_:: MSKtask,j_:: Int32,lower_:: Int32,finite_:: Int32,value_:: Float64)
  res = @msk_ccall( "chgvarbound",Int32,(Ptr{Void},Int32,Int32,Int32,Float64,),task_.task,j_-1,lower_,finite_,value_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

function commitchanges(task_:: MSKtask)
  res = @msk_ccall( "commitchanges",Int32,(Ptr{Void},),task_.task)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

function conetypetostr(task_:: MSKtask,conetype_:: Int32)
  str_ = zeros(Uint8,1024)
  res = @msk_ccall( "conetypetostr",Int32,(Ptr{Void},Int32,Ptr{Uint8},),task_.task,conetype_,str_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (bytestring(str_))
end

function deletesolution(task_:: MSKtask,whichsol_:: Int32)
  res = @msk_ccall( "deletesolution",Int32,(Ptr{Void},Int32,),task_.task,whichsol_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

dualsensitivity(task:: MSKtask,subj:: Array) = dualsensitivity(convert(MSKtask,task),convert(Array{Int32},subj))
function dualsensitivity(task_:: MSKtask,subj_:: Array{Int32})
  numj_ = minimum([ length(subj_) ])
  __tmp_var_0 = if (typeof(subj_) != Array{Int32}) convert(Array{Int32},subj_) else subj_ end
  __tmp_var_1 = (numj_)
  __tmp_var_2 = zeros(Float64,__tmp_var_1)
  __tmp_var_3 = (numj_)
  __tmp_var_4 = zeros(Float64,__tmp_var_3)
  __tmp_var_5 = (numj_)
  __tmp_var_6 = zeros(Float64,__tmp_var_5)
  __tmp_var_7 = (numj_)
  __tmp_var_8 = zeros(Float64,__tmp_var_7)
  res = @msk_ccall( "dualsensitivity",Int32,(Ptr{Void},Int32,Ptr{Int32},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},),task_.task,numj_,__tmp_var_0-1,__tmp_var_2,__tmp_var_4,__tmp_var_6,__tmp_var_8)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (__tmp_var_2,__tmp_var_4,__tmp_var_6,__tmp_var_8)
end

getacol(task:: MSKtask,j) = getacol(convert(MSKtask,task),convert(Int32,j))
function getacol(task_:: MSKtask,j_:: Int32)
  nzj_ = Array(Int32,(1,))
  __tmp_var_0 = getacolnumnz(task_,(j_))
  __tmp_var_1 = zeros(Int32,__tmp_var_0)
  __tmp_var_2 = getacolnumnz(task_,(j_))
  __tmp_var_3 = zeros(Float64,__tmp_var_2)
  res = @msk_ccall( "getacol",Int32,(Ptr{Void},Int32,Ptr{Int32},Ptr{Int32},Ptr{Float64},),task_.task,j_-1,nzj_,__tmp_var_1,__tmp_var_3)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,nzj_[1]),__tmp_var_1,__tmp_var_3)
end

getacolnumnz(task:: MSKtask,i) = getacolnumnz(convert(MSKtask,task),convert(Int32,i))
function getacolnumnz(task_:: MSKtask,i_:: Int32)
  nzj_ = Array(Int32,(1,))
  res = @msk_ccall( "getacolnumnz",Int32,(Ptr{Void},Int32,Ptr{Int32},),task_.task,i_-1,nzj_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,nzj_[1]))
end

getacolslicetrip(task:: MSKtask,first,last) = getacolslicetrip(convert(MSKtask,task),convert(Int32,first),convert(Int32,last))
function getacolslicetrip(task_:: MSKtask,first_:: Int32,last_:: Int32)
  maxnumnz_ = minimum([ length(subi_),length(subj_),length(val_) ])
  surp_ = convert(Int64,length(subi_))
  __tmp_var_0 = getaslicenumnz(task_,MSK_ACC_CON,(first_),(last_))
  __tmp_var_1 = zeros(Int32,__tmp_var_0)
  __tmp_var_2 = getaslicenumnz(task_,MSK_ACC_CON,(first_),(last_))
  __tmp_var_3 = zeros(Int32,__tmp_var_2)
  __tmp_var_4 = getaslicenumnz(task_,MSK_ACC_CON,(first_),(last_))
  __tmp_var_5 = zeros(Float64,__tmp_var_4)
  res = @msk_ccall( "getacolslicetrip",Int32,(Ptr{Void},Int32,Int32,Int64,Ptr{Int64},Ptr{Int32},Ptr{Int32},Ptr{Float64},),task_.task,first_-1,last_-1,maxnumnz_,&surp_,__tmp_var_1,__tmp_var_3,__tmp_var_5)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (__tmp_var_1,__tmp_var_3,__tmp_var_5)
end

getaij(task:: MSKtask,i,j) = getaij(convert(MSKtask,task),convert(Int32,i),convert(Int32,j))
function getaij(task_:: MSKtask,i_:: Int32,j_:: Int32)
  aij_ = Array(Float64,(1,))
  res = @msk_ccall( "getaij",Int32,(Ptr{Void},Int32,Int32,Ptr{Float64},),task_.task,i_-1,j_-1,aij_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Float64,aij_[1]))
end

getapiecenumnz(task:: MSKtask,firsti,lasti,firstj,lastj) = getapiecenumnz(convert(MSKtask,task),convert(Int32,firsti),convert(Int32,lasti),convert(Int32,firstj),convert(Int32,lastj))
function getapiecenumnz(task_:: MSKtask,firsti_:: Int32,lasti_:: Int32,firstj_:: Int32,lastj_:: Int32)
  numnz_ = Array(Int32,(1,))
  res = @msk_ccall( "getapiecenumnz",Int32,(Ptr{Void},Int32,Int32,Int32,Int32,Ptr{Int32},),task_.task,firsti_-1,lasti_-1,firstj_-1,lastj_-1,numnz_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,numnz_[1]))
end

getarow(task:: MSKtask,i) = getarow(convert(MSKtask,task),convert(Int32,i))
function getarow(task_:: MSKtask,i_:: Int32)
  nzi_ = Array(Int32,(1,))
  __tmp_var_0 = getarownumnz(task_,(i_))
  __tmp_var_1 = zeros(Int32,__tmp_var_0)
  __tmp_var_2 = getarownumnz(task_,(i_))
  __tmp_var_3 = zeros(Float64,__tmp_var_2)
  res = @msk_ccall( "getarow",Int32,(Ptr{Void},Int32,Ptr{Int32},Ptr{Int32},Ptr{Float64},),task_.task,i_-1,nzi_,__tmp_var_1,__tmp_var_3)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,nzi_[1]),__tmp_var_1,__tmp_var_3)
end

getarownumnz(task:: MSKtask,i) = getarownumnz(convert(MSKtask,task),convert(Int32,i))
function getarownumnz(task_:: MSKtask,i_:: Int32)
  nzi_ = Array(Int32,(1,))
  res = @msk_ccall( "getarownumnz",Int32,(Ptr{Void},Int32,Ptr{Int32},),task_.task,i_-1,nzi_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,nzi_[1]))
end

getarowslicetrip(task:: MSKtask,first,last) = getarowslicetrip(convert(MSKtask,task),convert(Int32,first),convert(Int32,last))
function getarowslicetrip(task_:: MSKtask,first_:: Int32,last_:: Int32)
  maxnumnz_ = minimum([ length(subi_),length(subj_),length(val_) ])
  surp_ = convert(Int64,length(subi_))
  __tmp_var_0 = getaslicenumnz(task_,MSK_ACC_CON,(first_),(last_))
  __tmp_var_1 = zeros(Int32,__tmp_var_0)
  __tmp_var_2 = getaslicenumnz(task_,MSK_ACC_CON,(first_),(last_))
  __tmp_var_3 = zeros(Int32,__tmp_var_2)
  __tmp_var_4 = getaslicenumnz(task_,MSK_ACC_CON,(first_),(last_))
  __tmp_var_5 = zeros(Float64,__tmp_var_4)
  res = @msk_ccall( "getarowslicetrip",Int32,(Ptr{Void},Int32,Int32,Int64,Ptr{Int64},Ptr{Int32},Ptr{Int32},Ptr{Float64},),task_.task,first_-1,last_-1,maxnumnz_,&surp_,__tmp_var_1,__tmp_var_3,__tmp_var_5)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (__tmp_var_1,__tmp_var_3,__tmp_var_5)
end

getaslice(task:: MSKtask,accmode:: Int32,first,last) = getaslice(convert(MSKtask,task),convert(Int32,accmode),convert(Int32,first),convert(Int32,last))
function getaslice(task_:: MSKtask,accmode_:: Int32,first_:: Int32,last_:: Int32)
  maxnumnz_ = getaslicenumnz(task_,(accmode_),(first_),(last_))
  surp_ = convert(Int64,length(sub_))
  __tmp_var_0 = ((last_) - (first_))
  __tmp_var_1 = zeros(Int64,__tmp_var_0)
  __tmp_var_2 = ((last_) - (first_))
  __tmp_var_3 = zeros(Int64,__tmp_var_2)
  __tmp_var_4 = (maxnumnz_)
  __tmp_var_5 = zeros(Int32,__tmp_var_4)
  __tmp_var_6 = (maxnumnz_)
  __tmp_var_7 = zeros(Float64,__tmp_var_6)
  res = @msk_ccall( "getaslice64",Int32,(Ptr{Void},Int32,Int32,Int32,Int64,Ptr{Int64},Ptr{Int64},Ptr{Int64},Ptr{Int32},Ptr{Float64},),task_.task,accmode_,first_-1,last_-1,maxnumnz_,&surp_,__tmp_var_1,__tmp_var_3,__tmp_var_5,__tmp_var_7)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (__tmp_var_1,__tmp_var_3,__tmp_var_5,__tmp_var_7)
end

getaslicenumnz(task:: MSKtask,accmode:: Int32,first,last) = getaslicenumnz(convert(MSKtask,task),convert(Int32,accmode),convert(Int32,first),convert(Int32,last))
function getaslicenumnz(task_:: MSKtask,accmode_:: Int32,first_:: Int32,last_:: Int32)
  numnz_ = Array(Int64,(1,))
  res = @msk_ccall( "getaslicenumnz64",Int32,(Ptr{Void},Int32,Int32,Int32,Ptr{Int64},),task_.task,accmode_,first_-1,last_-1,numnz_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int64,numnz_[1]))
end

function getbarablocktriplet(task_:: MSKtask)
  maxnum_ = getnumbarablocktriplets(task_)
  num_ = Array(Int64,(1,))
  __tmp_var_0 = (maxnum_)
  __tmp_var_1 = zeros(Int32,__tmp_var_0)
  __tmp_var_2 = (maxnum_)
  __tmp_var_3 = zeros(Int32,__tmp_var_2)
  __tmp_var_4 = (maxnum_)
  __tmp_var_5 = zeros(Int32,__tmp_var_4)
  __tmp_var_6 = (maxnum_)
  __tmp_var_7 = zeros(Int32,__tmp_var_6)
  __tmp_var_8 = (maxnum_)
  __tmp_var_9 = zeros(Float64,__tmp_var_8)
  res = @msk_ccall( "getbarablocktriplet",Int32,(Ptr{Void},Int64,Ptr{Int64},Ptr{Int32},Ptr{Int32},Ptr{Int32},Ptr{Int32},Ptr{Float64},),task_.task,maxnum_,num_,__tmp_var_1,__tmp_var_3,__tmp_var_5,__tmp_var_7,__tmp_var_9)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int64,num_[1]),__tmp_var_1,__tmp_var_3,__tmp_var_5,__tmp_var_7,__tmp_var_9)
end

getbaraidx(task:: MSKtask,idx) = getbaraidx(convert(MSKtask,task),convert(Int64,idx))
function getbaraidx(task_:: MSKtask,idx_:: Int64)
  maxnum_ = getbaraidxinfo(task_,(idx_))
  i_ = Array(Int32,(1,))
  j_ = Array(Int32,(1,))
  num_ = Array(Int64,(1,))
  __tmp_var_0 = (maxnum_)
  __tmp_var_1 = zeros(Int64,__tmp_var_0)
  __tmp_var_2 = (maxnum_)
  __tmp_var_3 = zeros(Float64,__tmp_var_2)
  res = @msk_ccall( "getbaraidx",Int32,(Ptr{Void},Int64,Int64,Ptr{Int32},Ptr{Int32},Ptr{Int64},Ptr{Int64},Ptr{Float64},),task_.task,idx_-1,maxnum_,i_,j_,num_,__tmp_var_1,__tmp_var_3)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,i_[1]+1),convert(Int32,j_[1]+1),convert(Int64,num_[1]),__tmp_var_1,__tmp_var_3)
end

getbaraidxij(task:: MSKtask,idx) = getbaraidxij(convert(MSKtask,task),convert(Int64,idx))
function getbaraidxij(task_:: MSKtask,idx_:: Int64)
  i_ = Array(Int32,(1,))
  j_ = Array(Int32,(1,))
  res = @msk_ccall( "getbaraidxij",Int32,(Ptr{Void},Int64,Ptr{Int32},Ptr{Int32},),task_.task,idx_-1,i_,j_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,i_[1]+1),convert(Int32,j_[1]+1))
end

getbaraidxinfo(task:: MSKtask,idx) = getbaraidxinfo(convert(MSKtask,task),convert(Int64,idx))
function getbaraidxinfo(task_:: MSKtask,idx_:: Int64)
  num_ = Array(Int64,(1,))
  res = @msk_ccall( "getbaraidxinfo",Int32,(Ptr{Void},Int64,Ptr{Int64},),task_.task,idx_-1,num_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int64,num_[1]))
end

function getbarasparsity(task_:: MSKtask)
  maxnumnz_ = getnumbaranz(task_)
  numnz_ = Array(Int64,(1,))
  __tmp_var_0 = (maxnumnz_)
  __tmp_var_1 = zeros(Int64,__tmp_var_0)
  res = @msk_ccall( "getbarasparsity",Int32,(Ptr{Void},Int64,Ptr{Int64},Ptr{Int64},),task_.task,maxnumnz_,numnz_,__tmp_var_1)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int64,numnz_[1]),__tmp_var_1)
end

function getbarcblocktriplet(task_:: MSKtask)
  maxnum_ = getnumbarcblocktriplets(task_)
  num_ = Array(Int64,(1,))
  __tmp_var_0 = (maxnum_)
  __tmp_var_1 = zeros(Int32,__tmp_var_0)
  __tmp_var_2 = (maxnum_)
  __tmp_var_3 = zeros(Int32,__tmp_var_2)
  __tmp_var_4 = (maxnum_)
  __tmp_var_5 = zeros(Int32,__tmp_var_4)
  __tmp_var_6 = (maxnum_)
  __tmp_var_7 = zeros(Float64,__tmp_var_6)
  res = @msk_ccall( "getbarcblocktriplet",Int32,(Ptr{Void},Int64,Ptr{Int64},Ptr{Int32},Ptr{Int32},Ptr{Int32},Ptr{Float64},),task_.task,maxnum_,num_,__tmp_var_1,__tmp_var_3,__tmp_var_5,__tmp_var_7)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int64,num_[1]),__tmp_var_1,__tmp_var_3,__tmp_var_5,__tmp_var_7)
end

getbarcidx(task:: MSKtask,idx) = getbarcidx(convert(MSKtask,task),convert(Int64,idx))
function getbarcidx(task_:: MSKtask,idx_:: Int64)
  maxnum_ = getbarcidxinfo(task_,(idx_))
  j_ = Array(Int32,(1,))
  num_ = Array(Int64,(1,))
  __tmp_var_0 = (maxnum_)
  __tmp_var_1 = zeros(Int64,__tmp_var_0)
  __tmp_var_2 = (maxnum_)
  __tmp_var_3 = zeros(Float64,__tmp_var_2)
  res = @msk_ccall( "getbarcidx",Int32,(Ptr{Void},Int64,Int64,Ptr{Int32},Ptr{Int64},Ptr{Int64},Ptr{Float64},),task_.task,idx_-1,maxnum_,j_,num_,__tmp_var_1,__tmp_var_3)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,j_[1]+1),convert(Int64,num_[1]),__tmp_var_1,__tmp_var_3)
end

getbarcidxinfo(task:: MSKtask,idx) = getbarcidxinfo(convert(MSKtask,task),convert(Int64,idx))
function getbarcidxinfo(task_:: MSKtask,idx_:: Int64)
  num_ = Array(Int64,(1,))
  res = @msk_ccall( "getbarcidxinfo",Int32,(Ptr{Void},Int64,Ptr{Int64},),task_.task,idx_-1,num_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int64,num_[1]))
end

getbarcidxj(task:: MSKtask,idx) = getbarcidxj(convert(MSKtask,task),convert(Int64,idx))
function getbarcidxj(task_:: MSKtask,idx_:: Int64)
  j_ = Array(Int32,(1,))
  res = @msk_ccall( "getbarcidxj",Int32,(Ptr{Void},Int64,Ptr{Int32},),task_.task,idx_-1,j_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,j_[1]+1))
end

function getbarcsparsity(task_:: MSKtask)
  maxnumnz_ = getnumbarcnz(task_)
  numnz_ = Array(Int64,(1,))
  __tmp_var_0 = (maxnumnz_)
  __tmp_var_1 = zeros(Int64,__tmp_var_0)
  res = @msk_ccall( "getbarcsparsity",Int32,(Ptr{Void},Int64,Ptr{Int64},Ptr{Int64},),task_.task,maxnumnz_,numnz_,__tmp_var_1)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int64,numnz_[1]),__tmp_var_1)
end

getbarsj(task:: MSKtask,whichsol:: Int32,j) = getbarsj(convert(MSKtask,task),convert(Int32,whichsol),convert(Int32,j))
function getbarsj(task_:: MSKtask,whichsol_:: Int32,j_:: Int32)
  __tmp_var_0 = getlenbarvarj(task_,(j_))
  __tmp_var_1 = zeros(Float64,__tmp_var_0)
  res = @msk_ccall( "getbarsj",Int32,(Ptr{Void},Int32,Int32,Ptr{Float64},),task_.task,whichsol_,j_-1,__tmp_var_1)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (__tmp_var_1)
end

getbarvarname(task:: MSKtask,i) = getbarvarname(convert(MSKtask,task),convert(Int32,i))
function getbarvarname(task_:: MSKtask,i_:: Int32)
  maxlen_ = (1 + getbarvarnamelen(task_,(i_)))
  name_ = zeros(Uint8,(maxlen_))
  res = @msk_ccall( "getbarvarname",Int32,(Ptr{Void},Int32,Int32,Ptr{Uint8},),task_.task,i_-1,maxlen_,name_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (bytestring(name_))
end

function getbarvarnameindex(task_:: MSKtask,somename_:: String)
  asgn_ = Array(Int32,(1,))
  index_ = Array(Int32,(1,))
  res = @msk_ccall( "getbarvarnameindex",Int32,(Ptr{Void},Ptr{Uint8},Ptr{Int32},Ptr{Int32},),task_.task,bytestring(somename_),asgn_,index_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,asgn_[1]),convert(Int32,index_[1]))
end

getbarvarnamelen(task:: MSKtask,i) = getbarvarnamelen(convert(MSKtask,task),convert(Int32,i))
function getbarvarnamelen(task_:: MSKtask,i_:: Int32)
  len_ = Array(Int32,(1,))
  res = @msk_ccall( "getbarvarnamelen",Int32,(Ptr{Void},Int32,Ptr{Int32},),task_.task,i_-1,len_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,len_[1]))
end

getbarxj(task:: MSKtask,whichsol:: Int32,j) = getbarxj(convert(MSKtask,task),convert(Int32,whichsol),convert(Int32,j))
function getbarxj(task_:: MSKtask,whichsol_:: Int32,j_:: Int32)
  __tmp_var_0 = getlenbarvarj(task_,(j_))
  __tmp_var_1 = zeros(Float64,__tmp_var_0)
  res = @msk_ccall( "getbarxj",Int32,(Ptr{Void},Int32,Int32,Ptr{Float64},),task_.task,whichsol_,j_-1,__tmp_var_1)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (__tmp_var_1)
end

getbound(task:: MSKtask,accmode:: Int32,i) = getbound(convert(MSKtask,task),convert(Int32,accmode),convert(Int32,i))
function getbound(task_:: MSKtask,accmode_:: Int32,i_:: Int32)
  bk_ = Array(Int32,(1,))
  bl_ = Array(Float64,(1,))
  bu_ = Array(Float64,(1,))
  res = @msk_ccall( "getbound",Int32,(Ptr{Void},Int32,Int32,Ptr{Int32},Ptr{Float64},Ptr{Float64},),task_.task,accmode_,i_-1,bk_,bl_,bu_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,bk_[1]),convert(Float64,bl_[1]),convert(Float64,bu_[1]))
end

getboundslice(task:: MSKtask,accmode:: Int32,first,last) = getboundslice(convert(MSKtask,task),convert(Int32,accmode),convert(Int32,first),convert(Int32,last))
function getboundslice(task_:: MSKtask,accmode_:: Int32,first_:: Int32,last_:: Int32)
  __tmp_var_0 = ((last_) - (first_))
  bk_ = zeros(Int32,__tmp_var_0)
  __tmp_var_1 = ((last_) - (first_))
  __tmp_var_2 = zeros(Float64,__tmp_var_1)
  __tmp_var_3 = ((last_) - (first_))
  __tmp_var_4 = zeros(Float64,__tmp_var_3)
  res = @msk_ccall( "getboundslice",Int32,(Ptr{Void},Int32,Int32,Int32,Ptr{Int32},Ptr{Float64},Ptr{Float64},),task_.task,accmode_,first_-1,last_-1,bk_,__tmp_var_2,__tmp_var_4)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (bk_,__tmp_var_2,__tmp_var_4)
end

function getc(task_:: MSKtask)
  __tmp_var_0 = getnumvar(task_)
  __tmp_var_1 = zeros(Float64,__tmp_var_0)
  res = @msk_ccall( "getc",Int32,(Ptr{Void},Ptr{Float64},),task_.task,__tmp_var_1)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (__tmp_var_1)
end

function getcfix(task_:: MSKtask)
  cfix_ = Array(Float64,(1,))
  res = @msk_ccall( "getcfix",Int32,(Ptr{Void},Ptr{Float64},),task_.task,cfix_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Float64,cfix_[1]))
end

getcj(task:: MSKtask,j) = getcj(convert(MSKtask,task),convert(Int32,j))
function getcj(task_:: MSKtask,j_:: Int32)
  cj_ = Array(Float64,(1,))
  res = @msk_ccall( "getcj",Int32,(Ptr{Void},Int32,Ptr{Float64},),task_.task,j_-1,cj_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Float64,cj_[1]))
end

getconbound(task:: MSKtask,i) = getconbound(convert(MSKtask,task),convert(Int32,i))
function getconbound(task_:: MSKtask,i_:: Int32)
  bk_ = Array(Int32,(1,))
  bl_ = Array(Float64,(1,))
  bu_ = Array(Float64,(1,))
  res = @msk_ccall( "getconbound",Int32,(Ptr{Void},Int32,Ptr{Int32},Ptr{Float64},Ptr{Float64},),task_.task,i_-1,bk_,bl_,bu_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,bk_[1]),convert(Float64,bl_[1]),convert(Float64,bu_[1]))
end

getconboundslice(task:: MSKtask,first,last) = getconboundslice(convert(MSKtask,task),convert(Int32,first),convert(Int32,last))
function getconboundslice(task_:: MSKtask,first_:: Int32,last_:: Int32)
  __tmp_var_0 = ((last_) - (first_))
  bk_ = zeros(Int32,__tmp_var_0)
  __tmp_var_1 = ((last_) - (first_))
  __tmp_var_2 = zeros(Float64,__tmp_var_1)
  __tmp_var_3 = ((last_) - (first_))
  __tmp_var_4 = zeros(Float64,__tmp_var_3)
  res = @msk_ccall( "getconboundslice",Int32,(Ptr{Void},Int32,Int32,Ptr{Int32},Ptr{Float64},Ptr{Float64},),task_.task,first_-1,last_-1,bk_,__tmp_var_2,__tmp_var_4)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (bk_,__tmp_var_2,__tmp_var_4)
end

getcone(task:: MSKtask,k) = getcone(convert(MSKtask,task),convert(Int32,k))
function getcone(task_:: MSKtask,k_:: Int32)
  conetype_ = Array(Int32,(1,))
  conepar_ = Array(Float64,(1,))
  nummem_ = Array(Int32,(1,))
  __tmp_var_0 = getconeinfo(task_,(k_))[3]
  __tmp_var_1 = zeros(Int32,__tmp_var_0)
  res = @msk_ccall( "getcone",Int32,(Ptr{Void},Int32,Ptr{Int32},Ptr{Float64},Ptr{Int32},Ptr{Int32},),task_.task,k_-1,conetype_,conepar_,nummem_,__tmp_var_1)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,conetype_[1]),convert(Float64,conepar_[1]),convert(Int32,nummem_[1]),__tmp_var_1)
end

getconeinfo(task:: MSKtask,k) = getconeinfo(convert(MSKtask,task),convert(Int32,k))
function getconeinfo(task_:: MSKtask,k_:: Int32)
  conetype_ = Array(Int32,(1,))
  conepar_ = Array(Float64,(1,))
  nummem_ = Array(Int32,(1,))
  res = @msk_ccall( "getconeinfo",Int32,(Ptr{Void},Int32,Ptr{Int32},Ptr{Float64},Ptr{Int32},),task_.task,k_-1,conetype_,conepar_,nummem_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,conetype_[1]),convert(Float64,conepar_[1]),convert(Int32,nummem_[1]))
end

getconename(task:: MSKtask,i) = getconename(convert(MSKtask,task),convert(Int32,i))
function getconename(task_:: MSKtask,i_:: Int32)
  maxlen_ = (1 + getconenamelen(task_,(i_)))
  name_ = zeros(Uint8,(maxlen_))
  res = @msk_ccall( "getconename",Int32,(Ptr{Void},Int32,Int32,Ptr{Uint8},),task_.task,i_-1,maxlen_,name_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (bytestring(name_))
end

function getconenameindex(task_:: MSKtask,somename_:: String)
  asgn_ = Array(Int32,(1,))
  index_ = Array(Int32,(1,))
  res = @msk_ccall( "getconenameindex",Int32,(Ptr{Void},Ptr{Uint8},Ptr{Int32},Ptr{Int32},),task_.task,bytestring(somename_),asgn_,index_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,asgn_[1]),convert(Int32,index_[1]))
end

getconenamelen(task:: MSKtask,i) = getconenamelen(convert(MSKtask,task),convert(Int32,i))
function getconenamelen(task_:: MSKtask,i_:: Int32)
  len_ = Array(Int32,(1,))
  res = @msk_ccall( "getconenamelen",Int32,(Ptr{Void},Int32,Ptr{Int32},),task_.task,i_-1,len_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,len_[1]))
end

getconname(task:: MSKtask,i) = getconname(convert(MSKtask,task),convert(Int32,i))
function getconname(task_:: MSKtask,i_:: Int32)
  maxlen_ = (1 + getconnamelen(task_,(i_)))
  name_ = zeros(Uint8,(maxlen_))
  res = @msk_ccall( "getconname",Int32,(Ptr{Void},Int32,Int32,Ptr{Uint8},),task_.task,i_-1,maxlen_,name_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (bytestring(name_))
end

function getconnameindex(task_:: MSKtask,somename_:: String)
  asgn_ = Array(Int32,(1,))
  index_ = Array(Int32,(1,))
  res = @msk_ccall( "getconnameindex",Int32,(Ptr{Void},Ptr{Uint8},Ptr{Int32},Ptr{Int32},),task_.task,bytestring(somename_),asgn_,index_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,asgn_[1]),convert(Int32,index_[1]))
end

getconnamelen(task:: MSKtask,i) = getconnamelen(convert(MSKtask,task),convert(Int32,i))
function getconnamelen(task_:: MSKtask,i_:: Int32)
  len_ = Array(Int32,(1,))
  res = @msk_ccall( "getconnamelen",Int32,(Ptr{Void},Int32,Ptr{Int32},),task_.task,i_-1,len_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,len_[1]))
end

getcslice(task:: MSKtask,first,last) = getcslice(convert(MSKtask,task),convert(Int32,first),convert(Int32,last))
function getcslice(task_:: MSKtask,first_:: Int32,last_:: Int32)
  __tmp_var_0 = ((last_) - (first_))
  __tmp_var_1 = zeros(Float64,__tmp_var_0)
  res = @msk_ccall( "getcslice",Int32,(Ptr{Void},Int32,Int32,Ptr{Float64},),task_.task,first_-1,last_-1,__tmp_var_1)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (__tmp_var_1)
end

getdbi(task:: MSKtask,whichsol:: Int32,accmode:: Int32,sub:: Array) = getdbi(convert(MSKtask,task),convert(Int32,whichsol),convert(Int32,accmode),convert(Array{Int32},sub))
function getdbi(task_:: MSKtask,whichsol_:: Int32,accmode_:: Int32,sub_:: Array{Int32})
  __tmp_var_0 = if (typeof(sub_) != Array{Int32}) convert(Array{Int32},sub_) else sub_ end
  len_ = minimum([ length(sub_) ])
  __tmp_var_1 = (len_)
  __tmp_var_2 = zeros(Float64,__tmp_var_1)
  res = @msk_ccall( "getdbi",Int32,(Ptr{Void},Int32,Int32,Ptr{Int32},Int32,Ptr{Float64},),task_.task,whichsol_,accmode_,__tmp_var_0-1,len_,__tmp_var_2)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (__tmp_var_2)
end

getdcni(task:: MSKtask,whichsol:: Int32,sub:: Array) = getdcni(convert(MSKtask,task),convert(Int32,whichsol),convert(Array{Int32},sub))
function getdcni(task_:: MSKtask,whichsol_:: Int32,sub_:: Array{Int32})
  __tmp_var_0 = if (typeof(sub_) != Array{Int32}) convert(Array{Int32},sub_) else sub_ end
  len_ = minimum([ length(sub_) ])
  __tmp_var_1 = (len_)
  __tmp_var_2 = zeros(Float64,__tmp_var_1)
  res = @msk_ccall( "getdcni",Int32,(Ptr{Void},Int32,Ptr{Int32},Int32,Ptr{Float64},),task_.task,whichsol_,__tmp_var_0-1,len_,__tmp_var_2)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (__tmp_var_2)
end

getdeqi(task:: MSKtask,whichsol:: Int32,accmode:: Int32,sub:: Array,normalize) = getdeqi(convert(MSKtask,task),convert(Int32,whichsol),convert(Int32,accmode),convert(Array{Int32},sub),convert(Int32,normalize))
function getdeqi(task_:: MSKtask,whichsol_:: Int32,accmode_:: Int32,sub_:: Array{Int32},normalize_:: Int32)
  __tmp_var_0 = if (typeof(sub_) != Array{Int32}) convert(Array{Int32},sub_) else sub_ end
  len_ = minimum([ length(sub_) ])
  __tmp_var_1 = (len_)
  __tmp_var_2 = zeros(Float64,__tmp_var_1)
  res = @msk_ccall( "getdeqi",Int32,(Ptr{Void},Int32,Int32,Ptr{Int32},Int32,Ptr{Float64},Int32,),task_.task,whichsol_,accmode_,__tmp_var_0-1,len_,__tmp_var_2,normalize_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (__tmp_var_2)
end

getdimbarvarj(task:: MSKtask,j) = getdimbarvarj(convert(MSKtask,task),convert(Int32,j))
function getdimbarvarj(task_:: MSKtask,j_:: Int32)
  dimbarvarj_ = Array(Int32,(1,))
  res = @msk_ccall( "getdimbarvarj",Int32,(Ptr{Void},Int32,Ptr{Int32},),task_.task,j_-1,dimbarvarj_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,dimbarvarj_[1]))
end

function getdouinf(task_:: MSKtask,whichdinf_:: Int32)
  dvalue_ = Array(Float64,(1,))
  res = @msk_ccall( "getdouinf",Int32,(Ptr{Void},Int32,Ptr{Float64},),task_.task,whichdinf_,dvalue_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Float64,dvalue_[1]))
end

function getdouparam(task_:: MSKtask,param_:: Int32)
  parvalue_ = Array(Float64,(1,))
  res = @msk_ccall( "getdouparam",Int32,(Ptr{Void},Int32,Ptr{Float64},),task_.task,param_,parvalue_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Float64,parvalue_[1]))
end

function getdualobj(task_:: MSKtask,whichsol_:: Int32)
  dualobj_ = Array(Float64,(1,))
  res = @msk_ccall( "getdualobj",Int32,(Ptr{Void},Int32,Ptr{Float64},),task_.task,whichsol_,dualobj_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Float64,dualobj_[1]))
end

getdviolbarvar(task:: MSKtask,whichsol:: Int32,sub:: Array) = getdviolbarvar(convert(MSKtask,task),convert(Int32,whichsol),convert(Array{Int32},sub))
function getdviolbarvar(task_:: MSKtask,whichsol_:: Int32,sub_:: Array{Int32})
  num_ = minimum([ length(sub_) ])
  __tmp_var_0 = if (typeof(sub_) != Array{Int32}) convert(Array{Int32},sub_) else sub_ end
  __tmp_var_1 = (num_)
  __tmp_var_2 = zeros(Float64,__tmp_var_1)
  res = @msk_ccall( "getdviolbarvar",Int32,(Ptr{Void},Int32,Int32,Ptr{Int32},Ptr{Float64},),task_.task,whichsol_,num_,__tmp_var_0-1,__tmp_var_2)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (__tmp_var_2)
end

getdviolcon(task:: MSKtask,whichsol:: Int32,sub:: Array) = getdviolcon(convert(MSKtask,task),convert(Int32,whichsol),convert(Array{Int32},sub))
function getdviolcon(task_:: MSKtask,whichsol_:: Int32,sub_:: Array{Int32})
  num_ = minimum([ length(sub_) ])
  __tmp_var_0 = if (typeof(sub_) != Array{Int32}) convert(Array{Int32},sub_) else sub_ end
  __tmp_var_1 = (num_)
  __tmp_var_2 = zeros(Float64,__tmp_var_1)
  res = @msk_ccall( "getdviolcon",Int32,(Ptr{Void},Int32,Int32,Ptr{Int32},Ptr{Float64},),task_.task,whichsol_,num_,__tmp_var_0-1,__tmp_var_2)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (__tmp_var_2)
end

getdviolcones(task:: MSKtask,whichsol:: Int32,sub:: Array) = getdviolcones(convert(MSKtask,task),convert(Int32,whichsol),convert(Array{Int32},sub))
function getdviolcones(task_:: MSKtask,whichsol_:: Int32,sub_:: Array{Int32})
  num_ = minimum([ length(sub_) ])
  __tmp_var_0 = if (typeof(sub_) != Array{Int32}) convert(Array{Int32},sub_) else sub_ end
  __tmp_var_1 = (num_)
  __tmp_var_2 = zeros(Float64,__tmp_var_1)
  res = @msk_ccall( "getdviolcones",Int32,(Ptr{Void},Int32,Int32,Ptr{Int32},Ptr{Float64},),task_.task,whichsol_,num_,__tmp_var_0-1,__tmp_var_2)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (__tmp_var_2)
end

getdviolvar(task:: MSKtask,whichsol:: Int32,sub:: Array) = getdviolvar(convert(MSKtask,task),convert(Int32,whichsol),convert(Array{Int32},sub))
function getdviolvar(task_:: MSKtask,whichsol_:: Int32,sub_:: Array{Int32})
  num_ = minimum([ length(sub_) ])
  __tmp_var_0 = if (typeof(sub_) != Array{Int32}) convert(Array{Int32},sub_) else sub_ end
  __tmp_var_1 = (num_)
  __tmp_var_2 = zeros(Float64,__tmp_var_1)
  res = @msk_ccall( "getdviolvar",Int32,(Ptr{Void},Int32,Int32,Ptr{Int32},Ptr{Float64},),task_.task,whichsol_,num_,__tmp_var_0-1,__tmp_var_2)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (__tmp_var_2)
end

function getinfeasiblesubproblem(task_:: MSKtask,whichsol_:: Int32)
  inftask_ = Array(Ptr{Void},(1,))
  res = @msk_ccall( "getinfeasiblesubproblem",Int32,(Ptr{Void},Int32,Ptr{Ptr{Void}},),task_.task,whichsol_,inftask_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(MSKtask,inftask_[1]))
end

getinfname(task:: MSKtask,inftype:: Int32,whichinf) = getinfname(convert(MSKtask,task),convert(Int32,inftype),convert(Int32,whichinf))
function getinfname(task_:: MSKtask,inftype_:: Int32,whichinf_:: Int32)
  infname_ = zeros(Uint8,MSK_MAX_STR_LEN)
  res = @msk_ccall( "getinfname",Int32,(Ptr{Void},Int32,Int32,Ptr{Uint8},),task_.task,inftype_,whichinf_,infname_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (bytestring(infname_))
end

getinti(task:: MSKtask,whichsol:: Int32,sub:: Array) = getinti(convert(MSKtask,task),convert(Int32,whichsol),convert(Array{Int32},sub))
function getinti(task_:: MSKtask,whichsol_:: Int32,sub_:: Array{Int32})
  __tmp_var_0 = if (typeof(sub_) != Array{Int32}) convert(Array{Int32},sub_) else sub_ end
  len_ = minimum([ length(sub_) ])
  __tmp_var_1 = (len_)
  __tmp_var_2 = zeros(Float64,__tmp_var_1)
  res = @msk_ccall( "getinti",Int32,(Ptr{Void},Int32,Ptr{Int32},Int32,Ptr{Float64},),task_.task,whichsol_,__tmp_var_0-1,len_,__tmp_var_2)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (__tmp_var_2)
end

function getintinf(task_:: MSKtask,whichiinf_:: Int32)
  ivalue_ = Array(Int32,(1,))
  res = @msk_ccall( "getintinf",Int32,(Ptr{Void},Int32,Ptr{Int32},),task_.task,whichiinf_,ivalue_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,ivalue_[1]))
end

function getintparam(task_:: MSKtask,param_:: Int32)
  parvalue_ = Array(Int32,(1,))
  res = @msk_ccall( "getintparam",Int32,(Ptr{Void},Int32,Ptr{Int32},),task_.task,param_,parvalue_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,parvalue_[1]))
end

getlenbarvarj(task:: MSKtask,j) = getlenbarvarj(convert(MSKtask,task),convert(Int32,j))
function getlenbarvarj(task_:: MSKtask,j_:: Int32)
  lenbarvarj_ = Array(Int64,(1,))
  res = @msk_ccall( "getlenbarvarj",Int32,(Ptr{Void},Int32,Ptr{Int64},),task_.task,j_-1,lenbarvarj_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int64,lenbarvarj_[1]))
end

function getlintinf(task_:: MSKtask,whichliinf_:: Int32)
  ivalue_ = Array(Int64,(1,))
  res = @msk_ccall( "getlintinf",Int32,(Ptr{Void},Int32,Ptr{Int64},),task_.task,whichliinf_,ivalue_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int64,ivalue_[1]))
end

function getmaxnumanz(task_:: MSKtask)
  maxnumanz_ = Array(Int64,(1,))
  res = @msk_ccall( "getmaxnumanz64",Int32,(Ptr{Void},Ptr{Int64},),task_.task,maxnumanz_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int64,maxnumanz_[1]))
end

function getmaxnumbarvar(task_:: MSKtask)
  maxnumbarvar_ = Array(Int32,(1,))
  res = @msk_ccall( "getmaxnumbarvar",Int32,(Ptr{Void},Ptr{Int32},),task_.task,maxnumbarvar_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,maxnumbarvar_[1]))
end

function getmaxnumcon(task_:: MSKtask)
  maxnumcon_ = Array(Int32,(1,))
  res = @msk_ccall( "getmaxnumcon",Int32,(Ptr{Void},Ptr{Int32},),task_.task,maxnumcon_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,maxnumcon_[1]))
end

function getmaxnumcone(task_:: MSKtask)
  maxnumcone_ = Array(Int32,(1,))
  res = @msk_ccall( "getmaxnumcone",Int32,(Ptr{Void},Ptr{Int32},),task_.task,maxnumcone_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,maxnumcone_[1]))
end

function getmaxnumqnz(task_:: MSKtask)
  maxnumqnz_ = Array(Int64,(1,))
  res = @msk_ccall( "getmaxnumqnz64",Int32,(Ptr{Void},Ptr{Int64},),task_.task,maxnumqnz_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int64,maxnumqnz_[1]))
end

function getmaxnumvar(task_:: MSKtask)
  maxnumvar_ = Array(Int32,(1,))
  res = @msk_ccall( "getmaxnumvar",Int32,(Ptr{Void},Ptr{Int32},),task_.task,maxnumvar_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,maxnumvar_[1]))
end

function getmemusage(task_:: MSKtask)
  meminuse_ = Array(Int64,(1,))
  maxmemuse_ = Array(Int64,(1,))
  res = @msk_ccall( "getmemusagetask",Int32,(Ptr{Void},Ptr{Int64},Ptr{Int64},),task_.task,meminuse_,maxmemuse_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int64,meminuse_[1]),convert(Int64,maxmemuse_[1]))
end

function getnadouinf(task_:: MSKtask,whichdinf_:: String)
  dvalue_ = Array(Float64,(1,))
  res = @msk_ccall( "getnadouinf",Int32,(Ptr{Void},Ptr{Uint8},Ptr{Float64},),task_.task,bytestring(whichdinf_),dvalue_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Float64,dvalue_[1]))
end

function getnadouparam(task_:: MSKtask,paramname_:: String)
  parvalue_ = Array(Float64,(1,))
  res = @msk_ccall( "getnadouparam",Int32,(Ptr{Void},Ptr{Uint8},Ptr{Float64},),task_.task,bytestring(paramname_),parvalue_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Float64,parvalue_[1]))
end

function getnaintinf(task_:: MSKtask,infitemname_:: String)
  ivalue_ = Array(Int32,(1,))
  res = @msk_ccall( "getnaintinf",Int32,(Ptr{Void},Ptr{Uint8},Ptr{Int32},),task_.task,bytestring(infitemname_),ivalue_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,ivalue_[1]))
end

function getnaintparam(task_:: MSKtask,paramname_:: String)
  parvalue_ = Array(Int32,(1,))
  res = @msk_ccall( "getnaintparam",Int32,(Ptr{Void},Ptr{Uint8},Ptr{Int32},),task_.task,bytestring(paramname_),parvalue_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,parvalue_[1]))
end

getnastrparam(task:: MSKtask,paramname:: String,maxlen) = getnastrparam(convert(MSKtask,task),convert(String,paramname),convert(Int32,maxlen))
function getnastrparam(task_:: MSKtask,paramname_:: String,maxlen_:: Int32)
  len_ = Array(Int32,(1,))
  parvalue_ = zeros(Uint8,(maxlen_))
  res = @msk_ccall( "getnastrparam",Int32,(Ptr{Void},Ptr{Uint8},Int32,Ptr{Int32},Ptr{Uint8},),task_.task,bytestring(paramname_),maxlen_,len_,parvalue_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,len_[1]),bytestring(parvalue_))
end

function getnumanz(task_:: MSKtask)
  numanz_ = Array(Int32,(1,))
  res = @msk_ccall( "getnumanz",Int32,(Ptr{Void},Ptr{Int32},),task_.task,numanz_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,numanz_[1]))
end

function getnumanz64(task_:: MSKtask)
  numanz_ = Array(Int64,(1,))
  res = @msk_ccall( "getnumanz64",Int32,(Ptr{Void},Ptr{Int64},),task_.task,numanz_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int64,numanz_[1]))
end

function getnumbarablocktriplets(task_:: MSKtask)
  num_ = Array(Int64,(1,))
  res = @msk_ccall( "getnumbarablocktriplets",Int32,(Ptr{Void},Ptr{Int64},),task_.task,num_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int64,num_[1]))
end

function getnumbaranz(task_:: MSKtask)
  nz_ = Array(Int64,(1,))
  res = @msk_ccall( "getnumbaranz",Int32,(Ptr{Void},Ptr{Int64},),task_.task,nz_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int64,nz_[1]))
end

function getnumbarcblocktriplets(task_:: MSKtask)
  num_ = Array(Int64,(1,))
  res = @msk_ccall( "getnumbarcblocktriplets",Int32,(Ptr{Void},Ptr{Int64},),task_.task,num_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int64,num_[1]))
end

function getnumbarcnz(task_:: MSKtask)
  nz_ = Array(Int64,(1,))
  res = @msk_ccall( "getnumbarcnz",Int32,(Ptr{Void},Ptr{Int64},),task_.task,nz_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int64,nz_[1]))
end

function getnumbarvar(task_:: MSKtask)
  numbarvar_ = Array(Int32,(1,))
  res = @msk_ccall( "getnumbarvar",Int32,(Ptr{Void},Ptr{Int32},),task_.task,numbarvar_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,numbarvar_[1]))
end

function getnumcon(task_:: MSKtask)
  numcon_ = Array(Int32,(1,))
  res = @msk_ccall( "getnumcon",Int32,(Ptr{Void},Ptr{Int32},),task_.task,numcon_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,numcon_[1]))
end

function getnumcone(task_:: MSKtask)
  numcone_ = Array(Int32,(1,))
  res = @msk_ccall( "getnumcone",Int32,(Ptr{Void},Ptr{Int32},),task_.task,numcone_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,numcone_[1]))
end

getnumconemem(task:: MSKtask,k) = getnumconemem(convert(MSKtask,task),convert(Int32,k))
function getnumconemem(task_:: MSKtask,k_:: Int32)
  nummem_ = Array(Int32,(1,))
  res = @msk_ccall( "getnumconemem",Int32,(Ptr{Void},Int32,Ptr{Int32},),task_.task,k_-1,nummem_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,nummem_[1]))
end

function getnumintvar(task_:: MSKtask)
  numintvar_ = Array(Int32,(1,))
  res = @msk_ccall( "getnumintvar",Int32,(Ptr{Void},Ptr{Int32},),task_.task,numintvar_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,numintvar_[1]))
end

function getnumparam(task_:: MSKtask,partype_:: Int32)
  numparam_ = Array(Int32,(1,))
  res = @msk_ccall( "getnumparam",Int32,(Ptr{Void},Int32,Ptr{Int32},),task_.task,partype_,numparam_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,numparam_[1]))
end

getnumqconknz(task:: MSKtask,k) = getnumqconknz(convert(MSKtask,task),convert(Int32,k))
function getnumqconknz(task_:: MSKtask,k_:: Int32)
  numqcnz_ = Array(Int32,(1,))
  res = @msk_ccall( "getnumqconknz",Int32,(Ptr{Void},Int32,Ptr{Int32},),task_.task,k_-1,numqcnz_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,numqcnz_[1]))
end

getnumqconknz64(task:: MSKtask,k) = getnumqconknz64(convert(MSKtask,task),convert(Int32,k))
function getnumqconknz64(task_:: MSKtask,k_:: Int32)
  numqcnz_ = Array(Int64,(1,))
  res = @msk_ccall( "getnumqconknz64",Int32,(Ptr{Void},Int32,Ptr{Int64},),task_.task,k_-1,numqcnz_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int64,numqcnz_[1]))
end

function getnumqobjnz(task_:: MSKtask)
  numqonz_ = Array(Int64,(1,))
  res = @msk_ccall( "getnumqobjnz64",Int32,(Ptr{Void},Ptr{Int64},),task_.task,numqonz_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int64,numqonz_[1]))
end

function getnumsymmat(task_:: MSKtask)
  num_ = Array(Int64,(1,))
  res = @msk_ccall( "getnumsymmat",Int32,(Ptr{Void},Ptr{Int64},),task_.task,num_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int64,num_[1]))
end

function getnumvar(task_:: MSKtask)
  numvar_ = Array(Int32,(1,))
  res = @msk_ccall( "getnumvar",Int32,(Ptr{Void},Ptr{Int32},),task_.task,numvar_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,numvar_[1]))
end

function getobjname(task_:: MSKtask)
  maxlen_ = (1 + getobjnamelen(task_))
  objname_ = zeros(Uint8,(maxlen_))
  res = @msk_ccall( "getobjname",Int32,(Ptr{Void},Int32,Ptr{Uint8},),task_.task,maxlen_,objname_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (bytestring(objname_))
end

function getobjnamelen(task_:: MSKtask)
  len_ = Array(Int32,(1,))
  res = @msk_ccall( "getobjnamelen",Int32,(Ptr{Void},Ptr{Int32},),task_.task,len_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,len_[1]))
end

function getobjsense(task_:: MSKtask)
  sense_ = Array(Int32,(1,))
  res = @msk_ccall( "getobjsense",Int32,(Ptr{Void},Ptr{Int32},),task_.task,sense_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,sense_[1]))
end

getparamname(task:: MSKtask,partype:: Int32,param) = getparamname(convert(MSKtask,task),convert(Int32,partype),convert(Int32,param))
function getparamname(task_:: MSKtask,partype_:: Int32,param_:: Int32)
  parname_ = zeros(Uint8,MSK_MAX_STR_LEN)
  res = @msk_ccall( "getparamname",Int32,(Ptr{Void},Int32,Int32,Ptr{Uint8},),task_.task,partype_,param_,parname_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (bytestring(parname_))
end

getpbi(task:: MSKtask,whichsol:: Int32,accmode:: Int32,sub:: Array,normalize) = getpbi(convert(MSKtask,task),convert(Int32,whichsol),convert(Int32,accmode),convert(Array{Int32},sub),convert(Int32,normalize))
function getpbi(task_:: MSKtask,whichsol_:: Int32,accmode_:: Int32,sub_:: Array{Int32},normalize_:: Int32)
  __tmp_var_0 = if (typeof(sub_) != Array{Int32}) convert(Array{Int32},sub_) else sub_ end
  len_ = minimum([ length(sub_) ])
  __tmp_var_1 = (len_)
  __tmp_var_2 = zeros(Float64,__tmp_var_1)
  res = @msk_ccall( "getpbi",Int32,(Ptr{Void},Int32,Int32,Ptr{Int32},Int32,Ptr{Float64},Int32,),task_.task,whichsol_,accmode_,__tmp_var_0-1,len_,__tmp_var_2,normalize_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (__tmp_var_2)
end

getpcni(task:: MSKtask,whichsol:: Int32,sub:: Array) = getpcni(convert(MSKtask,task),convert(Int32,whichsol),convert(Array{Int32},sub))
function getpcni(task_:: MSKtask,whichsol_:: Int32,sub_:: Array{Int32})
  __tmp_var_0 = if (typeof(sub_) != Array{Int32}) convert(Array{Int32},sub_) else sub_ end
  len_ = minimum([ length(sub_) ])
  __tmp_var_1 = (len_)
  __tmp_var_2 = zeros(Float64,__tmp_var_1)
  res = @msk_ccall( "getpcni",Int32,(Ptr{Void},Int32,Ptr{Int32},Int32,Ptr{Float64},),task_.task,whichsol_,__tmp_var_0-1,len_,__tmp_var_2)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (__tmp_var_2)
end

getpeqi(task:: MSKtask,whichsol:: Int32,sub:: Array,normalize) = getpeqi(convert(MSKtask,task),convert(Int32,whichsol),convert(Array{Int32},sub),convert(Int32,normalize))
function getpeqi(task_:: MSKtask,whichsol_:: Int32,sub_:: Array{Int32},normalize_:: Int32)
  __tmp_var_0 = if (typeof(sub_) != Array{Int32}) convert(Array{Int32},sub_) else sub_ end
  len_ = minimum([ length(sub_) ])
  __tmp_var_1 = (len_)
  __tmp_var_2 = zeros(Float64,__tmp_var_1)
  res = @msk_ccall( "getpeqi",Int32,(Ptr{Void},Int32,Ptr{Int32},Int32,Ptr{Float64},Int32,),task_.task,whichsol_,__tmp_var_0-1,len_,__tmp_var_2,normalize_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (__tmp_var_2)
end

function getprimalobj(task_:: MSKtask,whichsol_:: Int32)
  primalobj_ = Array(Float64,(1,))
  res = @msk_ccall( "getprimalobj",Int32,(Ptr{Void},Int32,Ptr{Float64},),task_.task,whichsol_,primalobj_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Float64,primalobj_[1]))
end

function getprobtype(task_:: MSKtask)
  probtype_ = Array(Int32,(1,))
  res = @msk_ccall( "getprobtype",Int32,(Ptr{Void},Ptr{Int32},),task_.task,probtype_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,probtype_[1]))
end

function getprosta(task_:: MSKtask,whichsol_:: Int32)
  prosta_ = Array(Int32,(1,))
  res = @msk_ccall( "getprosta",Int32,(Ptr{Void},Int32,Ptr{Int32},),task_.task,whichsol_,prosta_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,prosta_[1]))
end

getpviolbarvar(task:: MSKtask,whichsol:: Int32,sub:: Array) = getpviolbarvar(convert(MSKtask,task),convert(Int32,whichsol),convert(Array{Int32},sub))
function getpviolbarvar(task_:: MSKtask,whichsol_:: Int32,sub_:: Array{Int32})
  num_ = minimum([ length(sub_) ])
  __tmp_var_0 = if (typeof(sub_) != Array{Int32}) convert(Array{Int32},sub_) else sub_ end
  __tmp_var_1 = (num_)
  __tmp_var_2 = zeros(Float64,__tmp_var_1)
  res = @msk_ccall( "getpviolbarvar",Int32,(Ptr{Void},Int32,Int32,Ptr{Int32},Ptr{Float64},),task_.task,whichsol_,num_,__tmp_var_0-1,__tmp_var_2)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (__tmp_var_2)
end

getpviolcon(task:: MSKtask,whichsol:: Int32,sub:: Array) = getpviolcon(convert(MSKtask,task),convert(Int32,whichsol),convert(Array{Int32},sub))
function getpviolcon(task_:: MSKtask,whichsol_:: Int32,sub_:: Array{Int32})
  num_ = minimum([ length(sub_) ])
  __tmp_var_0 = if (typeof(sub_) != Array{Int32}) convert(Array{Int32},sub_) else sub_ end
  __tmp_var_1 = (num_)
  __tmp_var_2 = zeros(Float64,__tmp_var_1)
  res = @msk_ccall( "getpviolcon",Int32,(Ptr{Void},Int32,Int32,Ptr{Int32},Ptr{Float64},),task_.task,whichsol_,num_,__tmp_var_0-1,__tmp_var_2)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (__tmp_var_2)
end

getpviolcones(task:: MSKtask,whichsol:: Int32,sub:: Array) = getpviolcones(convert(MSKtask,task),convert(Int32,whichsol),convert(Array{Int32},sub))
function getpviolcones(task_:: MSKtask,whichsol_:: Int32,sub_:: Array{Int32})
  num_ = minimum([ length(sub_) ])
  __tmp_var_0 = if (typeof(sub_) != Array{Int32}) convert(Array{Int32},sub_) else sub_ end
  __tmp_var_1 = (num_)
  __tmp_var_2 = zeros(Float64,__tmp_var_1)
  res = @msk_ccall( "getpviolcones",Int32,(Ptr{Void},Int32,Int32,Ptr{Int32},Ptr{Float64},),task_.task,whichsol_,num_,__tmp_var_0-1,__tmp_var_2)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (__tmp_var_2)
end

getpviolvar(task:: MSKtask,whichsol:: Int32,sub:: Array) = getpviolvar(convert(MSKtask,task),convert(Int32,whichsol),convert(Array{Int32},sub))
function getpviolvar(task_:: MSKtask,whichsol_:: Int32,sub_:: Array{Int32})
  num_ = minimum([ length(sub_) ])
  __tmp_var_0 = if (typeof(sub_) != Array{Int32}) convert(Array{Int32},sub_) else sub_ end
  __tmp_var_1 = (num_)
  __tmp_var_2 = zeros(Float64,__tmp_var_1)
  res = @msk_ccall( "getpviolvar",Int32,(Ptr{Void},Int32,Int32,Ptr{Int32},Ptr{Float64},),task_.task,whichsol_,num_,__tmp_var_0-1,__tmp_var_2)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (__tmp_var_2)
end

getqconk(task:: MSKtask,k) = getqconk(convert(MSKtask,task),convert(Int32,k))
function getqconk(task_:: MSKtask,k_:: Int32)
  maxnumqcnz_ = getnumqconknz(task_,(k_))
  qcsurp_ = convert(Int64,length(qcsubi_))
  numqcnz_ = Array(Int64,(1,))
  __tmp_var_0 = (maxnumqcnz_)
  __tmp_var_1 = zeros(Int32,__tmp_var_0)
  __tmp_var_2 = (maxnumqcnz_)
  __tmp_var_3 = zeros(Int32,__tmp_var_2)
  __tmp_var_4 = (maxnumqcnz_)
  __tmp_var_5 = zeros(Float64,__tmp_var_4)
  res = @msk_ccall( "getqconk64",Int32,(Ptr{Void},Int32,Int64,Ptr{Int64},Ptr{Int64},Ptr{Int32},Ptr{Int32},Ptr{Float64},),task_.task,k_-1,maxnumqcnz_,&qcsurp_,numqcnz_,__tmp_var_1,__tmp_var_3,__tmp_var_5)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int64,numqcnz_[1]),__tmp_var_1,__tmp_var_3,__tmp_var_5)
end

function getqobj(task_:: MSKtask)
  maxnumqonz_ = getnumqobjnz(task_)
  qosurp_ = convert(Int32,length(qosubi_))
  numqonz_ = Array(Int32,(1,))
  __tmp_var_0 = (maxnumqonz_)
  __tmp_var_1 = zeros(Int32,__tmp_var_0)
  __tmp_var_2 = (maxnumqonz_)
  __tmp_var_3 = zeros(Int32,__tmp_var_2)
  __tmp_var_4 = (maxnumqonz_)
  __tmp_var_5 = zeros(Float64,__tmp_var_4)
  res = @msk_ccall( "getqobj",Int32,(Ptr{Void},Int32,Ptr{Int32},Ptr{Int32},Ptr{Int32},Ptr{Int32},Ptr{Float64},),task_.task,maxnumqonz_,&qosurp_,numqonz_,__tmp_var_1,__tmp_var_3,__tmp_var_5)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,numqonz_[1]),__tmp_var_1,__tmp_var_3,__tmp_var_5)
end

function getqobj64(task_:: MSKtask)
  maxnumqonz_ = getnumqobjnz(task_)
  qosurp_ = convert(Int64,length(qosubi_))
  numqonz_ = Array(Int64,(1,))
  __tmp_var_0 = (maxnumqonz_)
  __tmp_var_1 = zeros(Int32,__tmp_var_0)
  __tmp_var_2 = (maxnumqonz_)
  __tmp_var_3 = zeros(Int32,__tmp_var_2)
  __tmp_var_4 = (maxnumqonz_)
  __tmp_var_5 = zeros(Float64,__tmp_var_4)
  res = @msk_ccall( "getqobj64",Int32,(Ptr{Void},Int64,Ptr{Int64},Ptr{Int64},Ptr{Int32},Ptr{Int32},Ptr{Float64},),task_.task,maxnumqonz_,&qosurp_,numqonz_,__tmp_var_1,__tmp_var_3,__tmp_var_5)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int64,numqonz_[1]),__tmp_var_1,__tmp_var_3,__tmp_var_5)
end

getqobjij(task:: MSKtask,i,j) = getqobjij(convert(MSKtask,task),convert(Int32,i),convert(Int32,j))
function getqobjij(task_:: MSKtask,i_:: Int32,j_:: Int32)
  qoij_ = Array(Float64,(1,))
  res = @msk_ccall( "getqobjij",Int32,(Ptr{Void},Int32,Int32,Ptr{Float64},),task_.task,i_-1,j_-1,qoij_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Float64,qoij_[1]))
end

getreducedcosts(task:: MSKtask,whichsol:: Int32,first,last) = getreducedcosts(convert(MSKtask,task),convert(Int32,whichsol),convert(Int32,first),convert(Int32,last))
function getreducedcosts(task_:: MSKtask,whichsol_:: Int32,first_:: Int32,last_:: Int32)
  __tmp_var_0 = ((last_) - (first_))
  __tmp_var_1 = zeros(Float64,__tmp_var_0)
  res = @msk_ccall( "getreducedcosts",Int32,(Ptr{Void},Int32,Int32,Int32,Ptr{Float64},),task_.task,whichsol_,first_-1,last_-1,__tmp_var_1)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (__tmp_var_1)
end

function getskc(task_:: MSKtask,whichsol_:: Int32)
  __tmp_var_0 = getnumcon(task_)
  skc_ = zeros(Int32,__tmp_var_0)
  res = @msk_ccall( "getskc",Int32,(Ptr{Void},Int32,Ptr{Int32},),task_.task,whichsol_,skc_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (skc_)
end

getskcslice(task:: MSKtask,whichsol:: Int32,first,last) = getskcslice(convert(MSKtask,task),convert(Int32,whichsol),convert(Int32,first),convert(Int32,last))
function getskcslice(task_:: MSKtask,whichsol_:: Int32,first_:: Int32,last_:: Int32)
  __tmp_var_0 = ((last_) - (first_))
  skc_ = zeros(Int32,__tmp_var_0)
  res = @msk_ccall( "getskcslice",Int32,(Ptr{Void},Int32,Int32,Int32,Ptr{Int32},),task_.task,whichsol_,first_-1,last_-1,skc_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (skc_)
end

function getskx(task_:: MSKtask,whichsol_:: Int32)
  __tmp_var_0 = getnumvar(task_)
  skx_ = zeros(Int32,__tmp_var_0)
  res = @msk_ccall( "getskx",Int32,(Ptr{Void},Int32,Ptr{Int32},),task_.task,whichsol_,skx_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (skx_)
end

getskxslice(task:: MSKtask,whichsol:: Int32,first,last) = getskxslice(convert(MSKtask,task),convert(Int32,whichsol),convert(Int32,first),convert(Int32,last))
function getskxslice(task_:: MSKtask,whichsol_:: Int32,first_:: Int32,last_:: Int32)
  __tmp_var_0 = ((last_) - (first_))
  skx_ = zeros(Int32,__tmp_var_0)
  res = @msk_ccall( "getskxslice",Int32,(Ptr{Void},Int32,Int32,Int32,Ptr{Int32},),task_.task,whichsol_,first_-1,last_-1,skx_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (skx_)
end

function getslc(task_:: MSKtask,whichsol_:: Int32)
  __tmp_var_0 = getnumcon(task_)
  __tmp_var_1 = zeros(Float64,__tmp_var_0)
  res = @msk_ccall( "getslc",Int32,(Ptr{Void},Int32,Ptr{Float64},),task_.task,whichsol_,__tmp_var_1)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (__tmp_var_1)
end

getslcslice(task:: MSKtask,whichsol:: Int32,first,last) = getslcslice(convert(MSKtask,task),convert(Int32,whichsol),convert(Int32,first),convert(Int32,last))
function getslcslice(task_:: MSKtask,whichsol_:: Int32,first_:: Int32,last_:: Int32)
  __tmp_var_0 = ((last_) - (first_))
  __tmp_var_1 = zeros(Float64,__tmp_var_0)
  res = @msk_ccall( "getslcslice",Int32,(Ptr{Void},Int32,Int32,Int32,Ptr{Float64},),task_.task,whichsol_,first_-1,last_-1,__tmp_var_1)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (__tmp_var_1)
end

function getslx(task_:: MSKtask,whichsol_:: Int32)
  __tmp_var_0 = getnumvar(task_)
  __tmp_var_1 = zeros(Float64,__tmp_var_0)
  res = @msk_ccall( "getslx",Int32,(Ptr{Void},Int32,Ptr{Float64},),task_.task,whichsol_,__tmp_var_1)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (__tmp_var_1)
end

getslxslice(task:: MSKtask,whichsol:: Int32,first,last) = getslxslice(convert(MSKtask,task),convert(Int32,whichsol),convert(Int32,first),convert(Int32,last))
function getslxslice(task_:: MSKtask,whichsol_:: Int32,first_:: Int32,last_:: Int32)
  __tmp_var_0 = ((last_) - (first_))
  __tmp_var_1 = zeros(Float64,__tmp_var_0)
  res = @msk_ccall( "getslxslice",Int32,(Ptr{Void},Int32,Int32,Int32,Ptr{Float64},),task_.task,whichsol_,first_-1,last_-1,__tmp_var_1)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (__tmp_var_1)
end

function getsnx(task_:: MSKtask,whichsol_:: Int32)
  __tmp_var_0 = getnumvar(task_)
  __tmp_var_1 = zeros(Float64,__tmp_var_0)
  res = @msk_ccall( "getsnx",Int32,(Ptr{Void},Int32,Ptr{Float64},),task_.task,whichsol_,__tmp_var_1)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (__tmp_var_1)
end

getsnxslice(task:: MSKtask,whichsol:: Int32,first,last) = getsnxslice(convert(MSKtask,task),convert(Int32,whichsol),convert(Int32,first),convert(Int32,last))
function getsnxslice(task_:: MSKtask,whichsol_:: Int32,first_:: Int32,last_:: Int32)
  __tmp_var_0 = ((last_) - (first_))
  __tmp_var_1 = zeros(Float64,__tmp_var_0)
  res = @msk_ccall( "getsnxslice",Int32,(Ptr{Void},Int32,Int32,Int32,Ptr{Float64},),task_.task,whichsol_,first_-1,last_-1,__tmp_var_1)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (__tmp_var_1)
end

function getsolsta(task_:: MSKtask,whichsol_:: Int32)
  solsta_ = Array(Int32,(1,))
  res = @msk_ccall( "getsolsta",Int32,(Ptr{Void},Int32,Ptr{Int32},),task_.task,whichsol_,solsta_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,solsta_[1]))
end

function getsolution(task_:: MSKtask,whichsol_:: Int32)
  prosta_ = Array(Int32,(1,))
  solsta_ = Array(Int32,(1,))
  __tmp_var_0 = getnumcon(task_)
  skc_ = zeros(Int32,__tmp_var_0)
  __tmp_var_1 = getnumvar(task_)
  skx_ = zeros(Int32,__tmp_var_1)
  __tmp_var_2 = getnumcone(task_)
  skn_ = zeros(Int32,__tmp_var_2)
  __tmp_var_3 = getnumcon(task_)
  __tmp_var_4 = zeros(Float64,__tmp_var_3)
  __tmp_var_5 = getnumvar(task_)
  __tmp_var_6 = zeros(Float64,__tmp_var_5)
  __tmp_var_7 = getnumcon(task_)
  __tmp_var_8 = zeros(Float64,__tmp_var_7)
  __tmp_var_9 = getnumcon(task_)
  __tmp_var_10 = zeros(Float64,__tmp_var_9)
  __tmp_var_11 = getnumcon(task_)
  __tmp_var_12 = zeros(Float64,__tmp_var_11)
  __tmp_var_13 = getnumvar(task_)
  __tmp_var_14 = zeros(Float64,__tmp_var_13)
  __tmp_var_15 = getnumvar(task_)
  __tmp_var_16 = zeros(Float64,__tmp_var_15)
  __tmp_var_17 = getnumcone(task_)
  __tmp_var_18 = zeros(Float64,__tmp_var_17)
  res = @msk_ccall( "getsolution",Int32,(Ptr{Void},Int32,Ptr{Int32},Ptr{Int32},Ptr{Int32},Ptr{Int32},Ptr{Int32},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},),task_.task,whichsol_,prosta_,solsta_,skc_,skx_,skn_,__tmp_var_4,__tmp_var_6,__tmp_var_8,__tmp_var_10,__tmp_var_12,__tmp_var_14,__tmp_var_16,__tmp_var_18)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,prosta_[1]),convert(Int32,solsta_[1]),skc_,skx_,skn_,__tmp_var_4,__tmp_var_6,__tmp_var_8,__tmp_var_10,__tmp_var_12,__tmp_var_14,__tmp_var_16,__tmp_var_18)
end

getsolutioni(task:: MSKtask,accmode:: Int32,i,whichsol:: Int32) = getsolutioni(convert(MSKtask,task),convert(Int32,accmode),convert(Int32,i),convert(Int32,whichsol))
function getsolutioni(task_:: MSKtask,accmode_:: Int32,i_:: Int32,whichsol_:: Int32)
  sk_ = Array(Int32,(1,))
  x_ = Array(Float64,(1,))
  sl_ = Array(Float64,(1,))
  su_ = Array(Float64,(1,))
  sn_ = Array(Float64,(1,))
  res = @msk_ccall( "getsolutioni",Int32,(Ptr{Void},Int32,Int32,Int32,Ptr{Int32},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},),task_.task,accmode_,i_-1,whichsol_,sk_,x_,sl_,su_,sn_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,sk_[1]),convert(Float64,x_[1]),convert(Float64,sl_[1]),convert(Float64,su_[1]),convert(Float64,sn_[1]))
end

function getsolutioninf(task_:: MSKtask,whichsol_:: Int32)
  prosta_ = Array(Int32,(1,))
  solsta_ = Array(Int32,(1,))
  primalobj_ = Array(Float64,(1,))
  maxpbi_ = Array(Float64,(1,))
  maxpcni_ = Array(Float64,(1,))
  maxpeqi_ = Array(Float64,(1,))
  maxinti_ = Array(Float64,(1,))
  dualobj_ = Array(Float64,(1,))
  maxdbi_ = Array(Float64,(1,))
  maxdcni_ = Array(Float64,(1,))
  maxdeqi_ = Array(Float64,(1,))
  res = @msk_ccall( "getsolutioninf",Int32,(Ptr{Void},Int32,Ptr{Int32},Ptr{Int32},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},),task_.task,whichsol_,prosta_,solsta_,primalobj_,maxpbi_,maxpcni_,maxpeqi_,maxinti_,dualobj_,maxdbi_,maxdcni_,maxdeqi_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,prosta_[1]),convert(Int32,solsta_[1]),convert(Float64,primalobj_[1]),convert(Float64,maxpbi_[1]),convert(Float64,maxpcni_[1]),convert(Float64,maxpeqi_[1]),convert(Float64,maxinti_[1]),convert(Float64,dualobj_[1]),convert(Float64,maxdbi_[1]),convert(Float64,maxdcni_[1]),convert(Float64,maxdeqi_[1]))
end

function getsolutioninfo(task_:: MSKtask,whichsol_:: Int32)
  pobj_ = Array(Float64,(1,))
  pviolcon_ = Array(Float64,(1,))
  pviolvar_ = Array(Float64,(1,))
  pviolbarvar_ = Array(Float64,(1,))
  pviolcone_ = Array(Float64,(1,))
  pviolitg_ = Array(Float64,(1,))
  dobj_ = Array(Float64,(1,))
  dviolcon_ = Array(Float64,(1,))
  dviolvar_ = Array(Float64,(1,))
  dviolbarvar_ = Array(Float64,(1,))
  dviolcones_ = Array(Float64,(1,))
  res = @msk_ccall( "getsolutioninfo",Int32,(Ptr{Void},Int32,Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},),task_.task,whichsol_,pobj_,pviolcon_,pviolvar_,pviolbarvar_,pviolcone_,pviolitg_,dobj_,dviolcon_,dviolvar_,dviolbarvar_,dviolcones_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Float64,pobj_[1]),convert(Float64,pviolcon_[1]),convert(Float64,pviolvar_[1]),convert(Float64,pviolbarvar_[1]),convert(Float64,pviolcone_[1]),convert(Float64,pviolitg_[1]),convert(Float64,dobj_[1]),convert(Float64,dviolcon_[1]),convert(Float64,dviolvar_[1]),convert(Float64,dviolbarvar_[1]),convert(Float64,dviolcones_[1]))
end

getsolutionslice(task:: MSKtask,whichsol:: Int32,solitem:: Int32,first,last) = getsolutionslice(convert(MSKtask,task),convert(Int32,whichsol),convert(Int32,solitem),convert(Int32,first),convert(Int32,last))
function getsolutionslice(task_:: MSKtask,whichsol_:: Int32,solitem_:: Int32,first_:: Int32,last_:: Int32)
  __tmp_var_0 = ((last_) - (first_))
  __tmp_var_1 = zeros(Float64,__tmp_var_0)
  res = @msk_ccall( "getsolutionslice",Int32,(Ptr{Void},Int32,Int32,Int32,Int32,Ptr{Float64},),task_.task,whichsol_,solitem_,first_-1,last_-1,__tmp_var_1)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (__tmp_var_1)
end

getsparsesymmat(task:: MSKtask,idx) = getsparsesymmat(convert(MSKtask,task),convert(Int64,idx))
function getsparsesymmat(task_:: MSKtask,idx_:: Int64)
  maxlen_ = getsymmatinfo(task_,(idx_))[1]
  __tmp_var_0 = (maxlen_)
  __tmp_var_1 = zeros(Int32,__tmp_var_0)
  __tmp_var_2 = (maxlen_)
  __tmp_var_3 = zeros(Int32,__tmp_var_2)
  __tmp_var_4 = (maxlen_)
  __tmp_var_5 = zeros(Float64,__tmp_var_4)
  res = @msk_ccall( "getsparsesymmat",Int32,(Ptr{Void},Int64,Int64,Ptr{Int32},Ptr{Int32},Ptr{Float64},),task_.task,idx_-1,maxlen_,__tmp_var_1,__tmp_var_3,__tmp_var_5)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (__tmp_var_1,__tmp_var_3,__tmp_var_5)
end

function getstrparam(task_:: MSKtask,param_:: Int32)
  maxlen_ = (1 + getstrparamlen(task_,(param_)))
  len_ = Array(Int32,(1,))
  parvalue_ = zeros(Uint8,(maxlen_))
  res = @msk_ccall( "getstrparam",Int32,(Ptr{Void},Int32,Int32,Ptr{Int32},Ptr{Uint8},),task_.task,param_,maxlen_,len_,parvalue_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,len_[1]),bytestring(parvalue_))
end

function getstrparamlen(task_:: MSKtask,param_:: Int32)
  len_ = Array(Int32,(1,))
  res = @msk_ccall( "getstrparamlen",Int32,(Ptr{Void},Int32,Ptr{Int32},),task_.task,param_,len_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,len_[1]))
end

function getsuc(task_:: MSKtask,whichsol_:: Int32)
  __tmp_var_0 = getnumcon(task_)
  __tmp_var_1 = zeros(Float64,__tmp_var_0)
  res = @msk_ccall( "getsuc",Int32,(Ptr{Void},Int32,Ptr{Float64},),task_.task,whichsol_,__tmp_var_1)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (__tmp_var_1)
end

getsucslice(task:: MSKtask,whichsol:: Int32,first,last) = getsucslice(convert(MSKtask,task),convert(Int32,whichsol),convert(Int32,first),convert(Int32,last))
function getsucslice(task_:: MSKtask,whichsol_:: Int32,first_:: Int32,last_:: Int32)
  __tmp_var_0 = ((last_) - (first_))
  __tmp_var_1 = zeros(Float64,__tmp_var_0)
  res = @msk_ccall( "getsucslice",Int32,(Ptr{Void},Int32,Int32,Int32,Ptr{Float64},),task_.task,whichsol_,first_-1,last_-1,__tmp_var_1)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (__tmp_var_1)
end

function getsux(task_:: MSKtask,whichsol_:: Int32)
  __tmp_var_0 = getnumvar(task_)
  __tmp_var_1 = zeros(Float64,__tmp_var_0)
  res = @msk_ccall( "getsux",Int32,(Ptr{Void},Int32,Ptr{Float64},),task_.task,whichsol_,__tmp_var_1)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (__tmp_var_1)
end

getsuxslice(task:: MSKtask,whichsol:: Int32,first,last) = getsuxslice(convert(MSKtask,task),convert(Int32,whichsol),convert(Int32,first),convert(Int32,last))
function getsuxslice(task_:: MSKtask,whichsol_:: Int32,first_:: Int32,last_:: Int32)
  __tmp_var_0 = ((last_) - (first_))
  __tmp_var_1 = zeros(Float64,__tmp_var_0)
  res = @msk_ccall( "getsuxslice",Int32,(Ptr{Void},Int32,Int32,Int32,Ptr{Float64},),task_.task,whichsol_,first_-1,last_-1,__tmp_var_1)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (__tmp_var_1)
end

getsymmatinfo(task:: MSKtask,idx) = getsymmatinfo(convert(MSKtask,task),convert(Int64,idx))
function getsymmatinfo(task_:: MSKtask,idx_:: Int64)
  dim_ = Array(Int32,(1,))
  nz_ = Array(Int64,(1,))
  type_ = Array(Int32,(1,))
  res = @msk_ccall( "getsymmatinfo",Int32,(Ptr{Void},Int64,Ptr{Int32},Ptr{Int64},Ptr{Int32},),task_.task,idx_-1,dim_,nz_,type_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,dim_[1]),convert(Int64,nz_[1]),convert(Int32,type_[1]))
end

function gettaskname(task_:: MSKtask)
  maxlen_ = (1 + gettasknamelen(task_))
  taskname_ = zeros(Uint8,(maxlen_))
  res = @msk_ccall( "gettaskname",Int32,(Ptr{Void},Int32,Ptr{Uint8},),task_.task,maxlen_,taskname_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (bytestring(taskname_))
end

function gettasknamelen(task_:: MSKtask)
  len_ = Array(Int32,(1,))
  res = @msk_ccall( "gettasknamelen",Int32,(Ptr{Void},Ptr{Int32},),task_.task,len_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,len_[1]))
end

getvarbound(task:: MSKtask,i) = getvarbound(convert(MSKtask,task),convert(Int32,i))
function getvarbound(task_:: MSKtask,i_:: Int32)
  bk_ = Array(Int32,(1,))
  bl_ = Array(Float64,(1,))
  bu_ = Array(Float64,(1,))
  res = @msk_ccall( "getvarbound",Int32,(Ptr{Void},Int32,Ptr{Int32},Ptr{Float64},Ptr{Float64},),task_.task,i_-1,bk_,bl_,bu_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,bk_[1]),convert(Float64,bl_[1]),convert(Float64,bu_[1]))
end

getvarboundslice(task:: MSKtask,first,last) = getvarboundslice(convert(MSKtask,task),convert(Int32,first),convert(Int32,last))
function getvarboundslice(task_:: MSKtask,first_:: Int32,last_:: Int32)
  __tmp_var_0 = ((last_) - (first_))
  bk_ = zeros(Int32,__tmp_var_0)
  __tmp_var_1 = ((last_) - (first_))
  __tmp_var_2 = zeros(Float64,__tmp_var_1)
  __tmp_var_3 = ((last_) - (first_))
  __tmp_var_4 = zeros(Float64,__tmp_var_3)
  res = @msk_ccall( "getvarboundslice",Int32,(Ptr{Void},Int32,Int32,Ptr{Int32},Ptr{Float64},Ptr{Float64},),task_.task,first_-1,last_-1,bk_,__tmp_var_2,__tmp_var_4)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (bk_,__tmp_var_2,__tmp_var_4)
end

getvarbranchdir(task:: MSKtask,j) = getvarbranchdir(convert(MSKtask,task),convert(Int32,j))
function getvarbranchdir(task_:: MSKtask,j_:: Int32)
  direction_ = Array(Int32,(1,))
  res = @msk_ccall( "getvarbranchdir",Int32,(Ptr{Void},Int32,Ptr{Int32},),task_.task,j_-1,direction_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,direction_[1]))
end

getvarbranchpri(task:: MSKtask,j) = getvarbranchpri(convert(MSKtask,task),convert(Int32,j))
function getvarbranchpri(task_:: MSKtask,j_:: Int32)
  priority_ = Array(Int32,(1,))
  res = @msk_ccall( "getvarbranchpri",Int32,(Ptr{Void},Int32,Ptr{Int32},),task_.task,j_-1,priority_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,priority_[1]))
end

getvarname(task:: MSKtask,j) = getvarname(convert(MSKtask,task),convert(Int32,j))
function getvarname(task_:: MSKtask,j_:: Int32)
  maxlen_ = (1 + getvarnamelen(task_,(j_)))
  name_ = zeros(Uint8,(maxlen_))
  res = @msk_ccall( "getvarname",Int32,(Ptr{Void},Int32,Int32,Ptr{Uint8},),task_.task,j_-1,maxlen_,name_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (bytestring(name_))
end

function getvarnameindex(task_:: MSKtask,somename_:: String)
  asgn_ = Array(Int32,(1,))
  index_ = Array(Int32,(1,))
  res = @msk_ccall( "getvarnameindex",Int32,(Ptr{Void},Ptr{Uint8},Ptr{Int32},Ptr{Int32},),task_.task,bytestring(somename_),asgn_,index_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,asgn_[1]),convert(Int32,index_[1]))
end

getvarnamelen(task:: MSKtask,i) = getvarnamelen(convert(MSKtask,task),convert(Int32,i))
function getvarnamelen(task_:: MSKtask,i_:: Int32)
  len_ = Array(Int32,(1,))
  res = @msk_ccall( "getvarnamelen",Int32,(Ptr{Void},Int32,Ptr{Int32},),task_.task,i_-1,len_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,len_[1]))
end

getvartype(task:: MSKtask,j) = getvartype(convert(MSKtask,task),convert(Int32,j))
function getvartype(task_:: MSKtask,j_:: Int32)
  vartype_ = Array(Int32,(1,))
  res = @msk_ccall( "getvartype",Int32,(Ptr{Void},Int32,Ptr{Int32},),task_.task,j_-1,vartype_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,vartype_[1]))
end

getvartypelist(task:: MSKtask,subj:: Array) = getvartypelist(convert(MSKtask,task),convert(Array{Int32},subj))
function getvartypelist(task_:: MSKtask,subj_:: Array{Int32})
  num_ = minimum([ length(subj_) ])
  __tmp_var_0 = if (typeof(subj_) != Array{Int32}) convert(Array{Int32},subj_) else subj_ end
  __tmp_var_1 = (num_)
  vartype_ = zeros(Int32,__tmp_var_1)
  res = @msk_ccall( "getvartypelist",Int32,(Ptr{Void},Int32,Ptr{Int32},Ptr{Int32},),task_.task,num_,__tmp_var_0-1,vartype_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (vartype_)
end

function getxc(task_:: MSKtask,whichsol_:: Int32)
  __tmp_var_0 = getnumcon(task_)
  __tmp_var_1 = zeros(Float64,__tmp_var_0)
  res = @msk_ccall( "getxc",Int32,(Ptr{Void},Int32,Ptr{Float64},),task_.task,whichsol_,__tmp_var_1)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (__tmp_var_1)
end

getxcslice(task:: MSKtask,whichsol:: Int32,first,last) = getxcslice(convert(MSKtask,task),convert(Int32,whichsol),convert(Int32,first),convert(Int32,last))
function getxcslice(task_:: MSKtask,whichsol_:: Int32,first_:: Int32,last_:: Int32)
  __tmp_var_0 = ((last_) - (first_))
  __tmp_var_1 = zeros(Float64,__tmp_var_0)
  res = @msk_ccall( "getxcslice",Int32,(Ptr{Void},Int32,Int32,Int32,Ptr{Float64},),task_.task,whichsol_,first_-1,last_-1,__tmp_var_1)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (__tmp_var_1)
end

function getxx(task_:: MSKtask,whichsol_:: Int32)
  __tmp_var_0 = getnumvar(task_)
  __tmp_var_1 = zeros(Float64,__tmp_var_0)
  res = @msk_ccall( "getxx",Int32,(Ptr{Void},Int32,Ptr{Float64},),task_.task,whichsol_,__tmp_var_1)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (__tmp_var_1)
end

getxxslice(task:: MSKtask,whichsol:: Int32,first,last) = getxxslice(convert(MSKtask,task),convert(Int32,whichsol),convert(Int32,first),convert(Int32,last))
function getxxslice(task_:: MSKtask,whichsol_:: Int32,first_:: Int32,last_:: Int32)
  __tmp_var_0 = ((last_) - (first_))
  __tmp_var_1 = zeros(Float64,__tmp_var_0)
  res = @msk_ccall( "getxxslice",Int32,(Ptr{Void},Int32,Int32,Int32,Ptr{Float64},),task_.task,whichsol_,first_-1,last_-1,__tmp_var_1)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (__tmp_var_1)
end

function gety(task_:: MSKtask,whichsol_:: Int32)
  __tmp_var_0 = getnumcon(task_)
  __tmp_var_1 = zeros(Float64,__tmp_var_0)
  res = @msk_ccall( "gety",Int32,(Ptr{Void},Int32,Ptr{Float64},),task_.task,whichsol_,__tmp_var_1)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (__tmp_var_1)
end

getyslice(task:: MSKtask,whichsol:: Int32,first,last) = getyslice(convert(MSKtask,task),convert(Int32,whichsol),convert(Int32,first),convert(Int32,last))
function getyslice(task_:: MSKtask,whichsol_:: Int32,first_:: Int32,last_:: Int32)
  __tmp_var_0 = ((last_) - (first_))
  __tmp_var_1 = zeros(Float64,__tmp_var_0)
  res = @msk_ccall( "getyslice",Int32,(Ptr{Void},Int32,Int32,Int32,Ptr{Float64},),task_.task,whichsol_,first_-1,last_-1,__tmp_var_1)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (__tmp_var_1)
end

function initbasissolve(task_:: MSKtask)
  __tmp_var_0 = getnumcon(task_)
  __tmp_var_1 = zeros(Int32,__tmp_var_0)
  res = @msk_ccall( "initbasissolve",Int32,(Ptr{Void},Ptr{Int32},),task_.task,__tmp_var_1)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (__tmp_var_1)
end

inputdata(task:: MSKtask,maxnumcon,maxnumvar,c:: Array{Float64},cfix:: Float64,aptrb:: Array,aptre:: Array,asub:: Array,aval:: Array{Float64},bkc:: Array{Int32},blc:: Array{Float64},buc:: Array{Float64},bkx:: Array{Int32},blx:: Array{Float64},bux:: Array{Float64}) = inputdata(convert(MSKtask,task),convert(Int32,maxnumcon),convert(Int32,maxnumvar),convert(Array{Float64},c),convert(Float64,cfix),convert(Array{Int64},aptrb),convert(Array{Int64},aptre),convert(Array{Int32},asub),convert(Array{Float64},aval),convert(Array{Int32},bkc),convert(Array{Float64},blc),convert(Array{Float64},buc),convert(Array{Int32},bkx),convert(Array{Float64},blx),convert(Array{Float64},bux))
function inputdata(task:: MSKtask,maxnumcon,maxnumvar,c:: Array{Float64},cfix:: Float64,A:: SparseMatrixCSC{Float64},bkc:: Array{Int32},blc:: Array{Float64},buc:: Array{Float64},bkx:: Array{Int32},blx:: Array{Float64},bux:: Array{Float64})
  aptrb = A.colptr[1:size(A,2)-1]
  aptre = A.colptr[2:size(A,2)]
  asub = A.rowval
  aval = A.nzval
  inputdata(task,maxnumcon,maxnumvar,c,cfix,aptrb,aptre,asub,aval,bkc,blc,buc,bkx,blx,bux)
end
function inputdata(task_:: MSKtask,maxnumcon_:: Int32,maxnumvar_:: Int32,c_:: Array{Float64},cfix_:: Float64,aptrb_:: Array{Int64},aptre_:: Array{Int64},asub_:: Array{Int32},aval_:: Array{Float64},bkc_:: Array{Int32},blc_:: Array{Float64},buc_:: Array{Float64},bkx_:: Array{Int32},blx_:: Array{Float64},bux_:: Array{Float64})
  numcon_ = minimum([ length(buc_),length(blc_),length(bkc_) ])
  numvar_ = minimum([ length(c_),length(bux_),length(blx_),length(bkx_),length(aptrb_),length(aptre_) ])
  __tmp_var_0 = if (typeof(asub_) != Array{Int32}) convert(Array{Int32},asub_) else asub_ end
  res = @msk_ccall( "inputdata64",Int32,(Ptr{Void},Int32,Int32,Int32,Int32,Ptr{Float64},Float64,Ptr{Int64},Ptr{Int64},Ptr{Int32},Ptr{Float64},Ptr{Int32},Ptr{Float64},Ptr{Float64},Ptr{Int32},Ptr{Float64},Ptr{Float64},),task_.task,maxnumcon_,maxnumvar_,numcon_,numvar_,c_,cfix_,aptrb_-1,aptre_-1,__tmp_var_0-1,aval_,bkc_,blc_,buc_,bkx_,blx_,bux_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

function isdouparname(task_:: MSKtask,parname_:: String)
  param_ = Array(Int32,(1,))
  res = @msk_ccall( "isdouparname",Int32,(Ptr{Void},Ptr{Uint8},Ptr{Int32},),task_.task,bytestring(parname_),param_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,param_[1]))
end

function isintparname(task_:: MSKtask,parname_:: String)
  param_ = Array(Int32,(1,))
  res = @msk_ccall( "isintparname",Int32,(Ptr{Void},Ptr{Uint8},Ptr{Int32},),task_.task,bytestring(parname_),param_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,param_[1]))
end

function isstrparname(task_:: MSKtask,parname_:: String)
  param_ = Array(Int32,(1,))
  res = @msk_ccall( "isstrparname",Int32,(Ptr{Void},Ptr{Uint8},Ptr{Int32},),task_.task,bytestring(parname_),param_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,param_[1]))
end

linkfiletostream(task:: MSKtask,whichstream:: Int32,filename:: String,append) = linkfiletostream(convert(MSKtask,task),convert(Int32,whichstream),convert(String,filename),convert(Int32,append))
function linkfiletostream(task_:: MSKtask,whichstream_:: Int32,filename_:: String,append_:: Int32)
  res = @msk_ccall( "linkfiletotaskstream",Int32,(Ptr{Void},Int32,Ptr{Uint8},Int32,),task_.task,whichstream_,bytestring(filename_),append_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

function onesolutionsummary(task_:: MSKtask,whichstream_:: Int32,whichsol_:: Int32)
  res = @msk_ccall( "onesolutionsummary",Int32,(Ptr{Void},Int32,Int32,),task_.task,whichstream_,whichsol_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

function optimizeconcurrent(task_:: MSKtask,taskarray_:: Array{MSKtask})
  num_ = minimum([ length(taskarray_) ])
  __tmp_var_0 = (num_)
  if length(taskarray_) < __tmp_var_0
    println("Array argument taskarray is not long enough")
    throw(BoundsError())
  end
  _taskarray_tmp_:: Ptr{Void} = [ __tmp_var_1.task for __tmp_var_1 = taskarray_ ]
  res = @msk_ccall( "optimizeconcurrent2",Int32,(Ptr{Void},Int32,Ptr{Ptr{Void}},),task_.task,num_,_taskarray_tmp)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

function optimizersummary(task_:: MSKtask,whichstream_:: Int32)
  res = @msk_ccall( "optimizersummary",Int32,(Ptr{Void},Int32,),task_.task,whichstream_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

function optimize(task_:: MSKtask)
  trmcode_ = Array(Int32,(1,))
  res = @msk_ccall( "optimizetrm",Int32,(Ptr{Void},Ptr{Int32},),task_.task,trmcode_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,trmcode_[1]))
end

function primalrepair(task_:: MSKtask,wlc_:: Array{Float64},wuc_:: Array{Float64},wlx_:: Array{Float64},wux_:: Array{Float64})
  __tmp_var_0 = getnumcon(task_)
  if length(wlc_) < __tmp_var_0
    println("Array argument wlc is not long enough")
    throw(BoundsError())
  end
  __tmp_var_1 = getnumcon(task_)
  if length(wuc_) < __tmp_var_1
    println("Array argument wuc is not long enough")
    throw(BoundsError())
  end
  __tmp_var_2 = getnumvar(task_)
  if length(wlx_) < __tmp_var_2
    println("Array argument wlx is not long enough")
    throw(BoundsError())
  end
  __tmp_var_3 = getnumvar(task_)
  if length(wux_) < __tmp_var_3
    println("Array argument wux is not long enough")
    throw(BoundsError())
  end
  res = @msk_ccall( "primalrepair",Int32,(Ptr{Void},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},),task_.task,wlc_,wuc_,wlx_,wux_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

primalsensitivity(task:: MSKtask,subi:: Array,marki:: Array{Int32},subj:: Array,markj:: Array{Int32}) = primalsensitivity(convert(MSKtask,task),convert(Array{Int32},subi),convert(Array{Int32},marki),convert(Array{Int32},subj),convert(Array{Int32},markj))
function primalsensitivity(task_:: MSKtask,subi_:: Array{Int32},marki_:: Array{Int32},subj_:: Array{Int32},markj_:: Array{Int32})
  numi_ = minimum([ length(subi_),length(marki_) ])
  __tmp_var_0 = if (typeof(subi_) != Array{Int32}) convert(Array{Int32},subi_) else subi_ end
  numj_ = minimum([ length(subj_),length(markj_) ])
  __tmp_var_1 = if (typeof(subj_) != Array{Int32}) convert(Array{Int32},subj_) else subj_ end
  __tmp_var_2 = (numi_)
  __tmp_var_3 = zeros(Float64,__tmp_var_2)
  __tmp_var_4 = (numi_)
  __tmp_var_5 = zeros(Float64,__tmp_var_4)
  __tmp_var_6 = (numi_)
  __tmp_var_7 = zeros(Float64,__tmp_var_6)
  __tmp_var_8 = (numi_)
  __tmp_var_9 = zeros(Float64,__tmp_var_8)
  __tmp_var_10 = (numj_)
  __tmp_var_11 = zeros(Float64,__tmp_var_10)
  __tmp_var_12 = (numj_)
  __tmp_var_13 = zeros(Float64,__tmp_var_12)
  __tmp_var_14 = (numj_)
  __tmp_var_15 = zeros(Float64,__tmp_var_14)
  __tmp_var_16 = (numj_)
  __tmp_var_17 = zeros(Float64,__tmp_var_16)
  res = @msk_ccall( "primalsensitivity",Int32,(Ptr{Void},Int32,Ptr{Int32},Ptr{Int32},Int32,Ptr{Int32},Ptr{Int32},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},),task_.task,numi_,__tmp_var_0-1,marki_,numj_,__tmp_var_1-1,markj_,__tmp_var_3,__tmp_var_5,__tmp_var_7,__tmp_var_9,__tmp_var_11,__tmp_var_13,__tmp_var_15,__tmp_var_17)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (__tmp_var_3,__tmp_var_5,__tmp_var_7,__tmp_var_9,__tmp_var_11,__tmp_var_13,__tmp_var_15,__tmp_var_17)
end

printdata(task:: MSKtask,whichstream:: Int32,firsti,lasti,firstj,lastj,firstk,lastk,c,qo,a,qc,bc,bx,vartype,cones) = printdata(convert(MSKtask,task),convert(Int32,whichstream),convert(Int32,firsti),convert(Int32,lasti),convert(Int32,firstj),convert(Int32,lastj),convert(Int32,firstk),convert(Int32,lastk),convert(Int32,c),convert(Int32,qo),convert(Int32,a),convert(Int32,qc),convert(Int32,bc),convert(Int32,bx),convert(Int32,vartype),convert(Int32,cones))
function printdata(task_:: MSKtask,whichstream_:: Int32,firsti_:: Int32,lasti_:: Int32,firstj_:: Int32,lastj_:: Int32,firstk_:: Int32,lastk_:: Int32,c_:: Int32,qo_:: Int32,a_:: Int32,qc_:: Int32,bc_:: Int32,bx_:: Int32,vartype_:: Int32,cones_:: Int32)
  res = @msk_ccall( "printdata",Int32,(Ptr{Void},Int32,Int32,Int32,Int32,Int32,Int32,Int32,Int32,Int32,Int32,Int32,Int32,Int32,Int32,Int32,),task_.task,whichstream_,firsti_-1,lasti_-1,firstj_-1,lastj_-1,firstk_-1,lastk_-1,c_,qo_,a_,qc_,bc_,bx_,vartype_,cones_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

function printparam(task_:: MSKtask)
  res = @msk_ccall( "printparam",Int32,(Ptr{Void},),task_.task)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

putacol(task:: MSKtask,j,subj:: Array,valj:: Array{Float64}) = putacol(convert(MSKtask,task),convert(Int32,j),convert(Array{Int32},subj),convert(Array{Float64},valj))
function putacol(task_:: MSKtask,j_:: Int32,subj_:: Array{Int32},valj_:: Array{Float64})
  nzj_ = minimum([ length(subj_),length(valj_) ])
  __tmp_var_0 = if (typeof(subj_) != Array{Int32}) convert(Array{Int32},subj_) else subj_ end
  res = @msk_ccall( "putacol",Int32,(Ptr{Void},Int32,Int32,Ptr{Int32},Ptr{Float64},),task_.task,j_-1,nzj_,__tmp_var_0-1,valj_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

putacollist(task:: MSKtask,sub:: Array,ptrb:: Array,ptre:: Array,asub:: Array,aval:: Array{Float64}) = putacollist(convert(MSKtask,task),convert(Array{Int32},sub),convert(Array{Int64},ptrb),convert(Array{Int64},ptre),convert(Array{Int32},asub),convert(Array{Float64},aval))
function putacollist(task:: MSKtask,sub:: Array,A:: SparseMatrixCSC{Float64})
  ptrb = A.colptr[1:size(A,2)-1]
  ptre = A.colptr[2:size(A,2)]
  asub = A.rowval
  aval = A.nzval
  putacollist(task,sub,ptrb,ptre,asub,aval)
end
function putacollist(task_:: MSKtask,sub_:: Array{Int32},ptrb_:: Array{Int64},ptre_:: Array{Int64},asub_:: Array{Int32},aval_:: Array{Float64})
  num_ = minimum([ length(sub_),length(ptrb_),length(ptre_) ])
  __tmp_var_0 = if (typeof(sub_) != Array{Int32}) convert(Array{Int32},sub_) else sub_ end
  __tmp_var_1 = if (typeof(asub_) != Array{Int32}) convert(Array{Int32},asub_) else asub_ end
  res = @msk_ccall( "putacollist64",Int32,(Ptr{Void},Int32,Ptr{Int32},Ptr{Int64},Ptr{Int64},Ptr{Int32},Ptr{Float64},),task_.task,num_,__tmp_var_0-1,ptrb_-1,ptre_-1,__tmp_var_1-1,aval_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

putacolslice(task:: MSKtask,first,last,ptrb:: Array,ptre:: Array,asub:: Array,aval:: Array{Float64}) = putacolslice(convert(MSKtask,task),convert(Int32,first),convert(Int32,last),convert(Array{Int64},ptrb),convert(Array{Int64},ptre),convert(Array{Int32},asub),convert(Array{Float64},aval))
function putacolslice(task:: MSKtask,first,last,A:: SparseMatrixCSC{Float64})
  ptrb = A.colptr[1:size(A,2)-1]
  ptre = A.colptr[2:size(A,2)]
  asub = A.rowval
  aval = A.nzval
  putacolslice(task,first,last,ptrb,ptre,asub,aval)
end
function putacolslice(task_:: MSKtask,first_:: Int32,last_:: Int32,ptrb_:: Array{Int64},ptre_:: Array{Int64},asub_:: Array{Int32},aval_:: Array{Float64})
  __tmp_var_0 = ((last_) - (first_))
  if length(ptrb_) < __tmp_var_0
    println("Array argument ptrb is not long enough")
    throw(BoundsError())
  end
  __tmp_var_1 = ((last_) - (first_))
  if length(ptre_) < __tmp_var_1
    println("Array argument ptre is not long enough")
    throw(BoundsError())
  end
  __tmp_var_2 = if (typeof(asub_) != Array{Int32}) convert(Array{Int32},asub_) else asub_ end
  res = @msk_ccall( "putacolslice64",Int32,(Ptr{Void},Int32,Int32,Ptr{Int64},Ptr{Int64},Ptr{Int32},Ptr{Float64},),task_.task,first_-1,last_-1,ptrb_-1,ptre_-1,__tmp_var_2-1,aval_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

putaij(task:: MSKtask,i,j,aij:: Float64) = putaij(convert(MSKtask,task),convert(Int32,i),convert(Int32,j),convert(Float64,aij))
function putaij(task_:: MSKtask,i_:: Int32,j_:: Int32,aij_:: Float64)
  res = @msk_ccall( "putaij",Int32,(Ptr{Void},Int32,Int32,Float64,),task_.task,i_-1,j_-1,aij_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

putaijlist(task:: MSKtask,subi:: Array,subj:: Array,valij:: Array{Float64}) = putaijlist(convert(MSKtask,task),convert(Array{Int32},subi),convert(Array{Int32},subj),convert(Array{Float64},valij))
function putaijlist(task_:: MSKtask,subi_:: Array{Int32},subj_:: Array{Int32},valij_:: Array{Float64})
  num_ = minimum([ length(subi_),length(subj_),length(valij_) ])
  __tmp_var_0 = if (typeof(subi_) != Array{Int32}) convert(Array{Int32},subi_) else subi_ end
  __tmp_var_1 = if (typeof(subj_) != Array{Int32}) convert(Array{Int32},subj_) else subj_ end
  res = @msk_ccall( "putaijlist",Int32,(Ptr{Void},Int32,Ptr{Int32},Ptr{Int32},Ptr{Float64},),task_.task,num_,__tmp_var_0-1,__tmp_var_1-1,valij_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

putaijlist64(task:: MSKtask,subi:: Array,subj:: Array,valij:: Array{Float64}) = putaijlist64(convert(MSKtask,task),convert(Array{Int32},subi),convert(Array{Int32},subj),convert(Array{Float64},valij))
function putaijlist64(task_:: MSKtask,subi_:: Array{Int32},subj_:: Array{Int32},valij_:: Array{Float64})
  num_ = minimum([ length(subi_),length(subj_),length(valij_) ])
  __tmp_var_0 = if (typeof(subi_) != Array{Int32}) convert(Array{Int32},subi_) else subi_ end
  __tmp_var_1 = if (typeof(subj_) != Array{Int32}) convert(Array{Int32},subj_) else subj_ end
  res = @msk_ccall( "putaijlist64",Int32,(Ptr{Void},Int64,Ptr{Int32},Ptr{Int32},Ptr{Float64},),task_.task,num_,__tmp_var_0-1,__tmp_var_1-1,valij_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

putarow(task:: MSKtask,i,subi:: Array,vali:: Array{Float64}) = putarow(convert(MSKtask,task),convert(Int32,i),convert(Array{Int32},subi),convert(Array{Float64},vali))
function putarow(task_:: MSKtask,i_:: Int32,subi_:: Array{Int32},vali_:: Array{Float64})
  nzi_ = minimum([ length(subi_),length(vali_) ])
  __tmp_var_0 = if (typeof(subi_) != Array{Int32}) convert(Array{Int32},subi_) else subi_ end
  res = @msk_ccall( "putarow",Int32,(Ptr{Void},Int32,Int32,Ptr{Int32},Ptr{Float64},),task_.task,i_-1,nzi_,__tmp_var_0-1,vali_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

putarowlist(task:: MSKtask,sub:: Array,ptrb:: Array,ptre:: Array,asub:: Array,aval:: Array{Float64}) = putarowlist(convert(MSKtask,task),convert(Array{Int32},sub),convert(Array{Int64},ptrb),convert(Array{Int64},ptre),convert(Array{Int32},asub),convert(Array{Float64},aval))
function putarowlist(task:: MSKtask,sub:: Array,A:: SparseMatrixCSC{Float64})
  ptrb = A.colptr[1:size(A,2)-1]
  ptre = A.colptr[2:size(A,2)]
  asub = A.rowval
  aval = A.nzval
  putarowlist(task,sub,ptrb,ptre,asub,aval)
end
function putarowlist(task_:: MSKtask,sub_:: Array{Int32},ptrb_:: Array{Int64},ptre_:: Array{Int64},asub_:: Array{Int32},aval_:: Array{Float64})
  num_ = minimum([ length(sub_),length(ptrb_),length(ptre_) ])
  __tmp_var_0 = if (typeof(sub_) != Array{Int32}) convert(Array{Int32},sub_) else sub_ end
  __tmp_var_1 = if (typeof(asub_) != Array{Int32}) convert(Array{Int32},asub_) else asub_ end
  res = @msk_ccall( "putarowlist64",Int32,(Ptr{Void},Int32,Ptr{Int32},Ptr{Int64},Ptr{Int64},Ptr{Int32},Ptr{Float64},),task_.task,num_,__tmp_var_0-1,ptrb_-1,ptre_-1,__tmp_var_1-1,aval_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

putarowslice(task:: MSKtask,first,last,ptrb:: Array,ptre:: Array,asub:: Array,aval:: Array{Float64}) = putarowslice(convert(MSKtask,task),convert(Int32,first),convert(Int32,last),convert(Array{Int64},ptrb),convert(Array{Int64},ptre),convert(Array{Int32},asub),convert(Array{Float64},aval))
function putarowslice(task:: MSKtask,first,last,A:: SparseMatrixCSC{Float64})
  ptrb = A.colptr[1:size(A,2)-1]
  ptre = A.colptr[2:size(A,2)]
  asub = A.rowval
  aval = A.nzval
  putarowslice(task,first,last,ptrb,ptre,asub,aval)
end
function putarowslice(task_:: MSKtask,first_:: Int32,last_:: Int32,ptrb_:: Array{Int64},ptre_:: Array{Int64},asub_:: Array{Int32},aval_:: Array{Float64})
  __tmp_var_0 = ((last_) - (first_))
  if length(ptrb_) < __tmp_var_0
    println("Array argument ptrb is not long enough")
    throw(BoundsError())
  end
  __tmp_var_1 = ((last_) - (first_))
  if length(ptre_) < __tmp_var_1
    println("Array argument ptre is not long enough")
    throw(BoundsError())
  end
  __tmp_var_2 = if (typeof(asub_) != Array{Int32}) convert(Array{Int32},asub_) else asub_ end
  res = @msk_ccall( "putarowslice64",Int32,(Ptr{Void},Int32,Int32,Ptr{Int64},Ptr{Int64},Ptr{Int32},Ptr{Float64},),task_.task,first_-1,last_-1,ptrb_-1,ptre_-1,__tmp_var_2-1,aval_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

putbarablocktriplet(task:: MSKtask,num,subi:: Array,subj:: Array,subk:: Array,subl:: Array,valijkl:: Array{Float64}) = putbarablocktriplet(convert(MSKtask,task),convert(Int64,num),convert(Array{Int32},subi),convert(Array{Int32},subj),convert(Array{Int32},subk),convert(Array{Int32},subl),convert(Array{Float64},valijkl))
function putbarablocktriplet(task_:: MSKtask,num_:: Int64,subi_:: Array{Int32},subj_:: Array{Int32},subk_:: Array{Int32},subl_:: Array{Int32},valijkl_:: Array{Float64})
  __tmp_var_0 = (num_)
  if length(subi_) < __tmp_var_0
    println("Array argument subi is not long enough")
    throw(BoundsError())
  end
  __tmp_var_1 = if (typeof(subi_) != Array{Int32}) convert(Array{Int32},subi_) else subi_ end
  __tmp_var_2 = (num_)
  if length(subj_) < __tmp_var_2
    println("Array argument subj is not long enough")
    throw(BoundsError())
  end
  __tmp_var_3 = if (typeof(subj_) != Array{Int32}) convert(Array{Int32},subj_) else subj_ end
  __tmp_var_4 = (num_)
  if length(subk_) < __tmp_var_4
    println("Array argument subk is not long enough")
    throw(BoundsError())
  end
  __tmp_var_5 = if (typeof(subk_) != Array{Int32}) convert(Array{Int32},subk_) else subk_ end
  __tmp_var_6 = (num_)
  if length(subl_) < __tmp_var_6
    println("Array argument subl is not long enough")
    throw(BoundsError())
  end
  __tmp_var_7 = if (typeof(subl_) != Array{Int32}) convert(Array{Int32},subl_) else subl_ end
  __tmp_var_8 = (num_)
  if length(valijkl_) < __tmp_var_8
    println("Array argument valijkl is not long enough")
    throw(BoundsError())
  end
  res = @msk_ccall( "putbarablocktriplet",Int32,(Ptr{Void},Int64,Ptr{Int32},Ptr{Int32},Ptr{Int32},Ptr{Int32},Ptr{Float64},),task_.task,num_,__tmp_var_1-1,__tmp_var_3-1,__tmp_var_5-1,__tmp_var_7-1,valijkl_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

putbaraij(task:: MSKtask,i,j,sub:: Array,weights:: Array{Float64}) = putbaraij(convert(MSKtask,task),convert(Int32,i),convert(Int32,j),convert(Array{Int64},sub),convert(Array{Float64},weights))
function putbaraij(task_:: MSKtask,i_:: Int32,j_:: Int32,sub_:: Array{Int64},weights_:: Array{Float64})
  num_ = minimum([ length(sub_),length(weights_) ])
  res = @msk_ccall( "putbaraij",Int32,(Ptr{Void},Int32,Int32,Int64,Ptr{Int64},Ptr{Float64},),task_.task,i_-1,j_-1,num_,sub_-1,weights_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

putbarcblocktriplet(task:: MSKtask,num,subj:: Array,subk:: Array,subl:: Array,valjkl:: Array{Float64}) = putbarcblocktriplet(convert(MSKtask,task),convert(Int64,num),convert(Array{Int32},subj),convert(Array{Int32},subk),convert(Array{Int32},subl),convert(Array{Float64},valjkl))
function putbarcblocktriplet(task_:: MSKtask,num_:: Int64,subj_:: Array{Int32},subk_:: Array{Int32},subl_:: Array{Int32},valjkl_:: Array{Float64})
  __tmp_var_0 = (num_)
  if length(subj_) < __tmp_var_0
    println("Array argument subj is not long enough")
    throw(BoundsError())
  end
  __tmp_var_1 = if (typeof(subj_) != Array{Int32}) convert(Array{Int32},subj_) else subj_ end
  __tmp_var_2 = (num_)
  if length(subk_) < __tmp_var_2
    println("Array argument subk is not long enough")
    throw(BoundsError())
  end
  __tmp_var_3 = if (typeof(subk_) != Array{Int32}) convert(Array{Int32},subk_) else subk_ end
  __tmp_var_4 = (num_)
  if length(subl_) < __tmp_var_4
    println("Array argument subl is not long enough")
    throw(BoundsError())
  end
  __tmp_var_5 = if (typeof(subl_) != Array{Int32}) convert(Array{Int32},subl_) else subl_ end
  __tmp_var_6 = (num_)
  if length(valjkl_) < __tmp_var_6
    println("Array argument valjkl is not long enough")
    throw(BoundsError())
  end
  res = @msk_ccall( "putbarcblocktriplet",Int32,(Ptr{Void},Int64,Ptr{Int32},Ptr{Int32},Ptr{Int32},Ptr{Float64},),task_.task,num_,__tmp_var_1-1,__tmp_var_3-1,__tmp_var_5-1,valjkl_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

putbarcj(task:: MSKtask,j,sub:: Array,weights:: Array{Float64}) = putbarcj(convert(MSKtask,task),convert(Int32,j),convert(Array{Int64},sub),convert(Array{Float64},weights))
function putbarcj(task_:: MSKtask,j_:: Int32,sub_:: Array{Int64},weights_:: Array{Float64})
  num_ = minimum([ length(sub_),length(weights_) ])
  res = @msk_ccall( "putbarcj",Int32,(Ptr{Void},Int32,Int64,Ptr{Int64},Ptr{Float64},),task_.task,j_-1,num_,sub_-1,weights_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

putbarsj(task:: MSKtask,whichsol:: Int32,j,barsj:: Array{Float64}) = putbarsj(convert(MSKtask,task),convert(Int32,whichsol),convert(Int32,j),convert(Array{Float64},barsj))
function putbarsj(task_:: MSKtask,whichsol_:: Int32,j_:: Int32,barsj_:: Array{Float64})
  __tmp_var_0 = getlenbarvarj(task_,(j_))
  if length(barsj_) < __tmp_var_0
    println("Array argument barsj is not long enough")
    throw(BoundsError())
  end
  res = @msk_ccall( "putbarsj",Int32,(Ptr{Void},Int32,Int32,Ptr{Float64},),task_.task,whichsol_,j_-1,barsj_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

putbarvarname(task:: MSKtask,j,name:: String) = putbarvarname(convert(MSKtask,task),convert(Int32,j),convert(String,name))
function putbarvarname(task_:: MSKtask,j_:: Int32,name_:: String)
  res = @msk_ccall( "putbarvarname",Int32,(Ptr{Void},Int32,Ptr{Uint8},),task_.task,j_-1,bytestring(name_))
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

putbarxj(task:: MSKtask,whichsol:: Int32,j,barxj:: Array{Float64}) = putbarxj(convert(MSKtask,task),convert(Int32,whichsol),convert(Int32,j),convert(Array{Float64},barxj))
function putbarxj(task_:: MSKtask,whichsol_:: Int32,j_:: Int32,barxj_:: Array{Float64})
  __tmp_var_0 = getlenbarvarj(task_,(j_))
  if length(barxj_) < __tmp_var_0
    println("Array argument barxj is not long enough")
    throw(BoundsError())
  end
  res = @msk_ccall( "putbarxj",Int32,(Ptr{Void},Int32,Int32,Ptr{Float64},),task_.task,whichsol_,j_-1,barxj_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

putbound(task:: MSKtask,accmode:: Int32,i,bk:: Int32,bl:: Float64,bu:: Float64) = putbound(convert(MSKtask,task),convert(Int32,accmode),convert(Int32,i),convert(Int32,bk),convert(Float64,bl),convert(Float64,bu))
function putbound(task_:: MSKtask,accmode_:: Int32,i_:: Int32,bk_:: Int32,bl_:: Float64,bu_:: Float64)
  res = @msk_ccall( "putbound",Int32,(Ptr{Void},Int32,Int32,Int32,Float64,Float64,),task_.task,accmode_,i_-1,bk_,bl_,bu_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

putboundlist(task:: MSKtask,accmode:: Int32,sub:: Array,bk:: Array{Int32},bl:: Array{Float64},bu:: Array{Float64}) = putboundlist(convert(MSKtask,task),convert(Int32,accmode),convert(Array{Int32},sub),convert(Array{Int32},bk),convert(Array{Float64},bl),convert(Array{Float64},bu))
function putboundlist(task_:: MSKtask,accmode_:: Int32,sub_:: Array{Int32},bk_:: Array{Int32},bl_:: Array{Float64},bu_:: Array{Float64})
  num_ = minimum([ length(sub_),length(bk_),length(bl_),length(bu_) ])
  __tmp_var_0 = if (typeof(sub_) != Array{Int32}) convert(Array{Int32},sub_) else sub_ end
  res = @msk_ccall( "putboundlist",Int32,(Ptr{Void},Int32,Int32,Ptr{Int32},Ptr{Int32},Ptr{Float64},Ptr{Float64},),task_.task,accmode_,num_,__tmp_var_0-1,bk_,bl_,bu_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

putboundslice(task:: MSKtask,con:: Int32,first,last,bk:: Array{Int32},bl:: Array{Float64},bu:: Array{Float64}) = putboundslice(convert(MSKtask,task),convert(Int32,con),convert(Int32,first),convert(Int32,last),convert(Array{Int32},bk),convert(Array{Float64},bl),convert(Array{Float64},bu))
function putboundslice(task_:: MSKtask,con_:: Int32,first_:: Int32,last_:: Int32,bk_:: Array{Int32},bl_:: Array{Float64},bu_:: Array{Float64})
  __tmp_var_0 = ((last_) - (first_))
  if length(bk_) < __tmp_var_0
    println("Array argument bk is not long enough")
    throw(BoundsError())
  end
  __tmp_var_1 = ((last_) - (first_))
  if length(bl_) < __tmp_var_1
    println("Array argument bl is not long enough")
    throw(BoundsError())
  end
  __tmp_var_2 = ((last_) - (first_))
  if length(bu_) < __tmp_var_2
    println("Array argument bu is not long enough")
    throw(BoundsError())
  end
  res = @msk_ccall( "putboundslice",Int32,(Ptr{Void},Int32,Int32,Int32,Ptr{Int32},Ptr{Float64},Ptr{Float64},),task_.task,con_,first_-1,last_-1,bk_,bl_,bu_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

function putcfix(task_:: MSKtask,cfix_:: Float64)
  res = @msk_ccall( "putcfix",Int32,(Ptr{Void},Float64,),task_.task,cfix_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

putcj(task:: MSKtask,j,cj:: Float64) = putcj(convert(MSKtask,task),convert(Int32,j),convert(Float64,cj))
function putcj(task_:: MSKtask,j_:: Int32,cj_:: Float64)
  res = @msk_ccall( "putcj",Int32,(Ptr{Void},Int32,Float64,),task_.task,j_-1,cj_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

putclist(task:: MSKtask,subj:: Array,val:: Array{Float64}) = putclist(convert(MSKtask,task),convert(Array{Int32},subj),convert(Array{Float64},val))
function putclist(task_:: MSKtask,subj_:: Array{Int32},val_:: Array{Float64})
  num_ = minimum([ length(subj_),length(val_) ])
  __tmp_var_0 = if (typeof(subj_) != Array{Int32}) convert(Array{Int32},subj_) else subj_ end
  res = @msk_ccall( "putclist",Int32,(Ptr{Void},Int32,Ptr{Int32},Ptr{Float64},),task_.task,num_,__tmp_var_0-1,val_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

putconbound(task:: MSKtask,i,bk:: Int32,bl:: Float64,bu:: Float64) = putconbound(convert(MSKtask,task),convert(Int32,i),convert(Int32,bk),convert(Float64,bl),convert(Float64,bu))
function putconbound(task_:: MSKtask,i_:: Int32,bk_:: Int32,bl_:: Float64,bu_:: Float64)
  res = @msk_ccall( "putconbound",Int32,(Ptr{Void},Int32,Int32,Float64,Float64,),task_.task,i_-1,bk_,bl_,bu_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

putconboundlist(task:: MSKtask,sub:: Array,bkc:: Array{Int32},blc:: Array{Float64},buc:: Array{Float64}) = putconboundlist(convert(MSKtask,task),convert(Array{Int32},sub),convert(Array{Int32},bkc),convert(Array{Float64},blc),convert(Array{Float64},buc))
function putconboundlist(task_:: MSKtask,sub_:: Array{Int32},bkc_:: Array{Int32},blc_:: Array{Float64},buc_:: Array{Float64})
  num_ = minimum([ length(sub_),length(bkc_),length(blc_),length(buc_) ])
  __tmp_var_0 = if (typeof(sub_) != Array{Int32}) convert(Array{Int32},sub_) else sub_ end
  res = @msk_ccall( "putconboundlist",Int32,(Ptr{Void},Int32,Ptr{Int32},Ptr{Int32},Ptr{Float64},Ptr{Float64},),task_.task,num_,__tmp_var_0-1,bkc_,blc_,buc_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

putconboundslice(task:: MSKtask,first,last,bk:: Array{Int32},bl:: Array{Float64},bu:: Array{Float64}) = putconboundslice(convert(MSKtask,task),convert(Int32,first),convert(Int32,last),convert(Array{Int32},bk),convert(Array{Float64},bl),convert(Array{Float64},bu))
function putconboundslice(task_:: MSKtask,first_:: Int32,last_:: Int32,bk_:: Array{Int32},bl_:: Array{Float64},bu_:: Array{Float64})
  __tmp_var_0 = ((last_) - (first_))
  if length(bk_) < __tmp_var_0
    println("Array argument bk is not long enough")
    throw(BoundsError())
  end
  __tmp_var_1 = ((last_) - (first_))
  if length(bl_) < __tmp_var_1
    println("Array argument bl is not long enough")
    throw(BoundsError())
  end
  __tmp_var_2 = ((last_) - (first_))
  if length(bu_) < __tmp_var_2
    println("Array argument bu is not long enough")
    throw(BoundsError())
  end
  res = @msk_ccall( "putconboundslice",Int32,(Ptr{Void},Int32,Int32,Ptr{Int32},Ptr{Float64},Ptr{Float64},),task_.task,first_-1,last_-1,bk_,bl_,bu_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

putcone(task:: MSKtask,k,conetype:: Int32,conepar:: Float64,submem:: Array) = putcone(convert(MSKtask,task),convert(Int32,k),convert(Int32,conetype),convert(Float64,conepar),convert(Array{Int32},submem))
function putcone(task_:: MSKtask,k_:: Int32,conetype_:: Int32,conepar_:: Float64,submem_:: Array{Int32})
  nummem_ = minimum([ length(submem_) ])
  __tmp_var_0 = if (typeof(submem_) != Array{Int32}) convert(Array{Int32},submem_) else submem_ end
  res = @msk_ccall( "putcone",Int32,(Ptr{Void},Int32,Int32,Float64,Int32,Ptr{Int32},),task_.task,k_-1,conetype_,conepar_,nummem_,__tmp_var_0-1)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

putconename(task:: MSKtask,j,name:: String) = putconename(convert(MSKtask,task),convert(Int32,j),convert(String,name))
function putconename(task_:: MSKtask,j_:: Int32,name_:: String)
  res = @msk_ccall( "putconename",Int32,(Ptr{Void},Int32,Ptr{Uint8},),task_.task,j_-1,bytestring(name_))
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

putconname(task:: MSKtask,i,name:: String) = putconname(convert(MSKtask,task),convert(Int32,i),convert(String,name))
function putconname(task_:: MSKtask,i_:: Int32,name_:: String)
  res = @msk_ccall( "putconname",Int32,(Ptr{Void},Int32,Ptr{Uint8},),task_.task,i_-1,bytestring(name_))
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

putcslice(task:: MSKtask,first,last,slice:: Array{Float64}) = putcslice(convert(MSKtask,task),convert(Int32,first),convert(Int32,last),convert(Array{Float64},slice))
function putcslice(task_:: MSKtask,first_:: Int32,last_:: Int32,slice_:: Array{Float64})
  __tmp_var_0 = ((last_) - (first_))
  if length(slice_) < __tmp_var_0
    println("Array argument slice is not long enough")
    throw(BoundsError())
  end
  res = @msk_ccall( "putcslice",Int32,(Ptr{Void},Int32,Int32,Ptr{Float64},),task_.task,first_-1,last_-1,slice_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

function putdouparam(task_:: MSKtask,param_:: Int32,parvalue_:: Float64)
  res = @msk_ccall( "putdouparam",Int32,(Ptr{Void},Int32,Float64,),task_.task,param_,parvalue_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

putintparam(task:: MSKtask,param:: Int32,parvalue) = putintparam(convert(MSKtask,task),convert(Int32,param),convert(Int32,parvalue))
function putintparam(task_:: MSKtask,param_:: Int32,parvalue_:: Int32)
  res = @msk_ccall( "putintparam",Int32,(Ptr{Void},Int32,Int32,),task_.task,param_,parvalue_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

putmaxnumanz(task:: MSKtask,maxnumanz) = putmaxnumanz(convert(MSKtask,task),convert(Int64,maxnumanz))
function putmaxnumanz(task_:: MSKtask,maxnumanz_:: Int64)
  res = @msk_ccall( "putmaxnumanz",Int32,(Ptr{Void},Int64,),task_.task,maxnumanz_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

putmaxnumbarvar(task:: MSKtask,maxnumbarvar) = putmaxnumbarvar(convert(MSKtask,task),convert(Int32,maxnumbarvar))
function putmaxnumbarvar(task_:: MSKtask,maxnumbarvar_:: Int32)
  res = @msk_ccall( "putmaxnumbarvar",Int32,(Ptr{Void},Int32,),task_.task,maxnumbarvar_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

putmaxnumcon(task:: MSKtask,maxnumcon) = putmaxnumcon(convert(MSKtask,task),convert(Int32,maxnumcon))
function putmaxnumcon(task_:: MSKtask,maxnumcon_:: Int32)
  res = @msk_ccall( "putmaxnumcon",Int32,(Ptr{Void},Int32,),task_.task,maxnumcon_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

putmaxnumcone(task:: MSKtask,maxnumcone) = putmaxnumcone(convert(MSKtask,task),convert(Int32,maxnumcone))
function putmaxnumcone(task_:: MSKtask,maxnumcone_:: Int32)
  res = @msk_ccall( "putmaxnumcone",Int32,(Ptr{Void},Int32,),task_.task,maxnumcone_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

putmaxnumqnz(task:: MSKtask,maxnumqnz) = putmaxnumqnz(convert(MSKtask,task),convert(Int64,maxnumqnz))
function putmaxnumqnz(task_:: MSKtask,maxnumqnz_:: Int64)
  res = @msk_ccall( "putmaxnumqnz",Int32,(Ptr{Void},Int64,),task_.task,maxnumqnz_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

putmaxnumvar(task:: MSKtask,maxnumvar) = putmaxnumvar(convert(MSKtask,task),convert(Int32,maxnumvar))
function putmaxnumvar(task_:: MSKtask,maxnumvar_:: Int32)
  res = @msk_ccall( "putmaxnumvar",Int32,(Ptr{Void},Int32,),task_.task,maxnumvar_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

function putnadouparam(task_:: MSKtask,paramname_:: String,parvalue_:: Float64)
  res = @msk_ccall( "putnadouparam",Int32,(Ptr{Void},Ptr{Uint8},Float64,),task_.task,bytestring(paramname_),parvalue_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

putnaintparam(task:: MSKtask,paramname:: String,parvalue) = putnaintparam(convert(MSKtask,task),convert(String,paramname),convert(Int32,parvalue))
function putnaintparam(task_:: MSKtask,paramname_:: String,parvalue_:: Int32)
  res = @msk_ccall( "putnaintparam",Int32,(Ptr{Void},Ptr{Uint8},Int32,),task_.task,bytestring(paramname_),parvalue_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

function putnastrparam(task_:: MSKtask,paramname_:: String,parvalue_:: String)
  res = @msk_ccall( "putnastrparam",Int32,(Ptr{Void},Ptr{Uint8},Ptr{Uint8},),task_.task,bytestring(paramname_),bytestring(parvalue_))
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

function putobjname(task_:: MSKtask,objname_:: String)
  res = @msk_ccall( "putobjname",Int32,(Ptr{Void},Ptr{Uint8},),task_.task,bytestring(objname_))
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

function putobjsense(task_:: MSKtask,sense_:: Int32)
  res = @msk_ccall( "putobjsense",Int32,(Ptr{Void},Int32,),task_.task,sense_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

function putparam(task_:: MSKtask,parname_:: String,parvalue_:: String)
  res = @msk_ccall( "putparam",Int32,(Ptr{Void},Ptr{Uint8},Ptr{Uint8},),task_.task,bytestring(parname_),bytestring(parvalue_))
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

putqcon(task:: MSKtask,qcsubk:: Array,qcsubi:: Array,qcsubj:: Array,qcval:: Array{Float64}) = putqcon(convert(MSKtask,task),convert(Array{Int32},qcsubk),convert(Array{Int32},qcsubi),convert(Array{Int32},qcsubj),convert(Array{Float64},qcval))
function putqcon(task_:: MSKtask,qcsubk_:: Array{Int32},qcsubi_:: Array{Int32},qcsubj_:: Array{Int32},qcval_:: Array{Float64})
  numqcnz_ = minimum([ length(qcsubi_),length(qcsubj_),length(qcval_) ])
  __tmp_var_0 = if (typeof(qcsubk_) != Array{Int32}) convert(Array{Int32},qcsubk_) else qcsubk_ end
  __tmp_var_1 = if (typeof(qcsubi_) != Array{Int32}) convert(Array{Int32},qcsubi_) else qcsubi_ end
  __tmp_var_2 = if (typeof(qcsubj_) != Array{Int32}) convert(Array{Int32},qcsubj_) else qcsubj_ end
  res = @msk_ccall( "putqcon",Int32,(Ptr{Void},Int32,Ptr{Int32},Ptr{Int32},Ptr{Int32},Ptr{Float64},),task_.task,numqcnz_,__tmp_var_0-1,__tmp_var_1-1,__tmp_var_2-1,qcval_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

putqconk(task:: MSKtask,k,qcsubi:: Array,qcsubj:: Array,qcval:: Array{Float64}) = putqconk(convert(MSKtask,task),convert(Int32,k),convert(Array{Int32},qcsubi),convert(Array{Int32},qcsubj),convert(Array{Float64},qcval))
function putqconk(task_:: MSKtask,k_:: Int32,qcsubi_:: Array{Int32},qcsubj_:: Array{Int32},qcval_:: Array{Float64})
  numqcnz_ = minimum([ length(qcsubi_),length(qcsubj_),length(qcval_) ])
  __tmp_var_0 = if (typeof(qcsubi_) != Array{Int32}) convert(Array{Int32},qcsubi_) else qcsubi_ end
  __tmp_var_1 = if (typeof(qcsubj_) != Array{Int32}) convert(Array{Int32},qcsubj_) else qcsubj_ end
  res = @msk_ccall( "putqconk",Int32,(Ptr{Void},Int32,Int32,Ptr{Int32},Ptr{Int32},Ptr{Float64},),task_.task,k_-1,numqcnz_,__tmp_var_0-1,__tmp_var_1-1,qcval_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

putqobj(task:: MSKtask,qosubi:: Array,qosubj:: Array,qoval:: Array{Float64}) = putqobj(convert(MSKtask,task),convert(Array{Int32},qosubi),convert(Array{Int32},qosubj),convert(Array{Float64},qoval))
function putqobj(task_:: MSKtask,qosubi_:: Array{Int32},qosubj_:: Array{Int32},qoval_:: Array{Float64})
  numqonz_ = minimum([ length(qosubi_),length(qosubj_),length(qoval_) ])
  __tmp_var_0 = if (typeof(qosubi_) != Array{Int32}) convert(Array{Int32},qosubi_) else qosubi_ end
  __tmp_var_1 = if (typeof(qosubj_) != Array{Int32}) convert(Array{Int32},qosubj_) else qosubj_ end
  res = @msk_ccall( "putqobj",Int32,(Ptr{Void},Int32,Ptr{Int32},Ptr{Int32},Ptr{Float64},),task_.task,numqonz_,__tmp_var_0-1,__tmp_var_1-1,qoval_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

putqobjij(task:: MSKtask,i,j,qoij:: Float64) = putqobjij(convert(MSKtask,task),convert(Int32,i),convert(Int32,j),convert(Float64,qoij))
function putqobjij(task_:: MSKtask,i_:: Int32,j_:: Int32,qoij_:: Float64)
  res = @msk_ccall( "putqobjij",Int32,(Ptr{Void},Int32,Int32,Float64,),task_.task,i_-1,j_-1,qoij_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

function putskc(task_:: MSKtask,whichsol_:: Int32,skc_:: Array{Int32})
  __tmp_var_0 = getnumcon(task_)
  if length(skc_) < __tmp_var_0
    println("Array argument skc is not long enough")
    throw(BoundsError())
  end
  res = @msk_ccall( "putskc",Int32,(Ptr{Void},Int32,Ptr{Int32},),task_.task,whichsol_,skc_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

putskcslice(task:: MSKtask,whichsol:: Int32,first,last,skc:: Array{Int32}) = putskcslice(convert(MSKtask,task),convert(Int32,whichsol),convert(Int32,first),convert(Int32,last),convert(Array{Int32},skc))
function putskcslice(task_:: MSKtask,whichsol_:: Int32,first_:: Int32,last_:: Int32,skc_:: Array{Int32})
  __tmp_var_0 = ((last_) - (first_))
  if length(skc_) < __tmp_var_0
    println("Array argument skc is not long enough")
    throw(BoundsError())
  end
  res = @msk_ccall( "putskcslice",Int32,(Ptr{Void},Int32,Int32,Int32,Ptr{Int32},),task_.task,whichsol_,first_-1,last_-1,skc_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

function putskx(task_:: MSKtask,whichsol_:: Int32,skx_:: Array{Int32})
  __tmp_var_0 = getnumvar(task_)
  if length(skx_) < __tmp_var_0
    println("Array argument skx is not long enough")
    throw(BoundsError())
  end
  res = @msk_ccall( "putskx",Int32,(Ptr{Void},Int32,Ptr{Int32},),task_.task,whichsol_,skx_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

putskxslice(task:: MSKtask,whichsol:: Int32,first,last,skx:: Array{Int32}) = putskxslice(convert(MSKtask,task),convert(Int32,whichsol),convert(Int32,first),convert(Int32,last),convert(Array{Int32},skx))
function putskxslice(task_:: MSKtask,whichsol_:: Int32,first_:: Int32,last_:: Int32,skx_:: Array{Int32})
  __tmp_var_0 = ((last_) - (first_))
  if length(skx_) < __tmp_var_0
    println("Array argument skx is not long enough")
    throw(BoundsError())
  end
  res = @msk_ccall( "putskxslice",Int32,(Ptr{Void},Int32,Int32,Int32,Ptr{Int32},),task_.task,whichsol_,first_-1,last_-1,skx_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

function putslc(task_:: MSKtask,whichsol_:: Int32,slc_:: Array{Float64})
  __tmp_var_0 = getnumcon(task_)
  if length(slc_) < __tmp_var_0
    println("Array argument slc is not long enough")
    throw(BoundsError())
  end
  res = @msk_ccall( "putslc",Int32,(Ptr{Void},Int32,Ptr{Float64},),task_.task,whichsol_,slc_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

putslcslice(task:: MSKtask,whichsol:: Int32,first,last,slc:: Array{Float64}) = putslcslice(convert(MSKtask,task),convert(Int32,whichsol),convert(Int32,first),convert(Int32,last),convert(Array{Float64},slc))
function putslcslice(task_:: MSKtask,whichsol_:: Int32,first_:: Int32,last_:: Int32,slc_:: Array{Float64})
  __tmp_var_0 = ((last_) - (first_))
  if length(slc_) < __tmp_var_0
    println("Array argument slc is not long enough")
    throw(BoundsError())
  end
  res = @msk_ccall( "putslcslice",Int32,(Ptr{Void},Int32,Int32,Int32,Ptr{Float64},),task_.task,whichsol_,first_-1,last_-1,slc_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

function putslx(task_:: MSKtask,whichsol_:: Int32,slx_:: Array{Float64})
  __tmp_var_0 = getnumvar(task_)
  if length(slx_) < __tmp_var_0
    println("Array argument slx is not long enough")
    throw(BoundsError())
  end
  res = @msk_ccall( "putslx",Int32,(Ptr{Void},Int32,Ptr{Float64},),task_.task,whichsol_,slx_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

putslxslice(task:: MSKtask,whichsol:: Int32,first,last,slx:: Array{Float64}) = putslxslice(convert(MSKtask,task),convert(Int32,whichsol),convert(Int32,first),convert(Int32,last),convert(Array{Float64},slx))
function putslxslice(task_:: MSKtask,whichsol_:: Int32,first_:: Int32,last_:: Int32,slx_:: Array{Float64})
  __tmp_var_0 = ((last_) - (first_))
  if length(slx_) < __tmp_var_0
    println("Array argument slx is not long enough")
    throw(BoundsError())
  end
  res = @msk_ccall( "putslxslice",Int32,(Ptr{Void},Int32,Int32,Int32,Ptr{Float64},),task_.task,whichsol_,first_-1,last_-1,slx_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

function putsnx(task_:: MSKtask,whichsol_:: Int32,sux_:: Array{Float64})
  __tmp_var_0 = getnumvar(task_)
  if length(sux_) < __tmp_var_0
    println("Array argument sux is not long enough")
    throw(BoundsError())
  end
  res = @msk_ccall( "putsnx",Int32,(Ptr{Void},Int32,Ptr{Float64},),task_.task,whichsol_,sux_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

putsnxslice(task:: MSKtask,whichsol:: Int32,first,last,snx:: Array{Float64}) = putsnxslice(convert(MSKtask,task),convert(Int32,whichsol),convert(Int32,first),convert(Int32,last),convert(Array{Float64},snx))
function putsnxslice(task_:: MSKtask,whichsol_:: Int32,first_:: Int32,last_:: Int32,snx_:: Array{Float64})
  __tmp_var_0 = ((last_) - (first_))
  if length(snx_) < __tmp_var_0
    println("Array argument snx is not long enough")
    throw(BoundsError())
  end
  res = @msk_ccall( "putsnxslice",Int32,(Ptr{Void},Int32,Int32,Int32,Ptr{Float64},),task_.task,whichsol_,first_-1,last_-1,snx_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

function putsolution(task_:: MSKtask,whichsol_:: Int32,skc_:: Array{Int32},skx_:: Array{Int32},skn_:: Array{Int32},xc_:: Array{Float64},xx_:: Array{Float64},y_:: Array{Float64},slc_:: Array{Float64},suc_:: Array{Float64},slx_:: Array{Float64},sux_:: Array{Float64},snx_:: Array{Float64})
  res = @msk_ccall( "putsolution",Int32,(Ptr{Void},Int32,Ptr{Int32},Ptr{Int32},Ptr{Int32},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},),task_.task,whichsol_,skc_,skx_,skn_,xc_,xx_,y_,slc_,suc_,slx_,sux_,snx_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

putsolutioni(task:: MSKtask,accmode:: Int32,i,whichsol:: Int32,sk:: Int32,x:: Float64,sl:: Float64,su:: Float64,sn:: Float64) = putsolutioni(convert(MSKtask,task),convert(Int32,accmode),convert(Int32,i),convert(Int32,whichsol),convert(Int32,sk),convert(Float64,x),convert(Float64,sl),convert(Float64,su),convert(Float64,sn))
function putsolutioni(task_:: MSKtask,accmode_:: Int32,i_:: Int32,whichsol_:: Int32,sk_:: Int32,x_:: Float64,sl_:: Float64,su_:: Float64,sn_:: Float64)
  res = @msk_ccall( "putsolutioni",Int32,(Ptr{Void},Int32,Int32,Int32,Int32,Float64,Float64,Float64,Float64,),task_.task,accmode_,i_-1,whichsol_,sk_,x_,sl_,su_,sn_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

putsolutionyi(task:: MSKtask,i,whichsol:: Int32,y:: Float64) = putsolutionyi(convert(MSKtask,task),convert(Int32,i),convert(Int32,whichsol),convert(Float64,y))
function putsolutionyi(task_:: MSKtask,i_:: Int32,whichsol_:: Int32,y_:: Float64)
  res = @msk_ccall( "putsolutionyi",Int32,(Ptr{Void},Int32,Int32,Float64,),task_.task,i_-1,whichsol_,y_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

function putstrparam(task_:: MSKtask,param_:: Int32,parvalue_:: String)
  res = @msk_ccall( "putstrparam",Int32,(Ptr{Void},Int32,Ptr{Uint8},),task_.task,param_,bytestring(parvalue_))
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

function putsuc(task_:: MSKtask,whichsol_:: Int32,suc_:: Array{Float64})
  __tmp_var_0 = getnumcon(task_)
  if length(suc_) < __tmp_var_0
    println("Array argument suc is not long enough")
    throw(BoundsError())
  end
  res = @msk_ccall( "putsuc",Int32,(Ptr{Void},Int32,Ptr{Float64},),task_.task,whichsol_,suc_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

putsucslice(task:: MSKtask,whichsol:: Int32,first,last,suc:: Array{Float64}) = putsucslice(convert(MSKtask,task),convert(Int32,whichsol),convert(Int32,first),convert(Int32,last),convert(Array{Float64},suc))
function putsucslice(task_:: MSKtask,whichsol_:: Int32,first_:: Int32,last_:: Int32,suc_:: Array{Float64})
  __tmp_var_0 = ((last_) - (first_))
  if length(suc_) < __tmp_var_0
    println("Array argument suc is not long enough")
    throw(BoundsError())
  end
  res = @msk_ccall( "putsucslice",Int32,(Ptr{Void},Int32,Int32,Int32,Ptr{Float64},),task_.task,whichsol_,first_-1,last_-1,suc_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

function putsux(task_:: MSKtask,whichsol_:: Int32,sux_:: Array{Float64})
  __tmp_var_0 = getnumvar(task_)
  if length(sux_) < __tmp_var_0
    println("Array argument sux is not long enough")
    throw(BoundsError())
  end
  res = @msk_ccall( "putsux",Int32,(Ptr{Void},Int32,Ptr{Float64},),task_.task,whichsol_,sux_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

putsuxslice(task:: MSKtask,whichsol:: Int32,first,last,sux:: Array{Float64}) = putsuxslice(convert(MSKtask,task),convert(Int32,whichsol),convert(Int32,first),convert(Int32,last),convert(Array{Float64},sux))
function putsuxslice(task_:: MSKtask,whichsol_:: Int32,first_:: Int32,last_:: Int32,sux_:: Array{Float64})
  __tmp_var_0 = ((last_) - (first_))
  if length(sux_) < __tmp_var_0
    println("Array argument sux is not long enough")
    throw(BoundsError())
  end
  res = @msk_ccall( "putsuxslice",Int32,(Ptr{Void},Int32,Int32,Int32,Ptr{Float64},),task_.task,whichsol_,first_-1,last_-1,sux_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

function puttaskname(task_:: MSKtask,taskname_:: String)
  res = @msk_ccall( "puttaskname",Int32,(Ptr{Void},Ptr{Uint8},),task_.task,bytestring(taskname_))
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

putvarbound(task:: MSKtask,j,bk:: Int32,bl:: Float64,bu:: Float64) = putvarbound(convert(MSKtask,task),convert(Int32,j),convert(Int32,bk),convert(Float64,bl),convert(Float64,bu))
function putvarbound(task_:: MSKtask,j_:: Int32,bk_:: Int32,bl_:: Float64,bu_:: Float64)
  res = @msk_ccall( "putvarbound",Int32,(Ptr{Void},Int32,Int32,Float64,Float64,),task_.task,j_-1,bk_,bl_,bu_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

putvarboundlist(task:: MSKtask,sub:: Array,bkx:: Array{Int32},blx:: Array{Float64},bux:: Array{Float64}) = putvarboundlist(convert(MSKtask,task),convert(Array{Int32},sub),convert(Array{Int32},bkx),convert(Array{Float64},blx),convert(Array{Float64},bux))
function putvarboundlist(task_:: MSKtask,sub_:: Array{Int32},bkx_:: Array{Int32},blx_:: Array{Float64},bux_:: Array{Float64})
  num_ = minimum([ length(sub_),length(bkx_),length(blx_),length(bux_) ])
  __tmp_var_0 = if (typeof(sub_) != Array{Int32}) convert(Array{Int32},sub_) else sub_ end
  res = @msk_ccall( "putvarboundlist",Int32,(Ptr{Void},Int32,Ptr{Int32},Ptr{Int32},Ptr{Float64},Ptr{Float64},),task_.task,num_,__tmp_var_0-1,bkx_,blx_,bux_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

putvarboundslice(task:: MSKtask,first,last,bk:: Array{Int32},bl:: Array{Float64},bu:: Array{Float64}) = putvarboundslice(convert(MSKtask,task),convert(Int32,first),convert(Int32,last),convert(Array{Int32},bk),convert(Array{Float64},bl),convert(Array{Float64},bu))
function putvarboundslice(task_:: MSKtask,first_:: Int32,last_:: Int32,bk_:: Array{Int32},bl_:: Array{Float64},bu_:: Array{Float64})
  __tmp_var_0 = ((last_) - (first_))
  if length(bk_) < __tmp_var_0
    println("Array argument bk is not long enough")
    throw(BoundsError())
  end
  __tmp_var_1 = ((last_) - (first_))
  if length(bl_) < __tmp_var_1
    println("Array argument bl is not long enough")
    throw(BoundsError())
  end
  __tmp_var_2 = ((last_) - (first_))
  if length(bu_) < __tmp_var_2
    println("Array argument bu is not long enough")
    throw(BoundsError())
  end
  res = @msk_ccall( "putvarboundslice",Int32,(Ptr{Void},Int32,Int32,Ptr{Int32},Ptr{Float64},Ptr{Float64},),task_.task,first_-1,last_-1,bk_,bl_,bu_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

putvarbranchorder(task:: MSKtask,j,priority,direction:: Int32) = putvarbranchorder(convert(MSKtask,task),convert(Int32,j),convert(Int32,priority),convert(Int32,direction))
function putvarbranchorder(task_:: MSKtask,j_:: Int32,priority_:: Int32,direction_:: Int32)
  res = @msk_ccall( "putvarbranchorder",Int32,(Ptr{Void},Int32,Int32,Int32,),task_.task,j_-1,priority_,direction_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

putvarname(task:: MSKtask,j,name:: String) = putvarname(convert(MSKtask,task),convert(Int32,j),convert(String,name))
function putvarname(task_:: MSKtask,j_:: Int32,name_:: String)
  res = @msk_ccall( "putvarname",Int32,(Ptr{Void},Int32,Ptr{Uint8},),task_.task,j_-1,bytestring(name_))
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

putvartype(task:: MSKtask,j,vartype:: Int32) = putvartype(convert(MSKtask,task),convert(Int32,j),convert(Int32,vartype))
function putvartype(task_:: MSKtask,j_:: Int32,vartype_:: Int32)
  res = @msk_ccall( "putvartype",Int32,(Ptr{Void},Int32,Int32,),task_.task,j_-1,vartype_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

putvartypelist(task:: MSKtask,subj:: Array,vartype:: Array{Int32}) = putvartypelist(convert(MSKtask,task),convert(Array{Int32},subj),convert(Array{Int32},vartype))
function putvartypelist(task_:: MSKtask,subj_:: Array{Int32},vartype_:: Array{Int32})
  num_ = minimum([ length(subj_),length(vartype_) ])
  __tmp_var_0 = if (typeof(subj_) != Array{Int32}) convert(Array{Int32},subj_) else subj_ end
  res = @msk_ccall( "putvartypelist",Int32,(Ptr{Void},Int32,Ptr{Int32},Ptr{Int32},),task_.task,num_,__tmp_var_0-1,vartype_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

function putxc(task_:: MSKtask,whichsol_:: Int32)
  __tmp_var_0 = getnumcon(task_)
  __tmp_var_1 = zeros(Float64,__tmp_var_0)
  res = @msk_ccall( "putxc",Int32,(Ptr{Void},Int32,Ptr{Float64},),task_.task,whichsol_,__tmp_var_1)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (__tmp_var_1)
end

putxcslice(task:: MSKtask,whichsol:: Int32,first,last,xc:: Array{Float64}) = putxcslice(convert(MSKtask,task),convert(Int32,whichsol),convert(Int32,first),convert(Int32,last),convert(Array{Float64},xc))
function putxcslice(task_:: MSKtask,whichsol_:: Int32,first_:: Int32,last_:: Int32,xc_:: Array{Float64})
  __tmp_var_0 = ((last_) - (first_))
  if length(xc_) < __tmp_var_0
    println("Array argument xc is not long enough")
    throw(BoundsError())
  end
  res = @msk_ccall( "putxcslice",Int32,(Ptr{Void},Int32,Int32,Int32,Ptr{Float64},),task_.task,whichsol_,first_-1,last_-1,xc_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

function putxx(task_:: MSKtask,whichsol_:: Int32,xx_:: Array{Float64})
  __tmp_var_0 = getnumvar(task_)
  if length(xx_) < __tmp_var_0
    println("Array argument xx is not long enough")
    throw(BoundsError())
  end
  res = @msk_ccall( "putxx",Int32,(Ptr{Void},Int32,Ptr{Float64},),task_.task,whichsol_,xx_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

putxxslice(task:: MSKtask,whichsol:: Int32,first,last,xx:: Array{Float64}) = putxxslice(convert(MSKtask,task),convert(Int32,whichsol),convert(Int32,first),convert(Int32,last),convert(Array{Float64},xx))
function putxxslice(task_:: MSKtask,whichsol_:: Int32,first_:: Int32,last_:: Int32,xx_:: Array{Float64})
  __tmp_var_0 = ((last_) - (first_))
  if length(xx_) < __tmp_var_0
    println("Array argument xx is not long enough")
    throw(BoundsError())
  end
  res = @msk_ccall( "putxxslice",Int32,(Ptr{Void},Int32,Int32,Int32,Ptr{Float64},),task_.task,whichsol_,first_-1,last_-1,xx_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

function puty(task_:: MSKtask,whichsol_:: Int32,y_:: Array{Float64})
  __tmp_var_0 = getnumcon(task_)
  if length(y_) < __tmp_var_0
    println("Array argument y is not long enough")
    throw(BoundsError())
  end
  res = @msk_ccall( "puty",Int32,(Ptr{Void},Int32,Ptr{Float64},),task_.task,whichsol_,y_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

putyslice(task:: MSKtask,whichsol:: Int32,first,last,y:: Array{Float64}) = putyslice(convert(MSKtask,task),convert(Int32,whichsol),convert(Int32,first),convert(Int32,last),convert(Array{Float64},y))
function putyslice(task_:: MSKtask,whichsol_:: Int32,first_:: Int32,last_:: Int32,y_:: Array{Float64})
  __tmp_var_0 = ((last_) - (first_))
  if length(y_) < __tmp_var_0
    println("Array argument y is not long enough")
    throw(BoundsError())
  end
  res = @msk_ccall( "putyslice",Int32,(Ptr{Void},Int32,Int32,Int32,Ptr{Float64},),task_.task,whichsol_,first_-1,last_-1,y_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

function readbranchpriorities(task_:: MSKtask,filename_:: String)
  res = @msk_ccall( "readbranchpriorities",Int32,(Ptr{Void},Ptr{Uint8},),task_.task,bytestring(filename_))
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

function readdata(task_:: MSKtask,filename_:: String)
  res = @msk_ccall( "readdataautoformat",Int32,(Ptr{Void},Ptr{Uint8},),task_.task,bytestring(filename_))
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

function readdataformat(task_:: MSKtask,filename_:: String,format_:: Int32,compress_:: Int32)
  res = @msk_ccall( "readdataformat",Int32,(Ptr{Void},Ptr{Uint8},Int32,Int32,),task_.task,bytestring(filename_),format_,compress_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

function readparamfile(task_:: MSKtask)
  res = @msk_ccall( "readparamfile",Int32,(Ptr{Void},),task_.task)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

function readsolution(task_:: MSKtask,whichsol_:: Int32,filename_:: String)
  res = @msk_ccall( "readsolution",Int32,(Ptr{Void},Int32,Ptr{Uint8},),task_.task,whichsol_,bytestring(filename_))
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

function readsummary(task_:: MSKtask,whichstream_:: Int32)
  res = @msk_ccall( "readsummary",Int32,(Ptr{Void},Int32,),task_.task,whichstream_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

function readtask(task_:: MSKtask,filename_:: String)
  res = @msk_ccall( "readtask",Int32,(Ptr{Void},Ptr{Uint8},),task_.task,bytestring(filename_))
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

function relaxprimal(task_:: MSKtask,wlc_:: Array{Float64},wuc_:: Array{Float64},wlx_:: Array{Float64},wux_:: Array{Float64})
  relaxedtask_ = Array(Ptr{Void},(1,))
  __tmp_var_0 = getnumcon(task_)
  if length(wlc_) < __tmp_var_0
    println("Array argument wlc is not long enough")
    throw(BoundsError())
  end
  __tmp_var_1 = getnumcon(task_)
  if length(wuc_) < __tmp_var_1
    println("Array argument wuc is not long enough")
    throw(BoundsError())
  end
  __tmp_var_2 = getnumvar(task_)
  if length(wlx_) < __tmp_var_2
    println("Array argument wlx is not long enough")
    throw(BoundsError())
  end
  __tmp_var_3 = getnumvar(task_)
  if length(wux_) < __tmp_var_3
    println("Array argument wux is not long enough")
    throw(BoundsError())
  end
  res = @msk_ccall( "relaxprimal",Int32,(Ptr{Void},Ptr{Ptr{Void}},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},),task_.task,relaxedtask_,wlc_,wuc_,wlx_,wux_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(MSKtask,relaxedtask_[1]))
end

removebarvars(task:: MSKtask,subset:: Array) = removebarvars(convert(MSKtask,task),convert(Array{Int32},subset))
function removebarvars(task_:: MSKtask,subset_:: Array{Int32})
  num_ = minimum([ length(subset_) ])
  __tmp_var_0 = if (typeof(subset_) != Array{Int32}) convert(Array{Int32},subset_) else subset_ end
  res = @msk_ccall( "removebarvars",Int32,(Ptr{Void},Int32,Ptr{Int32},),task_.task,num_,__tmp_var_0)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

removecones(task:: MSKtask,subset:: Array) = removecones(convert(MSKtask,task),convert(Array{Int32},subset))
function removecones(task_:: MSKtask,subset_:: Array{Int32})
  num_ = minimum([ length(subset_) ])
  __tmp_var_0 = if (typeof(subset_) != Array{Int32}) convert(Array{Int32},subset_) else subset_ end
  res = @msk_ccall( "removecones",Int32,(Ptr{Void},Int32,Ptr{Int32},),task_.task,num_,__tmp_var_0)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

removecons(task:: MSKtask,subset:: Array) = removecons(convert(MSKtask,task),convert(Array{Int32},subset))
function removecons(task_:: MSKtask,subset_:: Array{Int32})
  num_ = minimum([ length(subset_) ])
  __tmp_var_0 = if (typeof(subset_) != Array{Int32}) convert(Array{Int32},subset_) else subset_ end
  res = @msk_ccall( "removecons",Int32,(Ptr{Void},Int32,Ptr{Int32},),task_.task,num_,__tmp_var_0)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

removevars(task:: MSKtask,subset:: Array) = removevars(convert(MSKtask,task),convert(Array{Int32},subset))
function removevars(task_:: MSKtask,subset_:: Array{Int32})
  num_ = minimum([ length(subset_) ])
  __tmp_var_0 = if (typeof(subset_) != Array{Int32}) convert(Array{Int32},subset_) else subset_ end
  res = @msk_ccall( "removevars",Int32,(Ptr{Void},Int32,Ptr{Int32},),task_.task,num_,__tmp_var_0)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

resizetask(task:: MSKtask,maxnumcon,maxnumvar,maxnumcone,maxnumanz,maxnumqnz) = resizetask(convert(MSKtask,task),convert(Int32,maxnumcon),convert(Int32,maxnumvar),convert(Int32,maxnumcone),convert(Int64,maxnumanz),convert(Int64,maxnumqnz))
function resizetask(task_:: MSKtask,maxnumcon_:: Int32,maxnumvar_:: Int32,maxnumcone_:: Int32,maxnumanz_:: Int64,maxnumqnz_:: Int64)
  res = @msk_ccall( "resizetask",Int32,(Ptr{Void},Int32,Int32,Int32,Int64,Int64,),task_.task,maxnumcon_,maxnumvar_,maxnumcone_,maxnumanz_,maxnumqnz_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

function sensitivityreport(task_:: MSKtask,whichstream_:: Int32)
  res = @msk_ccall( "sensitivityreport",Int32,(Ptr{Void},Int32,),task_.task,whichstream_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

function setdefaults(task_:: MSKtask)
  res = @msk_ccall( "setdefaults",Int32,(Ptr{Void},),task_.task)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

function solutiondef(task_:: MSKtask,whichsol_:: Int32)
  isdef_ = Array(Int32,(1,))
  res = @msk_ccall( "solutiondef",Int32,(Ptr{Void},Int32,Ptr{Int32},),task_.task,whichsol_,isdef_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Bool,isdef_[1]))
end

function solutionsummary(task_:: MSKtask,whichstream_:: Int32)
  res = @msk_ccall( "solutionsummary",Int32,(Ptr{Void},Int32,),task_.task,whichstream_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

solvewithbasis(task:: MSKtask,transp,numnz,sub:: Array,val:: Array{Float64}) = solvewithbasis(convert(MSKtask,task),convert(Int32,transp),convert(Int32,numnz),convert(Array{Int32},sub),convert(Array{Float64},val))
function solvewithbasis(task_:: MSKtask,transp_:: Int32,numnz_:: Int32,sub_:: Array{Int32},val_:: Array{Float64})
  __tmp_var_0 = [ numnz ]
  __tmp_var_1 = getnumcon(task_)
  if length(sub_) < __tmp_var_1
    println("Array argument sub is not long enough")
    throw(BoundsError())
  end
  __tmp_var_2 = if (typeof(sub_) != Array{Int32}) convert(Array{Int32},sub_) else sub_ end
  __tmp_var_3 = getnumcon(task_)
  if length(val_) < __tmp_var_3
    println("Array argument val is not long enough")
    throw(BoundsError())
  end
  res = @msk_ccall( "solvewithbasis",Int32,(Ptr{Void},Int32,Ptr{Int32},Ptr{Int32},Ptr{Float64},),task_.task,transp_,__tmp_var_0,__tmp_var_2-1,val_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (__tmp_var_0[1])
end

function startstat(task_:: MSKtask)
  res = @msk_ccall( "startstat",Int32,(Ptr{Void},),task_.task)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

function stopstat(task_:: MSKtask)
  res = @msk_ccall( "stopstat",Int32,(Ptr{Void},),task_.task)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

function strtoconetype(task_:: MSKtask,str_:: String)
  conetype_ = Array(Int32,(1,))
  res = @msk_ccall( "strtoconetype",Int32,(Ptr{Void},Ptr{Uint8},Ptr{Int32},),task_.task,bytestring(str_),conetype_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,conetype_[1]))
end

function strtosk(task_:: MSKtask,str_:: String)
  sk_ = Array(Int32,(1,))
  res = @msk_ccall( "strtosk",Int32,(Ptr{Void},Ptr{Uint8},Ptr{Int32},),task_.task,bytestring(str_),sk_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  (convert(Int32,sk_[1]))
end

function updatesolutioninfo(task_:: MSKtask,whichsol_:: Int32)
  res = @msk_ccall( "updatesolutioninfo",Int32,(Ptr{Void},Int32,),task_.task,whichsol_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

function writebranchpriorities(task_:: MSKtask,filename_:: String)
  res = @msk_ccall( "writebranchpriorities",Int32,(Ptr{Void},Ptr{Uint8},),task_.task,bytestring(filename_))
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

function writedata(task_:: MSKtask,filename_:: String)
  res = @msk_ccall( "writedata",Int32,(Ptr{Void},Ptr{Uint8},),task_.task,bytestring(filename_))
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

function writeparamfile(task_:: MSKtask,filename_:: String)
  res = @msk_ccall( "writeparamfile",Int32,(Ptr{Void},Ptr{Uint8},),task_.task,bytestring(filename_))
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

function writesolution(task_:: MSKtask,whichsol_:: Int32,filename_:: String)
  res = @msk_ccall( "writesolution",Int32,(Ptr{Void},Int32,Ptr{Uint8},),task_.task,whichsol_,bytestring(filename_))
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

function writetask(task_:: MSKtask,filename_:: String)
  res = @msk_ccall( "writetask",Int32,(Ptr{Void},Ptr{Uint8},),task_.task,bytestring(filename_))
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end

function checkinlicense(env_:: MSKenv,feature_:: Int32)
  res = @msk_ccall( "checkinlicense",Int32,(Ptr{Void},Int32,),env_.env,feature_)
  if res != 0
    throw (MosekError(res,""))
  end
end

function checkoutlicense(env_:: MSKenv,feature_:: Int32)
  res = @msk_ccall( "checkoutlicense",Int32,(Ptr{Void},Int32,),env_.env,feature_)
  if res != 0
    throw (MosekError(res,""))
  end
end

echointro(env:: MSKenv,longver) = echointro(convert(MSKenv,env),convert(Int32,longver))
function echointro(env_:: MSKenv,longver_:: Int32)
  res = @msk_ccall( "echointro",Int32,(Ptr{Void},Int32,),env_.env,longver_)
  if res != 0
    throw (MosekError(res,""))
  end
end

function getcodedesc(code_:: Int32)
  symname_ = zeros(Uint8,MSK_MAX_STR_LEN)
  str_ = zeros(Uint8,MSK_MAX_STR_LEN)
  res = @msk_ccall( "getcodedesc",Int32,(Int32,Ptr{Uint8},Ptr{Uint8},),code_,symname_,str_)
  if res != 0
    throw (MosekError(res,""))
  end
  (bytestring(symname_),bytestring(str_))
end

function getversion()
  major_ = Array(Int32,(1,))
  minor_ = Array(Int32,(1,))
  build_ = Array(Int32,(1,))
  revision_ = Array(Int32,(1,))
  res = @msk_ccall( "getversion",Int32,(Ptr{Int32},Ptr{Int32},Ptr{Int32},Ptr{Int32},),major_,minor_,build_,revision_)
  if res != 0
    throw (MosekError(res,""))
  end
  (convert(Int32,major_[1]),convert(Int32,minor_[1]),convert(Int32,build_[1]),convert(Int32,revision_[1]))
end

function licensecleanup()
  res = @msk_ccall( "licensecleanup",Int32,())
  if res != 0
    throw (MosekError(res,""))
  end
end

linkfiletostream(env:: MSKenv,whichstream:: Int32,filename:: String,append) = linkfiletostream(convert(MSKenv,env),convert(Int32,whichstream),convert(String,filename),convert(Int32,append))
function linkfiletostream(env_:: MSKenv,whichstream_:: Int32,filename_:: String,append_:: Int32)
  res = @msk_ccall( "linkfiletoenvstream",Int32,(Ptr{Void},Int32,Ptr{Uint8},Int32,),env_.env,whichstream_,bytestring(filename_),append_)
  if res != 0
    throw (MosekError(res,""))
  end
end

function putdllpath(env_:: MSKenv,dllpath_:: String)
  res = @msk_ccall( "putdllpath",Int32,(Ptr{Void},Ptr{Uint8},),env_.env,bytestring(dllpath_))
  if res != 0
    throw (MosekError(res,""))
  end
end

putkeepdlls(env:: MSKenv,keepdlls) = putkeepdlls(convert(MSKenv,env),convert(Int32,keepdlls))
function putkeepdlls(env_:: MSKenv,keepdlls_:: Int32)
  res = @msk_ccall( "putkeepdlls",Int32,(Ptr{Void},Int32,),env_.env,keepdlls_)
  if res != 0
    throw (MosekError(res,""))
  end
end

putlicensecode(env:: MSKenv,code:: Array) = putlicensecode(convert(MSKenv,env),convert(Array{Int32},code))
function putlicensecode(env_:: MSKenv,code_:: Array{Int32})
  __tmp_var_0 = MSK_LICENSE_BUFFER_LENGTH
  if length(code_) < __tmp_var_0
    println("Array argument code is not long enough")
    throw(BoundsError())
  end
  __tmp_var_1 = if (typeof(code_) != Array{Int32}) convert(Array{Int32},code_) else code_ end
  res = @msk_ccall( "putlicensecode",Int32,(Ptr{Void},Ptr{Int32},),env_.env,__tmp_var_1)
  if res != 0
    throw (MosekError(res,""))
  end
end

putlicensedebug(env:: MSKenv,licdebug) = putlicensedebug(convert(MSKenv,env),convert(Int32,licdebug))
function putlicensedebug(env_:: MSKenv,licdebug_:: Int32)
  res = @msk_ccall( "putlicensedebug",Int32,(Ptr{Void},Int32,),env_.env,licdebug_)
  if res != 0
    throw (MosekError(res,""))
  end
end

function putlicensepath(env_:: MSKenv,licensepath_:: String)
  res = @msk_ccall( "putlicensepath",Int32,(Ptr{Void},Ptr{Uint8},),env_.env,bytestring(licensepath_))
  if res != 0
    throw (MosekError(res,""))
  end
end

putlicensewait(env:: MSKenv,licwait) = putlicensewait(convert(MSKenv,env),convert(Int32,licwait))
function putlicensewait(env_:: MSKenv,licwait_:: Int32)
  res = @msk_ccall( "putlicensewait",Int32,(Ptr{Void},Int32,),env_.env,licwait_)
  if res != 0
    throw (MosekError(res,""))
  end
end

