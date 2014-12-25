# This file contains the implementation of the geometric programming interface

export
  MSK_OPR_ENT,
  MSK_OPR_EXP,
  MSK_OPR_LOG,
  MSK_OPR_POW,
  scbegin,
  scend

@checked_lib libmosekscopt "/Users/yeesian/mosek/7/tools/platform/osx64x86/bin/libmosekscopt7_1.dylib"

# A macro to make calling C API a little cleaner
macro scopt_ccall(func, args...)
  f = Base.Meta.quot(symbol("MSK_$(func)"))
  args = [esc(a) for a in args]
  quote
    ccall(($f,libmosekscopt), $(args...))
  end
end

# msk_enums
const MSK_OPR_ENT = convert(Int32,23)
const MSK_OPR_EXP = convert(Int32,1)
const MSK_OPR_LOG = convert(Int32,2)
const MSK_OPR_POW = convert(Int32,3)

# msk_functions
function scbegin(task_::MSKtask,
                 opro_::Vector{Cint},     # Defines the functions used for the objective.
                 oprjo_::Vector{Cint},    # Defines the variable indexes used in non-linear objective function.
                 oprfo_::Vector{Cdouble}, # Defines constants used in the objective.
                 oprgo_::Vector{Cdouble}, # Defines constants used in the objective.
                 oprho_::Vector{Cdouble}, # Defines constants used in the objective.
                 oprc_::Vector{Cint},     # Defines the functions used for the constraints.
                 opric_::Vector{Cint},    # Defines the variable indexes used in the non-linear constraint functions.
                 oprjc_::Vector{Cint},    # Defines the constraint indexes where non-linear functions appear.
                 oprfc_::Vector{Cdouble}, # Defines constants used in the non-linear constraints.
                 oprgc_::Vector{Cdouble}, # Defines constants used in the non-linear constraints.
                 oprhc_::Vector{Cdouble}) # Defines constants used in the non-linear constraints.
  numopro_ = int32(length(opro_))
  numoprc_ = int32(length(oprc_))
  sch = Array(Ptr{Void}, 1)
  res = @scopt_ccall(scbegin,Cint,(Ptr{Void},Cint,Ptr{Cint},Ptr{Cint},Ptr{Cdouble},Ptr{Cdouble},
                     Ptr{Cdouble},Cint,Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cdouble},Ptr{Cdouble},
                     Ptr{Cdouble},Ptr{Void}),task_.task,numopro_,opro_,oprjo_,oprfo_,oprgo_,
                     oprho_,numoprc_,oprc_,opric_,oprjc_,oprfc_,oprgc_,oprhc_,sch)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
  return sch
end

function scend(task_::MSKtask,
               sch_::Vector{Ptr{Void}})
  res = @scopt_ccall(scend,Cint,(Ptr{Void},Ptr{Void}),task_.task,sch_)
  if res != MSK_RES_OK
    msg = getlasterror(task_)
    throw (MosekError(res,msg))
  end
end