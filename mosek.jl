const msklibpath = "/home/idunning/mosek/7/tools/platform/linux64x86/bin/libmosek64.so"

macro msk_ccall(func, args...)
  f = "MSK_$(func)"
  quote
    ccall(($f,msklibpath), $(args...))
  end
end

# Environment
# typedef void * MSKenv_t;
type MSKenv
  env::Ptr{Void}
end
temp = Array(Ptr{Void}, 1)
res = @msk_ccall(makeenv, Int32, (Ptr{Ptr{Void}}, Ptr{Uint8}), temp, C_NULL)
println("MSK_makeenv ", res)
myenv = MSKenv(temp[1])

# Task
# typedef void * MSKtask_t;
type MSKtask
  task::Ptr{Void}
end
res = @msk_ccall(maketask, Int32, (Ptr{Void}, Int32, Int32, Ptr{Void}), myenv.env, 0, 0, temp)
println("MSK_maketask ", res)
mytask = MSKtask(temp[1])

# Read data
filename = "/home/idunning/mosek/7/tools/examples/data/25fv47.mps"
res = @msk_ccall(readdata, Int32, (Ptr{Void}, Ptr{Uint8}), mytask.task, filename)
println("MSK_readdata ", res)

# Solve
res = @msk_ccall(optimize, Int32, (Ptr{Void},), mytask.task)
println("MSK_optimize ", res)

# Write output
filename = "/home/idunning/out.txt"
res = @msk_ccall(writedata, Int32, (Ptr{Void}, Ptr{Uint8}), mytask.task, filename)
println("MSK_writedata ", res)
