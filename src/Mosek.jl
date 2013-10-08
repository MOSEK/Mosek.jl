module Mosek
  
  # Exported functions
  export makeenv, maketask, readdata, optimize, writedata

  # Temporary: use BINDEPS to find path
  const msklibpath = "/home/idunning/mosek/7/tools/platform/linux64x86/bin/libmosek64.so"

  # A macro to make calling C API a little cleaner
  macro msk_ccall(func, args...)
    f = "MSK_$(func)"
    quote
      ccall(($f,msklibpath), $(args...))
    end
  end

  # -----
  # Types
  # -----
  # Environment: typedef void * MSKenv_t;
  type MSKenv
    env::Ptr{Void}
  end
  # Task: typedef void * MSKtask_t;
  type MSKtask
    task::Ptr{Void}
  end

  # ------------
  # API wrappers
  # ------------
  # TODO: Support other argument
  function makeenv()
    temp = Array(Ptr{Void}, 1)
    res = @msk_ccall(makeenv, Int32, (Ptr{Ptr{Void}}, Ptr{Uint8}), temp, C_NULL)
    if res != 0
      # TODO: Actually use result code
      error("MOSEK: Error creating environment")
    end
    return( MSKenv(temp[1]) )
  end

  # TODO: Support other arguments
  function maketask(env::MSKenv)
    temp = Array(Ptr{Void}, 1)
    res = @msk_ccall(maketask, Int32, (Ptr{Void}, Int32, Int32, Ptr{Void}), env.env, 0, 0, temp)
    if res != 0
      # TODO: Actually use result code
      error("MOSEK: Error creating task")
    end
    return( MSKtask(temp[1]) )
  end

  function readdata(task::MSKtask, filename)
    res = @msk_ccall(readdata, Int32, (Ptr{Void}, Ptr{Uint8}), task.task, filename)
    if res != 0
      # TODO: Actually use result code
      error("MOSEK: Error reading file")
    end
  end

  function optimize(task::MSKtask)
    res = @msk_ccall(optimize, Int32, (Ptr{Void},), task.task)
    if res != 0
      # TODO: Actually use result code
      error("MOSEK: Error optimizing")
    end
  end

  function writedata(task::MSKtask, filename)
    res = @msk_ccall(writedata, Int32, (Ptr{Void}, Ptr{Uint8}), task.task, filename)
    if res != 0
      # TODO: Actually use result code
      error("MOSEK: Error optimizing")
    end
  end

end
