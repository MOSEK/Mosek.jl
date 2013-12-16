module Mosek
  using BinDeps
  @BinDeps.load_dependencies


  export 
    makeenv, maketask,
    MosekError
  

  #require(joinpath(Pkg.dir("MathProgBase"),"src","LinprogSolverInterface.jl"))
  
  # I am not entirely sure where this code belongs... In deps/build.jl? How to do that?
  # search paths should be something like this: 
  #   1) library path (DYLD_LIBRARY_PATH, LD_LIBRARY_PATH or PATH)
  #   2) $HOME/mosek/7/tools/platform/[platformname]/bin/[libmosekname]
  #   3) 
#  const libmosek =
#    begin
#      libname,pfname = 
#      if     OS_NAME == :Linux
#        if WORD_SIZE == 32 "libmosek","linux32x86"
#        else               "libmosek64","linux64x86"
#        end
#      elseif OS_NAME == :Darwin
#        "libmosek64","osx64x86"
#      elseif OS_NAME == :Windows
#        if WORD_SIZE == 32 "libmosek7_0","win32x86"
#        else               "libmosek64_7_0","win64x86"
#        end
#      else
#        error("Platform not supported")
#      end
#
#      dfltlibpath = joinpath(ENV["HOME"],"mosek","7","tools","platform",pfname,"bin")
#     
#      msklib = find_library([libname],[])
#      if msklib == ""
#        msklib = find_library([libname],[dfltlibpath])
#      end
#      msklib
#    end

  # Temporary: use BINDEPS to find path  

  # A macro to make calling C API a little cleaner
  macro msk_ccall(func, args...)
    f = "MSK_$(func)"
    quote
      ccall(($f,libmosek), $(args...))
    end
  end

  # -----
  # Types
  # -----
  type MosekError
    rcode :: Int32
    msg   :: ASCIIString
  end

  # Environment: typedef void * MSKenv_t;
  type MSKenv
    env::Ptr{Void}
    streamcallbackfunc::Any
  end






  
  # Task: typedef void * MSKtask_t;
  type MSKtask
    env::MSKenv
    task::Ptr{Void}
    
    streamcallbackfunc:: Any
    callbackfunc:: Any
    nlinfo:: Any

    function MSKtask(env::MSKenv)
      temp = Array(Ptr{Void}, 1)
      res = @msk_ccall(maketask, Int32, (Ptr{Void}, Int32, Int32, Ptr{Void}), env.env, 0, 0, temp)

      if res != MSK_RES_OK
        throw(MosekError(res,""))
      end     
      
      task = new(env,temp[1],nothing,nothing,nothing) 

      finalizer(task,deletetask)

      task
    end
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
    MSKenv(temp[1],nothing)
  end

  const msk_global_env = makeenv() :: MSKenv

  function maketask(env::MSKenv)
    MSKtask(env)
  end
  
  function maketask()
    MSKtask(msk_global_env)
  end

  function deletetask(t::MSKtask)
    if t.task != C_NULL
      temp = Array(Ptr{Void},1)
      temp[1] = t.task
      @msk_ccall(deletetask,Int32,(Ptr{Ptr{Void}},), temp)
      t.task = C_NULL
    end
  end
  
  function deleteenv(e::MSKenv)
    if e.env != C_NULL
      temp = Array(Ptr{Void},1)
      temp[1] = t.env
      @msk_ccall(deleteenv,Int32,(Ptr{Ptr{Void}},), temp)
      e.env = C_NULL
    end
  end

  function getlasterror(t::MSKtask)
    lasterrcode = Array(Cint,1)
    lastmsglen = Array(Cint,1)
    
    @msk_ccall(getlasterror,Cint,(Ptr{Void},Ptr{Cint},Cint,Ptr{Cint},Ptr{Uint8}),
               t.task, lasterrcode, 0, lastmsglen, C_NULL)
    lastmsg = Array(Uint8,lastmsglen[1])
    @msk_ccall(getlasterror,Cint,(Ptr{Void},Ptr{Cint},Cint,Ptr{Cint},Ptr{Uint8}), 
               t.task, lasterrcode, lastmsglen[1], lastmsglen, lastmsg)
    convert(ASCIIString,lastmsg[1:lastmsglen[1]-1])
  end


  #include("msk_callback.jl")
  # Generated content
  include("msk_enums.jl")
  include("msk_functions.jl")
  include("msk_callback.jl")
  include("msk_geco.jl")

  include("MosekSolverInterface.jl")
  using Mosek.MosekMathProgSolverInterface
  export MosekSolver
end
