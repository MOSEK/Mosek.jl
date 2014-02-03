module Mosek
  if isfile(joinpath(Pkg.dir("Mosek"),"deps","deps.jl"))
    include("../deps/deps.jl")
  else
    error("Mosek not properly installed. Please run Pkg.build(\"Mosek\")")
  end


  export
    makeenv, maketask,
    MosekError

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

    function MSKtask(t::MSKtask)
      temp = Array(Ptr{Void}, 1)
      res = @msk_ccall(clonetask, Int32, (Ptr{Void}, Ptr{Void}), t.task, temp)

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

  function maketask(task::MSKtask)
    MSKtask(task)
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
