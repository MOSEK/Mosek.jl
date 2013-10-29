module Mosek
  export MosekError

# I am not entirely sure where this code belongs... In deps/build.jl? How to do that?
  const libmosek =
    if     OS_NAME == :Linux
      msklib = find_library(["libmosek64"],[])
      if msklib == ""
        msklib = 
          try
            f = open(joinpath(ENV["HOME"],".mosek"),"rt")
            msklibdir = readall(f)
            close(f)
            find_library(["libmosek64"],[msklibdir])
          catch
            ""
          end
      end
      msklib
    elseif OS_NAME == :Darwin
      find_library(["libmosek64"],[])
    elseif OS_NAME == :Windows
      find_library(["libmosek64_7_0"],[])
    end



  # Exported functions
  export makeenv, maketask

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
  # Environment: typedef void * MSKenv_t;
  type MSKenv
    env::Ptr{Void}
  end
  
  # Task: typedef void * MSKtask_t;
  type MSKtask
    task::Ptr{Void}
  end

  type MosekError
    rcode :: Int32
    msg   :: ASCIIString
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
end
