using BinDeps

import BinDeps.libdir
import BinDeps.provides
import BinDeps.generate_steps

@BinDeps.setup

StringType = VERSION < v"0.4.0-dev" ? String : AbstractString

# define current version:
mskvmajor = "7"
mskvminor = "1"

libmosek      = library_dependency("libmosek",      aliases=["libmosek.so.7.0",     "libmosek.so.7.1",       # linux32x86
                                                             "libmosek64.so.7.0",   "libmosek64.so.7.1",     # linux64x86
                                                             "libmosek64.7.0.dylib","libmosek64.7.1.dylib",  # osx64x86
                                                             "mosek7_0.dll",        "mosek7_1.dll",          # win32x86
                                                             "mosek64_7_0.dll",     "mosek64_7_1.dll",       # win64x86
                                                             ])
libmosekscopt = library_dependency("libmosekscopt", aliases=["libmosekscopt7_0.so",    "libmosekscopt7_1.so",    # linux
                                                             "libmosekscopt7_0.dylib", "libmosekscopt7_1.dylib", # osx
                                                             "mosekscopt7_0.dll",      "mosekscopt7_1.dll",      # windows
                                                             ])
mskplatform,distroext =
  if WORD_SIZE == 32
    if     OS_NAME == :Linux "linux32x86",  ".tar.bz2"
    elseif OS_NAME == :Windows "win32x86",  ".zip"
    else   error("Platform not supported")
    end                                   
  else                                    
    if     OS_NAME == :Linux   "linux64x86",".tar.bz2"
    elseif OS_NAME == :Darwin  "osx64x86",  ".tar.bz2"
    elseif OS_NAME == :Windows "win64x86",  ".zip"
    else   error("Platform not supported")
    end
  end

# 1. Is MOSEKBINDIR set? If so this must point to the binaries dir in the MOSEK DISTRO
if haskey(ENV,"MOSEKBINDIR")
  provides(Binaries, ENV["MOSEKBINDIR"], libmosek)
  provides(Binaries, ENV["MOSEKBINDIR"], libmosekscopt)

elseif haskey(ENV,"MOSEK_7_1_BINDIR")
  provides(Binaries, ENV["MOSEK_7_1_BINDIR"], libmosek)

# 2a. Otherwise, use the default installation path (Linux)
elseif haskey(ENV,"HOME") && isdir(joinpath(ENV["HOME"],"mosek","7","tools","platform",mskplatform))
  provides(Binaries, joinpath(ENV["HOME"],"mosek","7","tools","platform",mskplatform,"bin"), libmosek)
  provides(Binaries, joinpath(ENV["HOME"],"mosek","7","tools","platform",mskplatform,"bin"), libmosekscopt)

# 2b. Windows default install path
elseif ( haskey(ENV,"HOMEDRIVE") && 
         haskey(ENV,"HOMEPATH") && 
         isdir(joinpath(string(ENV["HOMEDRIVE"],ENV["HOMEPATH"]),"mosek","7","tools","platform",mskplatform)) )
  home = string(ENV["HOMEDRIVE"],ENV["HOMEPATH"])
  provides(Binaries, joinpath(home,"mosek","7","tools","platform",mskplatform,"bin"), libmosek)
  provides(Binaries, joinpath(home,"mosek","7","tools","platform",mskplatform,"bin"), libmosekscopt)

# 3. Otherwise, fetch the MOSEK distro and unpack it
else

    # Custom BuildProcess: Basically the same as SimpleBuild, except
    # it accepts a 'path' argument that will be used as library path
    # (like in Binaries above).

    type LessSimpleBuild <: BuildProcess
        steps
        path
    end
    
    libdir(p::LessSimpleBuild,dep::BinDeps.LibraryDependency) = p.path
    provides(::Type{LessSimpleBuild},steps,path::StringType,dep; opts...) = provides(LessSimpleBuild(steps,path),dep; opts...)
    generate_steps(dep::BinDeps.LibraryDependency,h::LessSimpleBuild,opts) = h.steps


    srcdir  = joinpath(BinDeps.depsdir(libmosek),"src")
    tarname = string("mosektools",mskplatform, distroext)
    
    provides(Sources, URI(string("http://download.mosek.com/stable/7/",tarname)), libmosek, unpacked_dir="mosek")
    provides(LessSimpleBuild,
             GetSources(libmosek),
             "$srcdir/mosek/7/tools/platform/$mskplatform/bin", # path to binaries
             libmosek,
             os = :Unix)

    provides(LessSimpleBuild,
             GetSources(libmosek),
             "$srcdir/mosek/7/tools/platform/$mskplatform/bin", # path to binaries
             libmosekscopt,
             os = :Unix)
    # since it was never tested: disable windows install - should be
    # easy enough to add if requested.
end
@BinDeps.install [ :libmosek => :libmosek, 
                   :libmosekscopt => :libmosekscopt ]
