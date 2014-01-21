using BinDeps

@BinDeps.setup

libmosek = library_dependency("libmosek", aliases=["libmosek64","mosek7_0","mosek64_7_0"])

mskplatform,pfdlls,distroext,libmosekname =
  if WORD_SIZE == 32
    if     OS_NAME == :Linux "linux32x86",  ["libmosek.so",       "libmosek.so.7.0",     "libiomp5.so",    "libmosekglb.so.7.0"],   ".tar.bz2", "libmosek"
    elseif OS_NAME == :Windows "win32x86",  ["mosek_7_0.dll",     "mosek_7_0.dll",       "libiomp5md.dll", "mosekglb_7_0"        ], ".zip",     "mosek7_0"
    else   error("Platform not supported")                         
    end                                   
  else                                    
    if     OS_NAME == :Linux   "linux64x86",["libmosek64.so",     "libmosek64.so.7.0",   "libiomp5.so",    "libmosekglb64.so.7.0"],    ".tar.bz2", "libmosek64"
    elseif OS_NAME == :Darwin  "osx4x86",   ["libmosek64.dylib",  "libmosek64.7.0.dylib","libiomp5.dylib", "libmosekglb64.7.0.dylib"], ".tar.bz2", "libmosek64"
    elseif OS_NAME == :Windows "win64x86",  ["libmosek64_7_0.dll","libmosek64_7_0.dll",  "libiomp5md.dll", "mosekglb64_7_0.dll"],      ".zip",     "moske64_7_0"
    else   error("Platform not supported")
    end
  end

# 1. Is MOSEKBINDIR set? If so this must point to the binaries dir in the MOSEK DISTRO  
if haskey(ENV,"MOSEKBINDIR")
  provides(Binaries, ENV["MOSEKBINDIR"], libmosek)
# 2. Otherwise, use the default installation path
elseif haskey(ENV,"HOME") && isdir(joinpath(ENV["HOME"],"mosek","7","tools","platform",mskplatform))
  if WORD_SIZE == 32
    provides(Binaries, joinpath(ENV["HOME"],"mosek","7","tools","platform","linux32x86","bin"), libmosek, os = :Linux)
    provides(Binaries, joinpath(ENV["HOME"],"mosek","7","tools","platform","win32x86","bin"),   libmosek, os = :Windows)
  else
    provides(Binaries, joinpath(ENV["HOME"],"mosek","7","tools","platform","linux64x86","bin"), libmosek, os = :Linux)
    provides(Binaries, joinpath(ENV["HOME"],"mosek","7","tools","platform","osx64x86","bin"),   libmosek, os = :Darwin)
    provides(Binaries, joinpath(ENV["HOME"],"mosek","7","tools","platform","win64x86","bin"),   libmosek, os = :Windows)
  end
elseif haskey(ENV,"HOMEDRIVE") && haskey(ENV,"HOMEPATH") && isdir(joinpath(string(ENV["HOMEDRIVE"],ENV["HOMEPATH"]),"mosek","7","tools","platform",mskplatform))
  home = string(ENV["HOMEDRIVE"],ENV["HOMEPATH"])
  if WORD_SIZE == 32
    provides(Binaries, joinpath(home,"mosek","7","tools","platform","linux32x86","bin"), libmosek, os = :Linux)
    provides(Binaries, joinpath(home,"mosek","7","tools","platform","win32x86","bin"),   libmosek, os = :Windows)
  else
    provides(Binaries, joinpath(home,"mosek","7","tools","platform","linux64x86","bin"), libmosek, os = :Linux)
    provides(Binaries, joinpath(home,"mosek","7","tools","platform","osx64x86","bin"),   libmosek, os = :Darwin)
    provides(Binaries, joinpath(home,"mosek","7","tools","platform","win64x86","bin"),   libmosek, os = :Windows)
  end
# 3. Otherwise, fetch the MOSEK distro and unpack it
else
  copycmd = if OS_NAME == :Windows "copy" else "cp" end
  srcdir  = joinpath(BinDeps.depsdir(libmosek),"src")
  tarname = string("mosektools",mskplatform, distroext)
  prefix  = joinpath(BinDeps.depsdir(libmosek),"usr")

  provides(Sources, URI(string("http://download.mosek.com/stable/7/",tarname)), libmosek, unpacked_dir="mosek")
  provides(SimpleBuild
    ( @build_steps begin
        GetSources(libmosek)
        CreateDirectory(joinpath(BinDeps.depsdir(libmosek),"usr","lib"))
        ( @build_steps begin
          `$copycmd "$srcdir/mosek/7/tools/platform/$mskplatform/bin/$(pfdlls[2])" "$prefix/lib"`
          `$copycmd "$srcdir/mosek/7/tools/platform/$mskplatform/bin/$(pfdlls[3])" "$prefix/lib"`
          `$copycmd "$srcdir/mosek/7/tools/platform/$mskplatform/bin/$(pfdlls[4])" "$prefix/lib"`
          `$copycmd "$srcdir/mosek/7/tools/platform/$mskplatform/bin/$(pfdlls[1])" "$prefix/lib"`
          end
        )
      end),
      libmosek)
end
@BinDeps.install
