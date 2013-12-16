using BinDeps

@BinDeps.setup

libmosek = library_dependency("libmosek", aliases=["libmosek64","libmosek7_0","libmosek64_7_0"])

mskplatform,pfdlls,distroext =
 if WORD_SIZE == 32
   if     OS_NAME == :Linux "linux32x86",  "lib*.so*", ".tar.bz2"
   elseif OS_NAME == :Windows "win32x86",  "*.dll",    ".zip"
   else   error("Platform not supported")                         
   end                                   
 else                                    
   if     OS_NAME == :Linux   "linux64x86","lib*.so*", ".tar.bz2"
   elseif OS_NAME == :Darwin  "osx4x86",   "*.dylib",  ".tar.bz2"
   elseif OS_NAME == :Windows "win64x86",  "*.dll",    ".zip"
   else   error("Platform not supported")
   end
 end

# 1. Is MOSEKBINDIR set? If so this must point to the binaries dir in the MOSEK DISTRO  
if haskey(ENV,"MOSEKBINDIR")
  provides(Binaries, ENV["MOSEKBINDIR"], libmosek)
# 2. Otherwise, use the default installation path
else haskey(ENV,"HOME") && isdir(joinpath(ENV["HOME"],"mosek","7","tools","platform",mskplatform))
  if WORD_SIZE == 32
    provides(Binaries, joinpath(ENV["HOME"],"mosek","7","tools","platform","linux32x86","bin"), libmosek, os = :Linux)
    provides(Binaries, joinpath(ENV["HOME"],"mosek","7","tools","platform","win32x86","bin"),   libmosek, os = :Windows)
  else
    provides(Binaries, joinpath(ENV["HOME"],"mosek","7","tools","platform","linux64x86","bin"), libmosek, os = :Linux)
    provides(Binaries, joinpath(ENV["HOME"],"mosek","7","tools","platform","osx64x86","bin"),   libmosek, os = :Darwin)
    provides(Binaries, joinpath(ENV["HOME"],"mosek","7","tools","platform","win64x86","bin"),   libmosek, os = :Windows)
  end
# 3. Otherwise, fetch the MOSEK distro and unpack it
# else
#   mskplatform,pfdlls,distroext =
#     if WORD_SIZE == 32
#       if     OS_NAME == :Linux "linux32x86",  "lib*.so*", ".tar.bz2"
#       elseif OS_NAME == :Windows "win32x86",  "*.dll",    ".zip"
#       else   error("Platform not supported")                         
#       end                                   
#     else                                    
#       if     OS_NAME == :Linux   "linux64x86","lib*.so*", ".tar.bz2"
#       elseif OS_NAME == :Darwin  "osx4x86",   "*.dylib",  ".tar.bz2"
#       elseif OS_NAME == :Windows "win64x86",  "*.dll",    ".zip"
#       else   error("Platform not supported")
#       end
#     end
# 
#   copycmd = if OS_NAME == :Windows "copy" else "cp" end
#   srcdir  = joinpath(BinDeps.depsdir(libmosek),"src")
#   tarname = string("mosektools",mskplatform, distroext)
#   prefix  = joinpath(BinDeps.depsdir(libmosek),"usr")
# 
#   provides(Sources, URI(string("http://download.mosek.com/stable/7/",tarname)), libmosek, unpacked_dir="mosek")
# 
#   provides(SimpleBuild
#     ( @build_steps begin
#         GetSources(libmosek)
#         `$copycmd "$srcdir/mosek/7/tools/platform/$mskplatform/bin/$pfdlls" "$prefix/lib"`
#       end),
#       libmosek)
end

@BinDeps.install
