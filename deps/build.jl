include("liftedfrombindeps.jl")

# define current version:
mskvmajor = "7"
mskvminor = "1"

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

libalternatives = 
  if     mskplatform == "linux32x86" 
                                       [ ("libmosek.so.7.0",     "libmosekscopt7_0.so"),    ("libmosek.so.7.1",     "libmosekscopt7_1.so") ]
  elseif mskplatform == "linux64x86" 
                                       [ ("libmosek64.so.7.0",   "libmosekscopt7_0.so"),    ("libmosek64.so.7.1",   "libmosekscopt7_1.so") ]
  elseif mskplatform == "osx64x86"   
                                       [ ("libmosek64.7.0.dylib","libmosekscopt7_0.dylib"), ("libmosek64.7.1.dylib","libmosekscopt7_1.dylib") ]
  elseif mskplatform == "win32x86"   
                                       [ ("mosek7_0.dll",        "mosekscopt7_0.dll"),      ("mosek7_1.dll",        "mosekscopt7_1.dll") ]
  elseif mskplatform == "win64x86"   
                                       [ ("mosek64_7_0.dll",     "mosekscopt7_0.dll"),      ("mosek64_7_1.dll",     "mosekscopt7_1.dll") ]
  else   error("Platform not supported")
  end

bindepsdir = dirname(@__FILE__)
     
mskbindir = 
# 1. Is MOSEKBINDIR set? If so this must point to the binaries dir in the MOSEK DISTRO
    if haskey(ENV,"MOSEKBINDIR")
        ENV["MOSEKBINDIR"],idxs
    elseif haskey(ENV,"MOSEK_7_1_BINDIR")
        ENV["MOSEK_7_1_BINDIR"]
# 2a. Otherwise, use the default installation path (Linux)
    elseif haskey(ENV,"HOME") && isdir(joinpath(ENV["HOME"],"mosek","7","tools","platform",mskplatform))
        joinpath(ENV["HOME"],"mosek","7","tools","platform",mskplatform,"bin")
# 2b. Windows default install path
    elseif ( haskey(ENV,"HOMEDRIVE") && 
            haskey(ENV,"HOMEPATH") && 
            isdir(joinpath(string(ENV["HOMEDRIVE"],ENV["HOMEPATH"]),"mosek","7","tools","platform",mskplatform)) )
        home = string(ENV["HOMEDRIVE"],ENV["HOMEPATH"])
        joinpath(home,"mosek","7","tools","platform",mskplatform,"bin")
# 3. Otherwise, fetch the MOSEK distro and unpack it
    else
        srcdir  = joinpath(bindepsdir,"src")
        dldir   = joinpath(bindepsdir,"downloads")
        tarname = string("mosektools",mskplatform, distroext)
        
        basename,ext,sndext = splittarpath(tarname)

        mkpath(dldir)
        info("Download MOSEK distro (http://download.mosek.com/stable/7/$tarname)")
        success(download_cmd(string("http://download.mosek.com/stable/7/",tarname), joinpath(dldir,tarname))) || error("Failed to download MOSEK distro")
        mkpath(srcdir)
        info("Unpack MOSEK distro ($dldir/$tarname -> $srcdir)")
        success(unpack_cmd(joinpath(dldir,tarname),srcdir, ext, sndext)) || error("Failed to unpack MOSEK distro")

        joinpath(srcdir,"mosek","7","tools","platform",mskplatform,"bin")
    end

idxs = reverse(find(libs -> all(lib -> isfile(joinpath(mskbindir,lib)), libs), libalternatives))
if length(idxs) == 0 
    error("Failed to find any usable MOSEK libraries")
else
    libmosek,libmosekscopt = libalternatives[idxs[1]]

    f = open(joinpath(bindepsdir,"deps.jl"),"w")
    write(f,"""# This is an auto-generated file; do not edit
# Macro to load a library
macro checked_lib(libname, path)
    (dlopen_e(path) == C_NULL) && error("Unable to load \\n\\n\$libname (\$path)\\n\\nPlease re-run Pkg.build(package), and restart Julia.")
    quote const \$(esc(libname)) = \$path end
end

# Load dependencies
@checked_lib libmosek      "$mskbindir/$libmosek"
@checked_lib libmosekscopt "$mskbindir/$libmosekscopt"
""")
    close(f)
end
