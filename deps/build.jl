include("liftedfrombindeps.jl")

# define current version:
mskvmajor = "8"
mskvminor = "0"


mskplatform,distroext =
  if Sys.ARCH == :i386 || Sys.ARCH == :i686
    if     is_linux()   "linux32x86",  ".tar.bz2"
    elseif is_windows() "win32x86",  ".zip"
    else   error("Platform not supported")
    end
  elseif Sys.ARCH == :x86_64
    if     is_linux()   "linux64x86",".tar.bz2"
    elseif is_apple()   "osx64x86",  ".tar.bz2"
    elseif is_windows() "win64x86",  ".zip"
    else   error("Platform not supported")
    end
  else
    error("Platform not supported")
  end

libalternatives =
  if     mskplatform == "linux32x86"
                                       [ ("libmosek.so.7.0",     "libmosekscopt7_0.so"),    ("libmosek.so.7.1",     "libmosekscopt7_1.so"),    ("libmosek.so.8.0",     "libmosekscopt8_0.so"),    ]
  elseif mskplatform == "linux64x86"
                                       [ ("libmosek64.so.7.0",   "libmosekscopt7_0.so"),    ("libmosek64.so.7.1",   "libmosekscopt7_1.so"),    ("libmosek64.so.8.0",   "libmosekscopt8_0.so"),    ]
  elseif mskplatform == "osx64x86"
                                       [ ("libmosek64.7.0.dylib","libmosekscopt7_0.dylib"), ("libmosek64.7.1.dylib","libmosekscopt7_1.dylib"), ("libmosek64.8.0.dylib","libmosekscopt8_0.dylib"), ]
  elseif mskplatform == "win32x86"
                                       [ ("mosek7_0.dll",        "mosekscopt7_0.dll"),      ("mosek7_1.dll",        "mosekscopt7_1.dll"),      ("mosek8_0.dll",        "mosekscopt8_0.dll"),      ]
  elseif mskplatform == "win64x86"
                                       [ ("mosek64_7_0.dll",     "mosekscopt7_0.dll"),      ("mosek64_7_1.dll",     "mosekscopt7_1.dll"),      ("mosek64_8_0.dll",     "mosekscopt8_0.dll"),      ]
  else   error("Platform not supported")
  end

bindepsdir = dirname(@__FILE__)

usepreinstalled = ! haskey(ENV,"MOSEKJL_FORCE_DOWNLOAD")

mskbindir =
# 1. Is MOSEKBINDIR set? If so this must point to the binaries dir in the MOSEK DISTRO
    if  usepreinstalled && haskey(ENV,"MOSEKBINDIR")
        ENV["MOSEKBINDIR"]
    elseif ! usepreinstalled && haskey(ENV,"MOSEK_8_0_BINDIR")
        ENV["MOSEK_8_0_BINDIR"]
    elseif ! usepreinstalled && haskey(ENV,"MOSEK_7_1_BINDIR")
        ENV["MOSEK_7_1_BINDIR"]
# 2a. Otherwise, use the default installation path (Linux)
    elseif usepreinstalled && ( haskey(ENV,"HOME") &&
                                  isdir(joinpath(ENV["HOME"],"mosek","8","tools","platform",mskplatform)))
        joinpath(ENV["HOME"],"mosek","8","tools","platform",mskplatform,"bin")
    elseif usepreinstalled && ( haskey(ENV,"HOME") &&
                                  isdir(joinpath(ENV["HOME"],"mosek","7","tools","platform",mskplatform)))
        joinpath(ENV["HOME"],"mosek","7","tools","platform",mskplatform,"bin")
# 2b. Windows default install path
    elseif usepreinstalled && (haskey(ENV,"HOMEDRIVE") &&
                                 haskey(ENV,"HOMEPATH") &&
                                 isdir(joinpath(string(ENV["HOMEDRIVE"],ENV["HOMEPATH"]),"mosek","8","tools","platform",mskplatform)))
        home = string(ENV["HOMEDRIVE"],ENV["HOMEPATH"])
        joinpath(home,"mosek","8","tools","platform",mskplatform,"bin")
    elseif usepreinstalled && (haskey(ENV,"HOMEDRIVE") &&
                                 haskey(ENV,"HOMEPATH") &&
                                 isdir(joinpath(string(ENV["HOMEDRIVE"],ENV["HOMEPATH"]),"mosek","7","tools","platform",mskplatform)))
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
mskbindir = replace(mskbindir,"\\","/")

idxs = reverse(find(libs -> all(lib -> isfile(joinpath(mskbindir,lib)), libs), libalternatives))
if length(idxs) == 0
    error("Failed to find any usable MOSEK libraries")
else
    libmosek,libmosekscopt = libalternatives[idxs[1]]

    libmosekpath = escape_string(normpath("$mskbindir/$libmosek"))
    libscoptpath = escape_string(normpath("$mskbindir/$libmosekscopt"))

    open(joinpath(bindepsdir,"deps.jl"),"w") do f
        write(f,"""# This is an auto-generated file; do not edit\n""")
        write(f,"""# Macro to load a library\n""")
        write(f,"""macro checked_lib(libname, path)\n""")
        write(f,"""    (Libdl.dlopen_e(path) == C_NULL) && error("Unable to load \\n\\n\$libname (\$path)\\n\\nPlease re-run Pkg.build(package), and restart Julia.")\n""")
        write(f,"""    quote const \$(esc(libname)) = \$path end\n""")
        write(f,"""end\n""")

        write(f,"""\n""")
        write(f,"""# Load dependencies\n""")
        write(f,string("""@checked_lib libmosek      \"$libmosekpath\"\n"""))
        write(f,string("""@checked_lib libmosekscopt \"$libscoptpath\"\n"""))
    end
end
