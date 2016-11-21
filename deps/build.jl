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
                                       [ ("libmosek.so.8.0",     "libmosekscopt8_0.so"),    ]
  elseif mskplatform == "linux64x86"
                                       [ ("libmosek64.so.8.0",   "libmosekscopt8_0.so"),    ]
  elseif mskplatform == "osx64x86"
                                       [ ("libmosek64.8.0.dylib","libmosekscopt8_0.dylib"), ]
  elseif mskplatform == "win32x86"
                                       [ ("mosek8_0.dll",        "mosekscopt8_0.dll"),      ]
  elseif mskplatform == "win64x86"
                                       [ ("mosek64_8_0.dll",     "mosekscopt8_0.dll"),      ]
  else   error("Platform not supported")
  end

bindepsdir = dirname(@__FILE__)





function versionFromBindir(bindir ::String)
    try
        mosekbin = if is_windows() "mosek.exe" else "mosek" end
        txt = readstring(`$bindir/$mosekbin`)
        m = match(r"\s*MOSEK Version ([0-9]+\.[0-9]+\.[0-9]+\.[0-9]+)",txt)
        if m == nothing
            return nothing
        else
            return m.captures[1]
        end
    catch
        return nothing
    end
end

function bindirIsCurrentVersion(bindir)
    ver = versionFromBindir(bindir)
    return ver != nothing && ver[1] == mskvmajor && ver[2] == mskvminor
end


# Detect previous installation method:
#   "internal" -> the MOSEK distro was downloaded and installed by the installer
#   "external" -> the MOSEK distro was detected in some other location and not downloaded by the installer
# If the previous installation was internal, we will not attempt to locate a MOSEK installation
instmethod =
    try
        open(joinpath(bindepsdir,"inst_method"),"r") do f
            strip(readstring(f))
        end
    catch
        nothing
    end

usepreinstalled = ! haskey(ENV,"MOSEKJL_FORCE_DOWNLOAD") && instmethod != "internal"
    

mskbindir =
# 1. Is MOSEKBINDIR set? If so this must point to the binaries dir in the MOSEK DISTRO
    if  usepreinstalled && haskey(ENV,"MOSEKBINDIR")  && bindirIsCurrentVersion(ENV["MOSEKBINDIR"])
        ENV["MOSEKBINDIR"]
    elseif ! usepreinstalled && haskey(ENV,"MOSEK_8_0_BINDIR") && bindirIsCurrentVersion(ENV["MOSEKBINDIR"])
        ENV["MOSEK_8_0_BINDIR"]
# 2a. Otherwise, use the default installation path (Linux)
    elseif usepreinstalled &&
        haskey(ENV,"HOME") &&
        bindirIsCurrentVersion(joinpath(ENV["HOME"],"mosek","8","tools","platform",mskplatform,"bin"))
        
        joinpath(ENV["HOME"],"mosek","8","tools","platform",mskplatform,"bin")
# 2b. Windows default install path
    elseif usepreinstalled &&
        haskey(ENV,"HOMEDRIVE") &&
        haskey(ENV,"HOMEPATH") &&
        bindirIsCurrentVersion(joinpath(string(ENV["HOMEDRIVE"],ENV["HOMEPATH"]),"mosek","8","tools","platform",mskplatform,"bin"))
        
        home = string(ENV["HOMEDRIVE"],ENV["HOMEPATH"])
        joinpath(home,"mosek","8","tools","platform",mskplatform,"bin")
# 3. Otherwise, fetch the MOSEK distro and unpack it
    else
        srcdir   = joinpath(bindepsdir,"src")
        dldir    = joinpath(bindepsdir,"downloads")
        archname = "mosektools$(mskplatform)$(distroext)"
        verurl   = "http://download.mosek.com/stable/$mskvmajor/version"

        mkpath(dldir)

        cur_version =
            if isfile(joinpath(dldir,"version"))
                open(joinpath(dldir,"version"),"r") do f
                    cur_version = strip(readstring(f))
                end
            else
                nothing
            end

        info("Get latest MOSEK version (http://download.mosek.com/stable/version)")
        success(download_cmd(verurl, joinpath(dldir,"new_version"))) || error("Failed to get MOSEK version")

        new_version =         
            open(joinpath(dldir,"new_version"),"r") do f
                strip(readstring(f))
            end
        info("Latest MOSEK version = $new_version")
        
        if cur_version == nothing || cur_version != new_version
            archurl = "http://download.mosek.com/stable/$(new_version)/$archname"
            info("Download MOSEK distro ($archurl)")
            
            basename,ext,sndext = splittarpath(archname)

            success(download_cmd(archurl, joinpath(dldir,archname))) || error("Failed to download MOSEK distro")

            mkpath(srcdir)
            info("Unpack MOSEK distro ($dldir/$archname -> $srcdir)")
            success(unpack_cmd(joinpath(dldir,archname),srcdir, ext, sndext)) || error("Failed to unpack MOSEK distro")

            
            bindir = joinpath(srcdir,"mosek",mskvmajor,"tools","platform",mskplatform,"bin")
            dl_version = versionFromBindir(bindir)
            if new_version == nothing
                error("MOSEK package is broken")
            else
                open(joinpath(dldir,"version"),"w") do f
                    write(f,dl_version)
                end
            end
            
            bindir
        else
            joinpath(srcdir,"mosek",mskvmajor,"tools","platform",mskplatform,"bin")
        end
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
