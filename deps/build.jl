include("liftedfrombindeps.jl")

# define current version:
mskvmajor = "9"
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

bindepsdir = dirname(@__FILE__)

function findlibs(path::AbstractString,mskvmajor::AbstractString,mskvminor::AbstractString)
    moseklib,scoptlib =
        if Sys.ARCH == :i386 || Sys.ARCH == :i686
            if     is_windows() "mosek$(mskvmajor)_$(mskvminor).dll",        "mosekscopt$(mskvmajor)_$(mskvminor).dll"
            else   error("Platform not supported")
            end
        elseif Sys.ARCH == :x86_64
            if     is_linux()   "libmosek64.so.$(mskvmajor).$(mskvminor)",    "libmosekscopt$(mskvmajor)_$(mskvminor).so"
            elseif is_apple()   "libmosek64.$(mskvmajor).$(mskvminor).dylib", "libmosekscopt$(mskvmajor)_$(mskvminor).dylib"
            elseif is_windows() "mosek64_$(mskvmajor)_$(mskvminor).dll",      "mosekscopt$(mskvmajor)_$(mskvminor).dll"
            else   error("Platform not supported")
            end
        else
            error("Platform not supported")
        end

    let moseklibpath = joinpath(path,moseklib),
        scoptlibpath = joinpath(path,scoptlib)

        if ! isfile(moseklibpath)
            error("Library '$moseklib' not found ($path)")
        elseif ! isfile(scoptlibpath)
            error("Library '$scoptlib' not found ($path)")
        else
            moseklibpath,scoptlibpath
        end
    end
end

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
    if ver != nothing
        ver = split(ver,".")

        return ver[1] == mskvmajor && ver[2] == mskvminor
        #return ver[1] == mskvmajor && ver[2] == mskvminor
    else
        return false
    end
end


# Detect previous installation method:
#   "internal" -> the MOSEK distro was downloaded and installed by the installer
#   "external" -> the MOSEK distro was detected in some other location and not downloaded by the installer
# If the previous installation was internal, we will not attempt to locate a MOSEK installation
#
# This ensures that you can rebuild an internal installation without
# the installer checking for other installations.
instmethod =
    try
        open(joinpath(bindepsdir,"inst_method"),"r") do f
            strip(readstring(f))
        end
    catch
        nothing
    end

curmosekbindir =
    try
        open(joinpath(bindepsdir,"mosekbindir"),"r") do f
            strip(readstring(f))
        end
    catch
        nothing
    end
forcedownload = haskey(ENV,"MOSEKJL_FORCE_DOWNLOAD")

mskbindir =
# 1. If MOSEKBINDIR we use that path (and if it is not valid we produce an error message)
    if  ! forcedownload && haskey(ENV,"MOSEKBINDIR")
        mosekbindir = ENV["MOSEKBINDIR"]

        if ! bindirIsCurrentVersion(mosekbindir)
            error("MOSEKBINDIR ($mosekbindir) does not point to a MOSEK bin directory")
        end
        instmethod = "external"

        mosekbindir
# 2a. If last build used a user-specified MOSEK, we check if that is still valid and use that again
    elseif ! forcedownload &&
        instmethod == "external" &&
        curmosekbindir != nothing &&
        bindirIsCurrentVersion(curmosekbindir)

        curmosekbindir
# 2b. Otherwise, look in the UNIX default installation path
    elseif ! forcedownload &&
        haskey(ENV,"HOME") &&
        bindirIsCurrentVersion(joinpath(ENV["HOME"],"mosek","$mskvmajor","tools","platform",mskplatform,"bin"))

        instmethod = "external"

        joinpath(ENV["HOME"],"mosek","$mskvmajor","tools","platform",mskplatform,"bin")
# 2c. Windows default install path
    elseif ! forcedownload &&
        haskey(ENV,"HOMEDRIVE") &&
        haskey(ENV,"HOMEPATH") &&
        bindirIsCurrentVersion(joinpath(string(ENV["HOMEDRIVE"],ENV["HOMEPATH"]),"mosek","$mskvmajor","tools","platform",mskplatform,"bin"))

        home = string(ENV["HOMEDRIVE"],ENV["HOMEPATH"])

        instmethod = "external"

        joinpath(home,"mosek","$mskvmajor","tools","platform",mskplatform,"bin")
# 3. Otherwise, fetch the MOSEK distro and unpack it
    else
        srcdir   = joinpath(bindepsdir,"src")
        dldir    = joinpath(bindepsdir,"downloads")
        archname = "mosektools$(mskplatform)$(distroext)"
        hosturl  = "https://www.mosek.com/downloads/default_dns.txt"

        mkpath(dldir)
        success(download_cmd(hosturl, joinpath(dldir,"downloadhostname"))) || error("Failed to get MOSEK download host")
        downloadhost =
            open(joinpath(dldir,"downloadhostname"),"r") do f
                strip(readstring(f))
            end

        verurl   = "https://$downloadhost/stable/$mskvmajor/version"

        cur_version =
            if isfile(joinpath(bindepsdir,"version"))
                open(joinpath(bindepsdir,"version"),"r") do f
                    cur_version = strip(readstring(f))
                end
            else
                nothing
            end

        info("Get latest MOSEK version ($verurl)")
        success(download_cmd(verurl, joinpath(dldir,"new_version"))) || error("Failed to get MOSEK version")

        new_version =
            open(joinpath(dldir,"new_version"),"r") do f
                strip(readstring(f))
            end
        info("Latest MOSEK version = $new_version, currently installed = $cur_version")

        instmethod = "internal"





        if  cur_version == nothing ||
            cur_version != new_version ||
            !bindirIsCurrentVersion(joinpath(srcdir,"mosek",mskvmajor,"tools","platform",mskplatform,"bin"))

            verarr = split(new_version,'.')

            archurl = "https://$downloadhost/stable/$(new_version)/$archname"
            info("Download MOSEK distro ($archurl)")

            basename,ext,sndext = splittarpath(archname)

            success(download_cmd(archurl, joinpath(dldir,archname))) || error("Failed to download MOSEK distro")

            mkpath(srcdir)
            info("Unpack MOSEK distro ($dldir/$archname -> $srcdir)")
            success(unpack_cmd(joinpath(dldir,archname),srcdir, ext, sndext)) || error("Failed to unpack MOSEK distro")

            info("MOSEK installation complete.")
            joinpath(srcdir,"mosek",mskvmajor,"tools","platform",mskplatform,"bin")
        else
            info("Update not necessary")
            joinpath(srcdir,"mosek",mskvmajor,"tools","platform",mskplatform,"bin")
        end
    end

mskbindir = replace(mskbindir,"\\","/")


version = versionFromBindir(mskbindir)

if version == nothing
    error("MOSEK package is broken")
else
    open(joinpath(bindepsdir,"version"),"w") do f
        write(f,version)
    end
end

verarr = split(version,'.')

libmosekpath,libscoptpath = findlibs(mskbindir,verarr[1],verarr[2])

if instmethod == "external"
    info("""Found MOSEK $version at $mskbindir""")
end

open(joinpath(bindepsdir,"inst_method"),"w") do f
    write(f,instmethod)
end
open(joinpath(bindepsdir,"mosekbindir"),"w") do f
    write(f,mskbindir)
end

libmosekpath = escape_string(normpath(libmosekpath))
libscoptpath = escape_string(normpath(libscoptpath))

open(joinpath(bindepsdir,"deps.jl"),"w") do f
    write(f,"""
# This is an auto-generated file; do not edit
# Macro to load a library
macro checked_lib(libname, path)
    (Libdl.dlopen_e(path) == C_NULL) && error("Unable to load \\n\\n\$libname (\$path)\\n\\nPlease re-run Pkg.build(package), and restart Julia.")
    quote const \$(esc(libname)) = \$path end
end


# Load dependencies
@checked_lib libmosek      "$libmosekpath"
@checked_lib libmosekscopt "$libscoptpath"
""")


end # open
