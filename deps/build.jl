include("liftedfrombindeps.jl")


import WinReg

# define current version:
mosekver = open(joinpath(dirname(@__FILE__),"..","MOSEKVER")) do f readline(f) end
mskvmajor,mskvminor = split(mosekver,".")

mskplatform,distroext =
  if Sys.ARCH == :i386 || Sys.ARCH == :i686
    if     Sys.iswindows() "win32x86",  ".zip"
    else   error("Platform not supported in i686")
    end
  elseif Sys.ARCH == :x86_64
    if     Sys.islinux()   "linux64x86",".tar.bz2"
    elseif Sys.isapple()   "osx64x86",  ".tar.bz2"
    elseif Sys.iswindows() "win64x86",  ".zip"
    else   error("Platform not supported on AMD64")
    end
  elseif Sys.ARCH == :aarch64
    if     Sys.islinux()   "linuxaarch64",".tar.bz2"
    elseif Sys.isapple()   "osxaarch64",  ".tar.bz2"
    else   error("System not supported on AArch64")
    end
  else
    error("Platform not supported $(Sys.ARCH), $(Sys.isapple())")
  end

bindepsdir = dirname(@__FILE__)

function getregistrykey(base::UInt32,path::String,key::String)
    WinReg.querykey(base,path,key)
end

function hasregistrykey(base::UInt32,path::String,key::String)
    if ! Sys.iswindows()
        false
    else
        try
            getregistrykey(base,path,key)
	    true
        catch err
            false
        end
    end
end

function findlibs(path::AbstractString,mskvmajor::AbstractString,mskvminor::AbstractString)
    moseklib =
        if Sys.ARCH == :i386 || Sys.ARCH == :i686
            if     Sys.iswindows() "mosek$(mskvmajor)_$(mskvminor).dll"
            else   error("Platform not supported i686")
            end
        elseif Sys.ARCH == :x86_64
            if     Sys.islinux()   "libmosek64.so.$(mskvmajor).$(mskvminor)"
            elseif Sys.isapple()   "libmosek64.$(mskvmajor).$(mskvminor).dylib"
            elseif Sys.iswindows() "mosek64_$(mskvmajor)_$(mskvminor).dll"
            else   error("Unexpected platform for AMD64")
            end
        elseif Sys.ARCH == :aarch64
            if     Sys.islinux()   "libmosek64.so.$(mskvmajor).$(mskvminor)"
            elseif Sys.isapple()   "libmosek64.$(mskvmajor).$(mskvminor).dylib"
            else   error("Unexpected platform for Aarch64")
            end
        else
            error("Architecture not supported $(Sys.ARCH)")
        end

    let moseklibpath = joinpath(path,moseklib)

        if ! isfile(moseklibpath)
            error("Library '$moseklib' not found ($path)")
        else
            moseklibpath
        end
    end
end

function versionFromBindir(bindir ::AbstractString)
    try
        mosekbin = if Sys.iswindows() "mosek.exe" else "mosek" end
        txt = read(`$bindir/$mosekbin`,String)
        m = match(r"\s*MOSEK Version ([0-9]+\.[0-9]+\.[0-9])",txt)
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
        @info("Got version: $ver, expected version: $mskvmajor.$mskvminor")
        ver = split(ver,".")
        return ver[1] == mskvmajor && ver[2] == mskvminor
    else
        return false
    end
end

function collect_output(cmd::Cmd)
    out = Pipe()
    err = Pipe()

    process = run(pipeline(ignorestatus(cmd), stdout=out, stderr=err))

    stdout = @async (read(out))
    stderr = @async (read(err))

    wait(process)
    close(out.in)
    close(err.in)
    (
        process.exitcode,
        fetch(stdout),
        fetch(stderr)
    )
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
            strip(read(f,String))
        end
    catch
        nothing
    end

curmosekbindir =
    try
        open(joinpath(bindepsdir,"mosekbindir"),"r") do f
            strip(read(f,String))
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
            error("MOSEKBINDIR ($mosekbindir) does not point to a MOSEK $(mskvmajor).$(mskvminor) bin directory")
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
        bindirIsCurrentVersion(joinpath(ENV["HOME"],"mosek","$mskvmajor.$mskvminor","tools","platform",mskplatform,"bin"))

        instmethod = "external"

        joinpath(ENV["HOME"],"mosek","$mskvmajor.$mskvminor","tools","platform",mskplatform,"bin")
# 2c. Windows default home
    elseif ! forcedownload &&
        haskey(ENV,"HOMEDRIVE") &&
        haskey(ENV,"HOMEPATH") &&
        bindirIsCurrentVersion(joinpath(string(ENV["HOMEDRIVE"],ENV["HOMEPATH"]),"mosek","$mskvmajor.$mskvminor","tools","platform",mskplatform,"bin"))

        home = string(ENV["HOMEDRIVE"],ENV["HOMEPATH"])

        instmethod = "external"

        joinpath(home,"mosek","$mskvmajor.$mskvminor","tools","platform",mskplatform,"bin")
# 2d. Window global installation
    elseif ! forcedownload && hasregistrykey(WinReg.HKEY_LOCAL_MACHINE,"SOFTWARE\\mosek$mskvmajor$mskvminor","InstallDir")
        instmethod = "external"
        p = getregistrykey(WinReg.HKEY_LOCAL_MACHINE,"SOFTWARE\\MOSEK$mskvmajor$mskvminor","InstallDir")
	# Fix for broken  registry entry:
	p = replace(p,"\\\\" => "\\")
        joinpath(p,"tools","platform",mskplatform,"bin")
# 3. Otherwise, fetch the MOSEK distro and unpack it
    else
        srcdir   = joinpath(bindepsdir,"src")
        dldir    = joinpath(bindepsdir,"downloads")
        archname = "mosektools$(mskplatform)$(distroext)"
        hosturl  = "https://www.mosek.com/downloads/default_dns.txt"

        mkpath(dldir)


        downloadhost = "download.mosek.com"
        #dlcmd = download_cmd(hosturl, joinpath(dldir,"downloadhostname"))
        #@info("Download command: $dlcmd")
        #(res,stdout,stderr) = collect_output(dlcmd)
        #downloadhost =
        #    if res != 0
        #        @error(String(stderr))
        #        error("Failed to get MOSEK download host")
        #    else
        #        open(joinpath(dldir,"downloadhostname"),"r") do f
        #            strip(read(f,String))
        #        end
        #    end

        verurl   = "https://$downloadhost/stable/$mskvmajor.$mskvminor/version"

        cur_version =
            if isfile(joinpath(bindepsdir,"version"))
                open(joinpath(bindepsdir,"version"),"r") do f
                    cur_version = strip(read(f,String))
                end
            else
                nothing
            end

        @info("Get latest MOSEK version ($verurl)")
        success(download_cmd(verurl, joinpath(dldir,"new_version"))) || error("Failed to get MOSEK version")

        new_version =
            open(joinpath(dldir,"new_version"),"r") do f
                strip(read(f,String))
            end
        if cur_version === nothing
            @info("Latest MOSEK version = $new_version, nothing currently installed")
        else
            @info("Latest MOSEK version = $new_version, currently installed = $cur_version")
        end

        instmethod = "internal"





        if  cur_version == nothing ||
            cur_version != new_version ||
            !bindirIsCurrentVersion(joinpath(srcdir,"mosek","$mskvmajor.$mskvminor","tools","platform",mskplatform,"bin"))

            verarr = split(new_version,'.')

            archurl = "https://$downloadhost/stable/$(new_version)/$archname"
            @info("Download MOSEK distro ($archurl)")

            basename,ext,sndext = splittarpath(archname)

            success(download_cmd(archurl, joinpath(dldir,archname))) || error("Failed to download MOSEK distro")

            mkpath(srcdir)
            @info("Unpack MOSEK distro ($dldir/$archname -> $srcdir)")
            success(unpack_cmd(joinpath(dldir,archname),srcdir, ext, sndext)) || error("Failed to unpack MOSEK distro")

            @info("MOSEK installation complete.")
            joinpath(srcdir,"mosek","$mskvmajor.$mskvminor","tools","platform",mskplatform,"bin")
        else
            @info("Update not necessary")
            joinpath(srcdir,"mosek","$mskvmajor.$mskvminor","tools","platform",mskplatform,"bin")
        end
    end

mskbindir = replace(mskbindir,"\\" => "/")


version = versionFromBindir(mskbindir)
@info("mskbindir = $mskbindir")

if version == nothing
    error("MOSEK package is broken")
else
    open(joinpath(bindepsdir,"version"),"w") do f
        write(f,version)
    end
end

verarr = split(version,'.')

libmosekpath= findlibs(mskbindir,verarr[1],verarr[2])

if instmethod == "external"
    @info("""Found MOSEK $version at $mskbindir""")
end

open(joinpath(bindepsdir,"inst_method"),"w") do f
    write(f,instmethod)
end
open(joinpath(bindepsdir,"mosekbindir"),"w") do f
    write(f,mskbindir)
end

libmosekpath = escape_string(normpath(libmosekpath))

open(joinpath(bindepsdir,"deps.jl"),"w") do f
    write(f,"""
# This is an auto-generated file; do not edit
# Macro to load a library
import Libdl
macro checked_lib(libname, path)
    (Libdl.dlopen_e(path) == C_NULL) && error("Unable to load \\n\\n\$libname (\$path)\\n\\nPlease re-run Pkg.build(package), and restart Julia.")
    quote const \$(esc(libname)) = \$path end
end


# Load dependencies
@checked_lib libmosek      "$libmosekpath"
""")


end # open
