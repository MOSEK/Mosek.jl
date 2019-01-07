# Blatantly lifted from BinDeps.jl because for some reason the build
# script sometimes cannot find BinDeps.
# @windows_only unpack_cmd has been modified to use zip in addition to 7z.

function splittarpath(path)
    path,extension = splitext(path)
    base_filename,secondary_extension = splitext(path)
    if extension == ".tgz" || extension == ".tbz" || extension == ".zip" && !isempty(secondary_extension)
        base_filename *= secondary_extension
        secondary_extension = ""
    end
    (base_filename,extension,secondary_extension)
end


downloadcmd_candidates = @static if Sys.iswindows() 
  (:powershell, :curl, :wget, :fetch) 
else 
  (:curl, :wget, :fetch) 
end

downloadcmd = nothing
function mk_download_cmd()
    global downloadcmd
    whichcmd = @static if Sys.iswindows() "where" else "which" end
    if downloadcmd === nothing
        for checkcmd in downloadcmd_candidates
            try
                if success(`$whichcmd $checkcmd`)
                    downloadcmd = checkcmd
                    break
                end
            catch
                continue # don't bail if one of these fails
            end
        end
    end
end

function recmkdir(path :: String)    
    if ! isdir(path)
        pname = dirname(path)
        recmkdir(pname)
        if ! isdir(pname)
            error("Failed to recursively create directory")
        end
        mkdir(path)
    end        
end

function download_cmd(url::AbstractString, filename::AbstractString)
    if downloadcmd === nothing mk_download_cmd() end
    if downloadcmd == :wget
        return `wget -O $filename $url`
    elseif downloadcmd == :curl
        return `curl -o $filename -L $url`
    elseif downloadcmd == :fetch
        return `fetch -f $filename $url`
    elseif downloadcmd == :powershell
        recmkdir(dirname(filename))
        return `powershell -file $(joinpath(dirname(abspath(@__FILE__)),"winget.ps0")) $url $filename`
    else
        error("No download agent available; install curl, wget, fetch, or powershell.")
    end
end

@static if Sys.isunix()
    function unpack_cmd(file,directory,extension,secondary_extension)
        if (extension == ".gz" && secondary_extension == ".tar") || extension == ".tgz"
            return (`tar xzf $file --directory=$directory`)
        elseif (extension == ".bz2" && secondary_extension == ".tar") || extension == ".tbz"
            return (`tar xjf $file --directory=$directory`)
        elseif extension == ".xz" && secondary_extension == ".tar"
            return (`unxz -c $file `|>`tar xv --directory=$directory`)
        elseif extension == ".tar"
            return (`tar xf $file --directory=$directory`)
        elseif extension == ".zip"
            return (`unzip -x $file -d $directory`)
        elseif extension == ".gz"
            return (`mkdir $directory` |> `cp $file $directory` |> `gzip -d $directory/$file`)
        end
        error("I don't know how to unpack $file")
    end
end

@static if Sys.iswindows()
    has_7z  = nothing
    has_zip = nothing    
    function unpack_cmd(file,directory,extension,secondary_extension)
        global has_7z
        global has_zip

        if has_7z === nothing
            has_7z  = success(`where 7z`)
            has_zip = success(`where unzip`)
        end

        if has_7z
            if((extension == ".gz" || extension == ".xz" || extension == ".bz2") && secondary_extension == ".tar") ||
                extension == ".tgz" || extension == ".tbz"
                return (`7z x $file -y -so`|>`7z x -si -y -ttar -o$directory`)
            elseif extension == ".zip" || extension == ".7z"
                return (`7z x $file -y -o$directory`)
            else
                error("I don't know how to unpack $file")
            end
        elseif has_zip
            if extension == ".zip"
                return (`unzip $file -o -d$directory`)
            else
                error("I don't know how to unpack $file")
            end
        elseif extension == ".zip"
            rm(directory,recursive=true)
            return (`powershell -file $(joinpath(dirname(abspath(@__FILE__)),"winunzip.ps1")) $file $directory`)
        else
            error("I don't know how to unpack $file")
        end
    end
end
