using Compat
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

downloadcmd = nothing
function download_cmd(url::AbstractString, filename::AbstractString)
    global downloadcmd
    if downloadcmd === nothing
        for checkcmd in @windows? (:powershell, :curl, :wget, :fetch) : (:curl, :wget, :fetch)
            try
                if success(`$checkcmd --help`)
                    downloadcmd = checkcmd
                    break
                end
            catch
                continue # don't bail if one of these fails
            end
        end
    end
    if downloadcmd == :wget
        return `wget -O $filename $url`
    elseif downloadcmd == :curl
        return `curl -o $filename -L $url`
    elseif downloadcmd == :fetch
        return `fetch -f $filename $url`
    elseif downloadcmd == :powershell
        return `powershell -Command "(new-object net.webclient).DownloadFile(\"$url\", \"$filename\")"`
    else
        error("No download agent available; install curl, wget, or fetch.")
    end
end

@unix_only begin
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

@windows_only begin
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
        else
            error("I don't know how to unpack $file")
        end
    end
end
