using BinDeps

@BinDeps.setup

@unix_only begin
  if WORD_SIZE == 32
    libmosek = library_dependency("libmosek")
  else
    libmosek = library_dependency("libmosek64")
  end
end
@windows_only begin
  if WORD_SIZE == 32
    libmosek = library_dependency("mosek7_0")
  else
    libmosek = library_dependency("mosek64_7_0")
  end
end

@linux_only begin
  if WORD_SIZE == 32
    provides(Binaries, URI("http://download.mosek.com/stable/7/mosektoolslinux32x86.tar.bz2"), [libmosek])
  else
    provides(Binaries, URI("http://download.mosek.com/stable/7/mosektoolslinux64x86.tar.bz2"), [libmosek])
  end
end

@windows_only begin
  if WORD_SIZE == 32
    provides(Binaries, URI("http://download.mosek.com/stable/7/mosektoolswin32x86.zip"), [libmosek])
  else
    provides(Binaries, URI("http://download.mosek.com/stable/7/mosektoolswin64x86.zip"), [libmosek])
  end
end

@osx_only begin 
  provides(Binaries, URI("http://download.mosek.com/stable/7/mosektoolsosx64x86.tar.bz2"), [libmosek])
end

#provides(Binaries, joinpath("mosek","lib"), libmosek)

@BinDeps.install
