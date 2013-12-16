using BinDeps

@BinDeps.setup


if WORD_SIZE == 32
  libmosek = library_dependency("libmosek", os = :Unix)
  libmosek = library_dependency("mosek7_0", os = :Windows)
  provides(Binaries, URI("http://download.mosek.com/stable/7/mosektoolslinux32x86.tar.bz2"), [libmosek], unpacked_dir="mosek", os = :Linux)
  provides(Binaries, URI("http://download.mosek.com/stable/7/mosektoolswin32x86.zip"),       [libmosek], unpacked_dir="mosek", os = :Windows)
else
  libmosek = library_dependency("libmosek64", os = :Unix)
  libmosek = library_dependency("mosek64_7_0", os = :Windows)
  provides(Binaries, URI("http://download.mosek.com/stable/7/mosektoolslinux64x86.tar.bz2"), [libmosek], unpacked_dir="mosek", os = :Linux)
  provides(Binaries, URI("http://download.mosek.com/stable/7/mosektoolswin64x86.zip"),       [libmosek], unpacked_dir="mosek", os = :Windows)
  provides(Binaries, URI("http://download.mosek.com/stable/7/mosektoolsosx64x86.tar.bz2"),   [libmosek], unpacked_dir="mosek", os = :Darwin)
end
  
print("HELLO! BinDeps script\n")



@BinDeps.install
