## Note on versions and release

Since the `Mosek.jl` package is designed to match a specific MOSEK version (major+minor version), there are branches for the different MOSEK versions:
- Branch `b0.8` is compatible with MOSEK 8.0. Not actively updated. 
- Branch `b0.9` is compatible with MOSEK 8.1. Currently updated only for bugfixes.
- Branch `b1.1-msk9.1` is compatible with MOSEK 9.1. Not actively updated.
- Branch `b1.1-msk9.2` is compatible with MOSEK 9.2. Not actively updated.
- Branch `b1.1-msk9.3` is compatible with MOSEK 9.3. Currently updated only for bugfixes.
- Branch `b10.0` is compatible with MOSE 10.0, at the time of writing in beta. Since MOSEK 10.0 the branch names of Mosek.jl will follow Mosek, so the branch compatible with MOSEK 10.0 will be called b10.0. Release will be called `MAJORVER.MINORVER.N`, where `N` is incremented independently of MOSEK. This means that since MOSEK 10.0, Mosek.jl vX.Y.Z will require MOSEK version X.Y
- The `master` branch. This is more or less kept compatible with the latest MOSEK release, either latest stable release or, if available, the latest alpha or beta.

`Mosek.jl` releases are taken from the `b*.*` branches.

# Mosek.jl

Interface to the [MOSEK solver](https://www.mosek.com) in Julia.

Mosek.jl is a more or less complete mapping of the MOSEK functionality:
- Most MOSEK C API functions are available
- Callbacks for information retrival and log output during optimization
- Interface for the MOSEK general convex solver

MOSEK can solve LP (linear), Conic (second order conic, power, exponential
etc.), SDP (semi-definite), QP (quadratic objective, quadratic constraints),
and MIP (mixed-integer problems). These can be mixed as follows:
- LP+CONIC+SDP
- LP+CONIC+MIP

MOSEK is commercial software, but free licenses are available for
academic use. See [here](http://mosek.com/products/academic-licenses/)
for details.

## Installation

Use the Julia package manager to install Mosek.jl:

```julia
Pkg.add("Mosek")
```

The `Mosek.jl` package requires the MOSEK distribution binaries run. Upon
installation it will attempt to either local an installed MOSEK or download and
install from the MOSEK website (www.mosek.com):

1. If the environment variable `MOSEKBINDIR` is defined, the installer will assume that this directory contains the necessary libraries. If it does not, the installer will fail.
2. If the current `Mosek.jl` installation uses a user-defined MOSEK and this is a valid version, this will be used.
3. If MOSEK is installed in the default location in the users HOME directory, and this installation has the correct version, this will be used.
4. If no usable MOSEK installation is found here, the installer will
   attempt to download and unpack the latest distro. In this case doing
   `Pkg.build("Mosek")` will update the MOSEK distro if possible.`

If the MOSEK distro installation directory is moved it is necessary to rebuild the package using
```julia
Pkg.build("Mosek")
```

If you have previously installed `Mosek.jl` using a pre-installed
MOSEK distro, setting the `MOSEKJL_FORCE_DOWNLOAD=YES` will force the
installer to download MOSEK from the web instead of using the old
version.

Note that environment variables can be set temporarily from Julia as
```julia
ENV["MOSEKBINDIR"] = "/home/myname/lib"
```

Furthermore, a license file is required to use MOSEK (these are
free for academic use). MOSEK will look first for the enironment
variable `MOSEKLM_LICENSE_FILE` which, if defined, must point to the relevant
license file. If this is not defined, MOSEK will look for a file
called `mosek.lic` in the default install path, e.g.


```sh
$HOME/mosek/mosek.lic
```

### Updating the Mosek library

If the MOSEK distro was installed manually, it can be updated simply
by installing a newer distro in the same place. Otherwise, doing
`Pkg.build("Mosek")` will check the latest MOSEK distro and update if
possible.

You can see if the MOSEK distro was installed internally this way:

```julia
is_internal = open(joinpath(Pkg.dir("Mosek"),"deps","inst_method"),"r") do f readstring(f) == "internal" end
```

### When installation does not work

If you experience problems installing (in particular on Windows or OS X), you can try to pull the latest revision and see if that works
```julia
Pkg.add(Pkg.PackageSpec(url="https://github.com/mosek/Mosek.jl",rev="master"))
Pkg.build("Mosek")
```

If this also fails, please post an issue in Github.

## Documentation

All functions and constants in the Mosek.jl are briefly documented in docs strings, and an HTML reference can be built using `Documenter.jl`.

For a more complete description of functions, please refer to
[the MOSEK C API documentation](https://www.mosek.com/documentation/).

## Use with JuMP

The [MathOptInterface](https://github.com/jump-dev/MathOptInterface.jl) wrapper
for MOSEK is a separate package called [`MosekTools`](https://github.com/jump-dev/MosekTools.jl).
However, for consistency the optimizer is still named `Mosek.Optimizer`.

Use MOSEK with JuMP as follows:
```julia
using JuMP, MosekTools
model = Model(Mosek.Optimizer)
```

