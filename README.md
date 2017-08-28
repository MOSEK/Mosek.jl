Mosek.jl
========

Interface to the Mosek solver in Julia. 

Mosek.jl is a more or less complete mapping of the MOSEK functionality:
- Most MOSEK C API functions are available
- Callbacks for information retrival and log output during optimization
- Interface for the MOSEK general convex solver
- Implementation of the `LinprogSolver` interface and other interfaces for `JuMP` (https://github.com/JuliaOpt/JuMP.jl)

MOSEK can solve LP (linear), SOCP (second order conic), SDP (semi-definite), 
QP (quadratic objective, quadratic constraints), GECO (general
convex) and MIP (mixed-integer problems). These can be mixed as follows:
- LP+SOCP+SDP
- LP+SOCP+MIP
- LP+QP+MIP
- LP+QP+GECO

MOSEK is commercial software, but free licenses are available for academic use. See [here](http://mosek.com/resources/academic-license/) for details.

Installation
------------

Use the Julia package manager to install Mosek.jl:

```julia
Pkg.add("Mosek")
```
    
The `Mosek.jl` package requires the MOSEK distribution binaries run.Upon
installation it will attempt to either local an installed MOSEK or download and
install from the MOSEK website (www.mosek.com):

1. If the environment variable `MOSEKBINDIR` is defined, the installer will assume that this directory contains the necessary libraries. If it does not, the installer will fail.
2. If the current `MOSEK.lj` installation uses a user-defined MOSEK and this is a valid version, this will be used.
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
Pkg.checkout("Mosek","master")
Pkg.build("Mosek")
```

If this also fails, please post an issue in Github.


Documentation
-------------

All functions and constants in the Mosek.jl are briefly documented [here](doc/Mosek-Functions.rst).

For a more complete description of functions, please refer to 
[the MOSEK C API documentation](http://docs.mosek.com/7.0/capi/index.html).

The General Convex interface is not documented at all, but the example 
`nlo1.jl` show the general idea.

`MathProgBase` interface
------------------------

The [MathProgBase](https://github.com/JuliaOpt/MathProgBase.jl) for MOSEK is a separate package called `MathProgBaseMosek`.


`MathOptInterface`
------------------

The [MathOptInterface](https://github.com/JuliaOpt/MathOptInterface.jl) for MOSEK is a separate package called `MathOptInterfaceMosek`.
