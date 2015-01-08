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

    Pkg.add("Mosek")
    
When installing, the installer will look for the MOSEK library in following places:
- First, if the environment variable `MOSEKBINDIR` is defined, it will look for the MOSEK library in the directory it points to. I.e. it must point the the `bin/` directory in the MOSEK distro.
- Secondly, it will look for the MOSEK distribution in the user's home directory. This usually resolves to:
  - OS X: `/Users/username`
  - Linux: `/home/username`
  - Windows: `C:\Users\username`
- Thirdly, on Linux and OS X it will attempt to download MOSEK and install it in the Julia configuration directory.

If the MOSEK installation is moved it is necessary to rebuild the package using 

    Pkg.build("Mosek")

Furthermore, to run an optimization a license is required (these are free for academic use). Mosek will look first for the enironment variable `MOSEKLM_LICENSE_FILE` which must point to the relevant license file, and, if this is not present, it will look for a file called `mosek.lic` in the default install path.

### Updating the Mosek library
If the MOSEK distro was installed manually, it can be updated simply by installing a newer distro. 

If the Mosek.jl installer fetched and installed the MOSEK distro, it cannot automatically update, even when doing `Pkg.update()` or `Pkg.build("Mosek")`. Currently the simplest way to force the update is to do
```
Pkg.rm("Mosek") ; Pkg.add("Mosek") ; Pkg.build("Mosek")
```
If this is not desirable (e.g. because other packages depend on Mosek), you can remote the temporary files from the Mosek.jl installation. On Unix this is usually
```
rm -f ~/.julia/Mosek/downloads/*
```


Documentation
-------------

All functions and constants in the Mosek.jl are briefly documented [here](doc/Mosek-Functions.rst).

For a more complete description of functions, please refer to 
[the MOSEK C API documentation](http://docs.mosek.com/7.0/capi/index.html).

The General Convex interface is not documented at all, but the example 
<tt>nlo1.jl</tt> show the general idea.


MathProgBase interface
----------------------

Mosek implements the solver-independent [MathProgBase](https://github.com/JuliaOpt/MathProgBase.jl) interface,
and so can be used within modeling software like [JuMP](https://github.com/JuliaOpt/JuMP.jl).
The solver object is called ``MosekSolver``. Parameters are accepted and mirror the names in the Mosek documentation, except for the ``MSK_[I|D|S]PAR`` prefix.
For example, you can suppress output by saying ``MosekSolver(LOG=0)``, where ``LOG`` corresponds to the [MSK_IPAR_LOG](http://docs.mosek.com/7.0/capi/MSK_IPAR_LOG.html) parameter.
The type of the parameter is inferred by the provided value.

