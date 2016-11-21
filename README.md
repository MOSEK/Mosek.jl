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

```
Pkg.add("Mosek")
```
    
When installing, the installer will look for the MOSEK library in following places:
- First, if the environment variable `MOSEKBINDIR` is defined, it will look for the MOSEK library in the directory it points to. I.e. it must point the the `bin/` directory in the MOSEK distro.
- If `MOSEKBINDIR` is undefined or it points to an invalid directory,
  the installer will look for the MOSEK distribution in the user's
  home directory. This usually resolves to:
  - OS X: `/Users/username`
  - Linux: `/home/username`
  - Windows: `C:\Users\username`
- If no usable MOSEK installation is found here, the installer will
  attempt to download and unpack the latest distro. In this case doing
  `Pkg.build("Mosek")` will update the MOSEK distro if possible.`

If the MOSEK distro installation directory is moved it is necessary to rebuild the package using


```
Pkg.build("Mosek")
```

Furthermore, to run an optimization a license is required (these are
free for academic use). MOSEK will look first for the enironment
variable `MOSEKLM_LICENSE_FILE` which, if defined, must point to the relevant
license file. If this is not defined, MOSEK will look for a file
called `mosek.lic` in the default install path, e.g.


```
$HOME/mosek/mosek.lic
```

### Updating the Mosek library
If the MOSEK distro was installed manually, it can be updated simply
by installing a newer distro in the same place. Otherwise, doing
`Pkg.build("Mosek")` will check the latest MOSEK distro and update if
possible.

You can see if the MOSEK distro was installed internally this way:

```
is_internal = open(joinpath(Pkg.dir("Mosek"),"deps","inst_method"),"r") do f readstring(f) == "internal" end
```

### When installation does not work
If you experience problems installing (in particular on Windows or OS X), you can try to pull the latest revision and see if that works
```
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

MathProgBase interface
----------------------

Mosek implements the solver-independent
[MathProgBase](https://github.com/JuliaOpt/MathProgBase.jl) interface,
that allows it to works solver in e.g. [JuMP](https://github.com/JuliaOpt/JuMP.jl).

The solver object is called ``MosekSolver``. Parameters are accepted
and mirror the names in the Mosek documentation but the
``MSK_[IDS]PAR`` prefix is optional.

For example, you can suppress output by either saying
``MosekSolver(MSK_IPAR_LOG=0)`` or ``MosekSolver(LOG=0)``, where
``LOG`` corresponds to the `MSK_IPAR_LOG` parameter in the API. When
the prefix is excluded the type of the parameter is inferred by the
provided value.
