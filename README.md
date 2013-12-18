Mosek.jl
========

Interface to the Mosek solver in Julia. 

Mosek.jl is a more or less complete mapping of the MOSEK functionality:
- Most MOSEK C API functions are available
- Callbacks for information retrival and log output during optimization
- Interface for the MOSEK general convex solver (still quite buggy)
- Implementation of the LinprogSolver interface

MOSEK can solve LP (linear), SOCP (second order conic), SDP (semi-definite), 
QP (quadratic objective, quadratic constraints), GECO (general
convex) and MIP (mixed-integer problems). These can be mixed as follows:
- LP+SOCP+SDP
- LP+SOCP+MIP
- LP+QP+MIP
- LP+QP+GECO

MOSEK is commercial software, but free licenses are available for academic 
use. See [here](http://mosek.com/resources/academic-license/) for details.

Installation
------------
Mosek.jl is not yet in the official Julia repository. To install from the development repository use

    Pkg.clone("https://github.com/IainNZ/Mosek.jl.git")

When the Mosek module is loaded it attepts to locate the relevant MOSEK libraries. Mosek.jl will first look 
for an environment variable `MOSEKBINDIR`. If this is present it must point to the `bin/` directory
of the installed MOSEK distro. Otherwise it will look in the default install path:

    $HOME/mosek
    
Usually this resolves to 

    /Users/username
    
on OS X, 

    /home/username
    
on Linux, and

    C:\Users\username
    
for newer Windows versions.

Furthermore, to run an optimization a license is required (these are free for academic use). Mosek will look first for the enironment variable `MOSEKLM_LICENSE_FILE` which must point to the relevant license file, and, if this is not present, it will look for a file called `mosek.lic` in the default install path.


Documentation
-------------

All functions and constants in the Mosek.jl are briefly documented [here](https://github.com/IainNZ/Mosek.jl/wiki/Mosek-Functions).

For a more complete description of functions, please refer to 
[the MOSEK C API documentation](http://docs.mosek.com/7.0/capi/index.html).

The General Convex interface is not documented at all, but the example 
<tt>nlo1.jl</tt> show the general idea.

