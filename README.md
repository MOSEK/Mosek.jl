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

Documentation
-------------

There is not a lot of documentation yet, but Mosek.jl uses the same general
idea the MOSEK Python API which is documented
[here](http://docs.mosek.com/7.0/pythonapi/index.html). Functions have the
same names and arguments but all index arguments are 1-based in MOSEK.jl.

The General Convex interface is not documented at all, but the example 
<tt>nlo1.jl</tt>.

