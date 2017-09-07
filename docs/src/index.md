```@meta
CurrentModule = Mosek
```

# Overview

The `Mosek.jl` package defines a low-level interface to the MOSEK API,
mapping nearly all functions to Julia. The interface is array
oriented, and in most cases each function will change or get one
specific piece of the optimization problem.

## Examples

The package contains a set of examples demonstrating how to set up and
solve basic problems.

- `lo1.jl` Linear problem. This also introduces log printing callbacks and solution handling.
- `cqo1.jl` Second order cone problem.
- `qo1.jl` Problem with quadratic objective.
- `qcqo1.jl` Problem with quadratic objective and constraints.
- `sdo1.jl` Problem with semidefinite variables.
- `milo1.jl` Mixed integer linear problem.
- `nlo1.jl` Introduction to general non-linear interface.
- `callback.jl` Using progress callbacks to get integer solutions.
- `feasrepairex1.jl` Using `primalrepair()` to repair infeasibilities.
- `portfolio.jl` Second order conic problem implementing a portfolio model.
- `production.jl` Demonstrating how to modify and reoptimize a problem.

## Examining problems

Sometimes it is useful to be able to examine visually what is in a
task when debugging a model. `Mosek.jl` includes some functions to
format and examine an entire task or parts of it. This is defined in
the module `Mosek.Ext`.

### Formatting an entire task

The `show()`and `showall()` function can be used to format the task
```julia
using Mosek,Mosek.Ext

t = maketask(filename="myfile.task")
show(t)
```

### Examining parts of a task

The `show()` function can be used to examine certain parts of the task.

```julia
using Mosek,Mosek.Ext

t = maketask(filename="myfile.task")

# examine the objective
show(t[Obj()]

# examine first variable of the task
show(t[Var(1)])

# examine first constraint
show(t[Con(1)]

# examine first cone
show(t[Cone(1)]

# examine first semidefinite variable
show(t[Barvar(1)]

# examine first matrix from the symmetrix matrix store
show(t[Symmat(1)]

# examine first matrix from the symmetrix matrix store
show(t[Symmat(1)]
```

Finally, if there are solutions available, they can be printed as

```julia
sol = t[Sol(MSK_SOL_ITR)]

# examine entire solution (use showall instead to show everything)
show(sol)

# examine solution values for variable 1
show(sol[Var(1)]]

# examine solution values for constraint 1
show(sol[Con(1)]]

# examine solution values for semidefinite variable 1
show(sol[Barvar(1)]]
```

# MOSEK Solver API Reference

## Create and destroy MOSEK objects 
```@docs
makeenv
maketask
deleteenv
deletetask
```
## Callbacks
```@docs
putstreamfunc
putcallbackfunc
```
