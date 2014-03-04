
The MOSEK General Convex (GECO) interface allows for solving problems of the following type

.. math::

        \min\, & c^x + f_0(x)
        s.t.   & a_i^tx + f_i(x) < b_i ,~ i=1\ldots m
               b_l^x \leq x \leq b_u^x

where :math:`f_0` and all :math:`f_i` are convex on the whole of the domain of :math:`x`. 
Additionally, this can be mixed with linear constraints and quadratic terms (but not conic constraints).

The convex functions :math:`f_i` are provided to MOSEK as a set of call-back
functions that compute function values, gradient of the Lagrangian, gradients
of objective and constraints individually, and hessian of the Lagrangian.

There are two functions in the interface:

clearnlcallbacks
================

::

    function clearnlcallbacks(task::MSKtask)

Remove all non-linear callbacks from the problem.


putnlcallbacks
==============

::

    function putnlcallbacks
    ( task::MSKtask,
      grdobjsub :: Array{Int,1},
      grdconsub :: Array{Int,1},
      grdconptr :: Array{Int,1},
      hessubi   :: Array{Int,1},
      hessubj   :: Array{Int,1},
      evalobj   :: Function,
      evalconi  :: Function,
      grdlag    :: Function,
      grdobj    :: Function,
      grdconi   :: Function,
      heslag    :: Function)

This sets up the structure of the non-linear terms and sets the non-linear callback functions.

* ``grdobjsub :: Array{Int,1}`` The subscripts of the variables that appear in non-linear terms in the objective.
* ``grdconsub :: Array{Int,1}`` The subscripts of the variables that appear in constraints. 
* ``grdconptr :: Array{Int,1}`` Defines the positions in `grdconsub` where rows begin, so that `grdconsub[grdconptr[i]]` is the index of the first non-linear variable in constraint `i`. 
* ``hessubi   :: Array{Int,1}`` Row subscripts of the non-zero elements of the Hessian of the Lagrangian. This matrix is symmetrix and only elements in the lower triangular should be inputted.
* ``hessubj   :: Array{Int,1}`` Column subscripts of the non-zero elements of the Hessian of the Lagrangian. This matrix is symmetrix and only elements in the lower triangular should be inputted.
* ``evalobj   :: Function`` Function that evaluates the non-linear part of the objective at a given point. See below.
* ``evalconi  :: Function`` Function that evaluates the non-linear part of a constraint at a given point. See below.
* ``grdlag    :: Function`` Function that evaluates the gradient of the Lagrangian at a given point. See below.
* ``grdobj    :: Function`` Function that evluates the gradient of the objective at a given point. See below.
* ``grdconi   :: Function`` Function that evaluates the gradient of a constraint at a given point. See below.
* ``heslag    :: Function`` Function that evaluates the Hessian of the Lagrangian at a given point. See below.

The Lagrangian of the non-linear part of the problem has the following form:

.. math::

        \mathcal{L}(x,yo,yc) = yo\cdot f_0(x) + \sum^m_{i=1} yc_i\cdot f_i(x)

It is the first and second derivatives of this that should be computed. The following sections show the form of the callback functions. 

evalobj
-------

::

    function evalobj(x::Array{Float64,1}) -> Float64

Evaluate the non-linear part of the objective at the point `x`.


evalconi
--------

::

    function evalconi(x:: Array{Float64,1},i:: Int32) -> Float64

Evaluate the non-linear part of constraint `i` at the point `x`.

grdlag  
------

::

    function grdlag
    ( x   :: Array{Float64,1},
      yo  :: Float64,
      yc  :: Array{Float64,1},
      subi:: Array{Int32,1},
      val :: Array{Float64,1} )

Evaluate the gradient of the Lagrangian

.. math::
       \mbox{yo}\cdot f_0'(x) + \sum_i \mbox{yc}_i f_{\mbox{subi}[i]}'(x)

Evaluate the gradient of the Lagrangian at the point `x`. Here
* ``x`` is the point where the function should be evaluated.
* ``yo`` the multiplier for the objective gradient
* ``subi`` the indexes of the constraints that should be included in the Lagrangian
* ``yc`` the multipliers for the constraints that should be included in the Lagrangian
* ``val`` a vector of length `numvar` where the gradient values are returned. Only the non-zero places should be overwritten.

grdobj  
------

::

    function grdobj
    ( x  :: Array{Float64,1},
      sub:: Array{Int32,1}, 
      val:: Array{Float64,1} )

Evaluate the gradient of :math:`f_0`:

.. math::
      
        \mbox{val}[k] \leftarrow \frac{d}{dx_{\mbox{sub}[k]}} f_0'(x)

* ``x`` is the point where the function should be evaluated
* ``sub`` the variable subscripts corresponding to the non-zero places in the objective
* ``val`` the array that return the gradient values

grdconi 
-------

::

  function grdconi
  ( x  :: Array{Float64,1},
    i  :: Int32, 
    sub:: Array{Int32,1}, 
    val:: Array{Float64,1})

Evaluate the gradient of :math:`f_i`:

.. math::
      
        \mbox{val}[k] \leftarrow \frac{d}{dx_{\mbox{sub}[k]}} f_i'(x)

* ``x`` is the point where the function should be evaluated
* ``i`` the constraint index 
* ``sub`` the variable subscripts corresponding to the non-zero places in the constraint
* ``val`` the array that return the gradient values


heslag 
------ 

::

    function heslag
    ( x ::      Array{Float64,1},
      yo::      Float64,
      yc::      Array{Float64,1},
      subi::    Array{Int32,1},
      hessubi:: Array{Int32,1},
      hessubj:: Array{Int32,1},
      hesval::  Array{Float64,1})

Evaluate the Hessian of the Lagrangian:

.. math::
        \frac{d^2}{dx^2} \mathcal{L}(x,yo,yc)

Note that the Hessian is symmetric. Only elements from the lower triangular part should be inputted, i.e. all elements

.. math::
        
        \frac{d^2}{dx_idx_j} \mathcal{L}(x,yo,yc),~ j\leq i

* ``x`` is the point where the function should be evaluated.
* ``yo`` the multiplier for the objective gradient
* ``subi`` the indexes of the constraints that should be included in the Lagrangian
* ``yc`` the multipliers for the constraints that should be included in the Lagrangian
* ``hessubi`` row subscripts of the non-zeros
* ``hessubj`` column subscripts of the non-zeros
* ``hesval`` non-zero values of the Hessian



