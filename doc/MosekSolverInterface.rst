MosekSolverInterface
====================

Mosek.jl implements most of the MathProgSolver interface from ``MathProgBase.jl``<https://github.com/JuliaOpt/MathProgBase.jl>_. Specifically, linear, quadratic and SOCP problems can be solved, and integer variables are supported. 

Linear Problems
---------------
Linear problems should solve reliably and all result information should be as expected (primal and dual solutions, certificates etc.)

Quadratic Problems
------------------
Quadratic non-SOCP problems are by default passed to the solver as quadratic
programs. These are solved by specialized interior point method that can handle
quadratic terms in both objective and constraints. The dual information reported
should be fairly precise; this means that the primal-dual gap between the reported 
primal and dual solution conform to the tolerance defined by the solver parameters.

Second Order Conic Problems
---------------------------
For SOCP problems there are some gotchas:

* A variable can appear in *at most* one cone.
* The primal constraint solution values for conic constraints are meaningless,
  and the values should not be used.
* Dual variable information for variables is not the dual values for the quadratic program.

The latter point may need some explanation: MOSEK solves an SOCP of the form

.. math:: 

         \min\, & c^tx \\
         s.t.   & l^c \leq Ax \leq u^c \\
                & l^x \leq x \leq u^x \\
                & x\in \mathcal{C}

where :math:`\mathcal{C}` is a product of quadratic and rotated quadratic cones. The dual problem is

.. math:: 

         \max\, & (l^c)^ts^c_l - (u^c)^ts^c_u + (l^x)^ts^x_l - (u^x)^ts^x_u \\
         s.t.   & A^ty + s^x_l - s^x_u + s^n_x = c \\
                & -y + s^c_l - s^c_u = 0 \\
                & s^c_l,s^c_u,s^x_l,s^x_u \geq 0 \\
                & s^n_x\in\mathcal{C}^*

where :math:`\mathcal{C}^*` is the dual cone of :math:`\mathcal{C}` (which is self-dual, so they are identical).

In particular, we notice that there are three dual variables for :math:`x`; one for the
lower bound, one for the upper bound and one for the conic bound. When solving an SOCP problem, the dual 
solution reported is 

.. math:: 

        s^x_l - s^x_u + s^x_n

The relations between the primal and the dual problem is covered in more
details in ``The MOSEK
Manual``<http://docs/docs/manuals/dev/capi/Conic_quadratic_optimization__1.html>_.

In some cases MosekSolverInterface will add some extra variables when adding a
conic constraint. When the solution is reported these are omitted since they
are, in a sense, redundant. If a problem is written to a file and then read
from the same file, the information in these extra variables is lost, and their
solution values will be reported as for any other variable.

Mixed Integer Problems
----------------------
Integer variables can be used in linear and SOCP problems, but not with quadratic problems.
