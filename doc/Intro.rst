How to use Mosek.jl
===================

The overall structure of a program using Mosek is, somewhat simplified, this:

::

   using Mosek
   task = maketask()
   
   // set up task
   readdata(task, "somefile.opf")
   // solve
   optimize(task)
   // extract solution
   xx = getxx(task,MSK_SOL_ITR)
   
   println(xx)
   
   deletetask(task)

The constructor ``maketask()`` supports ``do``-block structure, so the above can be written:

::

   using Mosek
   maketask() do task
       // set up task
       readdata(task, "somefile.opf")
       // solve
       optimize(task)
       // extract solution
       xx = getxx(task,MSK_SOL_ITR)

       println(xx)
   end

