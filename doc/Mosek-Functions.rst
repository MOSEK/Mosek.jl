.. Contents of this file is generated. Do not edit by hand!
.. MOSEK 7.0.0.102



This page lists all MOSEK functions available from Julia. Please note that the documentation was generated from the documentation for the MOSEK C API, so in some cases it may be slightly invalid. Specifically,
* Index arguments may not be displayed or documented at starting at 1 (indexes in C start at 0). Nevertheless, in the MOSEK Julia API all indexes are 1-based.
* Values that are returned from functions may be documented as if they appeared in the argument list for functions.
* Probably a lot of other stuff. I will be trying to improve all this when I can.

For more verbose descriptions of the individual functions, it is a good idea to look at e.g. the Python documentation at http://docs.mosek.com.

      

Mosek.jl Functions
==================
analyzenames
------------
::

    function analyzenames
    ( task::        MSKtask,
      whichstream:: Int32,
      nametype::    Int32 )


* ``task:: MSKtask`` An optimization task.
* ``whichstream:: Int32`` (`Enum streamtype`_) Index of the stream.
* ``nametype:: Int32`` (`Enum nametype`_) The type of names e.g. valid in MPS or LP files.


Analyze the names and issue an error for the first invalid name.

analyzeproblem
--------------
::

    function analyzeproblem
    ( task::        MSKtask,
      whichstream:: Int32 )


* ``task:: MSKtask`` An optimization task.
* ``whichstream:: Int32`` (`Enum streamtype`_) Index of the stream.


Analyze the data of a task.

analyzesolution
---------------
::

    function analyzesolution
    ( task::        MSKtask,
      whichstream:: Int32,
      whichsol::    Int32 )


* ``task:: MSKtask`` An optimization task.
* ``whichstream:: Int32`` (`Enum streamtype`_) Index of the stream.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.


Print information related to the quality of the solution.

appendbarvars
-------------
::

    function appendbarvars
    ( task:: MSKtask,
      dim::  Array{Int32} )


* ``task:: MSKtask`` An optimization task.
* ``dim:: Array{Int32}`` Dimension of symmetric matrix variables to be added.


Appends a semidefinite  variable of dimension dim to the problem.

appendcone
----------
::

    function appendcone
    ( task::     MSKtask,
      conetype:: Int32,
      conepar::  Float64,
      submem::   Array{Int32} )


* ``task:: MSKtask`` An optimization task.
* ``conetype:: Int32`` (`Enum conetype`_) Specifies the type of the cone.
* ``conepar:: Float64`` This argument is currently not used. Can be set to 0.0.
* ``submem:: Array{Int32}`` Variable subscripts of the members in the cone.


Appends a new cone constraint to the problem.

appendconeseq
-------------
::

    function appendconeseq
    ( task::     MSKtask,
      conetype:: Int32,
      conepar::  Float64,
      nummem::   Int32,
      j::        Int32 )


* ``task:: MSKtask`` An optimization task.
* ``conetype:: Int32`` (`Enum conetype`_) Specifies the type of the cone.
* ``conepar:: Float64`` This argument is currently not used. Can be set to 0.0.
* ``nummem:: Int32`` Dimension of the conic constraint.
* ``j:: Int32`` Index of the first variable in the conic constraint.


Appends a new conic constraint to the problem.

appendconesseq
--------------
::

    function appendconesseq
    ( task::     MSKtask,
      conetype:: Array{Int32},
      conepar::  Array{Float64},
      nummem::   Array{Int32},
      j::        Int32 )


* ``task:: MSKtask`` An optimization task.
* ``conetype:: Array{Int32}`` (`Enum conetype`_) Specifies the type of the cone.
* ``conepar:: Array{Float64}`` This argument is currently not used. Can be set to 0.0.
* ``nummem:: Array{Int32}`` Number of member variables in the cone.
* ``j:: Int32`` Index of the first variable in the first cone to be appended.


Appends a multiple conic constraints to the problem.

appendcons
----------
::

    function appendcons
    ( task:: MSKtask,
      num::  Int32 )


* ``task:: MSKtask`` An optimization task.
* ``num:: Int32`` Number of constraints which should be appended.


Appends a number of constraints to the optimization task.

appendsparsesymmat
------------------
::

    function appendsparsesymmat
    ( task::  MSKtask,
      dim::   Int32,
      subi::  Array{Int32},
      subj::  Array{Int32},
      valij:: Array{Float64} )
    -> idx


* ``task:: MSKtask`` An optimization task.
* ``dim:: Int32`` Dimension of the symmetric matrix that is appended.
* ``subi:: Array{Int32}`` Row subscript in the triplets.
* ``subj:: Array{Int32}`` Column subscripts in the triplets.
* ``valij:: Array{Float64}`` Values of each triplet.
* ``idx:: Int64`` Unique index assigned to inputted matrix.


Appends a general sparse symmetric matrix to the vector E of symmetric matrixes.

appendstat
----------
::

    function appendstat(task:: MSKtask)


* ``task:: MSKtask`` An optimization task.


Appends a record the statistics file.

appendvars
----------
::

    function appendvars
    ( task:: MSKtask,
      num::  Int32 )


* ``task:: MSKtask`` An optimization task.
* ``num:: Int32`` Number of variables which should be appended.


Appends a number of variables to the optimization task.

basiscond
---------
::

    function basiscond(task:: MSKtask)
    -> nrmbasis,nrminvbasis


* ``task:: MSKtask`` An optimization task.
* ``nrmbasis:: Float64`` An estimate for the 1 norm of the basis.
* ``nrminvbasis:: Float64`` An estimate for the 1 norm of the inverse of the basis.


Computes conditioning information for the basis matrix.

bktostr
-------
::

    function bktostr
    ( task:: MSKtask,
      bk::   Int32 )
    -> str


* ``task:: MSKtask`` An optimization task.
* ``bk:: Int32`` (`Enum boundkey`_) Bound key.
* ``str:: String`` String corresponding to the bound key.


Obtains a bound key string identifier.

callbackcodetostr
-----------------
::

    function callbackcodetostr(code:: Int32)
    -> callbackcodestr


* ``code:: Int32`` (`Enum callbackcode`_) A call-back code.
* ``callbackcodestr:: String`` String corresponding to the call-back code.


Obtains a call-back code string identifier.

checkconvexity
--------------
::

    function checkconvexity(task:: MSKtask)


* ``task:: MSKtask`` An optimization task.


Checks if a quadratic optimization problem is convex.

checkmem
--------
::

    function checkmem
    ( task:: MSKtask,
      file:: String,
      line:: Int32 )


* ``task:: MSKtask`` An optimization task.
* ``file:: String`` File from which the function is called.
* ``line:: Int32`` Line in the file from which the function is called.


Checks the memory allocated by the task.

chgbound
--------
::

    function chgbound
    ( task::    MSKtask,
      accmode:: Int32,
      i::       Int32,
      lower::   Int32,
      finite::  Int32,
      value::   Float64 )


* ``task:: MSKtask`` An optimization task.
* ``accmode:: Int32`` (`Enum accmode`_) Defines if operations are performed row-wise (constraint-oriented) or column-wise (variable-oriented).
* ``i:: Int32`` Index of the constraint or variable for which the bounds should be changed.
* ``lower:: Int32`` If non-zero, then the lower bound is changed, otherwise
                            the upper bound is changed.
* ``finite:: Int32`` If non-zero, then the given value is assumed to be finite.
* ``value:: Float64`` New value for the bound.


Changes the bounds for one constraint or variable.

chgconbound
-----------
::

    function chgconbound
    ( task::   MSKtask,
      i::      Int32,
      lower::  Int32,
      finite:: Int32,
      value::  Float64 )


* ``task:: MSKtask`` An optimization task.
* ``i:: Int32`` Index of the constraint for which the bounds should be changed.
* ``lower:: Int32`` If non-zero, then the lower bound is changed, otherwise
                            the upper bound is changed.
* ``finite:: Int32`` If non-zero, then the given value is assumed to be finite.
* ``value:: Float64`` New value for the bound.


Changes the bounds for one constraint.

chgvarbound
-----------
::

    function chgvarbound
    ( task::   MSKtask,
      j::      Int32,
      lower::  Int32,
      finite:: Int32,
      value::  Float64 )


* ``task:: MSKtask`` An optimization task.
* ``j:: Int32`` Index of the variable for which the bounds should be changed.
* ``lower:: Int32`` If non-zero, then the lower bound is changed, otherwise
                            the upper bound is changed.
* ``finite:: Int32`` If non-zero, then the given value is assumed to be finite.
* ``value:: Float64`` New value for the bound.


Changes the bounds for one variable.

commitchanges
-------------
::

    function commitchanges(task:: MSKtask)


* ``task:: MSKtask`` An optimization task.


Commits all cached problem changes.

conetypetostr
-------------
::

    function conetypetostr
    ( task::     MSKtask,
      conetype:: Int32 )
    -> str


* ``task:: MSKtask`` An optimization task.
* ``conetype:: Int32`` (`Enum conetype`_) Specifies the type of the cone.
* ``str:: String`` String corresponding to the cone type.


Obtains a cone type string identifier.

deletesolution
--------------
::

    function deletesolution
    ( task::     MSKtask,
      whichsol:: Int32 )


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.


Undefine a solution and frees the memory it uses.

dualsensitivity
---------------
::

    function dualsensitivity
    ( task:: MSKtask,
      subj:: Array{Int32} )
    -> leftpricej,rightpricej,leftrangej,rightrangej


* ``task:: MSKtask`` An optimization task.
* ``subj:: Array{Int32}`` Index of objective coefficients to analyze.
* ``leftpricej:: Array{Float64}`` Left shadow prices for requested coefficients.
* ``rightpricej:: Array{Float64}`` Right shadow prices for requested coefficients.
* ``leftrangej:: Array{Float64}`` Left range for requested coefficients.
* ``rightrangej:: Array{Float64}`` Right range for requested coefficients.


Performs sensitivity analysis on objective coefficients.

getacol
-------
::

    function getacol
    ( task:: MSKtask,
      j::    Int32 )
    -> nzj,subj,valj


* ``task:: MSKtask`` An optimization task.
* ``j:: Int32`` Index of the column.
* ``nzj:: Int32`` Number of non-zeros in the column obtained.
* ``subj:: Array{Int32}`` Index of the non-zeros in the row obtained.
* ``valj:: Array{Float64}`` Numerical values of the column obtained.


Obtains one column of the linear constraint matrix.

getacolnumnz
------------
::

    function getacolnumnz
    ( task:: MSKtask,
      i::    Int32 )
    -> nzj


* ``task:: MSKtask`` An optimization task.
* ``i:: Int32`` Index of the column.
* ``nzj:: Int32`` Number of non-zeros in the j'th row or column of (A).


Obtains the number of non-zero elements in one column of the linear constraint matrix

getacolslicetrip
----------------
::

    function getacolslicetrip
    ( task::  MSKtask,
      first:: Int32,
      last::  Int32 )
    -> subi,subj,val


* ``task:: MSKtask`` An optimization task.
* ``first:: Int32`` Index of the first column in the sequence.
* ``last:: Int32`` Index of the last column in the sequence plus one.
* ``subi:: Array{Int32}`` Constraint subscripts.
* ``subj:: Array{Int32}`` Column subscripts.
* ``val:: Array{Float64}`` Values.


Obtains a sequence of columns from the coefficient matrix in triplet format.

getaij
------
::

    function getaij
    ( task:: MSKtask,
      i::    Int32,
      j::    Int32 )
    -> aij


* ``task:: MSKtask`` An optimization task.
* ``i:: Int32`` Row index of the coefficient to be returned.
* ``j:: Int32`` Column index of the coefficient to be returned.
* ``aij:: Float64`` Returns the requested coefficient.


Obtains a single coefficient in linear constraint matrix.

getapiecenumnz
--------------
::

    function getapiecenumnz
    ( task::   MSKtask,
      firsti:: Int32,
      lasti::  Int32,
      firstj:: Int32,
      lastj::  Int32 )
    -> numnz


* ``task:: MSKtask`` An optimization task.
* ``firsti:: Int32`` Index of the first row in the rectangular piece.
* ``lasti:: Int32`` Index of the last row plus one in the rectangular piece.
* ``firstj:: Int32`` Index of the first column in the rectangular piece.
* ``lastj:: Int32`` Index of the last column plus one in the rectangular piece.
* ``numnz:: Int32`` Number of non-zero elements in the rectangular piece of the linear constraint matrix.


Obtains the number non-zeros in a rectangular piece of the linear constraint matrix.

getarow
-------
::

    function getarow
    ( task:: MSKtask,
      i::    Int32 )
    -> nzi,subi,vali


* ``task:: MSKtask`` An optimization task.
* ``i:: Int32`` Index of the row or column.
* ``nzi:: Int32`` Number of non-zeros in the row obtained.
* ``subi:: Array{Int32}`` Index of the non-zeros in the row obtained.
* ``vali:: Array{Float64}`` Numerical values of the row obtained.


Obtains one row of the linear constraint matrix.

getarownumnz
------------
::

    function getarownumnz
    ( task:: MSKtask,
      i::    Int32 )
    -> nzi


* ``task:: MSKtask`` An optimization task.
* ``i:: Int32`` Index of the row or column.
* ``nzi:: Int32`` Number of non-zeros in the i'th row of (A).


Obtains the number of non-zero elements in one row of the linear constraint matrix

getarowslicetrip
----------------
::

    function getarowslicetrip
    ( task::  MSKtask,
      first:: Int32,
      last::  Int32 )
    -> subi,subj,val


* ``task:: MSKtask`` An optimization task.
* ``first:: Int32`` Index of the first row or column in the sequence.
* ``last:: Int32`` Index of the last row or column in the sequence plus one.
* ``subi:: Array{Int32}`` Constraint subscripts.
* ``subj:: Array{Int32}`` Column subscripts.
* ``val:: Array{Float64}`` Values.


Obtains a sequence of rows from the coefficient matrix in triplet format.

getaslice
---------
::

    function getaslice
    ( task::    MSKtask,
      accmode:: Int32,
      first::   Int32,
      last::    Int32 )
    -> ptrb,ptre,sub,val


* ``task:: MSKtask`` An optimization task.
* ``accmode:: Int32`` (`Enum accmode`_) Defines whether a column slice or a row slice is requested.
* ``first:: Int32`` Index of the first row or column in the sequence.
* ``last:: Int32`` Index of the last row or column in the sequence plus one.
* ``ptrb:: Array{Int64}`` Row or column start pointers.
* ``ptre:: Array{Int64}`` Row or column end pointers.
* ``sub:: Array{Int32}`` Contains the row or column subscripts.
* ``val:: Array{Float64}`` Contains the coefficient values.


Obtains a sequence of rows or columns from the coefficient matrix.

getaslicenumnz
--------------
::

    function getaslicenumnz
    ( task::    MSKtask,
      accmode:: Int32,
      first::   Int32,
      last::    Int32 )
    -> numnz


* ``task:: MSKtask`` An optimization task.
* ``accmode:: Int32`` (`Enum accmode`_) Defines whether non-zeros are counted in a column slice or a row slice.
* ``first:: Int32`` Index of the first row or column in the sequence.
* ``last:: Int32`` Index of the last row or column plus one in the sequence.
* ``numnz:: Int64`` Number of non-zeros in the slice.


Obtains the number of non-zeros in a slice of rows or columns of the coefficient matrix.

getbarablocktriplet
-------------------
::

    function getbarablocktriplet(task:: MSKtask)
    -> num,subi,subj,subk,subl,valijkl


* ``task:: MSKtask`` An optimization task.
* ``num:: Int64`` Number of elements in the block triplet form.
* ``subi:: Array{Int32}`` Constraint index.
* ``subj:: Array{Int32}`` Symmetric matrix variable index.
* ``subk:: Array{Int32}`` Block row index.
* ``subl:: Array{Int32}`` Block column index.
* ``valijkl:: Array{Float64}`` A list indexes   of the elements from symmetric matrix storage that appears in the weighted sum.


Obtains barA in block triplet form.

getbaraidx
----------
::

    function getbaraidx
    ( task:: MSKtask,
      idx::  Int64 )
    -> i,j,num,sub,weights


* ``task:: MSKtask`` An optimization task.
* ``idx:: Int64`` Position of the element in the vectorized form.
* ``i:: Int32`` Row index of the element at position idx.
* ``j:: Int32`` Column index of the element at position idx.
* ``num:: Int64`` Number of terms in weighted sum that forms the element.
* ``sub:: Array{Int64}`` A list indexes   of the elements from symmetric matrix storage that appears in the weighted sum.
* ``weights:: Array{Float64}`` The weights associated with each term in the weighted sum.


Obtains information about an element barA.

getbaraidxij
------------
::

    function getbaraidxij
    ( task:: MSKtask,
      idx::  Int64 )
    -> i,j


* ``task:: MSKtask`` An optimization task.
* ``idx:: Int64`` Position of the element in the vectorized form.
* ``i:: Int32`` Row index of the element at position idx.
* ``j:: Int32`` Column index of the element at position idx.


Obtains information about an element barA.

getbaraidxinfo
--------------
::

    function getbaraidxinfo
    ( task:: MSKtask,
      idx::  Int64 )
    -> num


* ``task:: MSKtask`` An optimization task.
* ``idx:: Int64`` The internal position of the element that should be obtained information for.
* ``num:: Int64`` Number of terms in the weighted sum that forms the specified element in barA.


Obtains the number terms in the weighted sum that forms a particular element in barA.

getbarasparsity
---------------
::

    function getbarasparsity(task:: MSKtask)
    -> numnz,idxij


* ``task:: MSKtask`` An optimization task.
* ``numnz:: Int64`` Number of nonzero elements in barA.
* ``idxij:: Array{Int64}`` Position of each nonzero element in the vector representation of barA.


Obtains the sparsity pattern of the barA matrix.

getbarcblocktriplet
-------------------
::

    function getbarcblocktriplet(task:: MSKtask)
    -> num,subj,subk,subl,valijkl


* ``task:: MSKtask`` An optimization task.
* ``num:: Int64`` Number of elements in the block triplet form.
* ``subj:: Array{Int32}`` Symmetric matrix variable index.
* ``subk:: Array{Int32}`` Block row index.
* ``subl:: Array{Int32}`` Block column index.
* ``valijkl:: Array{Float64}`` A list indexes   of the elements from symmetric matrix storage that appears in the weighted sum.


Obtains barc in block triplet form.

getbarcidx
----------
::

    function getbarcidx
    ( task:: MSKtask,
      idx::  Int64 )
    -> j,num,sub,weights


* ``task:: MSKtask`` An optimization task.
* ``idx:: Int64`` Index of the element that should be obtained information about.
* ``j:: Int32`` Row index in barc.
* ``num:: Int64`` Number of terms in the weighted sum.
* ``sub:: Array{Int64}`` Elements appearing the weighted sum.
* ``weights:: Array{Float64}`` Weights of terms in the weighted sum.


Obtains information about an element in barc.

getbarcidxinfo
--------------
::

    function getbarcidxinfo
    ( task:: MSKtask,
      idx::  Int64 )
    -> num


* ``task:: MSKtask`` An optimization task.
* ``idx:: Int64`` Index of element that should be obtained information about. The value is an index of a symmetric sparse variable.
* ``num:: Int64`` Number of terms that appears in weighted that forms the requested element.


Obtains information about an element in barc.

getbarcidxj
-----------
::

    function getbarcidxj
    ( task:: MSKtask,
      idx::  Int64 )
    -> j


* ``task:: MSKtask`` An optimization task.
* ``idx:: Int64`` Index of the element that should be obtained information about.
* ``j:: Int32`` Row index in barc.


Obtains the row index of an element in barc.

getbarcsparsity
---------------
::

    function getbarcsparsity(task:: MSKtask)
    -> numnz,idxj


* ``task:: MSKtask`` An optimization task.
* ``numnz:: Int64`` Number of nonzero elements in barc.
* ``idxj:: Array{Int64}`` Internal positions of the nonzeros elements in barc.


Get the positions of the nonzero elements in barc.

getbarsj
--------
::

    function getbarsj
    ( task::     MSKtask,
      whichsol:: Int32,
      j::        Int32 )
    -> barsj


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``j:: Int32`` Index of the semidefinite variable.
* ``barsj:: Array{Float64}`` Value of the j'th variable of barx.


Obtains the dual solution for a semidefinite variable.

getbarvarname
-------------
::

    function getbarvarname
    ( task:: MSKtask,
      i::    Int32 )
    -> name


* ``task:: MSKtask`` An optimization task.
* ``i:: Int32`` Index.
* ``name:: String`` The requested name is copied to this buffer.


Obtains a name of a semidefinite variable.

getbarvarnameindex
------------------
::

    function getbarvarnameindex
    ( task::     MSKtask,
      somename:: String )
    -> asgn,index


* ``task:: MSKtask`` An optimization task.
* ``somename:: String`` The requested name is copied to this buffer.
* ``asgn:: Int32`` Is non-zero if name somename is assigned to a semidefinite variable.
* ``index:: Int32`` If the name somename is assigned to a semidefinite variable, then index is the name of the constraint.


Obtains the index of name of semidefinite variable.

getbarvarnamelen
----------------
::

    function getbarvarnamelen
    ( task:: MSKtask,
      i::    Int32 )
    -> len


* ``task:: MSKtask`` An optimization task.
* ``i:: Int32`` Index.
* ``len:: Int32`` Returns the length of the indicated name.


Obtains the length of a name of a semidefinite variable.

getbarxj
--------
::

    function getbarxj
    ( task::     MSKtask,
      whichsol:: Int32,
      j::        Int32 )
    -> barxj


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``j:: Int32`` Index of the semidefinite variable.
* ``barxj:: Array{Float64}`` Value of the j'th variable of barx.


Obtains the primal solution for a semidefinite variable.

getbound
--------
::

    function getbound
    ( task::    MSKtask,
      accmode:: Int32,
      i::       Int32 )
    -> bk,bl,bu


* ``task:: MSKtask`` An optimization task.
* ``accmode:: Int32`` (`Enum accmode`_) Defines if operations are performed row-wise (constraint-oriented) or column-wise (variable-oriented).
* ``i:: Int32`` Index of the constraint or variable for which the bound information should be obtained.
* ``bk:: Int32`` Bound keys.
* ``bl:: Float64`` Values for lower bounds.
* ``bu:: Float64`` Values for upper bounds.


Obtains bound information for one constraint or variable.

getboundslice
-------------
::

    function getboundslice
    ( task::    MSKtask,
      accmode:: Int32,
      first::   Int32,
      last::    Int32 )
    -> bk,bl,bu


* ``task:: MSKtask`` An optimization task.
* ``accmode:: Int32`` (`Enum accmode`_) Defines if operations are performed row-wise (constraint-oriented) or column-wise (variable-oriented).
* ``first:: Int32`` First index in the sequence.
* ``last:: Int32`` Last index plus 1 in the sequence.
* ``bk:: Array{Int32}`` (`Enum boundkey`_) Bound keys.
* ``bl:: Array{Float64}`` Values for lower bounds.
* ``bu:: Array{Float64}`` Values for upper bounds.


Obtains bounds information for a sequence of variables or constraints.

getc
----
::

    function getc(task:: MSKtask)
    -> c


* ``task:: MSKtask`` An optimization task.
* ``c:: Array{Float64}`` Linear terms of the objective as a dense vector. The lengths is the number of variables.


Obtains all objective coefficients.

getcfix
-------
::

    function getcfix(task:: MSKtask)
    -> cfix


* ``task:: MSKtask`` An optimization task.
* ``cfix:: Float64`` Fixed term in the objective.


Obtains the fixed term in the objective.

getcj
-----
::

    function getcj
    ( task:: MSKtask,
      j::    Int32 )
    -> cj


* ``task:: MSKtask`` An optimization task.
* ``j:: Int32`` Index of the variable for which c coefficient should be obtained.
* ``cj:: Float64`` The c coefficient value.


Obtains one coefficient of c.

getconbound
-----------
::

    function getconbound
    ( task:: MSKtask,
      i::    Int32 )
    -> bk,bl,bu


* ``task:: MSKtask`` An optimization task.
* ``i:: Int32`` Index of the constraint for which the bound information should be obtained.
* ``bk:: Int32`` Bound keys.
* ``bl:: Float64`` Values for lower bounds.
* ``bu:: Float64`` Values for upper bounds.


Obtains bound information for one constraint.

getconboundslice
----------------
::

    function getconboundslice
    ( task::  MSKtask,
      first:: Int32,
      last::  Int32 )
    -> bk,bl,bu


* ``task:: MSKtask`` An optimization task.
* ``first:: Int32`` First index in the sequence.
* ``last:: Int32`` Last index plus 1 in the sequence.
* ``bk:: Array{Int32}`` (`Enum boundkey`_) Bound keys.
* ``bl:: Array{Float64}`` Values for lower bounds.
* ``bu:: Array{Float64}`` Values for upper bounds.


Obtains bounds information for a slice of the constraints.

getcone
-------
::

    function getcone
    ( task:: MSKtask,
      k::    Int32 )
    -> conetype,conepar,nummem,submem


* ``task:: MSKtask`` An optimization task.
* ``k:: Int32`` Index of the cone constraint.
* ``conetype:: Int32`` Specifies the type of the cone.
* ``conepar:: Float64`` This argument is currently not used. Can be set to 0.0.
* ``nummem:: Int32`` Number of member variables in the cone.
* ``submem:: Array{Int32}`` Variable subscripts of the members in the cone.


Obtains a conic constraint.

getconeinfo
-----------
::

    function getconeinfo
    ( task:: MSKtask,
      k::    Int32 )
    -> conetype,conepar,nummem


* ``task:: MSKtask`` An optimization task.
* ``k:: Int32`` Index of the conic constraint.
* ``conetype:: Int32`` Specifies the type of the cone.
* ``conepar:: Float64`` This argument is currently not used. Can be set to 0.0.
* ``nummem:: Int32`` Number of member variables in the cone.


Obtains information about a conic constraint.

getconename
-----------
::

    function getconename
    ( task:: MSKtask,
      i::    Int32 )
    -> name


* ``task:: MSKtask`` An optimization task.
* ``i:: Int32`` Index.
* ``name:: String`` Is assigned the required name.


Obtains a name of a cone.

getconenameindex
----------------
::

    function getconenameindex
    ( task::     MSKtask,
      somename:: String )
    -> asgn,index


* ``task:: MSKtask`` An optimization task.
* ``somename:: String`` The name which should be checked.
* ``asgn:: Int32`` Is non-zero if name somename is assigned to a cone.
* ``index:: Int32`` If the name somename is assigned to a cone, then index is the name of the cone.


Checks whether the name somename has been assigned  to any cone.

getconenamelen
--------------
::

    function getconenamelen
    ( task:: MSKtask,
      i::    Int32 )
    -> len


* ``task:: MSKtask`` An optimization task.
* ``i:: Int32`` Index.
* ``len:: Int32`` Returns the length of the indicated name.


Obtains the length of a name of a cone.

getconname
----------
::

    function getconname
    ( task:: MSKtask,
      i::    Int32 )
    -> name


* ``task:: MSKtask`` An optimization task.
* ``i:: Int32`` Index.
* ``name:: String`` Is assigned the required name.


Obtains a name of a constraint.

getconnameindex
---------------
::

    function getconnameindex
    ( task::     MSKtask,
      somename:: String )
    -> asgn,index


* ``task:: MSKtask`` An optimization task.
* ``somename:: String`` The name which should be checked.
* ``asgn:: Int32`` Is non-zero if name somename is assigned to a constraint.
* ``index:: Int32`` If the name somename is assigned to a constraint, then index is the name of the constraint.


Checks whether the name somename has been assigned  to any constraint.

getconnamelen
-------------
::

    function getconnamelen
    ( task:: MSKtask,
      i::    Int32 )
    -> len


* ``task:: MSKtask`` An optimization task.
* ``i:: Int32`` Index.
* ``len:: Int32`` Returns the length of the indicated name.


Obtains the length of a name of a constraint variable.

getcslice
---------
::

    function getcslice
    ( task::  MSKtask,
      first:: Int32,
      last::  Int32 )
    -> c


* ``task:: MSKtask`` An optimization task.
* ``first:: Int32`` First index in the sequence.
* ``last:: Int32`` Last index plus 1 in the sequence.
* ``c:: Array{Float64}`` Linear terms of the objective as a dense vector. The lengths is the number of variables.


Obtains a sequence of coefficients from the objective.

getdbi
------
::

    function getdbi
    ( task::     MSKtask,
      whichsol:: Int32,
      accmode::  Int32,
      sub::      Array{Int32} )
    -> dbi


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``accmode:: Int32`` (`Enum accmode`_) Defines whether sub contains constraint or variable indexes.
* ``sub:: Array{Int32}`` Indexes of constraints or variables.
* ``dbi:: Array{Float64}`` Dual bound infeasibility.


Deprecated.

getdcni
-------
::

    function getdcni
    ( task::     MSKtask,
      whichsol:: Int32,
      sub::      Array{Int32} )
    -> dcni


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``sub:: Array{Int32}`` Constraint indexes to calculate equation infeasibility for.
* ``dcni:: Array{Float64}`` Dual cone infeasibility.


Deprecated.

getdeqi
-------
::

    function getdeqi
    ( task::      MSKtask,
      whichsol::  Int32,
      accmode::   Int32,
      sub::       Array{Int32},
      normalize:: Int32 )
    -> deqi


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``accmode:: Int32`` (`Enum accmode`_) Defines whether equation infeasibilities for constraints or for variables are retrieved.
* ``sub:: Array{Int32}`` Indexes of constraints or variables.
* ``deqi:: Array{Float64}`` Dual equation infeasibilities corresponding to constraints or variables.
* ``normalize:: Int32`` If non-zero, normalize with largest absolute value of input data.


Deprecated.

getdimbarvarj
-------------
::

    function getdimbarvarj
    ( task:: MSKtask,
      j::    Int32 )
    -> dimbarvarj


* ``task:: MSKtask`` An optimization task.
* ``j:: Int32`` Index of the semidefinite variable whose dimension is requested.
* ``dimbarvarj:: Int32`` The dimension of the j'th semidefinite variable.


Obtains the dimension of a symmetric matrix variable.

getdouinf
---------
::

    function getdouinf
    ( task::      MSKtask,
      whichdinf:: Int32 )
    -> dvalue


* ``task:: MSKtask`` An optimization task.
* ``whichdinf:: Int32`` (`Enum dinfitem`_) A double float information item.
* ``dvalue:: Float64`` The value of the required double information item.


Obtains a double information item.

getdouparam
-----------
::

    function getdouparam
    ( task::  MSKtask,
      param:: Int32 )
    -> parvalue


* ``task:: MSKtask`` An optimization task.
* ``param:: Int32`` (`Enum dparam`_) Which parameter.
* ``parvalue:: Float64`` Parameter value.


Obtains a double parameter.

getdualobj
----------
::

    function getdualobj
    ( task::     MSKtask,
      whichsol:: Int32 )
    -> dualobj


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``dualobj:: Float64`` Objective value corresponding to the dual solution.


Computes the dual objective value associated with the solution.

getdviolbarvar
--------------
::

    function getdviolbarvar
    ( task::     MSKtask,
      whichsol:: Int32,
      sub::      Array{Int32} )
    -> viol


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``sub:: Array{Int32}`` An array of indexes of barx variables.
* ``viol:: Array{Float64}`` List of violations corresponding to sub.


Computes the violation of dual solution for a set of barx variables.

getdviolcon
-----------
::

    function getdviolcon
    ( task::     MSKtask,
      whichsol:: Int32,
      sub::      Array{Int32} )
    -> viol


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``sub:: Array{Int32}`` An array of indexes of constraints.
* ``viol:: Array{Float64}`` List of violations corresponding to sub.


Computes the violation of a dual solution associated with a set of constraints.

getdviolcones
-------------
::

    function getdviolcones
    ( task::     MSKtask,
      whichsol:: Int32,
      sub::      Array{Int32} )
    -> viol


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``sub:: Array{Int32}`` An array of indexes of barx variables.
* ``viol:: Array{Float64}`` List of violations corresponding to sub.


Computes the violation of a solution for set of dual conic constraints.

getdviolvar
-----------
::

    function getdviolvar
    ( task::     MSKtask,
      whichsol:: Int32,
      sub::      Array{Int32} )
    -> viol


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``sub:: Array{Int32}`` An array of indexes of x variables.
* ``viol:: Array{Float64}`` List of violations corresponding to sub.


Computes the violation of a dual solution associated with a set of x variables.

getinfeasiblesubproblem
-----------------------
::

    function getinfeasiblesubproblem
    ( task::     MSKtask,
      whichsol:: Int32 )
    -> inftask


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Which solution to use when determining the infeasible subproblem.
* ``inftask:: MSKtask`` A new task containing the infeasible subproblem.


Obtains an infeasible sub problem.

getinfname
----------
::

    function getinfname
    ( task::     MSKtask,
      inftype::  Int32,
      whichinf:: Int32 )
    -> infname


* ``task:: MSKtask`` An optimization task.
* ``inftype:: Int32`` (`Enum inftype`_) Type of the information item.
* ``whichinf:: Int32`` An information item.
* ``infname:: String`` Name of the information item.


Obtains the name of an information item.

getinti
-------
::

    function getinti
    ( task::     MSKtask,
      whichsol:: Int32,
      sub::      Array{Int32} )
    -> inti


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``sub:: Array{Int32}`` Variable indexes for which to calculate the integer infeasibility.
* ``inti:: Array{Float64}`` Integer infeasibility.


Deprecated.

getintinf
---------
::

    function getintinf
    ( task::      MSKtask,
      whichiinf:: Int32 )
    -> ivalue


* ``task:: MSKtask`` An optimization task.
* ``whichiinf:: Int32`` (`Enum iinfitem`_) Specifies an information item.
* ``ivalue:: Int32`` The value of the required integer information item.


Obtains an integer information item.

getintparam
-----------
::

    function getintparam
    ( task::  MSKtask,
      param:: Int32 )
    -> parvalue


* ``task:: MSKtask`` An optimization task.
* ``param:: Int32`` (`Enum iparam`_) Which parameter.
* ``parvalue:: Int32`` Parameter value.


Obtains an integer parameter.

getlenbarvarj
-------------
::

    function getlenbarvarj
    ( task:: MSKtask,
      j::    Int32 )
    -> lenbarvarj


* ``task:: MSKtask`` An optimization task.
* ``j:: Int32`` Index of the semidefinite variable whose length if requested.
* ``lenbarvarj:: Int64`` Number of scalar elements in the lower triangular part of the semidefinite variable.


Obtains the length if the j'th semidefinite variables.

getlintinf
----------
::

    function getlintinf
    ( task::       MSKtask,
      whichliinf:: Int32 )
    -> ivalue


* ``task:: MSKtask`` An optimization task.
* ``whichliinf:: Int32`` (`Enum liinfitem`_) Specifies an information item.
* ``ivalue:: Int64`` The value of the required integer information item.


Obtains an integer information item.

getmaxnumanz
------------
::

    function getmaxnumanz(task:: MSKtask)
    -> maxnumanz


* ``task:: MSKtask`` An optimization task.
* ``maxnumanz:: Int64`` Number of preallocated non-zero linear matrix elements.


Obtains number of preallocated non-zeros in the linear constraint matrix.

getmaxnumbarvar
---------------
::

    function getmaxnumbarvar(task:: MSKtask)
    -> maxnumbarvar


* ``task:: MSKtask`` An optimization task.
* ``maxnumbarvar:: Int32`` Obtains maximum number of semidefinite variable currently allowed.


Obtains the number of semidefinite variables.

getmaxnumcon
------------
::

    function getmaxnumcon(task:: MSKtask)
    -> maxnumcon


* ``task:: MSKtask`` An optimization task.
* ``maxnumcon:: Int32`` Number of preallocated constraints in the optimization task.


Obtains the number of preallocated constraints in the optimization task.

getmaxnumcone
-------------
::

    function getmaxnumcone(task:: MSKtask)
    -> maxnumcone


* ``task:: MSKtask`` An optimization task.
* ``maxnumcone:: Int32`` Number of preallocated conic constraints in the optimization task.


Obtains the number of preallocated cones in the optimization task.

getmaxnumqnz
------------
::

    function getmaxnumqnz(task:: MSKtask)
    -> maxnumqnz


* ``task:: MSKtask`` An optimization task.
* ``maxnumqnz:: Int64`` Number of non-zero elements preallocated in quadratic coefficient matrices.


Obtains the number of preallocated non-zeros for all quadratic terms in objective and constraints.

getmaxnumvar
------------
::

    function getmaxnumvar(task:: MSKtask)
    -> maxnumvar


* ``task:: MSKtask`` An optimization task.
* ``maxnumvar:: Int32`` Number of preallocated variables in the optimization task.


Obtains the maximum number variables allowed.

getmemusage
-----------
::

    function getmemusage(task:: MSKtask)
    -> meminuse,maxmemuse


* ``task:: MSKtask`` An optimization task.
* ``meminuse:: Int64`` Amount of memory currently used by the task.
* ``maxmemuse:: Int64`` Maximum amount of memory used by the task until now.


Obtains information about the amount of memory used by a task.

getnadouinf
-----------
::

    function getnadouinf
    ( task::      MSKtask,
      whichdinf:: String )
    -> dvalue


* ``task:: MSKtask`` An optimization task.
* ``whichdinf:: String`` A double float information item.
* ``dvalue:: Float64`` The value of the required double information item.


Obtains a double information item.

getnadouparam
-------------
::

    function getnadouparam
    ( task::      MSKtask,
      paramname:: String )
    -> parvalue


* ``task:: MSKtask`` An optimization task.
* ``paramname:: String`` Name of a parameter.
* ``parvalue:: Float64`` Parameter value.


Obtains a double parameter.

getnaintinf
-----------
::

    function getnaintinf
    ( task::        MSKtask,
      infitemname:: String )
    -> ivalue


* ``task:: MSKtask`` An optimization task.
* ``infitemname:: String`` <no description>
* ``ivalue:: Int32`` The value of the required integer information item.


Obtains an integer information item.

getnaintparam
-------------
::

    function getnaintparam
    ( task::      MSKtask,
      paramname:: String )
    -> parvalue


* ``task:: MSKtask`` An optimization task.
* ``paramname:: String`` Name of a parameter.
* ``parvalue:: Int32`` Parameter value.


Obtains an integer parameter.

getnastrparam
-------------
::

    function getnastrparam
    ( task::      MSKtask,
      paramname:: String,
      maxlen::    Int32 )
    -> len,parvalue


* ``task:: MSKtask`` An optimization task.
* ``paramname:: String`` Name of a parameter.
* ``maxlen:: Int32`` Length of the name parvalue buffer.
* ``len:: Int32`` Returns the length of the parameter value.
* ``parvalue:: String`` Parameter value.


Obtains a string parameter.

getnumanz
---------
::

    function getnumanz(task:: MSKtask)
    -> numanz


* ``task:: MSKtask`` An optimization task.
* ``numanz:: Int32`` Number of non-zero elements in the linear constraint matrix.


Obtains the number of non-zeros in the coefficient matrix.

getnumanz64
-----------
::

    function getnumanz64(task:: MSKtask)
    -> numanz


* ``task:: MSKtask`` An optimization task.
* ``numanz:: Int64`` Number of non-zero elements in the linear constraint matrix.


Obtains the number of non-zeros in the coefficient matrix.

getnumbarablocktriplets
-----------------------
::

    function getnumbarablocktriplets(task:: MSKtask)
    -> num


* ``task:: MSKtask`` An optimization task.
* ``num:: Int64`` Number elements in the block triplet form of bara.


Obtains an upper bound on the number of scalar elements in the block triplet form of bara.

getnumbaranz
------------
::

    function getnumbaranz(task:: MSKtask)
    -> nz


* ``task:: MSKtask`` An optimization task.
* ``nz:: Int64`` The number of nonzero block elements in barA.


Get the number of nonzero elements in barA.

getnumbarcblocktriplets
-----------------------
::

    function getnumbarcblocktriplets(task:: MSKtask)
    -> num


* ``task:: MSKtask`` An optimization task.
* ``num:: Int64`` An upper bound on the number elements in the block trip let form of barc.


Obtains an upper bound on the number of elements in the block triplet form of barc.

getnumbarcnz
------------
::

    function getnumbarcnz(task:: MSKtask)
    -> nz


* ``task:: MSKtask`` An optimization task.
* ``nz:: Int64`` The number of nonzero elements in barc.


Obtains the number of nonzero elements in barc.

getnumbarvar
------------
::

    function getnumbarvar(task:: MSKtask)
    -> numbarvar


* ``task:: MSKtask`` An optimization task.
* ``numbarvar:: Int32`` Number of semidefinite variable in the problem.


Obtains the number of semidefinite variables.

getnumcon
---------
::

    function getnumcon(task:: MSKtask)
    -> numcon


* ``task:: MSKtask`` An optimization task.
* ``numcon:: Int32`` Number of constraints.


Obtains the number of constraints.

getnumcone
----------
::

    function getnumcone(task:: MSKtask)
    -> numcone


* ``task:: MSKtask`` An optimization task.
* ``numcone:: Int32`` Number conic constraints.


Obtains the number of cones.

getnumconemem
-------------
::

    function getnumconemem
    ( task:: MSKtask,
      k::    Int32 )
    -> nummem


* ``task:: MSKtask`` An optimization task.
* ``k:: Int32`` Index of the cone.
* ``nummem:: Int32`` Number of member variables in the cone.


Obtains the number of members in a cone.

getnumintvar
------------
::

    function getnumintvar(task:: MSKtask)
    -> numintvar


* ``task:: MSKtask`` An optimization task.
* ``numintvar:: Int32`` Number of integer variables.


Obtains the number of integer-constrained variables.

getnumparam
-----------
::

    function getnumparam
    ( task::    MSKtask,
      partype:: Int32 )
    -> numparam


* ``task:: MSKtask`` An optimization task.
* ``partype:: Int32`` (`Enum parametertype`_) Parameter type.
* ``numparam:: Int32`` Returns the number of parameters of the requested type.


Obtains the number of parameters of a given type.

getnumqconknz
-------------
::

    function getnumqconknz
    ( task:: MSKtask,
      k::    Int32 )
    -> numqcnz


* ``task:: MSKtask`` An optimization task.
* ``k:: Int32`` Index of the constraint for which the number of non-zero quadratic terms should be obtained.
* ``numqcnz:: Int32`` Number of quadratic terms.


Obtains the number of non-zero quadratic terms in a constraint.

getnumqconknz64
---------------
::

    function getnumqconknz64
    ( task:: MSKtask,
      k::    Int32 )
    -> numqcnz


* ``task:: MSKtask`` An optimization task.
* ``k:: Int32`` Index of the constraint for which the number quadratic terms should be obtained.
* ``numqcnz:: Int64`` Number of quadratic terms.


Obtains the number of non-zero quadratic terms in a constraint.

getnumqobjnz
------------
::

    function getnumqobjnz(task:: MSKtask)
    -> numqonz


* ``task:: MSKtask`` An optimization task.
* ``numqonz:: Int64`` Number of non-zero elements in the quadratic objective terms.


Obtains the number of non-zero quadratic terms in the objective.

getnumsymmat
------------
::

    function getnumsymmat(task:: MSKtask)
    -> num


* ``task:: MSKtask`` An optimization task.
* ``num:: Int64`` Returns the number of symmetric sparse matrixes.


Get the number of symmetric matrixes stored.

getnumvar
---------
::

    function getnumvar(task:: MSKtask)
    -> numvar


* ``task:: MSKtask`` An optimization task.
* ``numvar:: Int32`` Number of variables.


Obtains the number of variables.

getobjname
----------
::

    function getobjname(task:: MSKtask)
    -> objname


* ``task:: MSKtask`` An optimization task.
* ``objname:: String`` Assigned the objective name.


Obtains the name assigned to the objective function.

getobjnamelen
-------------
::

    function getobjnamelen(task:: MSKtask)
    -> len


* ``task:: MSKtask`` An optimization task.
* ``len:: Int32`` Assigned the length of the objective name.


Obtains the length of the name assigned to the objective function.

getobjsense
-----------
::

    function getobjsense(task:: MSKtask)
    -> sense


* ``task:: MSKtask`` An optimization task.
* ``sense:: Int32`` The returned objective sense.


Gets the objective sense.

getparamname
------------
::

    function getparamname
    ( task::    MSKtask,
      partype:: Int32,
      param::   Int32 )
    -> parname


* ``task:: MSKtask`` An optimization task.
* ``partype:: Int32`` (`Enum parametertype`_) Parameter type.
* ``param:: Int32`` Which parameter.
* ``parname:: String`` Parameter name.


Obtains the name of a parameter.

getpbi
------
::

    function getpbi
    ( task::      MSKtask,
      whichsol::  Int32,
      accmode::   Int32,
      sub::       Array{Int32},
      normalize:: Int32 )
    -> pbi


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``accmode:: Int32`` (`Enum accmode`_) Defines whether bound infeasibilities related to constraints or variable are retrieved.
* ``sub:: Array{Int32}`` An array of constraint or variable indexes.
* ``pbi:: Array{Float64}`` Bound infeasibility.
* ``normalize:: Int32`` If non-zero, normalize with largest absolute value of input data.


Deprecated.

getpcni
-------
::

    function getpcni
    ( task::     MSKtask,
      whichsol:: Int32,
      sub::      Array{Int32} )
    -> pcni


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``sub:: Array{Int32}`` Constraint indexes for which to calculate the equation infeasibility.
* ``pcni:: Array{Float64}`` Primal cone infeasibility.


Deprecated.

getpeqi
-------
::

    function getpeqi
    ( task::      MSKtask,
      whichsol::  Int32,
      sub::       Array{Int32},
      normalize:: Int32 )
    -> peqi


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``sub:: Array{Int32}`` Constraint indexes for which to calculate the equation infeasibility.
* ``peqi:: Array{Float64}`` Equation infeasibility.
* ``normalize:: Int32`` If non-zero, normalize with largest absolute value of input data.


Deprecated.

getprimalobj
------------
::

    function getprimalobj
    ( task::     MSKtask,
      whichsol:: Int32 )
    -> primalobj


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``primalobj:: Float64`` Objective value corresponding to the primal solution.


Computes the primal objective value for the desired solution.

getprobtype
-----------
::

    function getprobtype(task:: MSKtask)
    -> probtype


* ``task:: MSKtask`` An optimization task.
* ``probtype:: Int32`` The problem type.


Obtains the problem type.

getprosta
---------
::

    function getprosta
    ( task::     MSKtask,
      whichsol:: Int32 )
    -> prosta


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``prosta:: Int32`` Problem status.


Obtains the problem status.

getpviolbarvar
--------------
::

    function getpviolbarvar
    ( task::     MSKtask,
      whichsol:: Int32,
      sub::      Array{Int32} )
    -> viol


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``sub:: Array{Int32}`` An array of indexes of barx variables.
* ``viol:: Array{Float64}`` List of violations corresponding to sub.


Computes the violation of a primal solution for a list of barx variables.

getpviolcon
-----------
::

    function getpviolcon
    ( task::     MSKtask,
      whichsol:: Int32,
      sub::      Array{Int32} )
    -> viol


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``sub:: Array{Int32}`` An array of indexes of constraints.
* ``viol:: Array{Float64}`` List of violations corresponding to sub.


Computes the violation of a primal solution for a list of xc variables.

getpviolcones
-------------
::

    function getpviolcones
    ( task::     MSKtask,
      whichsol:: Int32,
      sub::      Array{Int32} )
    -> viol


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``sub:: Array{Int32}`` An array of indexes of barx variables.
* ``viol:: Array{Float64}`` List of violations corresponding to sub.


Computes the violation of a solution for set of conic constraints.

getpviolvar
-----------
::

    function getpviolvar
    ( task::     MSKtask,
      whichsol:: Int32,
      sub::      Array{Int32} )
    -> viol


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``sub:: Array{Int32}`` An array of indexes of x variables.
* ``viol:: Array{Float64}`` List of violations corresponding to sub.


Computes the violation of a primal solution for a list of x variables.

getqconk
--------
::

    function getqconk
    ( task:: MSKtask,
      k::    Int32 )
    -> numqcnz,qcsubi,qcsubj,qcval


* ``task:: MSKtask`` An optimization task.
* ``k:: Int32`` Which constraint.
* ``numqcnz:: Int64`` Number of quadratic terms.
* ``qcsubi:: Array{Int32}`` Row subscripts for quadratic constraint matrix.
* ``qcsubj:: Array{Int32}`` Column subscripts for quadratic constraint matrix.
* ``qcval:: Array{Float64}`` Quadratic constraint coefficient values.


Obtains all the quadratic terms in a constraint.

getqobj
-------
::

    function getqobj(task:: MSKtask)
    -> numqonz,qosubi,qosubj,qoval


* ``task:: MSKtask`` An optimization task.
* ``numqonz:: Int32`` Number of non-zero elements in the quadratic objective terms.
* ``qosubi:: Array{Int32}`` Row subscripts for quadratic objective coefficients.
* ``qosubj:: Array{Int32}`` Column subscripts for quadratic objective coefficients.
* ``qoval:: Array{Float64}`` Quadratic objective coefficient values.


Obtains all the quadratic terms in the objective.

getqobj64
---------
::

    function getqobj64(task:: MSKtask)
    -> numqonz,qosubi,qosubj,qoval


* ``task:: MSKtask`` An optimization task.
* ``numqonz:: Int64`` Number of non-zero elements in the quadratic objective terms.
* ``qosubi:: Array{Int32}`` Row subscripts for quadratic objective coefficients.
* ``qosubj:: Array{Int32}`` Column subscripts for quadratic objective coefficients.
* ``qoval:: Array{Float64}`` Quadratic objective coefficient values.


Obtains all the quadratic terms in the objective.

getqobjij
---------
::

    function getqobjij
    ( task:: MSKtask,
      i::    Int32,
      j::    Int32 )
    -> qoij


* ``task:: MSKtask`` An optimization task.
* ``i:: Int32`` Row index of the coefficient.
* ``j:: Int32`` Column index of coefficient.
* ``qoij:: Float64`` The required coefficient.


Obtains one coefficient from the quadratic term of the objective

getreducedcosts
---------------
::

    function getreducedcosts
    ( task::     MSKtask,
      whichsol:: Int32,
      first::    Int32,
      last::     Int32 )
    -> redcosts


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``first:: Int32`` See the documentation for a full description.
* ``last:: Int32`` See the documentation for a full description.
* ``redcosts:: Array{Float64}`` Returns the requested reduced costs. See documentation for a full description.


Obtains the difference of (slx-sux) for a sequence of variables.

getskc
------
::

    function getskc
    ( task::     MSKtask,
      whichsol:: Int32 )
    -> skc


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``skc:: Array{Int32}`` (`Enum stakey`_) Status keys for the constraints.


Obtains the status keys for the constraints.

getskcslice
-----------
::

    function getskcslice
    ( task::     MSKtask,
      whichsol:: Int32,
      first::    Int32,
      last::     Int32 )
    -> skc


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``first:: Int32`` First index in the sequence.
* ``last:: Int32`` Last index plus 1 in the sequence.
* ``skc:: Array{Int32}`` (`Enum stakey`_) Status keys for the constraints.


Obtains the status keys for the constraints.

getskx
------
::

    function getskx
    ( task::     MSKtask,
      whichsol:: Int32 )
    -> skx


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``skx:: Array{Int32}`` (`Enum stakey`_) Status keys for the variables.


Obtains the status keys for the scalar variables.

getskxslice
-----------
::

    function getskxslice
    ( task::     MSKtask,
      whichsol:: Int32,
      first::    Int32,
      last::     Int32 )
    -> skx


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``first:: Int32`` First index in the sequence.
* ``last:: Int32`` Last index plus 1 in the sequence.
* ``skx:: Array{Int32}`` (`Enum stakey`_) Status keys for the variables.


Obtains the status keys for the variables.

getslc
------
::

    function getslc
    ( task::     MSKtask,
      whichsol:: Int32 )
    -> slc


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``slc:: Array{Float64}`` The slc vector.


Obtains the slc vector for a solution.

getslcslice
-----------
::

    function getslcslice
    ( task::     MSKtask,
      whichsol:: Int32,
      first::    Int32,
      last::     Int32 )
    -> slc


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``first:: Int32`` First index in the sequence.
* ``last:: Int32`` Last index plus 1 in the sequence.
* ``slc:: Array{Float64}`` Dual variables corresponding to the lower bounds on the constraints.


Obtains a slice of the slc vector for a solution.

getslx
------
::

    function getslx
    ( task::     MSKtask,
      whichsol:: Int32 )
    -> slx


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``slx:: Array{Float64}`` The slx vector.


Obtains the slx vector for a solution.

getslxslice
-----------
::

    function getslxslice
    ( task::     MSKtask,
      whichsol:: Int32,
      first::    Int32,
      last::     Int32 )
    -> slx


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``first:: Int32`` First index in the sequence.
* ``last:: Int32`` Last index plus 1 in the sequence.
* ``slx:: Array{Float64}`` Dual variables corresponding to the lower bounds on the variables.


Obtains a slice of the slx vector for a solution.

getsnx
------
::

    function getsnx
    ( task::     MSKtask,
      whichsol:: Int32 )
    -> snx


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``snx:: Array{Float64}`` The snx vector.


Obtains the snx vector for a solution.

getsnxslice
-----------
::

    function getsnxslice
    ( task::     MSKtask,
      whichsol:: Int32,
      first::    Int32,
      last::     Int32 )
    -> snx


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``first:: Int32`` First index in the sequence.
* ``last:: Int32`` Last index plus 1 in the sequence.
* ``snx:: Array{Float64}`` Dual variables corresponding to the conic constraints on the variables.


Obtains a slice of the snx vector for a solution.

getsolsta
---------
::

    function getsolsta
    ( task::     MSKtask,
      whichsol:: Int32 )
    -> solsta


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``solsta:: Int32`` Solution status.


Obtains the solution status.

getsolution
-----------
::

    function getsolution
    ( task::     MSKtask,
      whichsol:: Int32 )
    -> prosta,solsta,skc,skx,skn,xc,xx,y,slc,suc,slx,sux,snx


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``prosta:: Int32`` Problem status.
* ``solsta:: Int32`` Solution status.
* ``skc:: Array{Int32}`` (`Enum stakey`_) Status keys for the constraints.
* ``skx:: Array{Int32}`` (`Enum stakey`_) Status keys for the variables.
* ``skn:: Array{Int32}`` (`Enum stakey`_) Status keys for the conic constraints.
* ``xc:: Array{Float64}`` Primal constraint solution.
* ``xx:: Array{Float64}`` Primal variable solution.
* ``y:: Array{Float64}`` Vector of dual variables corresponding to the constraints.
* ``slc:: Array{Float64}`` Dual variables corresponding to the lower bounds on the constraints.
* ``suc:: Array{Float64}`` Dual variables corresponding to the upper bounds on the constraints.
* ``slx:: Array{Float64}`` Dual variables corresponding to the lower bounds on the variables.
* ``sux:: Array{Float64}`` Dual variables corresponding to the upper bounds on the variables.
* ``snx:: Array{Float64}`` Dual variables corresponding to the conic constraints on the variables.


Obtains the complete solution.

getsolutioni
------------
::

    function getsolutioni
    ( task::     MSKtask,
      accmode::  Int32,
      i::        Int32,
      whichsol:: Int32 )
    -> sk,x,sl,su,sn


* ``task:: MSKtask`` An optimization task.
* ``accmode:: Int32`` (`Enum accmode`_) Defines whether solution information for a constraint or for a variable is retrieved.
* ``i:: Int32`` Index of the constraint or variable.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``sk:: Int32`` Status key of the constraint of variable.
* ``x:: Float64`` Solution value of the primal variable.
* ``sl:: Float64`` Solution value of the dual variable associated with the lower bound.
* ``su:: Float64`` Solution value of the dual variable associated with the upper bound.
* ``sn:: Float64`` Solution value of the dual variable associated with the cone constraint.


Obtains the solution for a single constraint or variable.

getsolutioninf
--------------
::

    function getsolutioninf
    ( task::     MSKtask,
      whichsol:: Int32 )
    -> prosta,solsta,primalobj,maxpbi,maxpcni,maxpeqi,maxinti,dualobj,maxdbi,maxdcni,maxdeqi


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``prosta:: Int32`` Problem status.
* ``solsta:: Int32`` Solution status.
* ``primalobj:: Float64`` Value of the primal objective.
* ``maxpbi:: Float64`` Maximum infeasibility in primal bounds on variables.
* ``maxpcni:: Float64`` Maximum infeasibility in the primal conic constraints.
* ``maxpeqi:: Float64`` Maximum infeasibility in primal equality constraints.
* ``maxinti:: Float64`` Maximum infeasibility in primal equality constraints.
* ``dualobj:: Float64`` Value of the dual objective.
* ``maxdbi:: Float64`` Maximum infeasibility in bounds on dual variables.
* ``maxdcni:: Float64`` Maximum infeasibility in the dual conic constraints.
* ``maxdeqi:: Float64`` Maximum infeasibility in the dual equality constraints.


Deprecated

getsolutioninfo
---------------
::

    function getsolutioninfo
    ( task::     MSKtask,
      whichsol:: Int32 )
    -> pobj,pviolcon,pviolvar,pviolbarvar,pviolcone,pviolitg,dobj,dviolcon,dviolvar,dviolbarvar,dviolcones


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``pobj:: Float64`` The primal objective value.
* ``pviolcon:: Float64`` Maximal primal bound violation for a xc variable.
* ``pviolvar:: Float64`` Maximal primal bound violation for a xx variable.
* ``pviolbarvar:: Float64`` Maximal primal bound violation for a barx variable.
* ``pviolcone:: Float64`` Maximal primal violation of the solution with respect to the conic constraints.
* ``pviolitg:: Float64`` Maximal violation in the integer constraints.
* ``dobj:: Float64`` Dual objective value.
* ``dviolcon:: Float64`` Maximal dual bound violation a xc variable.
* ``dviolvar:: Float64`` Maximal dual bound violation xx variable.
* ``dviolbarvar:: Float64`` Maximal dual bound violation for a bars variable.
* ``dviolcones:: Float64`` Maximum violation of the dual solution in the dual conic constraints .


Obtains information about of a solution.

getsolutionslice
----------------
::

    function getsolutionslice
    ( task::     MSKtask,
      whichsol:: Int32,
      solitem::  Int32,
      first::    Int32,
      last::     Int32 )
    -> values


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``solitem:: Int32`` (`Enum solitem`_) Which part of the solution is required.
* ``first:: Int32`` Index of the first value in the slice.
* ``last:: Int32`` Value of the last index+1 in the slice.
* ``values:: Array{Float64}`` The values of the requested solution elements.


Obtains a slice of the solution.

getsparsesymmat
---------------
::

    function getsparsesymmat
    ( task:: MSKtask,
      idx::  Int64 )
    -> subi,subj,valij


* ``task:: MSKtask`` An optimization task.
* ``idx:: Int64`` Index of the matrix to get.
* ``subi:: Array{Int32}`` Row subscripts of the matrix non-zero elements.
* ``subj:: Array{Int32}`` Column subscripts of the matrix non-zero elements.
* ``valij:: Array{Float64}`` Coefficients of the matrix non-zero elements.


Gets a single symmetric matrix from the matrix store.

getstrparam
-----------
::

    function getstrparam
    ( task::  MSKtask,
      param:: Int32 )
    -> len,parvalue


* ``task:: MSKtask`` An optimization task.
* ``param:: Int32`` (`Enum sparam`_) Which parameter.
* ``len:: Int32`` The length of the parameter value.
* ``parvalue:: String`` If this is not NULL, the parameter value is stored here.


Obtains the value of a string parameter.

getstrparamlen
--------------
::

    function getstrparamlen
    ( task::  MSKtask,
      param:: Int32 )
    -> len


* ``task:: MSKtask`` An optimization task.
* ``param:: Int32`` (`Enum sparam`_) Which parameter.
* ``len:: Int32`` The length of the parameter value.


Obtains the length of a string parameter.

getsuc
------
::

    function getsuc
    ( task::     MSKtask,
      whichsol:: Int32 )
    -> suc


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``suc:: Array{Float64}`` The suc vector.


Obtains the suc vector for a solution.

getsucslice
-----------
::

    function getsucslice
    ( task::     MSKtask,
      whichsol:: Int32,
      first::    Int32,
      last::     Int32 )
    -> suc


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``first:: Int32`` First index in the sequence.
* ``last:: Int32`` Last index plus 1 in the sequence.
* ``suc:: Array{Float64}`` Dual variables corresponding to the upper bounds on the constraints.


Obtains a slice of the suc vector for a solution.

getsux
------
::

    function getsux
    ( task::     MSKtask,
      whichsol:: Int32 )
    -> sux


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``sux:: Array{Float64}`` The sux vector.


Obtains the sux vector for a solution.

getsuxslice
-----------
::

    function getsuxslice
    ( task::     MSKtask,
      whichsol:: Int32,
      first::    Int32,
      last::     Int32 )
    -> sux


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``first:: Int32`` First index in the sequence.
* ``last:: Int32`` Last index plus 1 in the sequence.
* ``sux:: Array{Float64}`` Dual variables corresponding to the upper bounds on the variables.


Obtains a slice of the sux vector for a solution.

getsymmatinfo
-------------
::

    function getsymmatinfo
    ( task:: MSKtask,
      idx::  Int64 )
    -> dim,nz,type


* ``task:: MSKtask`` An optimization task.
* ``idx:: Int64`` Index of the matrix that is requested information about.
* ``dim:: Int32`` Returns the dimension of the requested matrix.
* ``nz:: Int64`` Returns the number of non-zeros in the requested matrix.
* ``type:: Int32`` Returns the type of the requested matrix.


Obtains information of  a matrix from the symmetric matrix storage E.

gettaskname
-----------
::

    function gettaskname(task:: MSKtask)
    -> taskname


* ``task:: MSKtask`` An optimization task.
* ``taskname:: String`` Is assigned the task name.


Obtains the task name.

gettasknamelen
--------------
::

    function gettasknamelen(task:: MSKtask)
    -> len


* ``task:: MSKtask`` An optimization task.
* ``len:: Int32`` Returns the length of the task name.


Obtains the length the task name.

getvarbound
-----------
::

    function getvarbound
    ( task:: MSKtask,
      i::    Int32 )
    -> bk,bl,bu


* ``task:: MSKtask`` An optimization task.
* ``i:: Int32`` Index of the variable for which the bound information should be obtained.
* ``bk:: Int32`` Bound keys.
* ``bl:: Float64`` Values for lower bounds.
* ``bu:: Float64`` Values for upper bounds.


Obtains bound information for one variable.

getvarboundslice
----------------
::

    function getvarboundslice
    ( task::  MSKtask,
      first:: Int32,
      last::  Int32 )
    -> bk,bl,bu


* ``task:: MSKtask`` An optimization task.
* ``first:: Int32`` First index in the sequence.
* ``last:: Int32`` Last index plus 1 in the sequence.
* ``bk:: Array{Int32}`` (`Enum boundkey`_) Bound keys.
* ``bl:: Array{Float64}`` Values for lower bounds.
* ``bu:: Array{Float64}`` Values for upper bounds.


Obtains bounds information for a slice of the variables.

getvarbranchdir
---------------
::

    function getvarbranchdir
    ( task:: MSKtask,
      j::    Int32 )
    -> direction


* ``task:: MSKtask`` An optimization task.
* ``j:: Int32`` Index of the variable.
* ``direction:: Int32`` The branching direction assigned to the j'th variable.


Obtains the branching direction for a variable.

getvarbranchpri
---------------
::

    function getvarbranchpri
    ( task:: MSKtask,
      j::    Int32 )
    -> priority


* ``task:: MSKtask`` An optimization task.
* ``j:: Int32`` Index of the variable.
* ``priority:: Int32`` The branching priority assigned to the j'th variable.


Obtains the branching priority for a variable.

getvarname
----------
::

    function getvarname
    ( task:: MSKtask,
      j::    Int32 )
    -> name


* ``task:: MSKtask`` An optimization task.
* ``j:: Int32`` Index.
* ``name:: String`` Returns the required name.


Obtains a name of a variable.

getvarnameindex
---------------
::

    function getvarnameindex
    ( task::     MSKtask,
      somename:: String )
    -> asgn,index


* ``task:: MSKtask`` An optimization task.
* ``somename:: String`` The name which should be checked.
* ``asgn:: Int32`` Is non-zero if name somename is assigned to a variable.
* ``index:: Int32`` If the name somename is assigned to a variable, then index is the name of the variable.


Checks whether the name somename has been assigned  to any variable.

getvarnamelen
-------------
::

    function getvarnamelen
    ( task:: MSKtask,
      i::    Int32 )
    -> len


* ``task:: MSKtask`` An optimization task.
* ``i:: Int32`` Index.
* ``len:: Int32`` Returns the length of the indicated name.


Obtains the length of a name of a variable variable.

getvartype
----------
::

    function getvartype
    ( task:: MSKtask,
      j::    Int32 )
    -> vartype


* ``task:: MSKtask`` An optimization task.
* ``j:: Int32`` Index of the variable.
* ``vartype:: Int32`` Variable type of variable index j.


Gets the variable type of one variable.

getvartypelist
--------------
::

    function getvartypelist
    ( task:: MSKtask,
      subj:: Array{Int32} )
    -> vartype


* ``task:: MSKtask`` An optimization task.
* ``subj:: Array{Int32}`` A list of variable indexes.
* ``vartype:: Array{Int32}`` (`Enum variabletype`_) Returns the variables types corresponding the variable indexes requested.


Obtains the variable type for one or more variables.

getxc
-----
::

    function getxc
    ( task::     MSKtask,
      whichsol:: Int32 )
    -> xc


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``xc:: Array{Float64}`` The xc vector.


Obtains the xc vector for a solution.

getxcslice
----------
::

    function getxcslice
    ( task::     MSKtask,
      whichsol:: Int32,
      first::    Int32,
      last::     Int32 )
    -> xc


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``first:: Int32`` First index in the sequence.
* ``last:: Int32`` Last index plus 1 in the sequence.
* ``xc:: Array{Float64}`` Primal constraint solution.


Obtains a slice of the xc vector for a solution.

getxx
-----
::

    function getxx
    ( task::     MSKtask,
      whichsol:: Int32 )
    -> xx


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``xx:: Array{Float64}`` The xx vector.


Obtains the xx vector for a solution.

getxxslice
----------
::

    function getxxslice
    ( task::     MSKtask,
      whichsol:: Int32,
      first::    Int32,
      last::     Int32 )
    -> xx


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``first:: Int32`` First index in the sequence.
* ``last:: Int32`` Last index plus 1 in the sequence.
* ``xx:: Array{Float64}`` Primal variable solution.


Obtains a slice of the xx vector for a solution.

gety
----
::

    function gety
    ( task::     MSKtask,
      whichsol:: Int32 )
    -> y


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``y:: Array{Float64}`` The y vector.


Obtains the y vector for a solution.

getyslice
---------
::

    function getyslice
    ( task::     MSKtask,
      whichsol:: Int32,
      first::    Int32,
      last::     Int32 )
    -> y


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``first:: Int32`` First index in the sequence.
* ``last:: Int32`` Last index plus 1 in the sequence.
* ``y:: Array{Float64}`` Vector of dual variables corresponding to the constraints.


Obtains a slice of the y vector for a solution.

initbasissolve
--------------
::

    function initbasissolve(task:: MSKtask)
    -> basis


* ``task:: MSKtask`` An optimization task.
* ``basis:: Array{Int32}`` The array of basis indexes to use.


Prepare a task for basis solver.

inputdata
---------
::

    function inputdata
    ( task::      MSKtask,
      maxnumcon:: Int32,
      maxnumvar:: Int32,
      c::         Array{Float64},
      cfix::      Float64,
      aptrb::     Array{Int64},
      aptre::     Array{Int64},
      asub::      Array{Int32},
      aval::      Array{Float64},
      bkc::       Array{Int32},
      blc::       Array{Float64},
      buc::       Array{Float64},
      bkx::       Array{Int32},
      blx::       Array{Float64},
      bux::       Array{Float64} )

    function inputdata
    ( task::      MSKtask,
      maxnumcon,
      maxnumvar,
      c::         Array{Float64},
      cfix::      Float64,
      A::         SparseMatrixCSC{Float64},
      bkc::       Array{Int32},
      blc::       Array{Float64},
      buc::       Array{Float64},
      bkx::       Array{Int32},
      blx::       Array{Float64},
      bux::       Array{Float64} )


* ``task:: MSKtask`` An optimization task.
* ``maxnumcon:: Int32`` Number of preallocated constraints in the optimization task.
* ``maxnumvar:: Int32`` Number of preallocated variables in the optimization task.
* ``c:: Array{Float64}`` Linear terms of the objective as a dense vector. The lengths is the number of variables.
* ``cfix:: Float64`` Fixed term in the objective.
* ``aptrb:: Array{Int64}`` Row or column end pointers.
* ``aptre:: Array{Int64}`` Row or column start pointers.
* ``asub:: Array{Int32}`` Coefficient subscripts.
* ``aval:: Array{Float64}`` Coefficient values.
* ``bkc:: Array{Int32}`` (`Enum boundkey`_) Bound keys for the constraints.
* ``blc:: Array{Float64}`` Lower bounds for the constraints.
* ``buc:: Array{Float64}`` Upper bounds for the constraints.
* ``bkx:: Array{Int32}`` (`Enum boundkey`_) Bound keys for the variables.
* ``blx:: Array{Float64}`` Lower bounds for the variables.
* ``bux:: Array{Float64}`` Upper bounds for the variables.
* ``A:: SparseMatrixCSC{Float64}`` Sparse matrix defining the column values


Input the linear part of an optimization task in one function call.

isdouparname
------------
::

    function isdouparname
    ( task::    MSKtask,
      parname:: String )
    -> param


* ``task:: MSKtask`` An optimization task.
* ``parname:: String`` Parameter name.
* ``param:: Int32`` Which parameter.


Checks a double parameter name.

isintparname
------------
::

    function isintparname
    ( task::    MSKtask,
      parname:: String )
    -> param


* ``task:: MSKtask`` An optimization task.
* ``parname:: String`` Parameter name.
* ``param:: Int32`` Which parameter.


Checks an integer parameter name.

isstrparname
------------
::

    function isstrparname
    ( task::    MSKtask,
      parname:: String )
    -> param


* ``task:: MSKtask`` An optimization task.
* ``parname:: String`` Parameter name.
* ``param:: Int32`` Which parameter.


Checks a string parameter name.

linkfiletostream
----------------
::

    function linkfiletostream
    ( task::        MSKtask,
      whichstream:: Int32,
      filename::    String,
      append::      Int32 )


* ``task:: MSKtask`` An optimization task.
* ``whichstream:: Int32`` (`Enum streamtype`_) Index of the stream.
* ``filename:: String`` The name of the file where the stream is written.
* ``append:: Int32`` If this argument is 0 the output file will be overwritten, otherwise text is append to the output file.


Directs all output from a task stream to a file.

onesolutionsummary
------------------
::

    function onesolutionsummary
    ( task::        MSKtask,
      whichstream:: Int32,
      whichsol::    Int32 )


* ``task:: MSKtask`` An optimization task.
* ``whichstream:: Int32`` (`Enum streamtype`_) Index of the stream.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.


Prints a short summary for the specified solution.

optimizeconcurrent
------------------
::

    function optimizeconcurrent
    ( task::      MSKtask,
      taskarray:: Array{MSKtask} )


* ``task:: MSKtask`` An optimization task.
* ``taskarray:: Array{MSKtask}`` An array of tasks.


Optimize a given task with several optimizers concurrently.

optimizersummary
----------------
::

    function optimizersummary
    ( task::        MSKtask,
      whichstream:: Int32 )


* ``task:: MSKtask`` An optimization task.
* ``whichstream:: Int32`` (`Enum streamtype`_) Index of the stream.


Prints a short summary with optimizer statistics for last optimization.

optimize
--------
::

    function optimize(task:: MSKtask)
    -> trmcode


* ``task:: MSKtask`` An optimization task.
* ``trmcode:: Int32`` Is either OK or a termination response code.


Optimizes the problem.

primalrepair
------------
::

    function primalrepair
    ( task:: MSKtask,
      wlc::  Array{Float64},
      wuc::  Array{Float64},
      wlx::  Array{Float64},
      wux::  Array{Float64} )


* ``task:: MSKtask`` An optimization task.
* ``wlc:: Array{Float64}`` Weights associated with relaxing lower bounds on the constraints.
* ``wuc:: Array{Float64}`` Weights associated with relaxing the upper bound on the constraints.
* ``wlx:: Array{Float64}`` Weights associated with relaxing the lower bounds of the variables.
* ``wux:: Array{Float64}`` Weights associated with relaxing the upper bounds of variables.


The function repairs a primal infeasible optimization problem by adjusting the bounds on the constraints and variables.

primalsensitivity
-----------------
::

    function primalsensitivity
    ( task::  MSKtask,
      subi::  Array{Int32},
      marki:: Array{Int32},
      subj::  Array{Int32},
      markj:: Array{Int32} )
    -> leftpricei,rightpricei,leftrangei,rightrangei,leftpricej,rightpricej,leftrangej,rightrangej


* ``task:: MSKtask`` An optimization task.
* ``subi:: Array{Int32}`` Indexes of bounds on constraints to analyze.
* ``marki:: Array{Int32}`` (`Enum mark`_) Mark which constraint bounds to analyze.
* ``subj:: Array{Int32}`` Indexes of bounds on variables to analyze.
* ``markj:: Array{Int32}`` (`Enum mark`_) Mark which variable bounds to analyze.
* ``leftpricei:: Array{Float64}`` Left shadow price for constraints.
* ``rightpricei:: Array{Float64}`` Right shadow price for constraints.
* ``leftrangei:: Array{Float64}`` Left range for constraints.
* ``rightrangei:: Array{Float64}`` Right range for constraints.
* ``leftpricej:: Array{Float64}`` Left price for variables.
* ``rightpricej:: Array{Float64}`` Right price for variables.
* ``leftrangej:: Array{Float64}`` Left range for variables.
* ``rightrangej:: Array{Float64}`` Right range for variables.


Perform sensitivity analysis on bounds.

printdata
---------
::

    function printdata
    ( task::        MSKtask,
      whichstream:: Int32,
      firsti::      Int32,
      lasti::       Int32,
      firstj::      Int32,
      lastj::       Int32,
      firstk::      Int32,
      lastk::       Int32,
      c::           Int32,
      qo::          Int32,
      a::           Int32,
      qc::          Int32,
      bc::          Int32,
      bx::          Int32,
      vartype::     Int32,
      cones::       Int32 )


* ``task:: MSKtask`` An optimization task.
* ``whichstream:: Int32`` (`Enum streamtype`_) Index of the stream.
* ``firsti:: Int32`` Index of first constraint for which data should be printed.
* ``lasti:: Int32`` Index of last constraint plus 1 for which data should be printed.
* ``firstj:: Int32`` Index of first variable for which data should be printed.
* ``lastj:: Int32`` Index of last variable plus 1 for which data should be printed.
* ``firstk:: Int32`` Index of first cone for which data should be printed.
* ``lastk:: Int32`` Index of last cone plus 1 for which data should be printed.
* ``c:: Int32`` If non-zero the linear objective terms are printed.
* ``qo:: Int32`` If non-zero the quadratic objective terms are printed.
* ``a:: Int32`` If non-zero the linear constraint matrix is printed.
* ``qc:: Int32`` If non-zero q'th     quadratic constraint terms are printed for the relevant constraints.
* ``bc:: Int32`` If non-zero the constraints bounds are printed.
* ``bx:: Int32`` If non-zero the variable bounds are printed.
* ``vartype:: Int32`` If non-zero the variable types are printed.
* ``cones:: Int32`` If non-zero the  conic data is printed.


Prints a part of the problem data to a stream.

printparam
----------
::

    function printparam(task:: MSKtask)


* ``task:: MSKtask`` An optimization task.


Prints the current parameter settings.

putacol
-------
::

    function putacol
    ( task:: MSKtask,
      j::    Int32,
      subj:: Array{Int32},
      valj:: Array{Float64} )


* ``task:: MSKtask`` An optimization task.
* ``j:: Int32`` Column index.
* ``subj:: Array{Int32}`` Row indexes of non-zero values in column.
* ``valj:: Array{Float64}`` New non-zero values of column.


Replaces all elements in one column of A.

putacollist
-----------
::

    function putacollist
    ( task:: MSKtask,
      sub::  Array{Int32},
      ptrb:: Array{Int64},
      ptre:: Array{Int64},
      asub:: Array{Int32},
      aval:: Array{Float64} )

    function putacollist
    ( task:: MSKtask,
      sub::  Array,
      A::    SparseMatrixCSC{Float64} )


* ``task:: MSKtask`` An optimization task.
* ``sub:: Array{Int32}`` Indexes of columns that should be replaced.
* ``ptrb:: Array{Int64}`` Array of pointers to the first element in the columns.
* ``ptre:: Array{Int64}`` Array of pointers to the last element plus one in the columns.
* ``asub:: Array{Int32}`` Variable indexes.
* ``aval:: Array{Float64}`` Coefficient values.
* ``A:: SparseMatrixCSC{Float64}`` Sparse matrix defining the column values


Replaces all elements in several columns the linear constraint matrix by new values.

putacolslice
------------
::

    function putacolslice
    ( task::  MSKtask,
      first:: Int32,
      last::  Int32,
      ptrb::  Array{Int64},
      ptre::  Array{Int64},
      asub::  Array{Int32},
      aval::  Array{Float64} )

    function putacolslice
    ( task::  MSKtask,
      first,
      last,
      A::     SparseMatrixCSC{Float64} )


* ``task:: MSKtask`` An optimization task.
* ``first:: Int32`` First column in the slice.
* ``last:: Int32`` Last column plus one in the slice.
* ``ptrb:: Array{Int64}`` Array of pointers to the first element in the columns.
* ``ptre:: Array{Int64}`` Array of pointers to the last element plus one in the columns.
* ``asub:: Array{Int32}`` Variable indexes.
* ``aval:: Array{Float64}`` Coefficient values.
* ``A:: SparseMatrixCSC{Float64}`` Sparse matrix defining the column values


Replaces all elements in several columns the linear constraint matrix by new values.

putaij
------
::

    function putaij
    ( task:: MSKtask,
      i::    Int32,
      j::    Int32,
      aij::  Float64 )


* ``task:: MSKtask`` An optimization task.
* ``i:: Int32`` Index of the constraint in which the change should occur.
* ``j:: Int32`` Index of the variable in which the change should occur.
* ``aij:: Float64`` New coefficient.


Changes a single value in the linear coefficient matrix.

putaijlist
----------
::

    function putaijlist
    ( task::  MSKtask,
      subi::  Array{Int32},
      subj::  Array{Int32},
      valij:: Array{Float64} )


* ``task:: MSKtask`` An optimization task.
* ``subi:: Array{Int32}`` Constraint indexes in which the change should occur.
* ``subj:: Array{Int32}`` Variable indexes in which the change should occur.
* ``valij:: Array{Float64}`` New coefficient values.


Changes one or more coefficients in the linear constraint matrix.

putaijlist64
------------
::

    function putaijlist64
    ( task::  MSKtask,
      subi::  Array{Int32},
      subj::  Array{Int32},
      valij:: Array{Float64} )


* ``task:: MSKtask`` An optimization task.
* ``subi:: Array{Int32}`` Constraint indexes in which the change should occur.
* ``subj:: Array{Int32}`` Variable indexes in which the change should occur.
* ``valij:: Array{Float64}`` New coefficient values.


Changes one or more coefficients in the linear constraint matrix.

putarow
-------
::

    function putarow
    ( task:: MSKtask,
      i::    Int32,
      subi:: Array{Int32},
      vali:: Array{Float64} )


* ``task:: MSKtask`` An optimization task.
* ``i:: Int32`` row index.
* ``subi:: Array{Int32}`` Row indexes of non-zero values in row.
* ``vali:: Array{Float64}`` New non-zero values of row.


Replaces all elements in one row of A.

putarowlist
-----------
::

    function putarowlist
    ( task:: MSKtask,
      sub::  Array{Int32},
      ptrb:: Array{Int64},
      ptre:: Array{Int64},
      asub:: Array{Int32},
      aval:: Array{Float64} )

    function putarowlist
    ( task:: MSKtask,
      sub::  Array,
      At::   SparseMatrixCSC{Float64} )


* ``task:: MSKtask`` An optimization task.
* ``sub:: Array{Int32}`` Indexes of rows or columns that should be replaced.
* ``ptrb:: Array{Int64}`` Array of pointers to the first element in the rows or columns.
* ``ptre:: Array{Int64}`` Array of pointers to the last element plus one in the rows or columns.
* ``asub:: Array{Int32}`` Variable indexes.
* ``aval:: Array{Float64}`` Coefficient values.
* ``At:: SparseMatrixCSC{Float64}`` Transposed matrix defining the row values. Note that for efficiency reasons the *columns* of this matrix defines the *rows* to be replaced


Replaces all elements in several rows the linear constraint matrix by new values.

putarowslice
------------
::

    function putarowslice
    ( task::  MSKtask,
      first:: Int32,
      last::  Int32,
      ptrb::  Array{Int64},
      ptre::  Array{Int64},
      asub::  Array{Int32},
      aval::  Array{Float64} )

    function putarowslice
    ( task::  MSKtask,
      first,
      last,
      At::    SparseMatrixCSC{Float64} )


* ``task:: MSKtask`` An optimization task.
* ``first:: Int32`` First row in the slice.
* ``last:: Int32`` Last row plus one in the slice.
* ``ptrb:: Array{Int64}`` Array of pointers to the first element in the rows.
* ``ptre:: Array{Int64}`` Array of pointers to the last element plus one in the rows.
* ``asub:: Array{Int32}`` Variable indexes.
* ``aval:: Array{Float64}`` Coefficient values.
* ``At:: SparseMatrixCSC{Float64}`` Transposed matrix defining the row values. Note that for efficiency reasons the *columns* of this matrix defines the *rows* to be replaced


Replaces all elements in several rows the linear constraint matrix by new values.

putbarablocktriplet
-------------------
::

    function putbarablocktriplet
    ( task::    MSKtask,
      num::     Int64,
      subi::    Array{Int32},
      subj::    Array{Int32},
      subk::    Array{Int32},
      subl::    Array{Int32},
      valijkl:: Array{Float64} )


* ``task:: MSKtask`` An optimization task.
* ``num:: Int64`` Number of elements in the block triplet form.
* ``subi:: Array{Int32}`` Constraint index.
* ``subj:: Array{Int32}`` Symmetric matrix variable index.
* ``subk:: Array{Int32}`` Block row index.
* ``subl:: Array{Int32}`` Block column index.
* ``valijkl:: Array{Float64}`` The numerical value associated with the block triplet.


Inputs barA in block triplet form.

putbaraij
---------
::

    function putbaraij
    ( task::    MSKtask,
      i::       Int32,
      j::       Int32,
      sub::     Array{Int64},
      weights:: Array{Float64} )


* ``task:: MSKtask`` An optimization task.
* ``i:: Int32`` Row index of barA.
* ``j:: Int32`` Column index of barA.
* ``sub:: Array{Int64}`` See argument weights for an explanation.
* ``weights:: Array{Float64}`` Weights in the weighted sum.


Inputs an element of barA.

putbarcblocktriplet
-------------------
::

    function putbarcblocktriplet
    ( task::   MSKtask,
      num::    Int64,
      subj::   Array{Int32},
      subk::   Array{Int32},
      subl::   Array{Int32},
      valjkl:: Array{Float64} )


* ``task:: MSKtask`` An optimization task.
* ``num:: Int64`` Number of elements in the block triplet form.
* ``subj:: Array{Int32}`` Symmetric matrix variable index.
* ``subk:: Array{Int32}`` Block row index.
* ``subl:: Array{Int32}`` Block column index.
* ``valjkl:: Array{Float64}`` The numerical value associated with the block triplet.


Inputs barC in block triplet form.

putbarcj
--------
::

    function putbarcj
    ( task::    MSKtask,
      j::       Int32,
      sub::     Array{Int64},
      weights:: Array{Float64} )


* ``task:: MSKtask`` An optimization task.
* ``j:: Int32`` Index of the element in barc$ that should be changed.
* ``sub:: Array{Int64}`` sub is list of indexes of those symmetric matrices appearing in sum.
* ``weights:: Array{Float64}`` The weights of the terms in the weighted sum.


Changes one element in barc.

putbarsj
--------
::

    function putbarsj
    ( task::     MSKtask,
      whichsol:: Int32,
      j::        Int32,
      barsj::    Array{Float64} )


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``j:: Int32`` Index of the semidefinite variable.
* ``barsj:: Array{Float64}`` Value of the j'th variable of barx.


Sets the dual solution for a semidefinite variable.

putbarvarname
-------------
::

    function putbarvarname
    ( task:: MSKtask,
      j::    Int32,
      name:: String )


* ``task:: MSKtask`` An optimization task.
* ``j:: Int32`` Index of the variable.
* ``name:: String`` The variable name.


Puts the name of a semidefinite variable.

putbarxj
--------
::

    function putbarxj
    ( task::     MSKtask,
      whichsol:: Int32,
      j::        Int32,
      barxj::    Array{Float64} )


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``j:: Int32`` Index of the semidefinite variable.
* ``barxj:: Array{Float64}`` Value of the j'th variable of barx.


Sets the primal solution for a semidefinite variable.

putbound
--------
::

    function putbound
    ( task::    MSKtask,
      accmode:: Int32,
      i::       Int32,
      bk::      Int32,
      bl::      Float64,
      bu::      Float64 )


* ``task:: MSKtask`` An optimization task.
* ``accmode:: Int32`` (`Enum accmode`_) Defines whether the bound for a constraint or a variable is changed.
* ``i:: Int32`` Index of the constraint or variable.
* ``bk:: Int32`` (`Enum boundkey`_) New bound key.
* ``bl:: Float64`` New lower bound.
* ``bu:: Float64`` New upper bound.


Changes the bound for either one constraint or one variable.

putboundlist
------------
::

    function putboundlist
    ( task::    MSKtask,
      accmode:: Int32,
      sub::     Array{Int32},
      bk::      Array{Int32},
      bl::      Array{Float64},
      bu::      Array{Float64} )


* ``task:: MSKtask`` An optimization task.
* ``accmode:: Int32`` (`Enum accmode`_) Defines whether to access bounds on variables or constraints.
* ``sub:: Array{Int32}`` Subscripts of the bounds that should be changed.
* ``bk:: Array{Int32}`` (`Enum boundkey`_) Bound keys for variables or constraints.
* ``bl:: Array{Float64}`` Bound keys for variables or constraints.
* ``bu:: Array{Float64}`` Constraint or variable upper bounds.


Changes the bounds of constraints or variables.

putboundslice
-------------
::

    function putboundslice
    ( task::  MSKtask,
      con::   Int32,
      first:: Int32,
      last::  Int32,
      bk::    Array{Int32},
      bl::    Array{Float64},
      bu::    Array{Float64} )


* ``task:: MSKtask`` An optimization task.
* ``con:: Int32`` (`Enum accmode`_) Determines whether variables or constraints are modified.
* ``first:: Int32`` First index in the sequence.
* ``last:: Int32`` Last index plus 1 in the sequence.
* ``bk:: Array{Int32}`` (`Enum boundkey`_) Bound keys.
* ``bl:: Array{Float64}`` Values for lower bounds.
* ``bu:: Array{Float64}`` Values for upper bounds.


Modifies bounds.

putcfix
-------
::

    function putcfix
    ( task:: MSKtask,
      cfix:: Float64 )


* ``task:: MSKtask`` An optimization task.
* ``cfix:: Float64`` Fixed term in the objective.


Replaces the fixed term in the objective.

putcj
-----
::

    function putcj
    ( task:: MSKtask,
      j::    Int32,
      cj::   Float64 )


* ``task:: MSKtask`` An optimization task.
* ``j:: Int32`` Index of the variable whose objective coefficient should be changed.
* ``cj:: Float64`` New coefficient value.


Modifies one linear coefficient in the objective.

putclist
--------
::

    function putclist
    ( task:: MSKtask,
      subj:: Array{Int32},
      val::  Array{Float64} )


* ``task:: MSKtask`` An optimization task.
* ``subj:: Array{Int32}`` Index of variables for which objective coefficients should be changed.
* ``val:: Array{Float64}`` New numerical values for the objective coefficients that should be modified.


Modifies a part of the linear objective coefficients.

putconbound
-----------
::

    function putconbound
    ( task:: MSKtask,
      i::    Int32,
      bk::   Int32,
      bl::   Float64,
      bu::   Float64 )


* ``task:: MSKtask`` An optimization task.
* ``i:: Int32`` Index of the constraint.
* ``bk:: Int32`` (`Enum boundkey`_) New bound key.
* ``bl:: Float64`` New lower bound.
* ``bu:: Float64`` New upper bound.


Changes the bound for one constraint.

putconboundlist
---------------
::

    function putconboundlist
    ( task:: MSKtask,
      sub::  Array{Int32},
      bkc::  Array{Int32},
      blc::  Array{Float64},
      buc::  Array{Float64} )


* ``task:: MSKtask`` An optimization task.
* ``sub:: Array{Int32}`` List constraints indexes.
* ``bkc:: Array{Int32}`` (`Enum boundkey`_) New bound keys.
* ``blc:: Array{Float64}`` New lower bound values.
* ``buc:: Array{Float64}`` New upper bounds values.


Changes the bounds of a list of constraints.

putconboundslice
----------------
::

    function putconboundslice
    ( task::  MSKtask,
      first:: Int32,
      last::  Int32,
      bk::    Array{Int32},
      bl::    Array{Float64},
      bu::    Array{Float64} )


* ``task:: MSKtask`` An optimization task.
* ``first:: Int32`` Index of the first constraint in the slice.
* ``last:: Int32`` Index of the last constraint in the slice plus 1.
* ``bk:: Array{Int32}`` (`Enum boundkey`_) New bound keys.
* ``bl:: Array{Float64}`` New lower bounds.
* ``bu:: Array{Float64}`` New upper bounds.


Changes the bounds for a slice of the constraints.

putcone
-------
::

    function putcone
    ( task::     MSKtask,
      k::        Int32,
      conetype:: Int32,
      conepar::  Float64,
      submem::   Array{Int32} )


* ``task:: MSKtask`` An optimization task.
* ``k:: Int32`` Index of the cone.
* ``conetype:: Int32`` (`Enum conetype`_) Specifies the type of the cone.
* ``conepar:: Float64`` This argument is currently not used. Can be set to 0.0.
* ``submem:: Array{Int32}`` Variable subscripts of the members in the cone.


Replaces a conic constraint.

putconename
-----------
::

    function putconename
    ( task:: MSKtask,
      j::    Int32,
      name:: String )


* ``task:: MSKtask`` An optimization task.
* ``j:: Int32`` Index of the variable.
* ``name:: String`` The variable name.


Puts the name of a cone.

putconname
----------
::

    function putconname
    ( task:: MSKtask,
      i::    Int32,
      name:: String )


* ``task:: MSKtask`` An optimization task.
* ``i:: Int32`` Index of the variable.
* ``name:: String`` The variable name.


Puts the name of a constraint.

putcslice
---------
::

    function putcslice
    ( task::  MSKtask,
      first:: Int32,
      last::  Int32,
      slice:: Array{Float64} )


* ``task:: MSKtask`` An optimization task.
* ``first:: Int32`` First element in the slice of c.
* ``last:: Int32`` Last element plus 1 of the slice in c to be changed.
* ``slice:: Array{Float64}`` New numerical values for the objective coefficients that should be modified.


Modifies a slice of the linear objective coefficients.

putdouparam
-----------
::

    function putdouparam
    ( task::     MSKtask,
      param::    Int32,
      parvalue:: Float64 )


* ``task:: MSKtask`` An optimization task.
* ``param:: Int32`` (`Enum dparam`_) Which parameter.
* ``parvalue:: Float64`` Parameter value.


Sets a double parameter.

putintparam
-----------
::

    function putintparam
    ( task::     MSKtask,
      param::    Int32,
      parvalue:: Int32 )


* ``task:: MSKtask`` An optimization task.
* ``param:: Int32`` (`Enum iparam`_) Which parameter.
* ``parvalue:: Int32`` Parameter value.


Sets an integer parameter.

putmaxnumanz
------------
::

    function putmaxnumanz
    ( task::      MSKtask,
      maxnumanz:: Int64 )


* ``task:: MSKtask`` An optimization task.
* ``maxnumanz:: Int64`` New size of the storage reserved for storing the linear coefficient matrix.


The function changes the size of the preallocated storage for linear coefficients.

putmaxnumbarvar
---------------
::

    function putmaxnumbarvar
    ( task::         MSKtask,
      maxnumbarvar:: Int32 )


* ``task:: MSKtask`` An optimization task.
* ``maxnumbarvar:: Int32`` The maximum number of semidefinite variables.


Sets the number of preallocated symmetric matrix variables in the optimization task.

putmaxnumcon
------------
::

    function putmaxnumcon
    ( task::      MSKtask,
      maxnumcon:: Int32 )


* ``task:: MSKtask`` An optimization task.
* ``maxnumcon:: Int32`` Number of preallocated constraints in the optimization task.


Sets the number of preallocated constraints in the optimization task.

putmaxnumcone
-------------
::

    function putmaxnumcone
    ( task::       MSKtask,
      maxnumcone:: Int32 )


* ``task:: MSKtask`` An optimization task.
* ``maxnumcone:: Int32`` Number of preallocated conic constraints in the optimization task.


Sets the number of preallocated conic constraints in the optimization task.

putmaxnumqnz
------------
::

    function putmaxnumqnz
    ( task::      MSKtask,
      maxnumqnz:: Int64 )


* ``task:: MSKtask`` An optimization task.
* ``maxnumqnz:: Int64`` Number of non-zero elements preallocated in quadratic coefficient matrices.


Changes the size of the preallocated storage for quadratic terms.

putmaxnumvar
------------
::

    function putmaxnumvar
    ( task::      MSKtask,
      maxnumvar:: Int32 )


* ``task:: MSKtask`` An optimization task.
* ``maxnumvar:: Int32`` Number of preallocated variables in the optimization task.


Sets the number of preallocated variables in the optimization task.

putnadouparam
-------------
::

    function putnadouparam
    ( task::      MSKtask,
      paramname:: String,
      parvalue::  Float64 )


* ``task:: MSKtask`` An optimization task.
* ``paramname:: String`` Name of a parameter.
* ``parvalue:: Float64`` Parameter value.


Sets a double parameter.

putnaintparam
-------------
::

    function putnaintparam
    ( task::      MSKtask,
      paramname:: String,
      parvalue::  Int32 )


* ``task:: MSKtask`` An optimization task.
* ``paramname:: String`` Name of a parameter.
* ``parvalue:: Int32`` Parameter value.


Sets an integer parameter.

putnastrparam
-------------
::

    function putnastrparam
    ( task::      MSKtask,
      paramname:: String,
      parvalue::  String )


* ``task:: MSKtask`` An optimization task.
* ``paramname:: String`` Name of a parameter.
* ``parvalue:: String`` Parameter value.


Sets a string parameter.

putobjname
----------
::

    function putobjname
    ( task::    MSKtask,
      objname:: String )


* ``task:: MSKtask`` An optimization task.
* ``objname:: String`` Name of the objective.


Assigns a new name to the objective.

putobjsense
-----------
::

    function putobjsense
    ( task::  MSKtask,
      sense:: Int32 )


* ``task:: MSKtask`` An optimization task.
* ``sense:: Int32`` (`Enum objsense`_) The objective sense of the task


Sets the objective sense.

putparam
--------
::

    function putparam
    ( task::     MSKtask,
      parname::  String,
      parvalue:: String )


* ``task:: MSKtask`` An optimization task.
* ``parname:: String`` Parameter name.
* ``parvalue:: String`` Parameter value.


Modifies the value of parameter.

putqcon
-------
::

    function putqcon
    ( task::   MSKtask,
      qcsubk:: Array{Int32},
      qcsubi:: Array{Int32},
      qcsubj:: Array{Int32},
      qcval::  Array{Float64} )


* ``task:: MSKtask`` An optimization task.
* ``qcsubk:: Array{Int32}`` Constraint subscripts for quadratic coefficients.
* ``qcsubi:: Array{Int32}`` Row subscripts for quadratic constraint matrix.
* ``qcsubj:: Array{Int32}`` Column subscripts for quadratic constraint matrix.
* ``qcval:: Array{Float64}`` Quadratic constraint coefficient values.


Replaces all quadratic terms in constraints.

putqconk
--------
::

    function putqconk
    ( task::   MSKtask,
      k::      Int32,
      qcsubi:: Array{Int32},
      qcsubj:: Array{Int32},
      qcval::  Array{Float64} )

    function putqconk
    ( task::   MSKtask,
      k,
      Qk::     SparseMatrixCSC{Float64} )


* ``task:: MSKtask`` An optimization task.
* ``k:: Int32`` The constraint in which the new quadratic elements are inserted.
* ``qcsubi:: Array{Int32}`` Row subscripts for quadratic constraint matrix.
* ``qcsubj:: Array{Int32}`` Column subscripts for quadratic constraint matrix.
* ``qcval:: Array{Float64}`` Quadratic constraint coefficient values.
* ``Qk:: SparseMatrixCSC{Float64}`` The symmetric matrix 1/2 (Qk' + Qk) is used


Replaces all quadratic terms in a single constraint.

putqobj
-------
::

    function putqobj
    ( task::   MSKtask,
      qosubi:: Array{Int32},
      qosubj:: Array{Int32},
      qoval::  Array{Float64} )

    function putqobj
    ( task::   MSKtask,
      Qk::     SparseMatrixCSC{Float64} )


* ``task:: MSKtask`` An optimization task.
* ``qosubi:: Array{Int32}`` Row subscripts for quadratic objective coefficients.
* ``qosubj:: Array{Int32}`` Column subscripts for quadratic objective coefficients.
* ``qoval:: Array{Float64}`` Quadratic objective coefficient values.
* ``Qk:: SparseMatrixCSC{Float64}`` The symmetric matrix 1/2 (Qk' + Qk) is used


Replaces all quadratic terms in the objective.

putqobjij
---------
::

    function putqobjij
    ( task:: MSKtask,
      i::    Int32,
      j::    Int32,
      qoij:: Float64 )


* ``task:: MSKtask`` An optimization task.
* ``i:: Int32`` Row index for the coefficient to be replaced.
* ``j:: Int32`` Column index for the coefficient to be replaced.
* ``qoij:: Float64`` The new coefficient value.


Replaces one coefficient in the quadratic term in
the objective.

putskc
------
::

    function putskc
    ( task::     MSKtask,
      whichsol:: Int32,
      skc::      Array{Int32} )


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``skc:: Array{Int32}`` (`Enum stakey`_) Status keys for the constraints.


Sets the status keys for the constraints.

putskcslice
-----------
::

    function putskcslice
    ( task::     MSKtask,
      whichsol:: Int32,
      first::    Int32,
      last::     Int32,
      skc::      Array{Int32} )


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``first:: Int32`` First index in the sequence.
* ``last:: Int32`` Last index plus 1 in the sequence.
* ``skc:: Array{Int32}`` (`Enum stakey`_) Status keys for the constraints.


Sets the status keys for the constraints.

putskx
------
::

    function putskx
    ( task::     MSKtask,
      whichsol:: Int32,
      skx::      Array{Int32} )


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``skx:: Array{Int32}`` (`Enum stakey`_) Status keys for the variables.


Sets the status keys for the scalar variables.

putskxslice
-----------
::

    function putskxslice
    ( task::     MSKtask,
      whichsol:: Int32,
      first::    Int32,
      last::     Int32,
      skx::      Array{Int32} )


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``first:: Int32`` First index in the sequence.
* ``last:: Int32`` Last index plus 1 in the sequence.
* ``skx:: Array{Int32}`` (`Enum stakey`_) Status keys for the variables.


Sets the status keys for the variables.

putslc
------
::

    function putslc
    ( task::     MSKtask,
      whichsol:: Int32,
      slc::      Array{Float64} )


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``slc:: Array{Float64}`` The slc vector.


Sets the slc vector for a solution.

putslcslice
-----------
::

    function putslcslice
    ( task::     MSKtask,
      whichsol:: Int32,
      first::    Int32,
      last::     Int32,
      slc::      Array{Float64} )


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``first:: Int32`` First index in the sequence.
* ``last:: Int32`` Last index plus 1 in the sequence.
* ``slc:: Array{Float64}`` Dual variables corresponding to the lower bounds on the constraints.


Sets a slice of the slc vector for a solution.

putslx
------
::

    function putslx
    ( task::     MSKtask,
      whichsol:: Int32,
      slx::      Array{Float64} )


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``slx:: Array{Float64}`` The slx vector.


Sets the slx vector for a solution.

putslxslice
-----------
::

    function putslxslice
    ( task::     MSKtask,
      whichsol:: Int32,
      first::    Int32,
      last::     Int32,
      slx::      Array{Float64} )


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``first:: Int32`` First index in the sequence.
* ``last:: Int32`` Last index plus 1 in the sequence.
* ``slx:: Array{Float64}`` Dual variables corresponding to the lower bounds on the variables.


Sets a slice of the slx vector for a solution.

putsnx
------
::

    function putsnx
    ( task::     MSKtask,
      whichsol:: Int32,
      sux::      Array{Float64} )


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``sux:: Array{Float64}`` The snx vector.


Sets the snx vector for a solution.

putsnxslice
-----------
::

    function putsnxslice
    ( task::     MSKtask,
      whichsol:: Int32,
      first::    Int32,
      last::     Int32,
      snx::      Array{Float64} )


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``first:: Int32`` First index in the sequence.
* ``last:: Int32`` Last index plus 1 in the sequence.
* ``snx:: Array{Float64}`` Dual variables corresponding to the conic constraints on the variables.


Sets a slice of the snx vector for a solution.

putsolution
-----------
::

    function putsolution
    ( task::     MSKtask,
      whichsol:: Int32,
      skc::      Array{Int32},
      skx::      Array{Int32},
      skn::      Array{Int32},
      xc::       Array{Float64},
      xx::       Array{Float64},
      y::        Array{Float64},
      slc::      Array{Float64},
      suc::      Array{Float64},
      slx::      Array{Float64},
      sux::      Array{Float64},
      snx::      Array{Float64} )


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``skc:: Array{Int32}`` (`Enum stakey`_) Status keys for the constraints.
* ``skx:: Array{Int32}`` (`Enum stakey`_) Status keys for the variables.
* ``skn:: Array{Int32}`` (`Enum stakey`_) Status keys for the conic constraints.
* ``xc:: Array{Float64}`` Primal constraint solution.
* ``xx:: Array{Float64}`` Primal variable solution.
* ``y:: Array{Float64}`` Vector of dual variables corresponding to the constraints.
* ``slc:: Array{Float64}`` Dual variables corresponding to the lower bounds on the constraints.
* ``suc:: Array{Float64}`` Dual variables corresponding to the upper bounds on the constraints.
* ``slx:: Array{Float64}`` Dual variables corresponding to the lower bounds on the variables.
* ``sux:: Array{Float64}`` Dual variables corresponding to the upper bounds on the variables.
* ``snx:: Array{Float64}`` Dual variables corresponding to the conic constraints on the variables.


Inserts a solution.

putsolutioni
------------
::

    function putsolutioni
    ( task::     MSKtask,
      accmode::  Int32,
      i::        Int32,
      whichsol:: Int32,
      sk::       Int32,
      x::        Float64,
      sl::       Float64,
      su::       Float64,
      sn::       Float64 )


* ``task:: MSKtask`` An optimization task.
* ``accmode:: Int32`` (`Enum accmode`_) Defines whether solution information for a constraint or for a variable is modified.
* ``i:: Int32`` Index of the constraint or variable.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``sk:: Int32`` (`Enum stakey`_) Status key of the constraint or variable.
* ``x:: Float64`` Solution value of the primal constraint or variable.
* ``sl:: Float64`` Solution value of the dual variable associated with the lower bound.
* ``su:: Float64`` Solution value of the dual variable associated with the upper bound.
* ``sn:: Float64`` Solution value of the dual variable associated with the cone constraint.


Sets the primal and dual solution information for a single constraint or variable.

putsolutionyi
-------------
::

    function putsolutionyi
    ( task::     MSKtask,
      i::        Int32,
      whichsol:: Int32,
      y::        Float64 )


* ``task:: MSKtask`` An optimization task.
* ``i:: Int32`` Index of the dual variable.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``y:: Float64`` Solution value of the dual variable.


Inputs the dual variable of a solution.

putstrparam
-----------
::

    function putstrparam
    ( task::     MSKtask,
      param::    Int32,
      parvalue:: String )


* ``task:: MSKtask`` An optimization task.
* ``param:: Int32`` (`Enum sparam`_) Which parameter.
* ``parvalue:: String`` Parameter value.


Sets a string parameter.

putsuc
------
::

    function putsuc
    ( task::     MSKtask,
      whichsol:: Int32,
      suc::      Array{Float64} )


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``suc:: Array{Float64}`` The suc vector.


Sets the suc vector for a solution.

putsucslice
-----------
::

    function putsucslice
    ( task::     MSKtask,
      whichsol:: Int32,
      first::    Int32,
      last::     Int32,
      suc::      Array{Float64} )


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``first:: Int32`` First index in the sequence.
* ``last:: Int32`` Last index plus 1 in the sequence.
* ``suc:: Array{Float64}`` Dual variables corresponding to the upper bounds on the constraints.


Sets a slice of the suc vector for a solution.

putsux
------
::

    function putsux
    ( task::     MSKtask,
      whichsol:: Int32,
      sux::      Array{Float64} )


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``sux:: Array{Float64}`` The sux vector.


Sets the sux vector for a solution.

putsuxslice
-----------
::

    function putsuxslice
    ( task::     MSKtask,
      whichsol:: Int32,
      first::    Int32,
      last::     Int32,
      sux::      Array{Float64} )


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``first:: Int32`` First index in the sequence.
* ``last:: Int32`` Last index plus 1 in the sequence.
* ``sux:: Array{Float64}`` Dual variables corresponding to the upper bounds on the variables.


Sets a slice of the sux vector for a solution.

puttaskname
-----------
::

    function puttaskname
    ( task::     MSKtask,
      taskname:: String )


* ``task:: MSKtask`` An optimization task.
* ``taskname:: String`` Name assigned to the task.


Assigns a new name to the task.

putvarbound
-----------
::

    function putvarbound
    ( task:: MSKtask,
      j::    Int32,
      bk::   Int32,
      bl::   Float64,
      bu::   Float64 )


* ``task:: MSKtask`` An optimization task.
* ``j:: Int32`` Index of the variable.
* ``bk:: Int32`` (`Enum boundkey`_) New bound key.
* ``bl:: Float64`` New lower bound.
* ``bu:: Float64`` New upper bound.


Changes the bound for one variable.

putvarboundlist
---------------
::

    function putvarboundlist
    ( task:: MSKtask,
      sub::  Array{Int32},
      bkx::  Array{Int32},
      blx::  Array{Float64},
      bux::  Array{Float64} )


* ``task:: MSKtask`` An optimization task.
* ``sub:: Array{Int32}`` List of variable indexes.
* ``bkx:: Array{Int32}`` (`Enum boundkey`_) New bound keys.
* ``blx:: Array{Float64}`` New lower bound values.
* ``bux:: Array{Float64}`` New upper bounds values.


Changes the bounds of a list of variables.

putvarboundslice
----------------
::

    function putvarboundslice
    ( task::  MSKtask,
      first:: Int32,
      last::  Int32,
      bk::    Array{Int32},
      bl::    Array{Float64},
      bu::    Array{Float64} )


* ``task:: MSKtask`` An optimization task.
* ``first:: Int32`` Index of the first variable in the slice.
* ``last:: Int32`` Index of the last variable in the slice plus 1.
* ``bk:: Array{Int32}`` (`Enum boundkey`_) New bound keys.
* ``bl:: Array{Float64}`` New lower bounds.
* ``bu:: Array{Float64}`` New upper bounds.


Changes the bounds for a slice of the variables.

putvarbranchorder
-----------------
::

    function putvarbranchorder
    ( task::      MSKtask,
      j::         Int32,
      priority::  Int32,
      direction:: Int32 )


* ``task:: MSKtask`` An optimization task.
* ``j:: Int32`` Index of the variable.
* ``priority:: Int32`` The branching priority that should be assigned to the j'th variable.
* ``direction:: Int32`` (`Enum branchdir`_) Specifies the preferred branching direction for the j'th variable.


Assigns a branching priority and direction to a variable.

putvarname
----------
::

    function putvarname
    ( task:: MSKtask,
      j::    Int32,
      name:: String )


* ``task:: MSKtask`` An optimization task.
* ``j:: Int32`` Index of the variable.
* ``name:: String`` The variable name.


Puts the name of a variable.

putvartype
----------
::

    function putvartype
    ( task::    MSKtask,
      j::       Int32,
      vartype:: Int32 )


* ``task:: MSKtask`` An optimization task.
* ``j:: Int32`` Index of the variable.
* ``vartype:: Int32`` (`Enum variabletype`_) The new variable type.


Sets the variable type of one variable.

putvartypelist
--------------
::

    function putvartypelist
    ( task::    MSKtask,
      subj::    Array{Int32},
      vartype:: Array{Int32} )


* ``task:: MSKtask`` An optimization task.
* ``subj:: Array{Int32}`` A list of variable indexes for which the variable
                           type should be changed.
* ``vartype:: Array{Int32}`` (`Enum variabletype`_) A list of variable types.


Sets the variable type for one or more variables.

putxc
-----
::

    function putxc
    ( task::     MSKtask,
      whichsol:: Int32 )
    -> xc


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``xc:: Array{Float64}`` The xc vector.


Sets the xc vector for a solution.

putxcslice
----------
::

    function putxcslice
    ( task::     MSKtask,
      whichsol:: Int32,
      first::    Int32,
      last::     Int32,
      xc::       Array{Float64} )


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``first:: Int32`` First index in the sequence.
* ``last:: Int32`` Last index plus 1 in the sequence.
* ``xc:: Array{Float64}`` Primal constraint solution.


Sets a slice of the xc vector for a solution.

putxx
-----
::

    function putxx
    ( task::     MSKtask,
      whichsol:: Int32,
      xx::       Array{Float64} )


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``xx:: Array{Float64}`` The xx vector.


Sets the xx vector for a solution.

putxxslice
----------
::

    function putxxslice
    ( task::     MSKtask,
      whichsol:: Int32,
      first::    Int32,
      last::     Int32,
      xx::       Array{Float64} )


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``first:: Int32`` First index in the sequence.
* ``last:: Int32`` Last index plus 1 in the sequence.
* ``xx:: Array{Float64}`` Primal variable solution.


Obtains a slice of the xx vector for a solution.

puty
----
::

    function puty
    ( task::     MSKtask,
      whichsol:: Int32,
      y::        Array{Float64} )


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``y:: Array{Float64}`` The y vector.


Sets the y vector for a solution.

putyslice
---------
::

    function putyslice
    ( task::     MSKtask,
      whichsol:: Int32,
      first::    Int32,
      last::     Int32,
      y::        Array{Float64} )


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``first:: Int32`` First index in the sequence.
* ``last:: Int32`` Last index plus 1 in the sequence.
* ``y:: Array{Float64}`` Vector of dual variables corresponding to the constraints.


Sets a slice of the y vector for a solution.

readbranchpriorities
--------------------
::

    function readbranchpriorities
    ( task::     MSKtask,
      filename:: String )


* ``task:: MSKtask`` An optimization task.
* ``filename:: String`` Input file name.


Reads branching priority data from a file.

readdata
--------
::

    function readdata
    ( task::     MSKtask,
      filename:: String )


* ``task:: MSKtask`` An optimization task.
* ``filename:: String`` Input data file name.


Reads problem data from a file.

readdataformat
--------------
::

    function readdataformat
    ( task::     MSKtask,
      filename:: String,
      format::   Int32,
      compress:: Int32 )


* ``task:: MSKtask`` An optimization task.
* ``filename:: String`` Input data file name.
* ``format:: Int32`` (`Enum dataformat`_) File data format.
* ``compress:: Int32`` (`Enum compresstype`_) File compression type.


Reads problem data from a file.

readparamfile
-------------
::

    function readparamfile(task:: MSKtask)


* ``task:: MSKtask`` An optimization task.


Reads a parameter file.

readsolution
------------
::

    function readsolution
    ( task::     MSKtask,
      whichsol:: Int32,
      filename:: String )


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``filename:: String`` A valid file name.


Reads a solution from a file.

readsummary
-----------
::

    function readsummary
    ( task::        MSKtask,
      whichstream:: Int32 )


* ``task:: MSKtask`` An optimization task.
* ``whichstream:: Int32`` (`Enum streamtype`_) Index of the stream.


Prints information about last file read.

readtask
--------
::

    function readtask
    ( task::     MSKtask,
      filename:: String )


* ``task:: MSKtask`` An optimization task.
* ``filename:: String`` Input file name.


Load task data from a file.

relaxprimal
-----------
::

    function relaxprimal
    ( task:: MSKtask,
      wlc::  Array{Float64},
      wuc::  Array{Float64},
      wlx::  Array{Float64},
      wux::  Array{Float64} )
    -> relaxedtask


* ``task:: MSKtask`` An optimization task.
* ``relaxedtask:: MSKtask`` The returned task.
* ``wlc:: Array{Float64}`` Weights associated with lower bounds on the activity of constraints.
* ``wuc:: Array{Float64}`` Weights associated with upper bounds on
                            the activity of constraints.
* ``wlx:: Array{Float64}`` Weights associated with lower bounds on
                            the activity of variables.
* ``wux:: Array{Float64}`` Weights associated with upper bounds on
                            the activity of variables.


Deprecated.

removebarvars
-------------
::

    function removebarvars
    ( task::   MSKtask,
      subset:: Array{Int32} )


* ``task:: MSKtask`` An optimization task.
* ``subset:: Array{Int32}`` Indexes of symmetric matrix which should be removed.


The function removes a number of symmetric matrix.

removecones
-----------
::

    function removecones
    ( task::   MSKtask,
      subset:: Array{Int32} )


* ``task:: MSKtask`` An optimization task.
* ``subset:: Array{Int32}`` Indexes of cones which should be removed.


Removes a conic constraint from the problem.

removecons
----------
::

    function removecons
    ( task::   MSKtask,
      subset:: Array{Int32} )


* ``task:: MSKtask`` An optimization task.
* ``subset:: Array{Int32}`` Indexes of constraints which should be removed.


The function removes a number of constraints.

removevars
----------
::

    function removevars
    ( task::   MSKtask,
      subset:: Array{Int32} )


* ``task:: MSKtask`` An optimization task.
* ``subset:: Array{Int32}`` Indexes of variables which should be removed.


The function removes a number of variables.

resizetask
----------
::

    function resizetask
    ( task::       MSKtask,
      maxnumcon::  Int32,
      maxnumvar::  Int32,
      maxnumcone:: Int32,
      maxnumanz::  Int64,
      maxnumqnz::  Int64 )


* ``task:: MSKtask`` An optimization task.
* ``maxnumcon:: Int32`` New maximum number of constraints.
* ``maxnumvar:: Int32`` New maximum number of variables.
* ``maxnumcone:: Int32`` New maximum number of cones.
* ``maxnumanz:: Int64`` New maximum number of linear non-zero elements.
* ``maxnumqnz:: Int64`` New maximum number of quadratic non-zeros elements.


Resizes an optimization task.

sensitivityreport
-----------------
::

    function sensitivityreport
    ( task::        MSKtask,
      whichstream:: Int32 )


* ``task:: MSKtask`` An optimization task.
* ``whichstream:: Int32`` (`Enum streamtype`_) Index of the stream.


Creates a sensitivity report.

setdefaults
-----------
::

    function setdefaults(task:: MSKtask)


* ``task:: MSKtask`` An optimization task.


Resets all parameters values.

solutiondef
-----------
::

    function solutiondef
    ( task::     MSKtask,
      whichsol:: Int32 )
    -> isdef


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``isdef:: Bool`` Is non-zero if the requested solution is defined.


Checks whether a solution is defined.

solutionsummary
---------------
::

    function solutionsummary
    ( task::        MSKtask,
      whichstream:: Int32 )


* ``task:: MSKtask`` An optimization task.
* ``whichstream:: Int32`` (`Enum streamtype`_) Index of the stream.


Prints a short summary of the current solutions.

solvewithbasis
--------------
::

    function solvewithbasis
    ( task::   MSKtask,
      transp:: Int32,
      numnz::  Int32,
      sub::    Array{Int32},
      val::    Array{Float64} )
    -> numnz


* ``task:: MSKtask`` An optimization task.
* ``transp:: Int32`` Controls which problem formulation is solved.
* ``numnz:: Int32`` Input (number of non-zeros in right-hand side) and output (number of non-zeros in solution vector).
* ``sub:: Array{Int32}`` Input (indexes of non-zeros in right-hand side) and output (indexes of non-zeros in solution vector).
* ``val:: Array{Float64}`` Input (right-hand side values) and output (solution vector values).


Solve a linear equation system involving a basis matrix.

startstat
---------
::

    function startstat(task:: MSKtask)


* ``task:: MSKtask`` An optimization task.


Starts the statistics file.

stopstat
--------
::

    function stopstat(task:: MSKtask)


* ``task:: MSKtask`` An optimization task.


Stops the statistics file.

strtoconetype
-------------
::

    function strtoconetype
    ( task:: MSKtask,
      str::  String )
    -> conetype


* ``task:: MSKtask`` An optimization task.
* ``str:: String`` String corresponding to the cone type code.
* ``conetype:: Int32`` The cone type corresponding to str.


Obtains a cone type code.

strtosk
-------
::

    function strtosk
    ( task:: MSKtask,
      str::  String )
    -> sk


* ``task:: MSKtask`` An optimization task.
* ``str:: String`` Status key string.
* ``sk:: Int32`` Status key corresponding to the string.


Obtains a status key.

updatesolutioninfo
------------------
::

    function updatesolutioninfo
    ( task::     MSKtask,
      whichsol:: Int32 )


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.


Update the information items related to the solution.

writebranchpriorities
---------------------
::

    function writebranchpriorities
    ( task::     MSKtask,
      filename:: String )


* ``task:: MSKtask`` An optimization task.
* ``filename:: String`` Output file name.


Writes branching priority data to a file.

writedata
---------
::

    function writedata
    ( task::     MSKtask,
      filename:: String )


* ``task:: MSKtask`` An optimization task.
* ``filename:: String`` Output file name.


Writes problem data to a file.

writeparamfile
--------------
::

    function writeparamfile
    ( task::     MSKtask,
      filename:: String )


* ``task:: MSKtask`` An optimization task.
* ``filename:: String`` The name of parameter file.


Writes all the parameters to a parameter file.

writesolution
-------------
::

    function writesolution
    ( task::     MSKtask,
      whichsol:: Int32,
      filename:: String )


* ``task:: MSKtask`` An optimization task.
* ``whichsol:: Int32`` (`Enum soltype`_) Selects a solution.
* ``filename:: String`` A valid file name.


Write a solution to a file.

writetask
---------
::

    function writetask
    ( task::     MSKtask,
      filename:: String )


* ``task:: MSKtask`` An optimization task.
* ``filename:: String`` Output file name.


Write a complete binary dump of the task data.

checkinlicense
--------------
::

    function checkinlicense
    ( env::     MSKenv,
      feature:: Int32 )


* ``env:: MSKenv`` The MOSEK environment.
* ``feature:: Int32`` (`Enum feature`_) Feature to check in to the license system.


Check in a license feature from the license server ahead of time.

checkoutlicense
---------------
::

    function checkoutlicense
    ( env::     MSKenv,
      feature:: Int32 )


* ``env:: MSKenv`` The MOSEK environment.
* ``feature:: Int32`` (`Enum feature`_) Feature to check out from the license system.


Check out a license feature from the license server ahead of time.

echointro
---------
::

    function echointro
    ( env::     MSKenv,
      longver:: Int32 )


* ``env:: MSKenv`` The MOSEK environment.
* ``longver:: Int32`` If non-zero, then the intro is slightly longer.


Prints an intro to message stream.

getcodedesc
-----------
::

    function getcodedesc(code:: Int32)
    -> symname,str


* ``code:: Int32`` (`Enum rescode`_) A valid response code.
* ``symname:: String`` Symbolic name corresponding to the code.
* ``str:: String`` Obtains a short description of a response code.


Obtains a short description of a response code.

getversion
----------
::

    function getversion
    (  )
    -> major,minor,build,revision


* ``major:: Int32`` Major version number.
* ``minor:: Int32`` Minor version number.
* ``build:: Int32`` Build number.
* ``revision:: Int32`` Revision number.


Obtains MOSEK version information.

licensecleanup
--------------
::

    function licensecleanup
    (  )





Stops all threads and delete all handles used by the license system.

linkfiletostream
----------------
::

    function linkfiletostream
    ( env::         MSKenv,
      whichstream:: Int32,
      filename::    String,
      append::      Int32 )


* ``env:: MSKenv`` The MOSEK environment.
* ``whichstream:: Int32`` (`Enum streamtype`_) Index of the stream.
* ``filename:: String`` Name of the file to write stream data to.
* ``append:: Int32`` If this argument is non-zero, the output is appended to the file.


Directs all output from a stream to a file.

putdllpath
----------
::

    function putdllpath
    ( env::     MSKenv,
      dllpath:: String )


* ``env:: MSKenv`` The MOSEK environment.
* ``dllpath:: String`` A path to the dynamic MOSEK libraries.


Sets the path to the DLL/shared libraries that MOSEK is loading.

putkeepdlls
-----------
::

    function putkeepdlls
    ( env::      MSKenv,
      keepdlls:: Int32 )


* ``env:: MSKenv`` The MOSEK environment.
* ``keepdlls:: Int32`` Controls whether explicitly loaded DLLs should be kept.


Controls whether explicitly loaded DLLs should be kept.

putlicensecode
--------------
::

    function putlicensecode
    ( env::  MSKenv,
      code:: Array{Int32} )


* ``env:: MSKenv`` The MOSEK environment.
* ``code:: Array{Int32}`` A license key string.


The purpose of this function is to input a runtime license code.

putlicensedebug
---------------
::

    function putlicensedebug
    ( env::      MSKenv,
      licdebug:: Int32 )


* ``env:: MSKenv`` The MOSEK environment.
* ``licdebug:: Int32`` Enable output of license check-out debug information.


Enables debug information for the license system.

putlicensepath
--------------
::

    function putlicensepath
    ( env::         MSKenv,
      licensepath:: String )


* ``env:: MSKenv`` The MOSEK environment.
* ``licensepath:: String`` A path specifycing where to search for the license.


Set the path to the license file.

putlicensewait
--------------
::

    function putlicensewait
    ( env::     MSKenv,
      licwait:: Int32 )


* ``env:: MSKenv`` The MOSEK environment.
* ``licwait:: Int32`` Enable waiting for a license.


Control whether mosek should wait for an available license if no license is available.




Mosek.jl Constants
Enum solveform
--------------

MSK_SOLVE_DUAL = 2

    The optimizer should solve the dual problem.

MSK_SOLVE_FREE = 0

    The optimizer is free to solve either the primal or
    the dual problem.

MSK_SOLVE_PRIMAL = 1

    The optimizer should solve the primal problem.

Enum problemitem
----------------

MSK_PI_CON = 1

    Item is a constraint.

MSK_PI_CONE = 2

    Item is a cone.

MSK_PI_VAR = 0

    Item is a variable.

Enum accmode
------------

MSK_ACC_CON = 1

    Access data by rows (constraint oriented)

MSK_ACC_VAR = 0

    Access data by columns (variable oriented)

Enum blastranspose
------------------

MSK_BLAS_TRANSPOSE_NO = 0

    No transpose is applied.

MSK_BLAS_TRANSPOSE_YES = 1

    A transpose is applied.

Enum intpnthotstart
-------------------

MSK_INTPNT_HOTSTART_DUAL = 2

    The interior-point optimizer exploits the dual solution only.

MSK_INTPNT_HOTSTART_NONE = 0

    The interior-point optimizer performs a coldstart.

MSK_INTPNT_HOTSTART_PRIMAL = 1

    The interior-point optimizer exploits the primal solution only.

MSK_INTPNT_HOTSTART_PRIMAL_DUAL = 3

    The interior-point optimizer exploits both the primal and dual solution.

Enum sparam
-----------

MSK_SPAR_BAS_SOL_FILE_NAME = 0

    Name of the bas solution file.

MSK_SPAR_DATA_FILE_NAME = 1

    Data are read and written to this file.

MSK_SPAR_DEBUG_FILE_NAME = 2

    MOSEK debug file.

MSK_SPAR_FEASREPAIR_NAME_PREFIX = 3

    Feasibility repair name prefix.

MSK_SPAR_FEASREPAIR_NAME_SEPARATOR = 4

    Feasibility repair name separator.

MSK_SPAR_FEASREPAIR_NAME_WSUMVIOL = 5

    Feasibility repair name violation name.

MSK_SPAR_INT_SOL_FILE_NAME = 6

    Name of the int solution file.

MSK_SPAR_ITR_SOL_FILE_NAME = 7

    Name of the itr solution file.

MSK_SPAR_MIO_DEBUG_STRING = 8

    For internal use only.

MSK_SPAR_PARAM_COMMENT_SIGN = 9

    Solution file comment character.

MSK_SPAR_PARAM_READ_FILE_NAME = 10

    Modifications to the parameter
    database is read from this file.

MSK_SPAR_PARAM_WRITE_FILE_NAME = 11

    The parameter database is written to this file.

MSK_SPAR_READ_MPS_BOU_NAME = 12

    Name of the BOUNDS vector used.
    An empty name means that the first BOUNDS vector is used.

MSK_SPAR_READ_MPS_OBJ_NAME = 13

    Objective name in the MPS file.

MSK_SPAR_READ_MPS_RAN_NAME = 14

    Name of the RANGE vector  used.
    An empty name means that the first RANGE vector is used.

MSK_SPAR_READ_MPS_RHS_NAME = 15

    Name of the RHS used.
    An empty name means that the first RHS vector is used.

MSK_SPAR_SENSITIVITY_FILE_NAME = 16

    Sensitivity report file name.

MSK_SPAR_SENSITIVITY_RES_FILE_NAME = 17

    Name of the sensitivity report output file.

MSK_SPAR_SOL_FILTER_XC_LOW = 18

    Solution file filter.

MSK_SPAR_SOL_FILTER_XC_UPR = 19

    Solution file filter.

MSK_SPAR_SOL_FILTER_XX_LOW = 20

    Solution file filter.

MSK_SPAR_SOL_FILTER_XX_UPR = 21

    Solution file filter.

MSK_SPAR_STAT_FILE_NAME = 22

    Statistics file name.

MSK_SPAR_STAT_KEY = 23

    Key used when writing the summary file.

MSK_SPAR_STAT_NAME = 24

    Name used when writing the statistics file.

MSK_SPAR_WRITE_LP_GEN_VAR_NAME = 25

    Added variable names in the LP files.

Enum iparam
-----------

MSK_IPAR_ALLOC_ADD_QNZ = 0

    Controls how the quadratic matrixes are extended.

MSK_IPAR_ANA_SOL_BASIS = 1

    Controls whether the basis matrix is analyzed in solution analyzer.

MSK_IPAR_ANA_SOL_PRINT_VIOLATED = 2

    Controls whether a list of violated constraints is printed.

MSK_IPAR_AUTO_SORT_A_BEFORE_OPT = 3

    Controls whether the elements in each column of A are sorted before an optimization is performed.

MSK_IPAR_AUTO_UPDATE_SOL_INFO = 4

    Controls whether the solution information items are automatically updated after an optimization is performed.

MSK_IPAR_BASIS_SOLVE_USE_PLUS_ONE = 5

    Controls the sign of the columns in the basis matrix corresponding to slack variables.

MSK_IPAR_BI_CLEAN_OPTIMIZER = 6

    Controls which simplex optimizer is used in the clean-up phase.

MSK_IPAR_BI_IGNORE_MAX_ITER = 7

    Turns on basis identification in case the interior-point optimizer is terminated due to maximum number of iterations.

MSK_IPAR_BI_IGNORE_NUM_ERROR = 8

    Turns on basis identification in case the interior-point optimizer is terminated due to a numerical problem.

MSK_IPAR_BI_MAX_ITERATIONS = 9

    Maximum number of iterations after basis identification.

MSK_IPAR_CACHE_LICENSE = 10

    Control license caching.

MSK_IPAR_CHECK_CONVEXITY = 11

    Specify the level of convexity check on quadratic problems

MSK_IPAR_COMPRESS_STATFILE = 12

    Control compression of stat files.

MSK_IPAR_CONCURRENT_NUM_OPTIMIZERS = 13

    The maximum number of simultaneous optimizations that will be started
    by the concurrent optimizer.

MSK_IPAR_CONCURRENT_PRIORITY_DUAL_SIMPLEX = 14

    Priority of the dual simplex algorithm when selecting solvers for
    concurrent optimization.

MSK_IPAR_CONCURRENT_PRIORITY_FREE_SIMPLEX = 15

    Priority of the free simplex optimizer when selecting solvers for
    concurrent optimization.

MSK_IPAR_CONCURRENT_PRIORITY_INTPNT = 16

    Priority of the interior-point algorithm when selecting solvers for
    concurrent optimization.

MSK_IPAR_CONCURRENT_PRIORITY_PRIMAL_SIMPLEX = 17

    Priority of the primal simplex algorithm when selecting solvers for
    concurrent optimization.

MSK_IPAR_FEASREPAIR_OPTIMIZE = 18

    Controls which type of feasibility analysis is to be performed.

MSK_IPAR_INFEAS_GENERIC_NAMES = 19

    Controls the contents of the infeasibility report.

MSK_IPAR_INFEAS_PREFER_PRIMAL = 20

    Controls which certificate is used if both primal- and dual- certificate of infeasibility is available.

MSK_IPAR_INFEAS_REPORT_AUTO = 21

    Turns the feasibility report on or off.

MSK_IPAR_INFEAS_REPORT_LEVEL = 22

    Controls the contents of the infeasibility report.

MSK_IPAR_INTPNT_BASIS = 23

    Controls whether basis identification is performed.

MSK_IPAR_INTPNT_DIFF_STEP = 24

    Controls whether different step sizes
    are allowed in the primal and dual space.

MSK_IPAR_INTPNT_FACTOR_DEBUG_LVL = 25

    Controls factorization debug level.

MSK_IPAR_INTPNT_FACTOR_METHOD = 26

    Controls the method used to factor the Newton equation system.

MSK_IPAR_INTPNT_HOTSTART = 27

    Currently not in use.

MSK_IPAR_INTPNT_MAX_ITERATIONS = 28

    Controls the maximum number of iterations
    allowed in the interior-point optimizer.

MSK_IPAR_INTPNT_MAX_NUM_COR = 29

    Maximum number of correction steps.

MSK_IPAR_INTPNT_MAX_NUM_REFINEMENT_STEPS = 30

    Maximum number of steps to be used by the iterative
    search direction refinement.

MSK_IPAR_INTPNT_OFF_COL_TRH = 31

    Controls the aggressiveness of the offending column detection.

MSK_IPAR_INTPNT_ORDER_METHOD = 32

    Controls the ordering strategy.

MSK_IPAR_INTPNT_REGULARIZATION_USE = 33

    Controls whether regularization is allowed.

MSK_IPAR_INTPNT_SCALING = 34

    Controls how the problem is scaled
    before the interior-point optimizer
    is used.

MSK_IPAR_INTPNT_SOLVE_FORM = 35

    Controls whether the primal
    or the dual problem is solved.

MSK_IPAR_INTPNT_STARTING_POINT = 36

    Starting point used by the interior-point optimizer.

MSK_IPAR_LIC_TRH_EXPIRY_WRN = 37

    Controls when expiry warnings are issued.

MSK_IPAR_LICENSE_ALLOW_OVERUSE = 38

    Controls if license overuse is allowed when caching licenses

MSK_IPAR_LICENSE_DEBUG = 39

    Controls the license manager client debugging behavior.

MSK_IPAR_LICENSE_PAUSE_TIME = 40

    Controls license manager client behavior.

MSK_IPAR_LICENSE_SUPPRESS_EXPIRE_WRNS = 41

    Controls license manager client behavior.

MSK_IPAR_LICENSE_WAIT = 42

    Controls if MOSEK should queue for a license if none is available.

MSK_IPAR_LOG = 43

    Controls the amount of log information.

MSK_IPAR_LOG_BI = 44

    Controls the amount of output printed
    by the basis identification procedure. A higher level implies that more information is logged.

MSK_IPAR_LOG_BI_FREQ = 45

    Controls the logging frequency.

MSK_IPAR_LOG_CHECK_CONVEXITY = 46

    Controls logging in convexity check on quadratic problems.
    Set to a positive value to turn logging on.
    
    If a quadratic coefficient matrix is found to violate the requirement of PSD (NSD)
    then a list of negative (positive) pivot elements is printed. The absolute value of the pivot elements
    is also shown.

MSK_IPAR_LOG_CONCURRENT = 47

    Controls amount of output printed
    by the concurrent optimizer.

MSK_IPAR_LOG_CUT_SECOND_OPT = 48

    Controls the reduction in the log levels for the second and any subsequent optimizations.

MSK_IPAR_LOG_EXPAND = 49

    Controls the amount of logging when a data item such as the maximum number constrains is expanded.

MSK_IPAR_LOG_FACTOR = 50

    If turned on, then the factor log lines are added to the log.

MSK_IPAR_LOG_FEAS_REPAIR = 51

    Controls the amount of output printed when performing feasibility repair. A value higher than one means extensive logging.

MSK_IPAR_LOG_FILE = 52

    If turned on, then some log info is printed when a file is written or read.

MSK_IPAR_LOG_HEAD = 53

    If turned on, then a header line is added to the log.

MSK_IPAR_LOG_INFEAS_ANA = 54

    Controls log level for the infeasibility analyzer.

MSK_IPAR_LOG_INTPNT = 55

    Controls the amount of log information from the interior-point optimizers.

MSK_IPAR_LOG_MIO = 56

    Controls the amount of log information from the mixed-integer optimizers.

MSK_IPAR_LOG_MIO_FREQ = 57

    The mixed-integer solver logging frequency.

MSK_IPAR_LOG_NONCONVEX = 58

    Controls amount of output printed by the nonconvex optimizer.

MSK_IPAR_LOG_OPTIMIZER = 59

    Controls the amount of general optimizer information that is logged.

MSK_IPAR_LOG_ORDER = 60

    If turned on, then factor lines are added to the log.

MSK_IPAR_LOG_PARAM = 61

    Controls the amount of information printed out about parameter changes.

MSK_IPAR_LOG_PRESOLVE = 62

    Controls amount of output printed
    by the presolve procedure. A higher level implies that more information is logged.

MSK_IPAR_LOG_RESPONSE = 63

    Controls amount of output printed when response codes are reported. A higher level implies that more information is logged.

MSK_IPAR_LOG_SENSITIVITY = 64

    Control logging in sensitivity analyzer.

MSK_IPAR_LOG_SENSITIVITY_OPT = 65

    Control logging in sensitivity analyzer.

MSK_IPAR_LOG_SIM = 66

    Controls the amount of log information from the simplex optimizers.

MSK_IPAR_LOG_SIM_FREQ = 67

    Controls simplex logging frequency.

MSK_IPAR_LOG_SIM_MINOR = 68

    Currently not in use.

MSK_IPAR_LOG_SIM_NETWORK_FREQ = 69

    Controls the network simplex logging frequency.

MSK_IPAR_LOG_STORAGE = 70

    Controls the memory related log information.

MSK_IPAR_MAX_NUM_WARNINGS = 71

    Warning level. A higher value results in more warnings.

MSK_IPAR_MIO_BRANCH_DIR = 72

    Controls whether the mixed-integer optimizer is branching up or down by default.

MSK_IPAR_MIO_BRANCH_PRIORITIES_USE = 73

    Controls whether branching priorities are used by the mixed-integer optimizer.

MSK_IPAR_MIO_CONSTRUCT_SOL = 74

    Controls if an initial mixed integer solution should be constructed from the values of the integer variables.

MSK_IPAR_MIO_CONT_SOL = 75

    Controls the meaning of interior-point and basic solutions in mixed integer problems.

MSK_IPAR_MIO_CUT_LEVEL_ROOT = 76

    Controls the cut level employed by the mixed-integer optimizer at the root node.

MSK_IPAR_MIO_CUT_LEVEL_TREE = 77

    Controls the cut level employed by the mixed-integer optimizer in the tree.

MSK_IPAR_MIO_FEASPUMP_LEVEL = 78

    Controls the feasibility pump heuristic which is used to construct a good initial feasible solution.

MSK_IPAR_MIO_HEURISTIC_LEVEL = 79

    Controls the heuristic employed by the mixed-integer
    optimizer to locate an initial integer feasible
    solution.

MSK_IPAR_MIO_HOTSTART = 80

    Controls whether the integer optimizer is hot-started.

MSK_IPAR_MIO_KEEP_BASIS = 81

    Controls whether the integer presolve keeps bases in memory.

MSK_IPAR_MIO_LOCAL_BRANCH_NUMBER = 82

    Controls the size of the local search space when doing local branching.

MSK_IPAR_MIO_MAX_NUM_BRANCHES = 83

    Maximum number of branches allowed during the branch and bound search.

MSK_IPAR_MIO_MAX_NUM_RELAXS = 84

    Maximum number of relaxations in branch and bound search.

MSK_IPAR_MIO_MAX_NUM_SOLUTIONS = 85

    Controls how many feasible solutions the mixed-integer optimizer investigates.

MSK_IPAR_MIO_MODE = 86

    Turns on/off the mixed-integer mode.

MSK_IPAR_MIO_MT_USER_CB = 87

    It true user callbacks are called from each thread used by this optimizer. If false the user callback is only called from a single thread.

MSK_IPAR_MIO_NODE_OPTIMIZER = 88

    Controls which optimizer is employed at the non-root nodes in the mixed-integer optimizer.

MSK_IPAR_MIO_NODE_SELECTION = 89

    Controls the node selection strategy employed by the
    mixed-integer optimizer.

MSK_IPAR_MIO_OPTIMIZER_MODE = 90

    An experimental feature.

MSK_IPAR_MIO_PRESOLVE_AGGREGATE = 91

    Controls whether problem aggregation is performed in the mixed-integer presolve.

MSK_IPAR_MIO_PRESOLVE_PROBING = 92

    Controls whether probing is employed by the mixed-integer presolve.

MSK_IPAR_MIO_PRESOLVE_USE = 93

    Controls whether presolve is performed by the mixed-integer optimizer.

MSK_IPAR_MIO_RINS_MAX_NODES = 94

    Maximum number of nodes in each call to RINS.

MSK_IPAR_MIO_ROOT_OPTIMIZER = 95

    Controls which optimizer is employed at the root node in the mixed-integer optimizer.

MSK_IPAR_MIO_STRONG_BRANCH = 96

    The depth from the root in which strong branching is employed.

MSK_IPAR_MIO_USE_MULTITHREADED_OPTIMIZER = 97

    Controls whether the new multithreaded optimizer should be used for Mixed integer problems.

MSK_IPAR_MT_SPINCOUNT = 98

    Set the number of iterations to spin before sleeping.

MSK_IPAR_NONCONVEX_MAX_ITERATIONS = 99

    Maximum number of iterations that can be used by the nonconvex optimizer.

MSK_IPAR_NUM_THREADS = 100

    Controls the number of threads employed by the optimizer. If set to 0 the number of threads used will
    be equal to the number of cores detected on the machine.

MSK_IPAR_OPF_MAX_TERMS_PER_LINE = 101

    The maximum number of terms (linear and quadratic) per line when an OPF file is written.

MSK_IPAR_OPF_WRITE_HEADER = 102

    Write a text header with date and MOSEK version in an OPF file.

MSK_IPAR_OPF_WRITE_HINTS = 103

    Write a hint section with problem dimensions in the beginning of an OPF file.

MSK_IPAR_OPF_WRITE_PARAMETERS = 104

    Write a parameter section in an OPF file.

MSK_IPAR_OPF_WRITE_PROBLEM = 105

    Write objective, constraints, bounds etc. to an OPF file.

MSK_IPAR_OPF_WRITE_SOL_BAS = 106

    Controls what is written to the OPF files.

MSK_IPAR_OPF_WRITE_SOL_ITG = 107

    Controls what is written to the OPF files.

MSK_IPAR_OPF_WRITE_SOL_ITR = 108

    Controls what is written to the OPF files.

MSK_IPAR_OPF_WRITE_SOLUTIONS = 109

    Enable inclusion of solutions in the OPF files.

MSK_IPAR_OPTIMIZER = 110

    Controls which optimizer is used to optimize the task.

MSK_IPAR_PARAM_READ_CASE_NAME = 111

    If turned on, then names in the parameter file are case sensitive.

MSK_IPAR_PARAM_READ_IGN_ERROR = 112

    If turned on, then errors in parameter settings is ignored.

MSK_IPAR_PRESOLVE_ELIM_FILL = 113

    Maximum amount of fill-in in the elimination phase.

MSK_IPAR_PRESOLVE_ELIMINATOR_MAX_NUM_TRIES = 114

    Control the maximum number of times the eliminator is tried.

MSK_IPAR_PRESOLVE_ELIMINATOR_USE = 115

    Controls whether free or implied free
    variables are eliminated from the problem.

MSK_IPAR_PRESOLVE_LEVEL = 116

    Currently not used.

MSK_IPAR_PRESOLVE_LINDEP_ABS_WORK_TRH = 117

    Controls linear dependency check in presolve.

MSK_IPAR_PRESOLVE_LINDEP_REL_WORK_TRH = 118

    Controls linear dependency check in presolve.

MSK_IPAR_PRESOLVE_LINDEP_USE = 119

    Controls whether the linear constraints are checked for linear dependencies.

MSK_IPAR_PRESOLVE_MAX_NUM_REDUCTIONS = 120

    Controls the maximum number reductions performed by the presolve.

MSK_IPAR_PRESOLVE_USE = 121

    Controls whether the presolve is applied to a problem before it is optimized.

MSK_IPAR_PRIMAL_REPAIR_OPTIMIZER = 122

    Controls which optimizer that is used to find the optimal repair.

MSK_IPAR_QO_SEPARABLE_REFORMULATION = 123

    Determine if quadratic programing problems should be reformulated to separable form.

MSK_IPAR_READ_DATA_COMPRESSED = 124

    Controls the input file decompression.

MSK_IPAR_READ_DATA_FORMAT = 125

    Format of the data file to be read.

MSK_IPAR_READ_DEBUG = 126

    Turns on additional debugging information when reading files.

MSK_IPAR_READ_KEEP_FREE_CON = 127

    Controls whether the free constraints are included in
    the problem.

MSK_IPAR_READ_LP_DROP_NEW_VARS_IN_BOU = 128

    Controls how the LP files are interpreted.

MSK_IPAR_READ_LP_QUOTED_NAMES = 129

    If a name is in quotes when reading an LP file, the quotes will be removed.

MSK_IPAR_READ_MPS_FORMAT = 130

    Controls how strictly the MPS file reader interprets the MPS format.

MSK_IPAR_READ_MPS_KEEP_INT = 131

    Controls if integer constraints are read.

MSK_IPAR_READ_MPS_RELAX = 132

    Controls the meaning of integer constraints.

MSK_IPAR_READ_MPS_WIDTH = 133

    Controls the maximal number of characters allowed in one line of the MPS file.

MSK_IPAR_READ_TASK_IGNORE_PARAM = 134

    Controls what information is used from the task files.

MSK_IPAR_SENSITIVITY_ALL = 135

    Controls sensitivity report behavior.

MSK_IPAR_SENSITIVITY_OPTIMIZER = 136

    Controls which optimizer is used for optimal partition sensitivity analysis.

MSK_IPAR_SENSITIVITY_TYPE = 137

    Controls which type of sensitivity analysis is to be performed.

MSK_IPAR_SIM_BASIS_FACTOR_USE = 138

    Controls whether a (LU) factorization of the basis is used in a hot-start.
    Forcing a refactorization sometimes improves the stability of the simplex optimizers, but in most cases
    there is a performance penalty.

MSK_IPAR_SIM_DEGEN = 139

    Controls how aggressively degeneration is handled.

MSK_IPAR_SIM_DUAL_CRASH = 140

    Controls whether crashing is performed in the dual simplex optimizer.

MSK_IPAR_SIM_DUAL_PHASEONE_METHOD = 141

    An experimental feature.

MSK_IPAR_SIM_DUAL_RESTRICT_SELECTION = 142

    Controls how aggressively restricted selection is used.

MSK_IPAR_SIM_DUAL_SELECTION = 143

    Controls the dual simplex strategy.

MSK_IPAR_SIM_EXPLOIT_DUPVEC = 144

    Controls if the simplex optimizers are allowed to exploit duplicated columns.

MSK_IPAR_SIM_HOTSTART = 145

    Controls the type of hot-start that the simplex optimizer perform.

MSK_IPAR_SIM_HOTSTART_LU = 146

    Determines if the simplex optimizer should exploit the initial factorization.

MSK_IPAR_SIM_INTEGER = 147

    An experimental feature.

MSK_IPAR_SIM_MAX_ITERATIONS = 148

    Maximum number of iterations that can be used by a
    simplex optimizer.

MSK_IPAR_SIM_MAX_NUM_SETBACKS = 149

    Controls how many set-backs that are allowed within a
    simplex optimizer.

MSK_IPAR_SIM_NON_SINGULAR = 150

    Controls if the simplex optimizer ensures a non-singular basis, if possible.

MSK_IPAR_SIM_PRIMAL_CRASH = 151

    Controls the simplex crash.

MSK_IPAR_SIM_PRIMAL_PHASEONE_METHOD = 152

    An experimental feature.

MSK_IPAR_SIM_PRIMAL_RESTRICT_SELECTION = 153

    Controls how aggressively restricted selection is used.

MSK_IPAR_SIM_PRIMAL_SELECTION = 154

    Controls the primal simplex strategy.

MSK_IPAR_SIM_REFACTOR_FREQ = 155

    Controls the basis refactoring frequency.

MSK_IPAR_SIM_REFORMULATION = 156

    Controls if the simplex optimizers are allowed to reformulate the problem.

MSK_IPAR_SIM_SAVE_LU = 157

    Controls if the LU factorization stored should be replaced with the LU factorization
    corresponding to the initial basis.

MSK_IPAR_SIM_SCALING = 158

    Controls how much effort is used in scaling the problem
    before a simplex optimizer is used.

MSK_IPAR_SIM_SCALING_METHOD = 159

    Controls how the problem is scaled
    before a simplex optimizer is used.

MSK_IPAR_SIM_SOLVE_FORM = 160

    Controls whether the primal or the dual problem is solved by the primal-/dual- simplex optimizer.

MSK_IPAR_SIM_STABILITY_PRIORITY = 161

    Controls how high priority the numerical stability should be given.

MSK_IPAR_SIM_SWITCH_OPTIMIZER = 162

    Controls the simplex behavior.

MSK_IPAR_SOL_FILTER_KEEP_BASIC = 163

    Controls the license manager client behavior.

MSK_IPAR_SOL_FILTER_KEEP_RANGED = 164

    Control the contents of the solution files.

MSK_IPAR_SOL_READ_NAME_WIDTH = 165

    Controls the input solution file format.

MSK_IPAR_SOL_READ_WIDTH = 166

    Controls the input solution file format.

MSK_IPAR_SOLUTION_CALLBACK = 167

    Indicates whether solution call-backs will be
    performed during the optimization.

MSK_IPAR_TIMING_LEVEL = 168

    Controls the a amount of timing performed inside MOSEK.

MSK_IPAR_WARNING_LEVEL = 169

    Warning level.

MSK_IPAR_WRITE_BAS_CONSTRAINTS = 170

    Controls the basic solution file format.

MSK_IPAR_WRITE_BAS_HEAD = 171

    Controls the basic solution file format.

MSK_IPAR_WRITE_BAS_VARIABLES = 172

    Controls the basic solution file format.

MSK_IPAR_WRITE_DATA_COMPRESSED = 173

    Controls output file compression.

MSK_IPAR_WRITE_DATA_FORMAT = 174

    Controls the output file format.

MSK_IPAR_WRITE_DATA_PARAM = 175

    Controls output file data.

MSK_IPAR_WRITE_FREE_CON = 176

    Controls the output file data.

MSK_IPAR_WRITE_GENERIC_NAMES = 177

    Controls the output file data.

MSK_IPAR_WRITE_GENERIC_NAMES_IO = 178

    Index origin used in  generic names.

MSK_IPAR_WRITE_IGNORE_INCOMPATIBLE_CONIC_ITEMS = 179

    If the output format is not compatible with conic quadratic problems this parameter controls if the writer ignores the conic parts or produces an error.

MSK_IPAR_WRITE_IGNORE_INCOMPATIBLE_ITEMS = 180

    Controls if the writer ignores incompatible problem items when writing files.

MSK_IPAR_WRITE_IGNORE_INCOMPATIBLE_NL_ITEMS = 181

    Controls if the writer ignores general non-linear terms or produces an error.

MSK_IPAR_WRITE_IGNORE_INCOMPATIBLE_PSD_ITEMS = 182

    If the output format is not compatible with semidefinite problems this parameter controls if the writer ignores the conic parts or produces an error.

MSK_IPAR_WRITE_INT_CONSTRAINTS = 183

    Controls the integer solution file format.

MSK_IPAR_WRITE_INT_HEAD = 184

    Controls the integer solution file format.

MSK_IPAR_WRITE_INT_VARIABLES = 185

    Controls the integer solution file format.

MSK_IPAR_WRITE_LP_LINE_WIDTH = 186

    Controls the LP output file format.

MSK_IPAR_WRITE_LP_QUOTED_NAMES = 187

    Controls LP output file format.

MSK_IPAR_WRITE_LP_STRICT_FORMAT = 188

    Controls whether LP  output files satisfy the LP format strictly.

MSK_IPAR_WRITE_LP_TERMS_PER_LINE = 189

    Controls the LP output file format.

MSK_IPAR_WRITE_MPS_INT = 190

    Controls the output file data.

MSK_IPAR_WRITE_PRECISION = 191

    Controls data precision employed in when writing an MPS file.

MSK_IPAR_WRITE_SOL_BARVARIABLES = 192

    Controls the solution file format.

MSK_IPAR_WRITE_SOL_CONSTRAINTS = 193

    Controls the solution file format.

MSK_IPAR_WRITE_SOL_HEAD = 194

    Controls solution file format.

MSK_IPAR_WRITE_SOL_IGNORE_INVALID_NAMES = 195

    Controls whether the user specified names are employed even if they are invalid names.

MSK_IPAR_WRITE_SOL_VARIABLES = 196

    Controls the solution file format.

MSK_IPAR_WRITE_TASK_INC_SOL = 197

    Controls whether the solutions are  stored in the task file too.

MSK_IPAR_WRITE_XML_MODE = 198

    Controls if linear coefficients should be written by row or column when writing in the XML file format.

Enum solsta
-----------

MSK_SOL_STA_DUAL_FEAS = 3

    The solution is dual feasible.

MSK_SOL_STA_DUAL_INFEAS_CER = 6

    The solution is a certificate of dual infeasibility.

MSK_SOL_STA_INTEGER_OPTIMAL = 14

    The primal solution is integer optimal.

MSK_SOL_STA_NEAR_DUAL_FEAS = 10

    The solution is nearly dual feasible.

MSK_SOL_STA_NEAR_DUAL_INFEAS_CER = 13

    The solution is almost a certificate of dual infeasibility.

MSK_SOL_STA_NEAR_INTEGER_OPTIMAL = 15

    The primal solution is near integer optimal.

MSK_SOL_STA_NEAR_OPTIMAL = 8

    The solution is nearly optimal.

MSK_SOL_STA_NEAR_PRIM_AND_DUAL_FEAS = 11

    The solution is nearly both
    primal and dual feasible.

MSK_SOL_STA_NEAR_PRIM_FEAS = 9

    The solution is nearly primal feasible.

MSK_SOL_STA_NEAR_PRIM_INFEAS_CER = 12

    The solution is almost a certificate
    of primal infeasibility.

MSK_SOL_STA_OPTIMAL = 1

    The solution is optimal.

MSK_SOL_STA_PRIM_AND_DUAL_FEAS = 4

    The solution is both primal and dual feasible.

MSK_SOL_STA_PRIM_FEAS = 2

    The solution is primal feasible.

MSK_SOL_STA_PRIM_INFEAS_CER = 5

    The solution is a certificate
    of primal infeasibility.

MSK_SOL_STA_UNKNOWN = 0

    Status of the solution is unknown.

Enum objsense
-------------

MSK_OBJECTIVE_SENSE_MAXIMIZE = 1

    The problem should be maximized.

MSK_OBJECTIVE_SENSE_MINIMIZE = 0

    The problem should be minimized.

Enum solitem
------------

MSK_SOL_ITEM_SLC = 3

    Lagrange multipliers for lower
    bounds on the constraints.

MSK_SOL_ITEM_SLX = 5

    Lagrange multipliers for lower
    bounds on the variables.

MSK_SOL_ITEM_SNX = 7

    Lagrange multipliers corresponding to the conic constraints on the variables.

MSK_SOL_ITEM_SUC = 4

    Lagrange multipliers for upper
    bounds on the constraints.

MSK_SOL_ITEM_SUX = 6

    Lagrange multipliers for upper
    bounds on the variables.

MSK_SOL_ITEM_XC = 0

    Solution for the constraints.

MSK_SOL_ITEM_XX = 1

    Variable solution.

MSK_SOL_ITEM_Y = 2

    Lagrange multipliers for equations.

Enum boundkey
-------------

MSK_BK_FR = 3

    The constraint or variable is free.

MSK_BK_FX = 2

    The constraint or variable is fixed.

MSK_BK_LO = 0

    The constraint or variable has a finite
    lower bound and an infinite upper bound.

MSK_BK_RA = 4

    The constraint or variable is ranged.

MSK_BK_UP = 1

    The constraint or variable has an infinite
    lower bound and an finite upper bound.

Enum basindtype
---------------

MSK_BI_ALWAYS = 1

    Basis identification is always performed even if the interior-point optimizer terminates
    abnormally.

MSK_BI_IF_FEASIBLE = 3

    Basis identification is not performed if the interior-point optimizer terminates
    with a problem status saying that the problem is primal or dual infeasible.

MSK_BI_NEVER = 0

    Never do basis identification.

MSK_BI_NO_ERROR = 2

    Basis identification is performed if the interior-point optimizer terminates without an error.

MSK_BI_RESERVERED = 4

    Not currently in use.

Enum branchdir
--------------

MSK_BRANCH_DIR_DOWN = 2

    The mixed-integer optimizer always chooses the down branch first.

MSK_BRANCH_DIR_FREE = 0

    The mixed-integer optimizer decides which branch to choose.

MSK_BRANCH_DIR_UP = 1

    The mixed-integer optimizer always chooses the up branch first.

Enum sensitivitytype
--------------------

MSK_SENSITIVITY_TYPE_BASIS = 0

    Basis sensitivity analysis is performed.

MSK_SENSITIVITY_TYPE_OPTIMAL_PARTITION = 1

    Optimal partition sensitivity analysis is performed.

Enum liinfitem
--------------

MSK_LIINF_BI_CLEAN_DUAL_DEG_ITER = 0

    Number of dual degenerate clean iterations performed in the basis identification.

MSK_LIINF_BI_CLEAN_DUAL_ITER = 1

    Number of dual clean iterations performed in the basis identification.

MSK_LIINF_BI_CLEAN_PRIMAL_DEG_ITER = 2

    Number of primal degenerate clean iterations performed in the basis identification.

MSK_LIINF_BI_CLEAN_PRIMAL_DUAL_DEG_ITER = 3

    Number of primal-dual degenerate clean iterations performed in the basis identification.

MSK_LIINF_BI_CLEAN_PRIMAL_DUAL_ITER = 4

    Number of primal-dual clean iterations performed in the basis identification.

MSK_LIINF_BI_CLEAN_PRIMAL_DUAL_SUB_ITER = 5

    Number of primal-dual subproblem clean iterations performed in the basis identification.

MSK_LIINF_BI_CLEAN_PRIMAL_ITER = 6

    Number of primal clean iterations performed in the basis identification.

MSK_LIINF_BI_DUAL_ITER = 7

    Number of dual pivots performed in the basis identification.

MSK_LIINF_BI_PRIMAL_ITER = 8

    Number of primal pivots performed in the basis identification.

MSK_LIINF_INTPNT_FACTOR_NUM_NZ = 9

    Number of non-zeros in factorization.

MSK_LIINF_MIO_INTPNT_ITER = 10

    Number of interior-point iterations performed by the mixed-integer optimizer.

MSK_LIINF_MIO_SIMPLEX_ITER = 11

    Number of simplex iterations performed by the mixed-integer optimizer.

MSK_LIINF_RD_NUMANZ = 12

    Number of non-zeros in A that is read.

MSK_LIINF_RD_NUMQNZ = 13

    Number of Q non-zeros.

Enum streamtype
---------------

MSK_STREAM_ERR = 2

    Error stream. Error messages are written to this stream.

MSK_STREAM_LOG = 0

    Log stream. Contains the aggregated contents of all other streams. This means that a message written to any other stream will also be written to this stream.

MSK_STREAM_MSG = 1

    Message stream. Log information relating to performance and progress of the optimization is written to this stream.

MSK_STREAM_WRN = 3

    Warning stream. Warning messages are written to this stream.

Enum simhotstart
----------------

MSK_SIM_HOTSTART_FREE = 1

    The simplex optimize chooses the hot-start type.

MSK_SIM_HOTSTART_NONE = 0

    The simplex optimizer performs a coldstart.

MSK_SIM_HOTSTART_STATUS_KEYS = 2

    Only the status keys of the constraints and variables are used
    to choose the type of hot-start.

Enum callbackcode
-----------------

MSK_CALLBACK_BEGIN_BI = 0

    The basis identification procedure
    has been started.

MSK_CALLBACK_BEGIN_CONCURRENT = 1

    Concurrent optimizer is started.

MSK_CALLBACK_BEGIN_CONIC = 2

    The call-back function is called
    when the conic optimizer is started.

MSK_CALLBACK_BEGIN_DUAL_BI = 3

    The call-back function is called
    from within the basis identification procedure
    when the dual phase is started.

MSK_CALLBACK_BEGIN_DUAL_SENSITIVITY = 4

    Dual sensitivity analysis is started.

MSK_CALLBACK_BEGIN_DUAL_SETUP_BI = 5

    The call-back function is called when the dual BI phase is started.

MSK_CALLBACK_BEGIN_DUAL_SIMPLEX = 6

    The call-back function is called when the dual simplex optimizer started.

MSK_CALLBACK_BEGIN_DUAL_SIMPLEX_BI = 7

    The call-back function is called
    from within the basis identification procedure
    when the dual simplex clean-up phase is started.

MSK_CALLBACK_BEGIN_FULL_CONVEXITY_CHECK = 8

    Begin full convexity check.

MSK_CALLBACK_BEGIN_INFEAS_ANA = 9

    The call-back function is called when the infeasibility analyzer is started.

MSK_CALLBACK_BEGIN_INTPNT = 10

    The call-back function is called
    when the interior-point optimizer is started.

MSK_CALLBACK_BEGIN_LICENSE_WAIT = 11

    Begin waiting for license.

MSK_CALLBACK_BEGIN_MIO = 12

    The call-back function is called when the mixed-integer optimizer is started.

MSK_CALLBACK_BEGIN_NETWORK_DUAL_SIMPLEX = 13

    The call-back function is called when the dual network simplex optimizer is started.

MSK_CALLBACK_BEGIN_NETWORK_PRIMAL_SIMPLEX = 14

    The call-back function is called when the primal network simplex optimizer is started.

MSK_CALLBACK_BEGIN_NETWORK_SIMPLEX = 15

    The call-back function is called when the simplex network optimizer is started.

MSK_CALLBACK_BEGIN_NONCONVEX = 16

    The call-back function is called
    when the nonconvex optimizer is started.

MSK_CALLBACK_BEGIN_OPTIMIZER = 17

    The call-back function is called when the optimizer is started.

MSK_CALLBACK_BEGIN_PRESOLVE = 18

    The call-back function is called
    when the presolve is started.

MSK_CALLBACK_BEGIN_PRIMAL_BI = 19

    The call-back function is called
    from within the basis identification procedure
    when the primal phase is started.

MSK_CALLBACK_BEGIN_PRIMAL_DUAL_SIMPLEX = 20

    The call-back function is called when the primal-dual simplex optimizer is started.

MSK_CALLBACK_BEGIN_PRIMAL_DUAL_SIMPLEX_BI = 21

    The call-back function is called
    from within the basis identification procedure
    when the primal-dual simplex clean-up phase is started.

MSK_CALLBACK_BEGIN_PRIMAL_REPAIR = 22

    Begin primal feasibility repair.

MSK_CALLBACK_BEGIN_PRIMAL_SENSITIVITY = 23

    Primal sensitivity analysis is started.

MSK_CALLBACK_BEGIN_PRIMAL_SETUP_BI = 24

    The call-back function is called when the primal BI setup is started.

MSK_CALLBACK_BEGIN_PRIMAL_SIMPLEX = 25

    The call-back function is called when the primal simplex optimizer is started.

MSK_CALLBACK_BEGIN_PRIMAL_SIMPLEX_BI = 26

    The call-back function is called
    from within the basis identification procedure
    when the primal simplex clean-up phase is started.

MSK_CALLBACK_BEGIN_QCQO_REFORMULATE = 27

    Begin QCQO reformulation.

MSK_CALLBACK_BEGIN_READ = 28

    MOSEK has started reading a problem file.

MSK_CALLBACK_BEGIN_SIMPLEX = 29

    The call-back function is called when the simplex optimizer is started.

MSK_CALLBACK_BEGIN_SIMPLEX_BI = 30

    The call-back function is called
    from within the basis identification procedure
    when the simplex clean-up phase is started.

MSK_CALLBACK_BEGIN_SIMPLEX_NETWORK_DETECT = 31

    The call-back function is called when the network detection procedure is started.

MSK_CALLBACK_BEGIN_WRITE = 32

    MOSEK has started writing a problem file.

MSK_CALLBACK_CONIC = 33

    The call-back function is called from within the
    conic optimizer after the information database has been updated.

MSK_CALLBACK_DUAL_SIMPLEX = 34

    The call-back function is called
    from within the dual simplex optimizer.

MSK_CALLBACK_END_BI = 35

    The call-back function is called
    when the basis identification procedure
    is terminated.

MSK_CALLBACK_END_CONCURRENT = 36

    Concurrent optimizer is terminated.

MSK_CALLBACK_END_CONIC = 37

    The call-back function is called
    when the conic optimizer is terminated.

MSK_CALLBACK_END_DUAL_BI = 38

    The call-back function is called
    from within the basis identification procedure
    when the dual phase is terminated.

MSK_CALLBACK_END_DUAL_SENSITIVITY = 39

    Dual sensitivity analysis is terminated.

MSK_CALLBACK_END_DUAL_SETUP_BI = 40

    The call-back function is called when the dual BI phase is terminated.

MSK_CALLBACK_END_DUAL_SIMPLEX = 41

    The call-back function is called when the dual simplex optimizer is terminated.

MSK_CALLBACK_END_DUAL_SIMPLEX_BI = 42

    The call-back function is called
    from within the basis identification procedure
    when the dual clean-up phase is terminated.

MSK_CALLBACK_END_FULL_CONVEXITY_CHECK = 43

    End full convexity check.

MSK_CALLBACK_END_INFEAS_ANA = 44

    The call-back function is called when the infeasibility analyzer is terminated.

MSK_CALLBACK_END_INTPNT = 45

    The call-back function is called
    when the interior-point optimizer is terminated.

MSK_CALLBACK_END_LICENSE_WAIT = 46

    End waiting for license.

MSK_CALLBACK_END_MIO = 47

    The call-back function is called when the mixed-integer optimizer is terminated.

MSK_CALLBACK_END_NETWORK_DUAL_SIMPLEX = 48

    The call-back function is called when the dual network simplex optimizer is terminated.

MSK_CALLBACK_END_NETWORK_PRIMAL_SIMPLEX = 49

    The call-back function is called when the primal network simplex optimizer is terminated.

MSK_CALLBACK_END_NETWORK_SIMPLEX = 50

    The call-back function is called when the simplex network optimizer is terminated.

MSK_CALLBACK_END_NONCONVEX = 51

    The call-back function is called
    when the nonconvex optimizer is terminated.

MSK_CALLBACK_END_OPTIMIZER = 52

    The call-back function is called when the optimizer is terminated.

MSK_CALLBACK_END_PRESOLVE = 53

    The call-back function is called
    when the presolve is completed.

MSK_CALLBACK_END_PRIMAL_BI = 54

    The call-back function is called
    from within the basis identification procedure
    when the primal phase is terminated.

MSK_CALLBACK_END_PRIMAL_DUAL_SIMPLEX = 55

    The call-back function is called when the primal-dual simplex optimizer is terminated.

MSK_CALLBACK_END_PRIMAL_DUAL_SIMPLEX_BI = 56

    The call-back function is called
    from within the basis identification procedure
    when the primal-dual clean-up phase is terminated.

MSK_CALLBACK_END_PRIMAL_REPAIR = 57

    End primal feasibility repair.

MSK_CALLBACK_END_PRIMAL_SENSITIVITY = 58

    Primal sensitivity analysis is terminated.

MSK_CALLBACK_END_PRIMAL_SETUP_BI = 59

    The call-back function is called when the primal BI setup is terminated.

MSK_CALLBACK_END_PRIMAL_SIMPLEX = 60

    The call-back function is called when the primal simplex optimizer is terminated.

MSK_CALLBACK_END_PRIMAL_SIMPLEX_BI = 61

    The call-back function is called
    from within the basis identification procedure
    when the primal clean-up phase is terminated.

MSK_CALLBACK_END_QCQO_REFORMULATE = 62

    End QCQO reformulation.

MSK_CALLBACK_END_READ = 63

    MOSEK has finished reading a problem file.

MSK_CALLBACK_END_SIMPLEX = 64

    The call-back function is called when the simplex optimizer is terminated.

MSK_CALLBACK_END_SIMPLEX_BI = 65

    The call-back function is called
    from within the basis identification procedure
    when the simplex clean-up phase is terminated.

MSK_CALLBACK_END_SIMPLEX_NETWORK_DETECT = 66

    The call-back function is called when the network detection procedure is terminated.

MSK_CALLBACK_END_WRITE = 67

    MOSEK has finished writing a problem file.

MSK_CALLBACK_IM_BI = 68

    The call-back function is called
    from within the basis identification procedure
    at an intermediate point.

MSK_CALLBACK_IM_CONIC = 69

    The call-back function is called
    at an intermediate stage within the conic optimizer where
    the information database has not been updated.

MSK_CALLBACK_IM_DUAL_BI = 70

    The call-back function is called
    from within the basis identification procedure
    at an intermediate point in the dual phase.

MSK_CALLBACK_IM_DUAL_SENSIVITY = 71

    The call-back function is called at an intermediate stage of the dual sensitivity analysis.

MSK_CALLBACK_IM_DUAL_SIMPLEX = 72

    The call-back function is called at an intermediate point in the dual simplex optimizer.

MSK_CALLBACK_IM_FULL_CONVEXITY_CHECK = 73

    The call-back function is called at an intermediate stage of the full convexity check.

MSK_CALLBACK_IM_INTPNT = 74

    The call-back function is called
    at an intermediate stage within the interior-point optimizer where
    the information database has not been updated.

MSK_CALLBACK_IM_LICENSE_WAIT = 75

    MOSEK is waiting for a license.

MSK_CALLBACK_IM_LU = 76

    The call-back function is called
    from within the LU factorization procedure at an intermediate point.

MSK_CALLBACK_IM_MIO = 77

    The call-back function is called at an intermediate point in the mixed-integer optimizer.

MSK_CALLBACK_IM_MIO_DUAL_SIMPLEX = 78

    The call-back function is called at an intermediate point in the mixed-integer optimizer while running the
    dual simplex optimizer.

MSK_CALLBACK_IM_MIO_INTPNT = 79

    The call-back function is called at an intermediate point in the mixed-integer optimizer while running the
    interior-point optimizer.

MSK_CALLBACK_IM_MIO_PRESOLVE = 80

    The call-back function is called at an intermediate point in the mixed-integer optimizer while running the
    presolve.

MSK_CALLBACK_IM_MIO_PRIMAL_SIMPLEX = 81

    The call-back function is called at an intermediate point in the mixed-integer optimizer while running the
    primal simplex optimizer.

MSK_CALLBACK_IM_NETWORK_DUAL_SIMPLEX = 82

    The call-back function is called at an intermediate point in the dual network simplex optimizer.

MSK_CALLBACK_IM_NETWORK_PRIMAL_SIMPLEX = 83

    The call-back function is called at an intermediate point in the primal network simplex optimizer.

MSK_CALLBACK_IM_NONCONVEX = 84

    The call-back function is called
    at an intermediate stage within the nonconvex optimizer where
    the information database has not been updated.

MSK_CALLBACK_IM_ORDER = 85

    The call-back function is called
    from within the matrix ordering procedure at an intermediate point.

MSK_CALLBACK_IM_PRESOLVE = 86

    The call-back function is called
    from within the presolve procedure
    at an intermediate stage.

MSK_CALLBACK_IM_PRIMAL_BI = 87

    The call-back function is called
    from within the basis identification procedure
    at an intermediate point in the primal phase.

MSK_CALLBACK_IM_PRIMAL_DUAL_SIMPLEX = 88

    The call-back function is called at an intermediate point in the primal-dual simplex optimizer.

MSK_CALLBACK_IM_PRIMAL_SENSIVITY = 89

    The call-back function is called at an intermediate stage of the primal sensitivity analysis.

MSK_CALLBACK_IM_PRIMAL_SIMPLEX = 90

    The call-back function is called at an intermediate point in the primal simplex optimizer.

MSK_CALLBACK_IM_QO_REFORMULATE = 91

    The call-back function is called at an intermediate stage of the conic quadratic reformulation.

MSK_CALLBACK_IM_READ = 92

    Intermediate stage in reading.

MSK_CALLBACK_IM_SIMPLEX = 93

    The call-back function is called
    from within the simplex optimizer at an intermediate point.

MSK_CALLBACK_IM_SIMPLEX_BI = 94

    The call-back function is called
    from within the basis identification procedure
    at an intermediate point in the simplex clean-up phase.

MSK_CALLBACK_INTPNT = 95

    The call-back function is called from within the
    interior-point optimizer after the information database has been updated.

MSK_CALLBACK_NEW_INT_MIO = 96

    The call-back function is called after a new integer solution
    has been located by the mixed-integer optimizer.

MSK_CALLBACK_NONCOVEX = 97

    The call-back function is called from within the
    nonconvex optimizer after the information database has been updated.

MSK_CALLBACK_PRIMAL_SIMPLEX = 98

    The call-back function is called
    from within the primal simplex optimizer.

MSK_CALLBACK_READ_OPF = 99

    The call-back function is called
    from the OPF reader.

MSK_CALLBACK_READ_OPF_SECTION = 100

    A chunk of Q non-zeros has been read from a problem file.

MSK_CALLBACK_UPDATE_DUAL_BI = 101

    The call-back function is called
    from within the basis identification procedure
    at an intermediate point in the dual phase.

MSK_CALLBACK_UPDATE_DUAL_SIMPLEX = 102

    The call-back function is called in the dual simplex optimizer.

MSK_CALLBACK_UPDATE_DUAL_SIMPLEX_BI = 103

    The call-back function is called
    from within the basis identification procedure
    at an intermediate point in the dual simplex clean-up phase.

MSK_CALLBACK_UPDATE_NETWORK_DUAL_SIMPLEX = 104

    The call-back function is called in the dual network simplex optimizer.

MSK_CALLBACK_UPDATE_NETWORK_PRIMAL_SIMPLEX = 105

    The call-back function is called in the primal network simplex optimizer.

MSK_CALLBACK_UPDATE_NONCONVEX = 106

    The call-back function is called
    at an intermediate stage within the nonconvex optimizer where
    the information database has been updated.

MSK_CALLBACK_UPDATE_PRESOLVE = 107

    The call-back function is called
    from within the presolve procedure.

MSK_CALLBACK_UPDATE_PRIMAL_BI = 108

    The call-back function is called
    from within the basis identification procedure
    at an intermediate point in the primal phase.

MSK_CALLBACK_UPDATE_PRIMAL_DUAL_SIMPLEX = 109

    The call-back function is called  in the primal-dual simplex optimizer.

MSK_CALLBACK_UPDATE_PRIMAL_DUAL_SIMPLEX_BI = 110

    The call-back function is called
    from within the basis identification procedure
    at an intermediate point in the primal simplex clean-up phase.

MSK_CALLBACK_UPDATE_PRIMAL_SIMPLEX = 111

    The call-back function is called  in the primal simplex optimizer.

MSK_CALLBACK_UPDATE_PRIMAL_SIMPLEX_BI = 112

    The call-back function is called
    from within the basis identification procedure
    at an intermediate point in the primal simplex clean-up phase.

MSK_CALLBACK_WRITE_OPF = 113

    The call-back function is called
    from the OPF writer.

Enum symmattype
---------------

MSK_SYMMAT_TYPE_SPARSE = 0

    Sparse symmetric matrix.

Enum feature
------------

MSK_FEATURE_PTOM = 2

    Mixed-integer extension.

MSK_FEATURE_PTON = 1

    Nonlinear extension.

MSK_FEATURE_PTOX = 3

    Non-convex extension.

MSK_FEATURE_PTS = 0

    Base system.

Enum mark
---------

MSK_MARK_LO = 0

    The lower bound is selected for sensitivity analysis.

MSK_MARK_UP = 1

    The upper bound is selected for sensitivity analysis.

Enum conetype
-------------

MSK_CT_QUAD = 0

    The cone is a quadratic cone.

MSK_CT_RQUAD = 1

    The cone is a rotated quadratic cone.

Enum feasrepairtype
-------------------

MSK_FEASREPAIR_OPTIMIZE_COMBINED = 2

    Minimize with original objective subject to minimal weighted violation of bounds.

MSK_FEASREPAIR_OPTIMIZE_NONE = 0

    Do not optimize the feasibility repair problem.

MSK_FEASREPAIR_OPTIMIZE_PENALTY = 1

    Minimize weighted sum of violations.

Enum iomode
-----------

MSK_IOMODE_READ = 0

    The file is read-only.

MSK_IOMODE_READWRITE = 2

    The file is to read and written.

MSK_IOMODE_WRITE = 1

    The file is write-only. If the file exists then it is
    truncated when it is opened. Otherwise it is created when it is opened.

Enum simseltype
---------------

MSK_SIM_SELECTION_ASE = 2

    The optimizer uses approximate steepest-edge
    pricing.

MSK_SIM_SELECTION_DEVEX = 3

    The optimizer uses devex steepest-edge pricing (or if it is not available an
    approximate steep-edge selection).

MSK_SIM_SELECTION_FREE = 0

    The optimizer chooses the pricing strategy.

MSK_SIM_SELECTION_FULL = 1

    The optimizer uses full pricing.

MSK_SIM_SELECTION_PARTIAL = 5

    The optimizer uses a partial selection approach. The approach is usually
    beneficial if the number of variables is much larger than  the number of constraints.

MSK_SIM_SELECTION_SE = 4

    The optimizer uses steepest-edge selection (or if it is not available an
    approximate steep-edge selection).

Enum msgkey
-----------

MSK_MSG_MPS_SELECTED = 1100

    

MSK_MSG_READING_FILE = 1000

    

MSK_MSG_WRITING_FILE = 1001

    

Enum miomode
------------

MSK_MIO_MODE_IGNORED = 0

    The integer constraints are ignored and the problem is solved as a continuous problem.

MSK_MIO_MODE_LAZY = 2

    Integer restrictions should be satisfied if
    an optimizer is available for the problem.

MSK_MIO_MODE_SATISFIED = 1

    Integer restrictions should be satisfied.

Enum dinfitem
-------------

MSK_DINF_BI_CLEAN_DUAL_TIME = 0

    Time  spent within the dual clean-up optimizer of the basis identification
    procedure since its invocation.

MSK_DINF_BI_CLEAN_PRIMAL_DUAL_TIME = 1

    Time spent within the primal-dual clean-up optimizer of the basis identification
    procedure since its invocation.

MSK_DINF_BI_CLEAN_PRIMAL_TIME = 2

    Time spent within the primal clean-up optimizer of the basis identification
    procedure since its invocation.

MSK_DINF_BI_CLEAN_TIME = 3

    Time spent within the clean-up phase of the basis identification
    procedure since its invocation.

MSK_DINF_BI_DUAL_TIME = 4

    Time spent within the dual phase basis identification
    procedure since its invocation.

MSK_DINF_BI_PRIMAL_TIME = 5

    Time  spent within the primal phase of the basis identification
    procedure since its invocation.

MSK_DINF_BI_TIME = 6

    Time spent within the basis identification
    procedure since its invocation.

MSK_DINF_CONCURRENT_TIME = 7

    Time spent within the concurrent optimizer since its invocation.

MSK_DINF_INTPNT_DUAL_FEAS = 8

    Dual feasibility measure reported by the
    interior-point optimizer. (For the
    interior-point optimizer this measure does not
    directly related to the original problem because
    a homogeneous model is employed.)

MSK_DINF_INTPNT_DUAL_OBJ = 9

    Dual objective value reported by the
    interior-point optimizer.

MSK_DINF_INTPNT_FACTOR_NUM_FLOPS = 10

    An estimate of the number of flops used in the factorization.

MSK_DINF_INTPNT_OPT_STATUS = 11

    This measure should converge to +1 if the problem
    has a primal-dual optimal solution, and converge to -1
    if problem is (strictly) primal or dual infeasible. Furthermore, if the measure converges to 0
    the problem is usually ill-posed.

MSK_DINF_INTPNT_ORDER_TIME = 12

    Order time (in seconds).

MSK_DINF_INTPNT_PRIMAL_FEAS = 13

    Primal feasibility measure reported by the
    interior-point optimizers. (For the
    interior-point optimizer this measure does not
    directly related to the original problem because
    a homogeneous model is employed).

MSK_DINF_INTPNT_PRIMAL_OBJ = 14

    Primal objective value reported by the interior-point optimizer.

MSK_DINF_INTPNT_TIME = 15

    Time spent within the interior-point optimizer
    since its invocation.

MSK_DINF_MIO_CONSTRUCT_SOLUTION_OBJ = 16

    If MOSEK has successfully constructed an integer feasible solution, then this item
    contains the optimal objective value corresponding to the feasible solution.

MSK_DINF_MIO_HEURISTIC_TIME = 17

    Time spent in the optimizer while solving the relaxtions.

MSK_DINF_MIO_OBJ_ABS_GAP = 18

    If the mixed-integer optimizer has computed a feasible solution and a bound, this contains the absolute gap.

MSK_DINF_MIO_OBJ_BOUND = 19

    The best bound on the objective value known.

MSK_DINF_MIO_OBJ_INT = 20

    The primal objective value corresponding to the best integer feasible
    solution.

MSK_DINF_MIO_OBJ_REL_GAP = 21

    If the mixed-integer optimizer has computed a feasible solution and a bound, this contains the relative gap.

MSK_DINF_MIO_OPTIMIZER_TIME = 22

    Time spent in the optimizer while solving the relaxtions.

MSK_DINF_MIO_ROOT_OPTIMIZER_TIME = 23

    Time spent in the optimizer while solving the root relaxation.

MSK_DINF_MIO_ROOT_PRESOLVE_TIME = 24

    Time spent in while presolving the root relaxation.

MSK_DINF_MIO_TIME = 25

    Time spent in the mixed-integer optimizer.

MSK_DINF_MIO_USER_OBJ_CUT = 26

    If the objective cut is used, then this information item has the value of the cut.

MSK_DINF_OPTIMIZER_TIME = 27

    Total time spent in the optimizer since it was invoked.

MSK_DINF_PRESOLVE_ELI_TIME = 28

    Total time spent in the eliminator
    since the presolve was invoked.

MSK_DINF_PRESOLVE_LINDEP_TIME = 29

    Total time spent  in the linear dependency checker
    since the presolve was invoked.

MSK_DINF_PRESOLVE_TIME = 30

    Total time (in seconds) spent in the presolve
    since it was invoked.

MSK_DINF_PRIMAL_REPAIR_PENALTY_OBJ = 31

    The optimal objective value of the penalty function.

MSK_DINF_QCQO_REFORMULATE_TIME = 32

    Time spent with conic quadratic reformulation.

MSK_DINF_RD_TIME = 33

    Time spent reading the data file.

MSK_DINF_SIM_DUAL_TIME = 34

    Time spent in the dual simplex
    optimizer since invoking it.

MSK_DINF_SIM_FEAS = 35

    Feasibility measure reported by the
    simplex optimizer.

MSK_DINF_SIM_NETWORK_DUAL_TIME = 36

    Time spent in the dual network simplex
    optimizer since invoking it.

MSK_DINF_SIM_NETWORK_PRIMAL_TIME = 37

    Time spent in the primal network simplex
    optimizer since invoking it.

MSK_DINF_SIM_NETWORK_TIME = 38

    Time spent in the network simplex
    optimizer since invoking it.

MSK_DINF_SIM_OBJ = 39

    Objective value reported by the
    simplex optimizer.

MSK_DINF_SIM_PRIMAL_DUAL_TIME = 40

    Time spent in the primal-dual simplex optimizer
    since invoking it.

MSK_DINF_SIM_PRIMAL_TIME = 41

    Time spent in the primal simplex
    optimizer since invoking it.

MSK_DINF_SIM_TIME = 42

    Time spent in the simplex
    optimizer since invoking it.

MSK_DINF_SOL_BAS_DUAL_OBJ = 43

    Dual objective value of the basic solution. Updated by the function updatesolutioninfo.

MSK_DINF_SOL_BAS_DVIOLCON = 44

    Maximal dual bound violation for xx in the basic solution. Updated by the function updatesolutioninfo.

MSK_DINF_SOL_BAS_DVIOLVAR = 45

    Maximal dual bound violation for xx in the basic solution. Updated by the function updatesolutioninfo.

MSK_DINF_SOL_BAS_PRIMAL_OBJ = 46

    Primal objective value of the basic solution. Updated by the function updatesolutioninfo.

MSK_DINF_SOL_BAS_PVIOLCON = 47

    Maximal primal bound violation for xx in the basic solution. Updated by the function updatesolutioninfo.

MSK_DINF_SOL_BAS_PVIOLVAR = 48

    Maximal primal bound violation for xx in the basic solution. Updated by the function updatesolutioninfo.

MSK_DINF_SOL_ITG_PRIMAL_OBJ = 49

    Primal objective value of the integer solution. Updated by the function updatesolutioninfo.

MSK_DINF_SOL_ITG_PVIOLBARVAR = 50

    Maximal primal bound violation for barx in the integer solution. Updated by the function updatesolutioninfo.

MSK_DINF_SOL_ITG_PVIOLCON = 51

    Maximal primal bound violation for xx in the integer solution. Updated by the function updatesolutioninfo.

MSK_DINF_SOL_ITG_PVIOLCONES = 52

    Maximal primal violation for primal conic constraints in the integer solution. Updated by the function updatesolutioninfo.

MSK_DINF_SOL_ITG_PVIOLITG = 53

    Maximal violation for the integer constraints in the integer solution. Updated by the function updatesolutioninfo.

MSK_DINF_SOL_ITG_PVIOLVAR = 54

    Maximal primal bound violation for xx in the integer solution. Updated by the function updatesolutioninfo.

MSK_DINF_SOL_ITR_DUAL_OBJ = 55

    Dual objective value of the interior-point solution. Updated by the function updatesolutioninfo.

MSK_DINF_SOL_ITR_DVIOLBARVAR = 56

    Maximal dual bound violation for barx in the interior-point solution. Updated by the function updatesolutioninfo.

MSK_DINF_SOL_ITR_DVIOLCON = 57

    Maximal dual bound violation for xx in the interior-point solution. Updated by the function updatesolutioninfo.

MSK_DINF_SOL_ITR_DVIOLCONES = 58

    Maximal dual violation for dual conic constraints in the interior-point solution. Updated by the function updatesolutioninfo.

MSK_DINF_SOL_ITR_DVIOLVAR = 59

    Maximal dual bound violation for xx in the interior-point solution. Updated by the function updatesolutioninfo.

MSK_DINF_SOL_ITR_PRIMAL_OBJ = 60

    Primal objective value of the interior-point solution. Updated by the function updatesolutioninfo.

MSK_DINF_SOL_ITR_PVIOLBARVAR = 61

    Maximal primal bound violation for barx in the interior-point solution. Updated by the function updatesolutioninfo.

MSK_DINF_SOL_ITR_PVIOLCON = 62

    Maximal primal bound violation for xx in the interior-point solution. Updated by the function updatesolutioninfo.

MSK_DINF_SOL_ITR_PVIOLCONES = 63

    Maximal primal violation for primal conic constraints in the interior-point solution. Updated by the function updatesolutioninfo.

MSK_DINF_SOL_ITR_PVIOLVAR = 64

    Maximal primal bound violation for xx in the interior-point solution. Updated by the function updatesolutioninfo.

Enum parametertype
------------------

MSK_PAR_DOU_TYPE = 1

    Is a double parameter.

MSK_PAR_INT_TYPE = 2

    Is an integer parameter.

MSK_PAR_INVALID_TYPE = 0

    Not a valid parameter.

MSK_PAR_STR_TYPE = 3

    Is a string parameter.

Enum rescodetype
----------------

MSK_RESPONSE_ERR = 3

    The response code is an error.

MSK_RESPONSE_OK = 0

    The response code is OK.

MSK_RESPONSE_TRM = 2

    The response code is an optimizer termination status.

MSK_RESPONSE_UNK = 4

    The response code does not belong to any class.

MSK_RESPONSE_WRN = 1

    The response code is a warning.

Enum prosta
-----------

MSK_PRO_STA_DUAL_FEAS = 3

    The problem is dual feasible.

MSK_PRO_STA_DUAL_INFEAS = 5

    The problem is dual infeasible.

MSK_PRO_STA_ILL_POSED = 7

    The problem is ill-posed. For example,
    it may be primal and dual feasible but
    have a positive duality gap.

MSK_PRO_STA_NEAR_DUAL_FEAS = 10

    The problem is at least nearly dual feasible.

MSK_PRO_STA_NEAR_PRIM_AND_DUAL_FEAS = 8

    The problem is at least nearly primal and dual feasible.

MSK_PRO_STA_NEAR_PRIM_FEAS = 9

    The problem is at least nearly primal feasible.

MSK_PRO_STA_PRIM_AND_DUAL_FEAS = 1

    The problem is primal and dual feasible.

MSK_PRO_STA_PRIM_AND_DUAL_INFEAS = 6

    The problem is primal and dual infeasible.

MSK_PRO_STA_PRIM_FEAS = 2

    The problem is primal feasible.

MSK_PRO_STA_PRIM_INFEAS = 4

    The problem is primal infeasible.

MSK_PRO_STA_PRIM_INFEAS_OR_UNBOUNDED = 11

    The problem is either primal infeasible or unbounded. This may occur for
    mixed-integer problems.

MSK_PRO_STA_UNKNOWN = 0

    Unknown problem status.

Enum scalingtype
----------------

MSK_SCALING_AGGRESSIVE = 3

    A very aggressive scaling is performed.

MSK_SCALING_FREE = 0

    The optimizer chooses the scaling heuristic.

MSK_SCALING_MODERATE = 2

    A conservative scaling is performed.

MSK_SCALING_NONE = 1

    No scaling is performed.

Enum rescode
------------

MSK_RES_ERR_AD_INVALID_CODELIST = 3102

    The code list data was invalid.

MSK_RES_ERR_AD_INVALID_OPERAND = 3104

    The code list data was invalid.

MSK_RES_ERR_AD_INVALID_OPERATOR = 3103

    The code list data was invalid.

MSK_RES_ERR_AD_MISSING_OPERAND = 3105

    The code list data was invalid.

MSK_RES_ERR_AD_MISSING_RETURN = 3106

    The code list data was invalid.

MSK_RES_ERR_API_ARRAY_TOO_SMALL = 3001

    An input array was too short.

MSK_RES_ERR_API_CB_CONNECT = 3002

    Failed to connect a callback object.

MSK_RES_ERR_API_FATAL_ERROR = 3005

    An internal error occurred in the API. Please report this problem.

MSK_RES_ERR_API_INTERNAL = 3999

    An internal fatal error occurred in an interface function.:w

MSK_RES_ERR_ARG_IS_TOO_LARGE = 1227

    The value of a argument is too small.

MSK_RES_ERR_ARG_IS_TOO_SMALL = 1226

    The value of a argument is too small.

MSK_RES_ERR_ARGUMENT_DIMENSION = 1201

    A function argument is of incorrect dimension.

MSK_RES_ERR_ARGUMENT_IS_TOO_LARGE = 5005

    The value of a function argument is too large.

MSK_RES_ERR_ARGUMENT_LENNEQ = 1197

    Incorrect length of arguments.

MSK_RES_ERR_ARGUMENT_PERM_ARRAY = 1299

    An invalid permutation array is specified.

MSK_RES_ERR_ARGUMENT_TYPE = 1198

    Incorrect argument type.

MSK_RES_ERR_BAR_VAR_DIM = 3920

    The dimension of a symmetric matrix variable has to greater than 0.

MSK_RES_ERR_BASIS = 1266

    Invalid basis is specified.

MSK_RES_ERR_BASIS_FACTOR = 1610

    The factorization of the basis is invalid.

MSK_RES_ERR_BASIS_SINGULAR = 1615

    The basis is singular.

MSK_RES_ERR_BLANK_NAME = 1070

    An all blank name has been specified.

MSK_RES_ERR_BLAS_ARG_K = 7004

    Invalid argument k.

MSK_RES_ERR_BLAS_ARG_M = 7002

    Invalid argument m.

MSK_RES_ERR_BLAS_ARG_N = 7003

    Invalid argument n.

MSK_RES_ERR_BLAS_ARG_TRANSA = 7005

    Invalid argument transa.

MSK_RES_ERR_BLAS_ARG_TRANSB = 7006

    Invalid argument transb.

MSK_RES_ERR_BLAS_SINGULAR_MATRIX = 7000

    A matrix is singular.

MSK_RES_ERR_BLAS_UNKNOWN = 7001

    An unknown error.

MSK_RES_ERR_CANNOT_CLONE_NL = 2505

    A task with a nonlinear function call-back cannot be cloned.

MSK_RES_ERR_CANNOT_HANDLE_NL = 2506

    A function cannot handle a task with nonlinear function call-backs.

MSK_RES_ERR_CBF_DUP_ACOORD = 7116

    Duplicate index in ACOORD.

MSK_RES_ERR_CBF_DUP_BCOORD = 7115

    Duplicate index in BCOORD.

MSK_RES_ERR_CBF_DUP_OBJACOORD = 7114

    Duplicate index in OBJCOORD.

MSK_RES_ERR_CBF_DUPLICATE_CON = 7108

    Duplicate CON keyword.

MSK_RES_ERR_CBF_DUPLICATE_INT = 7110

    Duplicate INT keyword.

MSK_RES_ERR_CBF_DUPLICATE_OBJ = 7107

    Duplicate OBJ keyword.

MSK_RES_ERR_CBF_DUPLICATE_VAR = 7109

    Duplicate VAR keyword.

MSK_RES_ERR_CBF_INVALID_CON_TYPE = 7112

    Invalid constraint type.

MSK_RES_ERR_CBF_INVALID_DOMAIN_DIMENSION = 7113

    Invalid domain dimension.

MSK_RES_ERR_CBF_INVALID_INT_INDEX = 7121

    Invalid INT index.

MSK_RES_ERR_CBF_INVALID_VAR_TYPE = 7111

    Invalid variable type.

MSK_RES_ERR_CBF_NO_VARIABLES = 7102

    An invalid objective sense is specified.

MSK_RES_ERR_CBF_NO_VERSION_SPECIFIED = 7105

    No version specified.

MSK_RES_ERR_CBF_OBJ_SENSE = 7101

    An invalid objective sense is specified.

MSK_RES_ERR_CBF_PARSE = 7100

    An error occurred while parsing an CBF file.

MSK_RES_ERR_CBF_SYNTAX = 7106

    Invalid syntax.

MSK_RES_ERR_CBF_TOO_FEW_CONSTRAINTS = 7118

    Too few constraints defined.

MSK_RES_ERR_CBF_TOO_FEW_INTS = 7119

    Too ints specified.

MSK_RES_ERR_CBF_TOO_FEW_VARIABLES = 7117

    Too few variables defined.

MSK_RES_ERR_CBF_TOO_MANY_CONSTRAINTS = 7103

    Too many constraints specified.

MSK_RES_ERR_CBF_TOO_MANY_INTS = 7120

    Too ints specified.

MSK_RES_ERR_CBF_TOO_MANY_VARIABLES = 7104

    Too many variables specified.

MSK_RES_ERR_CBF_UNSUPPORTED = 7122

    Unsupported feature is present.

MSK_RES_ERR_CON_Q_NOT_NSD = 1294

    The quadratic constraint matrix is not NSD.

MSK_RES_ERR_CON_Q_NOT_PSD = 1293

    The quadratic constraint matrix is not PSD.

MSK_RES_ERR_CONCURRENT_OPTIMIZER = 3059

    An unsupported optimizer was chosen for use with the concurrent optimizer.

MSK_RES_ERR_CONE_INDEX = 1300

    An index of a non-existing cone has been specified.

MSK_RES_ERR_CONE_OVERLAP = 1302

    A new cone which variables overlap with an existing cone has been specified.

MSK_RES_ERR_CONE_OVERLAP_APPEND = 1307

    The cone to be appended has one variable which is already member of another cone.

MSK_RES_ERR_CONE_REP_VAR = 1303

    A variable is included multiple times in the cone.

MSK_RES_ERR_CONE_SIZE = 1301

    A cone with too few members is specified.

MSK_RES_ERR_CONE_TYPE = 1305

    Invalid cone type specified.

MSK_RES_ERR_CONE_TYPE_STR = 1306

    Invalid cone type specified.

MSK_RES_ERR_DATA_FILE_EXT = 1055

    The data file format cannot be determined from the file name.

MSK_RES_ERR_DUP_NAME = 1071

    Duplicate names specified.

MSK_RES_ERR_DUPLICATE_BARVARIABLE_NAMES = 4502

    Two barvariable names are identical.

MSK_RES_ERR_DUPLICATE_CONE_NAMES = 4503

    Two cone names are identical.

MSK_RES_ERR_DUPLICATE_CONSTRAINT_NAMES = 4500

    Two constraint names are identical.

MSK_RES_ERR_DUPLICATE_VARIABLE_NAMES = 4501

    Two variable names are identical.

MSK_RES_ERR_END_OF_FILE = 1059

    End of file reached.

MSK_RES_ERR_FACTOR = 1650

    An error occurred while factorizing a matrix.

MSK_RES_ERR_FEASREPAIR_CANNOT_RELAX = 1700

    An optimization problem cannot be relaxed.

MSK_RES_ERR_FEASREPAIR_INCONSISTENT_BOUND = 1702

    The upper bound is less than the lower bound for a variable or a constraint.

MSK_RES_ERR_FEASREPAIR_SOLVING_RELAXED = 1701

    The relaxed problem could not be solved to optimality.

MSK_RES_ERR_FILE_LICENSE = 1007

    Invalid license file.

MSK_RES_ERR_FILE_OPEN = 1052

    An error occurred while opening a file.

MSK_RES_ERR_FILE_READ = 1053

    An error occurred while reading file.

MSK_RES_ERR_FILE_WRITE = 1054

    An error occurred while writing to a file.

MSK_RES_ERR_FIRST = 1261

    Invalid first.

MSK_RES_ERR_FIRSTI = 1285

    Invalid firsti.

MSK_RES_ERR_FIRSTJ = 1287

    Invalid firstj.

MSK_RES_ERR_FIXED_BOUND_VALUES = 1425

    A fixed constraint/variable has been specified using the bound keys but the numerical bounds are different.

MSK_RES_ERR_FLEXLM = 1014

    The FLEXlm license manager reported an error.

MSK_RES_ERR_GLOBAL_INV_CONIC_PROBLEM = 1503

    The global optimizer can only be applied to problems without semidefinite variables.

MSK_RES_ERR_HUGE_AIJ = 1380

    A numerically huge value is specified for an element in A.

MSK_RES_ERR_HUGE_C = 1375

    A huge value in absolute size is specified for one an objective coefficient.

MSK_RES_ERR_IDENTICAL_TASKS = 3101

    Some tasks related to this function call were identical. Unique tasks were expected.

MSK_RES_ERR_IN_ARGUMENT = 1200

    A function argument is incorrect.

MSK_RES_ERR_INDEX = 1235

    An index is out of range.

MSK_RES_ERR_INDEX_ARR_IS_TOO_LARGE = 1222

    An index in an array argument is too large.

MSK_RES_ERR_INDEX_ARR_IS_TOO_SMALL = 1221

    An index in an array argument is too small.

MSK_RES_ERR_INDEX_IS_TOO_LARGE = 1204

    An index in an argument is too large.

MSK_RES_ERR_INDEX_IS_TOO_SMALL = 1203

    An index in an argument is too small.

MSK_RES_ERR_INF_DOU_INDEX = 1219

    A double information index is out of range for the specified type.

MSK_RES_ERR_INF_DOU_NAME = 1230

    A double information name is invalid.

MSK_RES_ERR_INF_INT_INDEX = 1220

    An integer information index is out of range for the specified type.

MSK_RES_ERR_INF_INT_NAME = 1231

    An integer information name is invalid.

MSK_RES_ERR_INF_LINT_INDEX = 1225

    A long integer information index is out of range for the specified type.

MSK_RES_ERR_INF_LINT_NAME = 1234

    A long integer information name is invalid.

MSK_RES_ERR_INF_TYPE = 1232

    The information type is invalid.

MSK_RES_ERR_INFEAS_UNDEFINED = 3910

    The requested value is not defined for this solution type.

MSK_RES_ERR_INFINITE_BOUND = 1400

    A numerically huge bound value is specified.

MSK_RES_ERR_INT64_TO_INT32_CAST = 3800

    An 32 bit integer could not cast to a 64 bit integer.

MSK_RES_ERR_INTERNAL = 3000

    An internal error occurred.

MSK_RES_ERR_INTERNAL_TEST_FAILED = 3500

    An internal unit test function failed.

MSK_RES_ERR_INV_APTRE = 1253

    aptre[j] is strictly smaller than aptrb[j] for some j.

MSK_RES_ERR_INV_BK = 1255

    Invalid bound key.

MSK_RES_ERR_INV_BKC = 1256

    Invalid bound key is specified for a constraint.

MSK_RES_ERR_INV_BKX = 1257

    An invalid bound key is specified for a variable.

MSK_RES_ERR_INV_CONE_TYPE = 1272

    Invalid cone type code encountered.

MSK_RES_ERR_INV_CONE_TYPE_STR = 1271

    Invalid cone type string encountered.

MSK_RES_ERR_INV_CONIC_PROBLEM = 1502

    The conic optimizer can only be applied to problems with linear objective and constraints.

MSK_RES_ERR_INV_MARKI = 2501

    Invalid value in marki.

MSK_RES_ERR_INV_MARKJ = 2502

    Invalid value in markj.

MSK_RES_ERR_INV_NAME_ITEM = 1280

    An invalid name item code is used.

MSK_RES_ERR_INV_NUMI = 2503

    Invalid numi.

MSK_RES_ERR_INV_NUMJ = 2504

    Invalid numj.

MSK_RES_ERR_INV_OPTIMIZER = 1550

    An invalid optimizer has been chosen for the problem.

MSK_RES_ERR_INV_PROBLEM = 1500

    Invalid problem type.

MSK_RES_ERR_INV_QCON_SUBI = 1405

    Invalid value in qcsubi.

MSK_RES_ERR_INV_QCON_SUBJ = 1406

    Invalid value in qcsubj.

MSK_RES_ERR_INV_QCON_SUBK = 1404

    Invalid value in qcsubk.

MSK_RES_ERR_INV_QCON_VAL = 1407

    Invalid value in qcval.

MSK_RES_ERR_INV_QOBJ_SUBI = 1401

    Invalid value %d at qosubi.

MSK_RES_ERR_INV_QOBJ_SUBJ = 1402

    Invalid value in qosubj.

MSK_RES_ERR_INV_QOBJ_VAL = 1403

    Invalid value in qoval.

MSK_RES_ERR_INV_SK = 1270

    Invalid status key code encountered.

MSK_RES_ERR_INV_SK_STR = 1269

    Invalid status key string encountered.

MSK_RES_ERR_INV_SKC = 1267

    Invalid value in skc encountered.

MSK_RES_ERR_INV_SKN = 1274

    Invalid value in skn encountered.

MSK_RES_ERR_INV_SKX = 1268

    Invalid value in skx encountered.

MSK_RES_ERR_INV_VAR_TYPE = 1258

    An invalid variable type is specified for a variable.

MSK_RES_ERR_INVALID_ACCMODE = 2520

    An invalid access mode is specified.

MSK_RES_ERR_INVALID_AMPL_STUB = 3700

    Invalid AMPL stub.

MSK_RES_ERR_INVALID_BARVAR_NAME = 1079

    An invalid symmetric matrix variable name is used.

MSK_RES_ERR_INVALID_BRANCH_DIRECTION = 3200

    An invalid branching direction is specified.

MSK_RES_ERR_INVALID_BRANCH_PRIORITY = 3201

    An invalid branching priority is specified.

MSK_RES_ERR_INVALID_COMPRESSION = 1800

    Invalid compression type.

MSK_RES_ERR_INVALID_CON_NAME = 1076

    An invalid constraint name is used.

MSK_RES_ERR_INVALID_CONE_NAME = 1078

    An invalid cone name is used.

MSK_RES_ERR_INVALID_FILE_FORMAT_FOR_CONES = 4005

    The file format does not support a problem with conic constraints.

MSK_RES_ERR_INVALID_FILE_FORMAT_FOR_GENERAL_NL = 4010

    The file format does not support a problem with general nonlinear terms.

MSK_RES_ERR_INVALID_FILE_FORMAT_FOR_SYM_MAT = 4000

    The file format does not support a problem with symmetric matrix variables.

MSK_RES_ERR_INVALID_FILE_NAME = 1056

    An invalid file name has been specified.

MSK_RES_ERR_INVALID_FORMAT_TYPE = 1283

    Invalid format type.

MSK_RES_ERR_INVALID_IDX = 1246

    A specified index is invalid.

MSK_RES_ERR_INVALID_IOMODE = 1801

    Invalid io mode.

MSK_RES_ERR_INVALID_MAX_NUM = 1247

    A specified index is invalid.

MSK_RES_ERR_INVALID_NAME_IN_SOL_FILE = 1170

    An invalid name occurred in a solution file.

MSK_RES_ERR_INVALID_NETWORK_PROBLEM = 1504

    The problem is not a network problem as expected.

MSK_RES_ERR_INVALID_OBJ_NAME = 1075

    An invalid objective name is specified.

MSK_RES_ERR_INVALID_OBJECTIVE_SENSE = 1445

    An invalid objective sense is specified.

MSK_RES_ERR_INVALID_PROBLEM_TYPE = 6000

    An invalid problem type.

MSK_RES_ERR_INVALID_SOL_FILE_NAME = 1057

    An invalid file name has been specified.

MSK_RES_ERR_INVALID_STREAM = 1062

    An invalid stream is referenced.

MSK_RES_ERR_INVALID_SURPLUS = 1275

    Invalid surplus.

MSK_RES_ERR_INVALID_SYM_MAT_DIM = 3950

    A sparse symmetric matrix of invalid dimension is specified.

MSK_RES_ERR_INVALID_TASK = 1064

    The task is invalid.

MSK_RES_ERR_INVALID_UTF8 = 2900

    An invalid UTF8 string is encountered.

MSK_RES_ERR_INVALID_VAR_NAME = 1077

    An invalid variable name is used.

MSK_RES_ERR_INVALID_WCHAR = 2901

    An invalid wchar string is encountered.

MSK_RES_ERR_INVALID_WHICHSOL = 1228

    whichsol is invalid.

MSK_RES_ERR_LAST = 1262

    Invalid last.

MSK_RES_ERR_LASTI = 1286

    Invalid lasti.

MSK_RES_ERR_LASTJ = 1288

    Invalid lastj.

MSK_RES_ERR_LICENSE = 1000

    Invalid license.

MSK_RES_ERR_LICENSE_CANNOT_ALLOCATE = 1020

    The license system cannot allocate the memory required.

MSK_RES_ERR_LICENSE_CANNOT_CONNECT = 1021

    MOSEK cannot connect to the license server.

MSK_RES_ERR_LICENSE_EXPIRED = 1001

    The license has expired.

MSK_RES_ERR_LICENSE_FEATURE = 1018

    A requested feature is not available in the license file(s).

MSK_RES_ERR_LICENSE_INVALID_HOSTID = 1025

    The host ID specified in the license file does not match the host ID of the computer.

MSK_RES_ERR_LICENSE_MAX = 1016

    Maximum number of licenses is reached.

MSK_RES_ERR_LICENSE_MOSEKLM_DAEMON = 1017

    The MOSEKLM license manager daemon is not up and running.

MSK_RES_ERR_LICENSE_NO_SERVER_LINE = 1028

    No SERVER lines in license file.

MSK_RES_ERR_LICENSE_NO_SERVER_SUPPORT = 1027

    The license server does not support the requested feature.

MSK_RES_ERR_LICENSE_SERVER = 1015

    The license server is not responding.

MSK_RES_ERR_LICENSE_SERVER_VERSION = 1026

    The version specified in the checkout request is greater than the highest version number the daemon supports.

MSK_RES_ERR_LICENSE_VERSION = 1002

    Invalid license version.

MSK_RES_ERR_LINK_FILE_DLL = 1040

    A file cannot be linked to a stream in the DLL version.

MSK_RES_ERR_LIVING_TASKS = 1066

    Not all tasks associated with the environment have been deleted.

MSK_RES_ERR_LOWER_BOUND_IS_A_NAN = 1390

    The lower bound specified is not a number (nan).

MSK_RES_ERR_LP_DUP_SLACK_NAME = 1152

    The name of the slack variable added to a ranged constraint already exists.

MSK_RES_ERR_LP_EMPTY = 1151

    The problem cannot be written to an LP formatted file.

MSK_RES_ERR_LP_FILE_FORMAT = 1157

    Syntax error in an LP file.

MSK_RES_ERR_LP_FORMAT = 1160

    Syntax error in an LP file.

MSK_RES_ERR_LP_FREE_CONSTRAINT = 1155

    Free constraints cannot be written in LP file format.

MSK_RES_ERR_LP_INCOMPATIBLE = 1150

    The problem cannot be written to an LP formatted file.

MSK_RES_ERR_LP_INVALID_CON_NAME = 1171

    A constraint name is invalid when used in an LP formatted file.

MSK_RES_ERR_LP_INVALID_VAR_NAME = 1154

    A variable name is invalid when used in an LP formatted file.

MSK_RES_ERR_LP_WRITE_CONIC_PROBLEM = 1163

    The problem contains cones that cannot be written to an LP formatted file.

MSK_RES_ERR_LP_WRITE_GECO_PROBLEM = 1164

    The problem contains general convex terms that cannot be written to an LP formatted file.

MSK_RES_ERR_LU_MAX_NUM_TRIES = 2800

    Could not compute the LU factors of the matrix within the maximum number of allowed tries.

MSK_RES_ERR_MAX_LEN_IS_TOO_SMALL = 1289

    An maximum length that is too small has been specified.

MSK_RES_ERR_MAXNUMBARVAR = 1242

    The maximum number of semidefinite variables limit is too small.

MSK_RES_ERR_MAXNUMCON = 1240

    Invalid maximum number of constraints specified.

MSK_RES_ERR_MAXNUMCONE = 1304

    The value specified for maxnumcone is too small.

MSK_RES_ERR_MAXNUMQNZ = 1243

    Too small maximum number of non-zeros for the Q matrices is specified.

MSK_RES_ERR_MAXNUMVAR = 1241

    The maximum number of variables limit is too small.

MSK_RES_ERR_MIO_INTERNAL = 5010

    A fatal error occurred in the mixed integer optimizer.  Please contact MOSEK support.

MSK_RES_ERR_MIO_NO_OPTIMIZER = 1551

    No optimizer is available for the current class of integer optimization problems.

MSK_RES_ERR_MIO_NOT_LOADED = 1553

    The mixed-integer optimizer is not loaded.

MSK_RES_ERR_MISSING_LICENSE_FILE = 1008

    MOSEK cannot locate a license.

MSK_RES_ERR_MIXED_PROBLEM = 1501

    The problem contains both conic and nonlinear constraints.

MSK_RES_ERR_MPS_CONE_OVERLAP = 1118

    A variable is specified to be a member of several cones.

MSK_RES_ERR_MPS_CONE_REPEAT = 1119

    A variable is repeated within the CSECTION.

MSK_RES_ERR_MPS_CONE_TYPE = 1117

    Invalid cone type specified in a  CSECTION.

MSK_RES_ERR_MPS_FILE = 1100

    An error occurred while reading an MPS file.

MSK_RES_ERR_MPS_INV_BOUND_KEY = 1108

    An invalid bound key occurred in an MPS file.

MSK_RES_ERR_MPS_INV_CON_KEY = 1107

    An invalid constraint key occurred in an MPS file.

MSK_RES_ERR_MPS_INV_FIELD = 1101

    Invalid field occurred while reading an MPS file.

MSK_RES_ERR_MPS_INV_MARKER = 1102

    An invalid marker has been specified in the MPS file.

MSK_RES_ERR_MPS_INV_SEC_NAME = 1109

    An invalid section name occurred in an MPS file.

MSK_RES_ERR_MPS_INV_SEC_ORDER = 1115

    The sections in an MPS file is not in the correct order.

MSK_RES_ERR_MPS_INVALID_OBJ_NAME = 1128

    An invalid objective name is specified.

MSK_RES_ERR_MPS_INVALID_OBJSENSE = 1122

    An invalid objective sense is specified.

MSK_RES_ERR_MPS_MUL_CON_NAME = 1112

    A constraint name is specified multiple times in the ROWS section in an MPS file.

MSK_RES_ERR_MPS_MUL_CSEC = 1116

    Multiple CSECTIONs are given the same name.

MSK_RES_ERR_MPS_MUL_QOBJ = 1114

    The Q term in the objective is specified multiple times.

MSK_RES_ERR_MPS_MUL_QSEC = 1113

    Multiple QSECTIONs are specified for a constraint.

MSK_RES_ERR_MPS_NO_OBJECTIVE = 1110

    No objective is defined in an MPS file.

MSK_RES_ERR_MPS_NULL_CON_NAME = 1103

    An empty constraint name is used in an MPS file.

MSK_RES_ERR_MPS_NULL_VAR_NAME = 1104

    An empty variable name is used in an MPS file.

MSK_RES_ERR_MPS_SPLITTED_VAR = 1111

    The non-zero elements in A corresponding to a variable in an MPS file must be specified consecutively.

MSK_RES_ERR_MPS_TAB_IN_FIELD2 = 1125

    A tab char occurred in field 2.

MSK_RES_ERR_MPS_TAB_IN_FIELD3 = 1126

    A tab char occurred in field 3.

MSK_RES_ERR_MPS_TAB_IN_FIELD5 = 1127

    A tab char occurred in field 5.

MSK_RES_ERR_MPS_UNDEF_CON_NAME = 1105

    An undefined constraint name occurred in an MPS file.

MSK_RES_ERR_MPS_UNDEF_VAR_NAME = 1106

    An undefined variable name occurred in an MPS file.

MSK_RES_ERR_MUL_A_ELEMENT = 1254

    An element in A is defined multiple times.

MSK_RES_ERR_NAME_IS_NULL = 1760

    The name buffer is a NULL pointer.

MSK_RES_ERR_NAME_MAX_LEN = 1750

    A name is longer than the buffer that is supposed to hold it.

MSK_RES_ERR_NAN_IN_AIJ = 1473

    a[i,j] contains an invalid floating point value, i.e. a NaN.

MSK_RES_ERR_NAN_IN_BLC = 1461

    blc contains an invalid floating point value, i.e. a NaN.

MSK_RES_ERR_NAN_IN_BLX = 1471

    blx contains an invalid floating point value, i.e. a NaN.

MSK_RES_ERR_NAN_IN_BUC = 1462

    buc contains an invalid floating point value, i.e. a NaN.

MSK_RES_ERR_NAN_IN_BUX = 1472

    bux contains an invalid floating point value, i.e. a NaN.

MSK_RES_ERR_NAN_IN_C = 1470

    c contains an invalid floating point value, i.e. a NaN.

MSK_RES_ERR_NAN_IN_DOUBLE_DATA = 1450

    An invalid floating value was used in some double data.

MSK_RES_ERR_NEGATIVE_APPEND = 1264

    Cannot append a negative number.

MSK_RES_ERR_NEGATIVE_SURPLUS = 1263

    Negative surplus.

MSK_RES_ERR_NEWER_DLL = 1036

    The dynamic link library is newer than the specified version.

MSK_RES_ERR_NO_BARS_FOR_SOLUTION = 3916

    There is no bars available for the solution specified.

MSK_RES_ERR_NO_BARX_FOR_SOLUTION = 3915

    There is no barx available for the solution specified.

MSK_RES_ERR_NO_BASIS_SOL = 1600

    No basic solution is defined.

MSK_RES_ERR_NO_DUAL_FOR_ITG_SOL = 2950

    No dual information is available for the integer solution.

MSK_RES_ERR_NO_DUAL_INFEAS_CER = 2001

    A certificate of dual infeasibility is not available.

MSK_RES_ERR_NO_INIT_ENV = 1063

    Environment is not initialized.

MSK_RES_ERR_NO_OPTIMIZER_VAR_TYPE = 1552

    No optimizer is available for this class of optimization problems.

MSK_RES_ERR_NO_PRIMAL_INFEAS_CER = 2000

    A certificate of primal infeasibility is not available.

MSK_RES_ERR_NO_SNX_FOR_BAS_SOL = 2953

    snx is not available for the basis solution.

MSK_RES_ERR_NO_SOLUTION_IN_CALLBACK = 2500

    The required solution is not available.

MSK_RES_ERR_NON_UNIQUE_ARRAY = 5000

    An array does not contain unique elements.

MSK_RES_ERR_NONCONVEX = 1291

    The optimization problem is nonconvex.

MSK_RES_ERR_NONLINEAR_EQUALITY = 1290

    The model contains a nonlinear equality.

MSK_RES_ERR_NONLINEAR_FUNCTIONS_NOT_ALLOWED = 1428

    An operation that is invalid for problems with nonlinear functions defined has been attempted.

MSK_RES_ERR_NONLINEAR_RANGED = 1292

    The model contains a nonlinear ranged constraint.

MSK_RES_ERR_NR_ARGUMENTS = 1199

    Incorrect number of function arguments.

MSK_RES_ERR_NULL_ENV = 1060

    env is a NULL pointer.

MSK_RES_ERR_NULL_POINTER = 1065

    An argument to a function is unexpectedly a NULL pointer.

MSK_RES_ERR_NULL_TASK = 1061

    task is a NULL pointer.

MSK_RES_ERR_NUMCONLIM = 1250

    Maximum number of constraints limit is exceeded.

MSK_RES_ERR_NUMVARLIM = 1251

    Maximum number of variables limit is exceeded.

MSK_RES_ERR_OBJ_Q_NOT_NSD = 1296

    The quadratic coefficient matrix in the objective is not NSD.

MSK_RES_ERR_OBJ_Q_NOT_PSD = 1295

    The quadratic coefficient matrix in the objective is not PSD.

MSK_RES_ERR_OBJECTIVE_RANGE = 1260

    Empty objective range.

MSK_RES_ERR_OLDER_DLL = 1035

    The dynamic link library is older than the specified version.

MSK_RES_ERR_OPEN_DL = 1030

    A dynamic link library could not be opened.

MSK_RES_ERR_OPF_FORMAT = 1168

    Syntax error in an OPF file

MSK_RES_ERR_OPF_NEW_VARIABLE = 1169

    Variable not previously defined.

MSK_RES_ERR_OPF_PREMATURE_EOF = 1172

    Premature end of file in an OPF file.

MSK_RES_ERR_OPTIMIZER_LICENSE = 1013

    The optimizer required is not licensed.

MSK_RES_ERR_ORD_INVALID = 1131

    Invalid content in branch ordering file.

MSK_RES_ERR_ORD_INVALID_BRANCH_DIR = 1130

    An invalid branch direction key is specified.

MSK_RES_ERR_OVERFLOW = 1590

    A computation produced an overflow.

MSK_RES_ERR_PARAM_INDEX = 1210

    Parameter index is out of range.

MSK_RES_ERR_PARAM_IS_TOO_LARGE = 1215

    A parameter value is too large.

MSK_RES_ERR_PARAM_IS_TOO_SMALL = 1216

    A parameter value is too small.

MSK_RES_ERR_PARAM_NAME = 1205

    A parameter name is not correct.

MSK_RES_ERR_PARAM_NAME_DOU = 1206

    A parameter name is not correct.

MSK_RES_ERR_PARAM_NAME_INT = 1207

    A parameter name is not correct.

MSK_RES_ERR_PARAM_NAME_STR = 1208

    A parameter name is not correct.

MSK_RES_ERR_PARAM_TYPE = 1218

    A parameter type is invalid.

MSK_RES_ERR_PARAM_VALUE_STR = 1217

    A parameter value string is incorrect.

MSK_RES_ERR_PLATFORM_NOT_LICENSED = 1019

    A requested license feature is not available for the required platform.

MSK_RES_ERR_POSTSOLVE = 1580

    An error occurred during the postsolve.

MSK_RES_ERR_PRO_ITEM = 1281

    An invalid problem item is used.

MSK_RES_ERR_PROB_LICENSE = 1006

    The software is not licensed to solve the problem.

MSK_RES_ERR_QCON_SUBI_TOO_LARGE = 1409

    Invalid value in qcsubi.

MSK_RES_ERR_QCON_SUBI_TOO_SMALL = 1408

    Invalid value in qcsubi.

MSK_RES_ERR_QCON_UPPER_TRIANGLE = 1417

    An element in the upper triangle of the quadratic term in a constraint.

MSK_RES_ERR_QOBJ_UPPER_TRIANGLE = 1415

    An element in the upper triangle of the quadratic term in the objective is specified.

MSK_RES_ERR_READ_FORMAT = 1090

    The specified format cannot be read.

MSK_RES_ERR_READ_LP_MISSING_END_TAG = 1159

    Missing End tag in LP file.

MSK_RES_ERR_READ_LP_NONEXISTING_NAME = 1162

    A variable never occurred in objective or constraints.

MSK_RES_ERR_REMOVE_CONE_VARIABLE = 1310

    A variable cannot be removed because it will make a cone invalid.

MSK_RES_ERR_REPAIR_INVALID_PROBLEM = 1710

    The feasibility repair does not support the specified problem type.

MSK_RES_ERR_REPAIR_OPTIMIZATION_FAILED = 1711

    Computation the optimal relaxation failed.

MSK_RES_ERR_SEN_BOUND_INVALID_LO = 3054

    Analysis of lower bound requested for an index, where no lower bound exists.

MSK_RES_ERR_SEN_BOUND_INVALID_UP = 3053

    Analysis of upper bound requested for an index, where no upper bound exists.

MSK_RES_ERR_SEN_FORMAT = 3050

    Syntax error in sensitivity analysis file.

MSK_RES_ERR_SEN_INDEX_INVALID = 3055

    Invalid range given in the sensitivity file.

MSK_RES_ERR_SEN_INDEX_RANGE = 3052

    Index out of range in the sensitivity analysis file.

MSK_RES_ERR_SEN_INVALID_REGEXP = 3056

    Syntax error in regexp or regexp longer than 1024.

MSK_RES_ERR_SEN_NUMERICAL = 3058

    Numerical difficulties encountered performing the sensitivity analysis.

MSK_RES_ERR_SEN_SOLUTION_STATUS = 3057

    No optimal solution found to the original problem given for sensitivity analysis.

MSK_RES_ERR_SEN_UNDEF_NAME = 3051

    An undefined name was encountered in the sensitivity analysis file.

MSK_RES_ERR_SEN_UNHANDLED_PROBLEM_TYPE = 3080

    Sensitivity analysis cannot be performed for the specified problem.

MSK_RES_ERR_SIZE_LICENSE = 1005

    The problem is bigger than the license.

MSK_RES_ERR_SIZE_LICENSE_CON = 1010

    The problem has too many constraints.

MSK_RES_ERR_SIZE_LICENSE_INTVAR = 1012

    The problem contains too many integer variables.

MSK_RES_ERR_SIZE_LICENSE_NUMCORES = 3900

    The computer contains more cpu cores than the license allows for.

MSK_RES_ERR_SIZE_LICENSE_VAR = 1011

    The problem has too many variables.

MSK_RES_ERR_SOL_FILE_INVALID_NUMBER = 1350

    An invalid number is specified in a solution file.

MSK_RES_ERR_SOLITEM = 1237

    The solution number  solemn does not exists.

MSK_RES_ERR_SOLVER_PROBTYPE = 1259

    Problem type does not match the chosen optimizer.

MSK_RES_ERR_SPACE = 1051

    Out of space.

MSK_RES_ERR_SPACE_LEAKING = 1080

    MOSEK is leaking memory.

MSK_RES_ERR_SPACE_NO_INFO = 1081

    No available information about the space usage.

MSK_RES_ERR_SYM_MAT_DUPLICATE = 3944

    A value in a symmetric matric as been specified more than once.

MSK_RES_ERR_SYM_MAT_INVALID_COL_INDEX = 3941

    A column index specified for sparse symmetric matrix is invalid.

MSK_RES_ERR_SYM_MAT_INVALID_ROW_INDEX = 3940

    A row index specified for sparse symmetric matrix is invalid.

MSK_RES_ERR_SYM_MAT_INVALID_VALUE = 3943

    The numerical value specified in a sparse symmetric matrix is not a value floating value.

MSK_RES_ERR_SYM_MAT_NOT_LOWER_TRINGULAR = 3942

    Only the lower triangular part of sparse symmetric matrix should be specified.

MSK_RES_ERR_TASK_INCOMPATIBLE = 2560

    The Task file is incompatible with  this platform.

MSK_RES_ERR_TASK_INVALID = 2561

    The Task file is invalid.

MSK_RES_ERR_THREAD_COND_INIT = 1049

    Could not initialize a condition.

MSK_RES_ERR_THREAD_CREATE = 1048

    Could not create a thread.

MSK_RES_ERR_THREAD_MUTEX_INIT = 1045

    Could not initialize a mutex.

MSK_RES_ERR_THREAD_MUTEX_LOCK = 1046

    Could not lock a mutex.

MSK_RES_ERR_THREAD_MUTEX_UNLOCK = 1047

    Could not unlock a mutex.

MSK_RES_ERR_TOO_MANY_CONCURRENT_TASKS = 3090

    Too many concurrent tasks specified.

MSK_RES_ERR_TOO_SMALL_MAX_NUM_NZ = 1245

    The maximum number of non-zeros specified is too small.

MSK_RES_ERR_TOO_SMALL_MAXNUMANZ = 1252

    Too small maximum number of non-zeros in A specified.

MSK_RES_ERR_UNB_STEP_SIZE = 3100

    A step-size in an optimizer was unexpectedly unbounded.

MSK_RES_ERR_UNDEF_SOLUTION = 1265

    The required solution is not defined.

MSK_RES_ERR_UNDEFINED_OBJECTIVE_SENSE = 1446

    The objective sense has not been specified before the optimization.

MSK_RES_ERR_UNHANDLED_SOLUTION_STATUS = 6010

    Unhandled solution status.

MSK_RES_ERR_UNKNOWN = 1050

    Unknown error.

MSK_RES_ERR_UPPER_BOUND_IS_A_NAN = 1391

    The upper bound specified is not a number (nan).

MSK_RES_ERR_UPPER_TRIANGLE = 6020

    An element in the upper triangle of a lower triangular matrix is specified.

MSK_RES_ERR_USER_FUNC_RET = 1430

    An user function reported an error.

MSK_RES_ERR_USER_FUNC_RET_DATA = 1431

    An user function returned invalid data.

MSK_RES_ERR_USER_NLO_EVAL = 1433

    The user-defined nonlinear function reported an error.

MSK_RES_ERR_USER_NLO_EVAL_HESSUBI = 1440

    The user-defined nonlinear function reported an Hessian an invalid subscript.

MSK_RES_ERR_USER_NLO_EVAL_HESSUBJ = 1441

    The user-defined nonlinear function reported an invalid subscript in the Hessian.

MSK_RES_ERR_USER_NLO_FUNC = 1432

    The user-defined nonlinear function reported an error.

MSK_RES_ERR_WHICHITEM_NOT_ALLOWED = 1238

    whichitem is unacceptable.

MSK_RES_ERR_WHICHSOL = 1236

    The solution defined by whichsol does not exists.

MSK_RES_ERR_WRITE_LP_FORMAT = 1158

    Problem cannot be written as an LP file.

MSK_RES_ERR_WRITE_LP_NON_UNIQUE_NAME = 1161

    An auto-generated name is not unique.

MSK_RES_ERR_WRITE_MPS_INVALID_NAME = 1153

    An invalid name is created while writing an MPS file.

MSK_RES_ERR_WRITE_OPF_INVALID_VAR_NAME = 1156

    Empty variable names cannot be written to OPF files.

MSK_RES_ERR_WRITING_FILE = 1166

    An error occurred while writing file

MSK_RES_ERR_XML_INVALID_PROBLEM_TYPE = 3600

    The problem type is not supported by the XML format.

MSK_RES_ERR_Y_IS_UNDEFINED = 1449

    The solution item y is undefined.

MSK_RES_OK = 0

    No error occurred.

MSK_RES_TRM_INTERNAL = 10030

    The optimizer terminated due to some internal reason.

MSK_RES_TRM_INTERNAL_STOP = 10031

    The optimizer terminated for internal reasons.

MSK_RES_TRM_MAX_ITERATIONS = 10000

    The optimizer terminated at the maximum number of iterations.

MSK_RES_TRM_MAX_NUM_SETBACKS = 10020

    The optimizer terminated as the maximum number of set-backs was reached.

MSK_RES_TRM_MAX_TIME = 10001

    The optimizer terminated at the maximum amount of time.

MSK_RES_TRM_MIO_NEAR_ABS_GAP = 10004

    The mixed-integer optimizer terminated because the near optimal absolute gap tolerance was satisfied.

MSK_RES_TRM_MIO_NEAR_REL_GAP = 10003

    The mixed-integer optimizer terminated because the near optimal relative gap tolerance was satisfied.

MSK_RES_TRM_MIO_NUM_BRANCHES = 10009

    The mixed-integer optimizer terminated as to the maximum number of branches was reached.

MSK_RES_TRM_MIO_NUM_RELAXS = 10008

    The mixed-integer optimizer terminated as the maximum number of relaxations was reached.

MSK_RES_TRM_NUM_MAX_NUM_INT_SOLUTIONS = 10015

    The mixed-integer optimizer terminated as the maximum number of feasible solutions was reached.

MSK_RES_TRM_NUMERICAL_PROBLEM = 10025

    The optimizer terminated due to a numerical problem.

MSK_RES_TRM_OBJECTIVE_RANGE = 10002

    The optimizer terminated on the bound of the objective range.

MSK_RES_TRM_STALL = 10006

    The optimizer is terminated due to slow progress.

MSK_RES_TRM_USER_CALLBACK = 10007

    The user-defined progress call-back function terminated the optimization.

MSK_RES_WRN_ANA_ALMOST_INT_BOUNDS = 904

    Warn against almost integral bounds.

MSK_RES_WRN_ANA_C_ZERO = 901

    Warn against all objective coefficients being zero.

MSK_RES_WRN_ANA_CLOSE_BOUNDS = 903

    Warn against close bounds.

MSK_RES_WRN_ANA_EMPTY_COLS = 902

    Warn against empty columns.

MSK_RES_WRN_ANA_LARGE_BOUNDS = 900

    Warn against very large bounds.

MSK_RES_WRN_CONSTRUCT_INVALID_SOL_ITG = 807

    The initial value for one or more  of the integer variables is not feasible.

MSK_RES_WRN_CONSTRUCT_NO_SOL_ITG = 810

    The construct solution requires an integer solution.

MSK_RES_WRN_CONSTRUCT_SOLUTION_INFEAS = 805

    After fixing the integer variables at the suggested values then the problem is infeasible.

MSK_RES_WRN_DROPPED_NZ_QOBJ = 201

    One or more non-zero elements were dropped in the Q matrix in the objective.

MSK_RES_WRN_DUPLICATE_BARVARIABLE_NAMES = 852

    Two barvariable names are identical.

MSK_RES_WRN_DUPLICATE_CONE_NAMES = 853

    Two cone names are identical.

MSK_RES_WRN_DUPLICATE_CONSTRAINT_NAMES = 850

    Two constraint names are identical.

MSK_RES_WRN_DUPLICATE_VARIABLE_NAMES = 851

    Two variable names are identical.

MSK_RES_WRN_ELIMINATOR_SPACE = 801

    The eliminator is skipped at least once due to lack of space.

MSK_RES_WRN_EMPTY_NAME = 502

    A variable or constraint name is empty. The output file may be invalid.

MSK_RES_WRN_IGNORE_INTEGER = 250

    Ignored integer constraints.

MSK_RES_WRN_INCOMPLETE_LINEAR_DEPENDENCY_CHECK = 800

    The linear dependency check(s) is incomplete.

MSK_RES_WRN_LARGE_AIJ = 62

    A numerically large value is specified for an element in A.

MSK_RES_WRN_LARGE_BOUND = 51

    A numerically large bound value is specified.

MSK_RES_WRN_LARGE_CJ = 57

    A numerically large value is specified for one element in A.

MSK_RES_WRN_LARGE_CON_FX = 54

    A equality constraint is fixed to numerically large value.

MSK_RES_WRN_LARGE_LO_BOUND = 52

    A numerically large lower bound value is specified.

MSK_RES_WRN_LARGE_UP_BOUND = 53

    A numerically large upper bound value is specified.

MSK_RES_WRN_LICENSE_EXPIRE = 500

    The license expires.

MSK_RES_WRN_LICENSE_FEATURE_EXPIRE = 505

    The license expires.

MSK_RES_WRN_LICENSE_SERVER = 501

    The license server is not responding.

MSK_RES_WRN_LP_DROP_VARIABLE = 85

    Ignore a variable because the variable was not previously defined.

MSK_RES_WRN_LP_OLD_QUAD_FORMAT = 80

    Missing '/2' after quadratic expressions in bound or objective.

MSK_RES_WRN_MIO_INFEASIBLE_FINAL = 270

    The final mixed-integer problem with all the integer variables fixed at their optimal values is infeasible.

MSK_RES_WRN_MPS_SPLIT_BOU_VECTOR = 72

    A BOUNDS vector is split into several nonadjacent parts in an MPS file.

MSK_RES_WRN_MPS_SPLIT_RAN_VECTOR = 71

    A RANGE vector is split into several nonadjacent parts in an MPS file.

MSK_RES_WRN_MPS_SPLIT_RHS_VECTOR = 70

    An RHS vector is split into several nonadjacent parts.

MSK_RES_WRN_NAME_MAX_LEN = 65

    A name is longer than the buffer that is supposed to hold it.

MSK_RES_WRN_NO_DUALIZER = 950

    No automatic dualizer is available for the specified problem.

MSK_RES_WRN_NO_GLOBAL_OPTIMIZER = 251

    No global optimizer is available.

MSK_RES_WRN_NO_NONLINEAR_FUNCTION_WRITE = 450

    The problem contains a general nonlinear function that cannot be written to a disk file.

MSK_RES_WRN_NZ_IN_UPR_TRI = 200

    Non-zero elements specified in the upper triangle of a matrix were ignored.

MSK_RES_WRN_OPEN_PARAM_FILE = 50

    The parameter file could not be opened.

MSK_RES_WRN_PARAM_IGNORED_CMIO = 516

    A parameter was ignored by the conic mixed integer optimizer.

MSK_RES_WRN_PARAM_NAME_DOU = 510

    Parameter name not recognized.

MSK_RES_WRN_PARAM_NAME_INT = 511

    Parameter name not recognized.

MSK_RES_WRN_PARAM_NAME_STR = 512

    Parameter name not recognized.

MSK_RES_WRN_PARAM_STR_VALUE = 515

    A parameter value is not correct.

MSK_RES_WRN_PRESOLVE_OUTOFSPACE = 802

    The presolve is incomplete due to lack of space.

MSK_RES_WRN_QUAD_CONES_WITH_ROOT_FIXED_AT_ZERO = 930

    For at least one quadratic cone the root is fixed at (nearly) zero.

MSK_RES_WRN_RQUAD_CONES_WITH_ROOT_FIXED_AT_ZERO = 931

    For at least one rotated quadratic cone the root is fixed at (nearly) zero.

MSK_RES_WRN_SOL_FILE_IGNORED_CON = 351

    One or more lines in the constraint section were ignored when reading a solution file.

MSK_RES_WRN_SOL_FILE_IGNORED_VAR = 352

    One or more lines in the variable section were ignored when reading a solution file.

MSK_RES_WRN_SOL_FILTER = 300

    Invalid solution filter is specified.

MSK_RES_WRN_SPAR_MAX_LEN = 66

    A value for a string parameter is longer than the buffer that is supposed to hold it.

MSK_RES_WRN_TOO_FEW_BASIS_VARS = 400

    An incomplete basis is specified.

MSK_RES_WRN_TOO_MANY_BASIS_VARS = 405

    A basis with too many variables is specified.

MSK_RES_WRN_TOO_MANY_THREADS_CONCURRENT = 750

    The concurrent optimizer employs more threads than available.

MSK_RES_WRN_UNDEF_SOL_FILE_NAME = 350

    Undefined name occurred in a solution.

MSK_RES_WRN_USING_GENERIC_NAMES = 503

    Generic names are used because a name is not valid.

MSK_RES_WRN_WRITE_CHANGED_NAMES = 803

    Some names were changed because they were invalid for the output file format.

MSK_RES_WRN_WRITE_DISCARDED_CFIX = 804

    The fixed objective term was discarded in the output file.

MSK_RES_WRN_ZERO_AIJ = 63

    One or more zero elements are specified in A.

MSK_RES_WRN_ZEROS_IN_SPARSE_COL = 710

    One or more (near) zero elements are specified in a sparse column of a matrix.

MSK_RES_WRN_ZEROS_IN_SPARSE_ROW = 705

    One or more (near) zero elements are specified in a sparse row of a matrix.

Enum mionodeseltype
-------------------

MSK_MIO_NODE_SELECTION_BEST = 2

    The optimizer employs a best bound node selection strategy.

MSK_MIO_NODE_SELECTION_FIRST = 1

    The optimizer employs a depth first node selection strategy.

MSK_MIO_NODE_SELECTION_FREE = 0

    The optimizer decides the node selection strategy.

MSK_MIO_NODE_SELECTION_HYBRID = 4

    The optimizer employs a hybrid strategy.

MSK_MIO_NODE_SELECTION_PSEUDO = 5

    The optimizer employs selects the node based on a pseudo cost estimate.

MSK_MIO_NODE_SELECTION_WORST = 3

    The optimizer employs a worst bound node selection strategy.

Enum onoffkey
-------------

MSK_OFF = 0

    Switch the option off.

MSK_ON = 1

    Switch the option on.

Enum simdegen
-------------

MSK_SIM_DEGEN_AGGRESSIVE = 2

    The simplex optimizer should use an aggressive degeneration strategy.

MSK_SIM_DEGEN_FREE = 1

    The simplex optimizer chooses the degeneration strategy.

MSK_SIM_DEGEN_MINIMUM = 4

    The simplex optimizer should use a minimum degeneration strategy.

MSK_SIM_DEGEN_MODERATE = 3

    The simplex optimizer should use a moderate degeneration strategy.

MSK_SIM_DEGEN_NONE = 0

    The simplex optimizer should use no degeneration strategy.

Enum dataformat
---------------

MSK_DATA_FORMAT_CB = 7

    Conic benchmark format,

MSK_DATA_FORMAT_EXTENSION = 0

    The file extension is used to determine the data file format.

MSK_DATA_FORMAT_FREE_MPS = 5

    The data a free MPS formatted file.

MSK_DATA_FORMAT_LP = 2

    The data file is LP formatted.

MSK_DATA_FORMAT_MPS = 1

    The data file is MPS formatted.

MSK_DATA_FORMAT_OP = 3

    The data file is an optimization problem formatted file.

MSK_DATA_FORMAT_TASK = 6

    Generic task dump file.

MSK_DATA_FORMAT_XML = 4

    The data file is an XML formatted file.

Enum orderingtype
-----------------

MSK_ORDER_METHOD_APPMINLOC = 1

    Approximate minimum local fill-in ordering is employed.

MSK_ORDER_METHOD_EXPERIMENTAL = 2

    This option should not be used.

MSK_ORDER_METHOD_FORCE_GRAPHPAR = 4

    Always use the graph partitioning based ordering even if it is worse that the approximate minimum local fill ordering.

MSK_ORDER_METHOD_FREE = 0

    The ordering method is chosen automatically.

MSK_ORDER_METHOD_NONE = 5

    No ordering is used.

MSK_ORDER_METHOD_TRY_GRAPHPAR = 3

    Always try the graph partitioning based ordering.

Enum problemtype
----------------

MSK_PROBTYPE_CONIC = 4

    A conic optimization.

MSK_PROBTYPE_GECO = 3

    General convex optimization.

MSK_PROBTYPE_LO = 0

    The problem is
    a linear optimization problem.

MSK_PROBTYPE_MIXED = 5

    General nonlinear constraints and conic constraints. This combination can not be solved by MOSEK.

MSK_PROBTYPE_QCQO = 2

    The problem is a quadratically constrained optimization problem.

MSK_PROBTYPE_QO = 1

    The problem is a quadratic optimization problem.

Enum inftype
------------

MSK_INF_DOU_TYPE = 0

    Is a double information type.

MSK_INF_INT_TYPE = 1

    Is an integer.

MSK_INF_LINT_TYPE = 2

    Is a long integer.

Enum dparam
-----------

MSK_DPAR_ANA_SOL_INFEAS_TOL = 0

    If a constraint violates its bound with an amount larger than this value,
    the constraint name, index and violation will be printed by the solution analyzer.

MSK_DPAR_BASIS_REL_TOL_S = 1

    Maximum relative dual bound violation allowed in an optimal
    basic solution.

MSK_DPAR_BASIS_TOL_S = 2

    Maximum absolute dual bound violation in
    an optimal basic solution.

MSK_DPAR_BASIS_TOL_X = 3

    Maximum absolute primal bound violation allowed
    in an optimal basic solution.

MSK_DPAR_CHECK_CONVEXITY_REL_TOL = 4

    Convexity check tolerance.

MSK_DPAR_DATA_TOL_AIJ = 5

    Data tolerance threshold.

MSK_DPAR_DATA_TOL_AIJ_HUGE = 6

    Data tolerance threshold.

MSK_DPAR_DATA_TOL_AIJ_LARGE = 7

    Data tolerance threshold.

MSK_DPAR_DATA_TOL_BOUND_INF = 8

    Data tolerance threshold.

MSK_DPAR_DATA_TOL_BOUND_WRN = 9

    Data tolerance threshold.

MSK_DPAR_DATA_TOL_C_HUGE = 10

    Data tolerance threshold.

MSK_DPAR_DATA_TOL_CJ_LARGE = 11

    Data tolerance threshold.

MSK_DPAR_DATA_TOL_QIJ = 12

    Data tolerance threshold.

MSK_DPAR_DATA_TOL_X = 13

    Data tolerance threshold.

MSK_DPAR_FEASREPAIR_TOL = 14

    Tolerance for constraint enforcing upper bound on
    sum of weighted violations in feasibility repair.

MSK_DPAR_INTPNT_CO_TOL_DFEAS = 15

    Dual feasibility tolerance used by the conic interior-point optimizer.

MSK_DPAR_INTPNT_CO_TOL_INFEAS = 16

    Infeasibility tolerance for the conic solver.

MSK_DPAR_INTPNT_CO_TOL_MU_RED = 17

    Optimality tolerance for the conic solver.

MSK_DPAR_INTPNT_CO_TOL_NEAR_REL = 18

    Optimality tolerance for the conic solver.

MSK_DPAR_INTPNT_CO_TOL_PFEAS = 19

    Primal feasibility tolerance used by the conic interior-point optimizer.

MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 20

    Relative gap termination tolerance used by the
    conic interior-point optimizer.

MSK_DPAR_INTPNT_NL_MERIT_BAL = 21

    Controls if the complementarity and infeasibility is converging to zero
    at about equal rates.

MSK_DPAR_INTPNT_NL_TOL_DFEAS = 22

    Dual feasibility tolerance used when a nonlinear
    model is solved.

MSK_DPAR_INTPNT_NL_TOL_MU_RED = 23

    Relative complementarity gap tolerance.

MSK_DPAR_INTPNT_NL_TOL_NEAR_REL = 24

    Nonlinear solver optimality tolerance parameter.

MSK_DPAR_INTPNT_NL_TOL_PFEAS = 25

    Primal feasibility tolerance used when a nonlinear
    model is solved.

MSK_DPAR_INTPNT_NL_TOL_REL_GAP = 26

    Relative gap termination tolerance for nonlinear problems.

MSK_DPAR_INTPNT_NL_TOL_REL_STEP = 27

    Relative step size to the boundary
    for general nonlinear optimization problems.

MSK_DPAR_INTPNT_TOL_DFEAS = 28

    Dual feasibility tolerance used for
    linear and quadratic optimization problems.

MSK_DPAR_INTPNT_TOL_DSAFE = 29

    Controls the interior-point dual starting point.

MSK_DPAR_INTPNT_TOL_INFEAS = 30

    Nonlinear solver infeasibility tolerance parameter.

MSK_DPAR_INTPNT_TOL_MU_RED = 31

    Relative complementarity gap tolerance.

MSK_DPAR_INTPNT_TOL_PATH = 32

    interior-point centering aggressiveness.

MSK_DPAR_INTPNT_TOL_PFEAS = 33

    Primal feasibility tolerance used for
    linear and quadratic optimization problems.

MSK_DPAR_INTPNT_TOL_PSAFE = 34

    Controls the interior-point primal starting point.

MSK_DPAR_INTPNT_TOL_REL_GAP = 35

    Relative gap termination tolerance.

MSK_DPAR_INTPNT_TOL_REL_STEP = 36

    Relative step size to the boundary
    for linear and quadratic optimization problems.

MSK_DPAR_INTPNT_TOL_STEP_SIZE = 37

    If the step size falls below the value of this parameter, then
    the interior-point optimizer assumes that it is stalled. In other words the interior-point optimizer does not make any progress and therefore it is better stop.

MSK_DPAR_LOWER_OBJ_CUT = 38

    Objective bound.

MSK_DPAR_LOWER_OBJ_CUT_FINITE_TRH = 39

    Objective bound.

MSK_DPAR_MIO_DISABLE_TERM_TIME = 40

    Certain termination criteria is disabled within the mixed-integer optimizer for period time specified by the parameter.

MSK_DPAR_MIO_HEURISTIC_TIME = 41

    Time limit for the mixed-integer heuristic.

MSK_DPAR_MIO_MAX_TIME = 42

    Time limit for the mixed-integer optimizer.

MSK_DPAR_MIO_MAX_TIME_APRX_OPT = 43

    Time limit for the mixed-integer optimizer.

MSK_DPAR_MIO_NEAR_TOL_ABS_GAP = 44

    Relaxed absolute optimality tolerance employed by the mixed-integer optimizer.

MSK_DPAR_MIO_NEAR_TOL_REL_GAP = 45

    The mixed-integer optimizer is terminated when this tolerance is satisfied.

MSK_DPAR_MIO_REL_ADD_CUT_LIMITED = 46

    Controls cut generation for mixed-integer optimizer.

MSK_DPAR_MIO_REL_GAP_CONST = 47

    This value is used to compute the relative gap for the solution to an integer optimization problem.

MSK_DPAR_MIO_TOL_ABS_GAP = 48

    Absolute optimality tolerance employed by the mixed-integer optimizer.

MSK_DPAR_MIO_TOL_ABS_RELAX_INT = 49

    Integer constraint tolerance.

MSK_DPAR_MIO_TOL_FEAS = 50

    Feasibility tolerance for mixed integer solver. Any solution with
    maximum infeasibility below this value will be considered feasible.

MSK_DPAR_MIO_TOL_REL_GAP = 51

    Relative optimality tolerance employed by the mixed-integer optimizer.

MSK_DPAR_MIO_TOL_REL_RELAX_INT = 52

    Integer constraint tolerance.

MSK_DPAR_MIO_TOL_X = 53

    Absolute solution tolerance used in mixed-integer optimizer.

MSK_DPAR_NONCONVEX_TOL_FEAS = 54

    Feasibility tolerance used by the nonconvex optimizer.

MSK_DPAR_NONCONVEX_TOL_OPT = 55

    Optimality tolerance used by the nonconvex optimizer.

MSK_DPAR_OPTIMIZER_MAX_TIME = 56

    Solver time limit.

MSK_DPAR_PRESOLVE_TOL_ABS_LINDEP = 57

    Absolute tolerance employed by the
    linear dependency checker.

MSK_DPAR_PRESOLVE_TOL_AIJ = 58

    Absolute zero tolerance employed for constraint coefficients in the presolve.

MSK_DPAR_PRESOLVE_TOL_REL_LINDEP = 59

    Relative tolerance employed by the
    linear dependency checker.

MSK_DPAR_PRESOLVE_TOL_S = 60

    Absolute zero tolerance employed for slack variables in the presolve.

MSK_DPAR_PRESOLVE_TOL_X = 61

    Absolute zero tolerance employed for variables in the presolve.

MSK_DPAR_QCQO_REFORMULATE_REL_DROP_TOL = 62

    This parameter determines when columns are dropped in incomplete Cholesky factorization doing reformulation of quadratic problems.

MSK_DPAR_SIM_LU_TOL_REL_PIV = 63

    Relative pivot tolerance employed when computing the LU factorization of the basis matrix.

MSK_DPAR_SIMPLEX_ABS_TOL_PIV = 64

    Absolute pivot tolerance employed by the simplex optimizers.

MSK_DPAR_UPPER_OBJ_CUT = 65

    Objective bound.

MSK_DPAR_UPPER_OBJ_CUT_FINITE_TRH = 66

    Objective bound.

Enum simdupvec
--------------

MSK_SIM_EXPLOIT_DUPVEC_FREE = 2

    The simplex optimizer can choose freely.

MSK_SIM_EXPLOIT_DUPVEC_OFF = 0

    Disallow the simplex optimizer to exploit duplicated columns.

MSK_SIM_EXPLOIT_DUPVEC_ON = 1

    Allow the simplex optimizer to exploit duplicated columns.

Enum compresstype
-----------------

MSK_COMPRESS_FREE = 1

    The type of compression used is chosen automatically.

MSK_COMPRESS_GZIP = 2

    The type of compression used is gzip compatible.

MSK_COMPRESS_NONE = 0

    No compression is used.

Enum nametype
-------------

MSK_NAME_TYPE_GEN = 0

    General names. However, no duplicate and blank names are allowed.

MSK_NAME_TYPE_LP = 2

    LP type names.

MSK_NAME_TYPE_MPS = 1

    MPS type names.

Enum mpsformat
--------------

MSK_MPS_FORMAT_FREE = 2

    It is assumed that the input file satisfies the free
    MPS format. This implies that spaces
    are not allowed in names. Otherwise
    the format is free.

MSK_MPS_FORMAT_RELAXED = 1

    It is assumed that the input file satisfies
    a slightly relaxed version of the MPS format.

MSK_MPS_FORMAT_STRICT = 0

    It is assumed that the input file satisfies
    the MPS format strictly.

Enum variabletype
-----------------

MSK_VAR_TYPE_CONT = 0

    Is a continuous variable.

MSK_VAR_TYPE_INT = 1

    Is an integer variable.

Enum checkconvexitytype
-----------------------

MSK_CHECK_CONVEXITY_FULL = 2

    Perform a full convexity check.

MSK_CHECK_CONVEXITY_NONE = 0

    No convexity check.

MSK_CHECK_CONVEXITY_SIMPLE = 1

    Perform simple and fast convexity check.

Enum language
-------------

MSK_LANG_DAN = 1

    Danish language selection

MSK_LANG_ENG = 0

    English language selection

Enum startpointtype
-------------------

MSK_STARTING_POINT_CONSTANT = 2

    The optimizer constructs a starting point by assigning a constant value to all primal and dual variables.
    This starting point is normally robust.

MSK_STARTING_POINT_FREE = 0

    The starting point is chosen automatically.

MSK_STARTING_POINT_GUESS = 1

    The optimizer guesses a starting point.

MSK_STARTING_POINT_SATISFY_BOUNDS = 3

    The starting point is chosen to satisfy all the simple bounds on nonlinear variables. If this starting point is employed,
    then more care than usual should employed when choosing the bounds on the nonlinear variables. In particular very tight bounds
    should be avoided.

Enum soltype
------------

MSK_SOL_BAS = 1

    The basic solution.

MSK_SOL_ITG = 2

    The integer solution.

MSK_SOL_ITR = 0

    The interior solution.

Enum scalingmethod
------------------

MSK_SCALING_METHOD_FREE = 1

    The optimizer chooses the scaling heuristic.

MSK_SCALING_METHOD_POW2 = 0

    Scales only with power of 2 leaving the mantissa untouched.

Enum value
----------

MSK_LICENSE_BUFFER_LENGTH = 20

    The length of a license key buffer.

MSK_MAX_STR_LEN = 1024

    Maximum string length allowed in MOSEK.

Enum stakey
-----------

MSK_SK_BAS = 1

    The constraint or variable is in the basis.

MSK_SK_FIX = 5

    The constraint or variable is fixed.

MSK_SK_INF = 6

    The constraint or variable is infeasible in the bounds.

MSK_SK_LOW = 3

    The constraint or variable is at its lower bound.

MSK_SK_SUPBAS = 2

    The constraint or variable is super basic.

MSK_SK_UNK = 0

    The status for the constraint or variable is unknown.

MSK_SK_UPR = 4

    The constraint or variable is at its upper bound.

Enum simreform
--------------

MSK_SIM_REFORMULATION_AGGRESSIVE = 3

    The simplex optimizer should use an aggressive reformulation strategy.

MSK_SIM_REFORMULATION_FREE = 2

    The simplex optimizer can choose freely.

MSK_SIM_REFORMULATION_OFF = 0

    Disallow the simplex optimizer to reformulate the problem.

MSK_SIM_REFORMULATION_ON = 1

    Allow the simplex optimizer to reformulate the problem.

Enum iinfitem
-------------

MSK_IINF_ANA_PRO_NUM_CON = 0

    Number of constraints in the problem.

MSK_IINF_ANA_PRO_NUM_CON_EQ = 1

    Number of equality constraints.

MSK_IINF_ANA_PRO_NUM_CON_FR = 2

    Number of unbounded constraints.

MSK_IINF_ANA_PRO_NUM_CON_LO = 3

    Number of constraints with a lower bound and an
    infinite upper bound.

MSK_IINF_ANA_PRO_NUM_CON_RA = 4

    Number of constraints with finite lower and upper bounds.

MSK_IINF_ANA_PRO_NUM_CON_UP = 5

    Number of constraints with an upper bound and an infinite lower bound.

MSK_IINF_ANA_PRO_NUM_VAR = 6

    Number of variables in the problem.

MSK_IINF_ANA_PRO_NUM_VAR_BIN = 7

    Number of binary variables.

MSK_IINF_ANA_PRO_NUM_VAR_CONT = 8

    Number of continuous variables.

MSK_IINF_ANA_PRO_NUM_VAR_EQ = 9

    Number of fixed variables.

MSK_IINF_ANA_PRO_NUM_VAR_FR = 10

    Number of unbounded constraints.

MSK_IINF_ANA_PRO_NUM_VAR_INT = 11

    Number of general integer variables.

MSK_IINF_ANA_PRO_NUM_VAR_LO = 12

    Number of variables with a lower bound and an
    infinite upper bound.

MSK_IINF_ANA_PRO_NUM_VAR_RA = 13

    Number of variables with finite lower and upper bounds.

MSK_IINF_ANA_PRO_NUM_VAR_UP = 14

    Number of variables with an upper bound and an infinite lower bound.

MSK_IINF_CONCURRENT_FASTEST_OPTIMIZER = 15

    The type of the optimizer that finished first in a concurrent optimization.

MSK_IINF_INTPNT_FACTOR_DIM_DENSE = 16

    Dimension of the dense sub system in factorization.

MSK_IINF_INTPNT_ITER = 17

    Number of interior-point iterations
    since invoking the interior-point optimizer.

MSK_IINF_INTPNT_NUM_THREADS = 18

    Number of threads that the interior-point optimizer is using.

MSK_IINF_INTPNT_SOLVE_DUAL = 19

    Non-zero if the interior-point optimizer is solving the dual problem.

MSK_IINF_MIO_CONSTRUCT_NUM_ROUNDINGS = 20

    Number of values in the integer solution that is rounded to an integer value.

MSK_IINF_MIO_CONSTRUCT_SOLUTION = 21

    If this item is positive, then MOSEK successfully constructed an initial integer feasible solution.

MSK_IINF_MIO_INITIAL_SOLUTION = 22

    Is non-zero if an initial integer solution is specified.

MSK_IINF_MIO_NUM_ACTIVE_NODES = 23

    Number of active branch bound nodes.

MSK_IINF_MIO_NUM_BASIS_CUTS = 24

    Number of basis cuts.

MSK_IINF_MIO_NUM_BRANCH = 25

    Number of branches performed during the optimization.

MSK_IINF_MIO_NUM_CARDGUB_CUTS = 26

    Number of cardgub cuts.

MSK_IINF_MIO_NUM_CLIQUE_CUTS = 27

    Number of clique cuts.

MSK_IINF_MIO_NUM_COEF_REDC_CUTS = 28

    Number of coef. redc. cuts.

MSK_IINF_MIO_NUM_CONTRA_CUTS = 29

    Number of contra cuts.

MSK_IINF_MIO_NUM_DISAGG_CUTS = 30

    Number of diasagg cuts.

MSK_IINF_MIO_NUM_FLOW_COVER_CUTS = 31

    Number of flow cover cuts.

MSK_IINF_MIO_NUM_GCD_CUTS = 32

    Number of gcd cuts.

MSK_IINF_MIO_NUM_GOMORY_CUTS = 33

    Number of Gomory cuts.

MSK_IINF_MIO_NUM_GUB_COVER_CUTS = 34

    Number of GUB cover cuts.

MSK_IINF_MIO_NUM_INT_SOLUTIONS = 35

    Number of integer feasible solutions that has been found.

MSK_IINF_MIO_NUM_KNAPSUR_COVER_CUTS = 36

    Number of knapsack cover cuts.

MSK_IINF_MIO_NUM_LATTICE_CUTS = 37

    Number of lattice cuts.

MSK_IINF_MIO_NUM_LIFT_CUTS = 38

    Number of lift cuts.

MSK_IINF_MIO_NUM_OBJ_CUTS = 39

    Number of obj cuts.

MSK_IINF_MIO_NUM_PLAN_LOC_CUTS = 40

    Number of loc cuts.

MSK_IINF_MIO_NUM_RELAX = 41

    Number of relaxations solved during the optimization.

MSK_IINF_MIO_NUMCON = 42

    Number of constraints in the problem solved be the mixed-integer optimizer.

MSK_IINF_MIO_NUMINT = 43

    Number of integer variables in the problem solved be the mixed-integer optimizer.

MSK_IINF_MIO_NUMVAR = 44

    Number of variables in the problem solved be the mixed-integer optimizer.

MSK_IINF_MIO_OBJ_BOUND_DEFINED = 45

    Non-zero if a valid objective bound has been found, otherwise zero.

MSK_IINF_MIO_TOTAL_NUM_CUTS = 46

    Total number of cuts generated by the mixed-integer optimizer.

MSK_IINF_MIO_USER_OBJ_CUT = 47

    If it is non-zero, then the objective cut is used.

MSK_IINF_OPT_NUMCON = 48

    Number of constraints in the problem solved when the optimizer is called.

MSK_IINF_OPT_NUMVAR = 49

    Number of variables in the problem solved when the optimizer is called

MSK_IINF_OPTIMIZE_RESPONSE = 50

    The response code returned by optimize.

MSK_IINF_RD_NUMBARVAR = 51

    Number of variables read.

MSK_IINF_RD_NUMCON = 52

    Number of constraints read.

MSK_IINF_RD_NUMCONE = 53

    Number of conic constraints read.

MSK_IINF_RD_NUMINTVAR = 54

    Number of integer-constrained variables read.

MSK_IINF_RD_NUMQ = 55

    Number of nonempty Q matrices read.

MSK_IINF_RD_NUMVAR = 56

    Number of variables read.

MSK_IINF_RD_PROTYPE = 57

    Problem type.

MSK_IINF_SIM_DUAL_DEG_ITER = 58

    The number of dual degenerate iterations.

MSK_IINF_SIM_DUAL_HOTSTART = 59

    If 1 then the dual simplex algorithm is solving from an advanced basis.

MSK_IINF_SIM_DUAL_HOTSTART_LU = 60

    If 1 then a valid basis factorization of full rank was located and used by the dual simplex algorithm.

MSK_IINF_SIM_DUAL_INF_ITER = 61

    The number of iterations taken with dual infeasibility.

MSK_IINF_SIM_DUAL_ITER = 62

    Number of dual simplex iterations during the last optimization.

MSK_IINF_SIM_NETWORK_DUAL_DEG_ITER = 63

    The number of dual network degenerate iterations.

MSK_IINF_SIM_NETWORK_DUAL_HOTSTART = 64

    If 1 then the dual network simplex algorithm is solving from an advanced basis.

MSK_IINF_SIM_NETWORK_DUAL_HOTSTART_LU = 65

    If 1 then a valid basis factorization of full rank was located and used by the dual network simplex algorithm.

MSK_IINF_SIM_NETWORK_DUAL_INF_ITER = 66

    The number of iterations taken with dual infeasibility in the network optimizer.

MSK_IINF_SIM_NETWORK_DUAL_ITER = 67

    Number of dual network simplex iterations during the last optimization.

MSK_IINF_SIM_NETWORK_PRIMAL_DEG_ITER = 68

    The number of primal network degenerate iterations.

MSK_IINF_SIM_NETWORK_PRIMAL_HOTSTART = 69

    If 1 then the primal network simplex algorithm is solving from an advanced basis.

MSK_IINF_SIM_NETWORK_PRIMAL_HOTSTART_LU = 70

    If 1 then a valid basis factorization of full rank was located and used by the primal network simplex algorithm.

MSK_IINF_SIM_NETWORK_PRIMAL_INF_ITER = 71

    The number of iterations taken with primal infeasibility in the network optimizer.

MSK_IINF_SIM_NETWORK_PRIMAL_ITER = 72

    Number of primal network simplex iterations during the last optimization.

MSK_IINF_SIM_NUMCON = 73

    Number of constraints in the problem solved by the simplex optimizer.

MSK_IINF_SIM_NUMVAR = 74

    Number of variables in the problem solved by the simplex optimizer.

MSK_IINF_SIM_PRIMAL_DEG_ITER = 75

    The number of primal degenerate iterations.

MSK_IINF_SIM_PRIMAL_DUAL_DEG_ITER = 76

    The number of degenerate major iterations taken by the primal dual simplex algorithm.

MSK_IINF_SIM_PRIMAL_DUAL_HOTSTART = 77

    If 1 then the primal dual simplex algorithm is solving from an advanced basis.

MSK_IINF_SIM_PRIMAL_DUAL_HOTSTART_LU = 78

    If 1 then a valid basis factorization of full rank was located and used by the primal dual simplex algorithm.

MSK_IINF_SIM_PRIMAL_DUAL_INF_ITER = 79

    The number of master iterations with dual infeasibility taken by the primal dual simplex algorithm.

MSK_IINF_SIM_PRIMAL_DUAL_ITER = 80

    Number of primal dual simplex iterations during the last optimization.

MSK_IINF_SIM_PRIMAL_HOTSTART = 81

    If 1 then the primal simplex algorithm is solving from an advanced basis.

MSK_IINF_SIM_PRIMAL_HOTSTART_LU = 82

    If 1 then a valid basis factorization of full rank was located and used by the primal simplex algorithm.

MSK_IINF_SIM_PRIMAL_INF_ITER = 83

    The number of iterations taken with primal infeasibility.

MSK_IINF_SIM_PRIMAL_ITER = 84

    Number of primal simplex iterations during the last optimization.

MSK_IINF_SIM_SOLVE_DUAL = 85

    Is non-zero if dual problem is solved.

MSK_IINF_SOL_BAS_PROSTA = 86

    Problem status of the basic solution. Updated after each optimization.

MSK_IINF_SOL_BAS_SOLSTA = 87

    Solution status of the basic solution. Updated after each optimization.

MSK_IINF_SOL_INT_PROSTA = 88

    Deprecated.

MSK_IINF_SOL_INT_SOLSTA = 89

    Deprecated.

MSK_IINF_SOL_ITG_PROSTA = 90

    Problem status of the integer solution. Updated after each optimization.

MSK_IINF_SOL_ITG_SOLSTA = 91

    Solution status of the integer solution. Updated after each optimization.

MSK_IINF_SOL_ITR_PROSTA = 92

    Problem status of the interior-point solution. Updated after each optimization.

MSK_IINF_SOL_ITR_SOLSTA = 93

    Solution status of the interior-point solution. Updated after each optimization.

MSK_IINF_STO_NUM_A_CACHE_FLUSHES = 94

    Number of times the cache of coefficient elements has been flushed.

MSK_IINF_STO_NUM_A_REALLOC = 95

    Number of times the storage for storing the linear coefficient matrix has been changed.

MSK_IINF_STO_NUM_A_TRANSPOSES = 96

    Number of times the coefficient matrix has been transposed.

Enum xmlwriteroutputtype
------------------------

MSK_WRITE_XML_MODE_COL = 1

    Write in column order.

MSK_WRITE_XML_MODE_ROW = 0

    Write in row order.

Enum optimizertype
------------------

MSK_OPTIMIZER_CONCURRENT = 10

    The optimizer for nonconvex nonlinear problems.

MSK_OPTIMIZER_CONIC = 2

    The optimizer for problems having conic constraints.

MSK_OPTIMIZER_DUAL_SIMPLEX = 4

    The dual simplex optimizer is used.

MSK_OPTIMIZER_FREE = 0

    The optimizer is chosen automatically.

MSK_OPTIMIZER_FREE_SIMPLEX = 6

    One of the simplex optimizers is used.

MSK_OPTIMIZER_INTPNT = 1

    The interior-point optimizer is used.

MSK_OPTIMIZER_MIXED_INT = 9

    The mixed-integer optimizer.

MSK_OPTIMIZER_MIXED_INT_CONIC = 8

    The mixed-integer optimizer for conic and linear problems.

MSK_OPTIMIZER_NETWORK_PRIMAL_SIMPLEX = 7

    The network primal simplex optimizer is used. It is only applicable to pure network problems.

MSK_OPTIMIZER_NONCONVEX = 11

    The optimizer for nonconvex nonlinear problems.

MSK_OPTIMIZER_PRIMAL_DUAL_SIMPLEX = 5

    The primal dual simplex optimizer is used.

MSK_OPTIMIZER_PRIMAL_SIMPLEX = 3

    The primal simplex optimizer is used.

Enum presolvemode
-----------------

MSK_PRESOLVE_MODE_FREE = 2

    It is decided automatically whether to presolve before the problem is optimized.

MSK_PRESOLVE_MODE_OFF = 0

    The problem is not presolved before it is optimized.

MSK_PRESOLVE_MODE_ON = 1

    The problem is presolved before it is optimized.

Enum miocontsoltype
-------------------

MSK_MIO_CONT_SOL_ITG = 2

    The reported interior-point and basic solutions are
    a solution to the problem with all integer variables
    fixed at the value they have in the integer solution.
    A solution is only reported in case the
    problem has a primal feasible solution.

MSK_MIO_CONT_SOL_ITG_REL = 3

    In case the problem is primal feasible
    then the reported interior-point and basic solutions
    are a solution to the problem with all integer variables
    fixed at the value they have in the integer solution.
    If the problem is primal infeasible, then the solution to the root node problem is reported.

MSK_MIO_CONT_SOL_NONE = 0

    No interior-point or basic solution are reported when the mixed-integer optimizer is used.

MSK_MIO_CONT_SOL_ROOT = 1

    The reported interior-point and basic solutions are a solution to the root node problem
    when mixed-integer optimizer is used.

