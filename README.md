README
======

Repository for prototype implementations of dynamic Buchberger algorithms based
on the implementation of
[Caboara and Perry](http://www.math.usm.edu/perry/Research/dynamic_gb.pyx).
Includes additional *unrestricted* dynamic algorithms.

Everything requires SageMath (tested with SageMath 8.8).

Calling from SageMath
---------------------

Load the code with `load("dynamicgb.pyx")` and call the `dynamic_gb` function.
See the function docstring for additional information on its arguments.

Running experiments
-------------------

Run `sage test.sage`.

Making tables
-------------

Requires R/Rstudio. See analysis.Rmd.
