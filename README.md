README
======

Repository for prototype implementations of dynamic Buchberger algorithms based
on the implementation of
[Caboara and Perry](http://www.math.usm.edu/perry/Research/dynamic_gb.pyx).
Includes additional *unrestricted* dynamic algorithms.

Everything requires SageMath (tested with SageMath 8.8 and 9.0 on Ubuntu 18.04,
Fedora 30 and Linux Mint 19.3). The code requires Python 3 support from
SageMath. SageMath 9.0 comes with Python 3 support by default, but in previous
versions it was necessary to recompile SageMath from scratch. See 
[this](https://wiki.sagemath.org/Python3-compatible%20code)
for instructions

Calling from SageMath
---------------------

Load the code with `load("dynamicgb.pyx")` and call the `dynamic_gb` function.
See the function docstring for additional information on its arguments.

Running experiments
-------------------

Run `sage test.sage`.

Making tables and exploring the data
------------------------------------

Requires R/Rstudio. The recommended way to explore the data is to load
analysis.Rmd in Rstudio, as it organizes the raw data into tables. It is
then easy to restrict the tables to the desired algorithm or input ideal
using e.g. dplyr. There are many examples of this use in analysis.Rmd.
