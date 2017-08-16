# Angular Distribution Calculator

Maintainer: Jonathan Williams

Contributors: K. Starosta, J. Williams


## Description

Calculates the directional distribution for gamma radiation emitted from an axially symmetric oriented source.  In part based on `g77` code by K. Starosta for gamma-gamma angular Installation.

## Installation

Use 'make' to compile.  Tested and seems to work with gfortran on Ubuntu 14.04/16.04 and Scientific Linux 6.  Requires `libmathlib` (in the `cernlib` package on Ubuntu 14.04/16.04 and SL6).
If on SL6, you may need to replace the LIB line in the `Makefile` with:

```
LIB = /usr/lib64/cernlib/2006/lib/libmathlib.so.2_gfortran
```

### Notes from porting the code from `g77` to `gfortran`:

* used -std=legacy flag to supress warnings
* used -ffixed-line-length-none to avoid errors caused by lines that go over the default 72 character limit
