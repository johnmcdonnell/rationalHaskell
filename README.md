The rational model of categorization (RMC)
=========================================

To build/ run sims
------------------

Just type make at the command line inside this directory:

    make

Be aware you will need to install the required packages listed below. After
installing the [Haskell Platform](http://www.haskell.org/platform/), you will
be able to install them by typing

    cabal install package-name

at the command line (where package-name is the name of the package you wanted
to install). If you still get import errors after installing all the required
packages, that probably means I forgot one so email me and I will update the
list. To get usage instructions for the binary, type:

    ./testanderson --help

Running the sims
----------------

Code in the fit\_glc.R file wraps the model. Sourcing this file from inside R
will run a bunch of sims, plot them, and leave them available to you in the
environment. Be aware the sims may take some time (20m+) to run; this can be
adjusted by changing the parameters to run and the number of reps to run, near
the end of the file.

Personally I run the fits from within [vim](http://vim.org) with the
[vim-r-plugin](http://www.vim.org/scripts/script.php?script_id=2628), using
folds to keep organized. I source all the functions, then adjust the code that
actually runs the models. To get appropriate folding in vim, type:

    :set foldmethod=marker

If you don't use vim, I'd recommend using the excellent 
[RStudio IDE](http://www.rstudio.com/ide/).

Required Packages
-----------------

- [hmatrix](http://hackage.haskell.org/package/hmatrix)
- [vector](http://hackage.haskell.org/package/vector)
- [vector-algorithms](http://hackage.haskell.org/package/vector-algorithms)
- [statistics](http://hackage.haskell.org/package/statistics)
- [random-shuffle](http://hackage.haskell.org/package/random-shuffle)
- [csv](http://hackage.haskell.org/package/csv)
- [cmdargs](http://hackage.haskell.org/package/cmdargs)

