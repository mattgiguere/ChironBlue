# ChironBlue

Repository for the CHIRON Reduction code extended to the blue.

###Quick Start
- Download the [ChironBlue Repository](https://github.com/mattgiguere/ChironBlue).
- Download the [idlutils Repository](https://github.com/mattgiguere/idlutils)
- Download [IDLAstro](https://github.com/wlandsman/IDLAstro)
- Download [MPFIT](http://www.physics.wisc.edu/~craigm/idl/fitting.html)
- `cd` to the `ChironBlue` directory
- either move `ChironBlue`,  `idlutil`, `IDLAstro`, and `MPFIT` to `~/projects/` on your machine, or modify `$ChironBlue/chiblue` as needed
- Open `$ChironBlue/REDUCTION/ctio.par` in a text editor and modify the directories to match your setup accordingly. Line by line comments are included for help.
- at the command line, type `./chiblue`
- IDL should start up without error
- type `chi_reduce_all, date='yymmdd'` to reduce data

###Dependencies
- IDL 8.1 or later. `ChironBlue` will work with older versions of IDL, but you will need to modify `chiblue` to specify the path to IDL.
- [`idlutil`](https://github.com/mattgiguere/idlutils): a repository of handy routines

- [IDLAstro](https://github.com/wlandsman/IDLAstro): A repository of Astronomy related IDL routines
- [MPFIT](http://www.physics.wisc.edu/~craigm/idl/fitting.html): A Levenberg-Marquardt Least-Squares fitting package made by Craig Markwardt.
