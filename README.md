# GeoBarS
A Geostrophic Barotropic Vorticity Solver


## Getting Started
Compiling the source code currently depends on the GNU Fortran compiler.
OpenMP can be enabled for the driver program by setting `OMP=yes` in `build/envfile`

Go to the `examples/LauBasinModel/` directory. Run 
`make geobars`

This makes the geobars executable.
`./geobars`
will run the model using the `runtime.params` namelist file and the included `box.mesh` file.

