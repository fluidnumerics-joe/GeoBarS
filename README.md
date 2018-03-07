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

## Examining Model Output
The model outputs a `residual.curve` and `GeoBarS.tec`.

`residual.curve` is a 2-D curve file with two columns. The first column is the iteration number of the GMRES-m algorithm. The second shows the L-2 Residual.

`GeoBarS.tec` is a tecplot file that contains the model parameters and output.

