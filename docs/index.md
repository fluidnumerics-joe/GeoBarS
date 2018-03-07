## GeoBarS
![useful image](geobars_logo.png)
GeoBarS is a numerical model that solves the geostrophic and barotropic vorticity equation. These equations are discretized using the Continuous Galerkin Spectral Element Method. The discrete system is inverted using GMRES. 

The permits domain tesselation as an unstructured mesh of quadrilateral elements. This is particularly useful for this equation set, given the natural tendency for westward intensification on a beta-plane. The unstructured mesh allows for local mesh refinement in areas where higher resolution is needed.

## Classic "Stommel" Gyre
This example can be found under the `examples/ClassicGyre/` directory in the GeoBarS repository. As in Stommel's 1948 paper on westward intensification, a sinusoidal wind stress drives the circulation on a beta-plane. Bottom drag provides the necessary dissipation for the model to have a well defined steady state. The combination of a varying coriolis parameter, wind stress, bottom drag, and no-normal flow boundary conditions result in the gyre circulation with a western boundary current as shown in the plot below.
![Westward Intensification Stream Function](classic_gyre.png)

## Lau Basin Idealized Abyssal Circulation

![useful image](streamfunction_3dview.0000.png)
![useful image](streamfunction_3dview.0008.png)
