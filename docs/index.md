## GeoBarS
GeoBarS is a numerical model that solves the geostrophic and barotropic vorticity equation. These equations are discretized using the Continuous Galerkin Spectral Element Method. The discrete system is inverted using GMRES. 

The permits domain tesselation as an unstructured mesh of quadrilateral elements. This is particularly useful for this equation set, given the natural tendency for westward intensification on a beta-plane. The unstructured mesh allows for local mesh refinement in areas where higher resolution is needed.
