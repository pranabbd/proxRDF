This is a python code to calculate the proximal radial distribution function (pRDF or pG(r)) of a surface. 
Here, the code is used to analyze the molecular dynamics hydration data of a zwitterionic surface (polymer-brush). The surface
is made with a polymer chain of TMAO (pTMAO). The generated pG(r) for the pTMAO surface is published in the following JACS article.

See Figure 4a in Huang; Zhang et al. J. Am. Chem. Soc. 2021, 143, 40, 16786–16795
DOI: https://doi.org/10.1021/jacs.1c08280

1. The surace is attached to SiO2 substrate, the atoms of which are not of our interest. We will consider atoms of pTMAO, water, and ions.
2. We have to choose the distance (zdist = zmax - zmin)) in Z-direction (perpendicular distance to the surface)
