NOTE: This is a python code to calculate the proximal radial distribution function (pRDF or pG(r)) of a surface. 
Here, the code is used to analyze the molecular dynamics hydration data of a zwitterionic surface (polymer-brush). The surface
is made with a polymer chain of TMAO (pTMAO). The generated pG(r) for the pTMAO surface is published in the following JACS article.

Huang; Zhang et al. J. Am. Chem. Soc. 2021, 143, 40, 16786–16795 (See Figure 4a)
DOI: https://doi.org/10.1021/jacs.1c08280

How it works:
1. The surface is attached to SiO2 substrate via -(CH2)10 connectors, 
   the atoms of which are not of our interest. We will consider atoms of pTMAO, water, and ions.

2. This code will separate TMAO, water, ion-atoms from the SiO2 and connector atoms upon providing a reference index.
   The reference index is the begining of non-substrate/connector atoms assuming that the substrate/connector atoms 
   grouply indexed first, followed by TMAO/water/ion-atoms. For example, if the first 100 atoms in the trajectory file
   are substrate/connector atoms, the reference index would be 101. All atoms with index > 100 will be treated as
   surface/water/ion atoms.

3. We have to choose the distance (zdist = zmax - zmin)) in Z-direction (perpendicular distance to the surface)
