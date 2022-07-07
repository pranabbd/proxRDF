# Surface Hydration Calculation 
This is a python code to calculate the proximal radial distribution function (pRDF or pG(r)) of a surface. 
Here, the code is used to analyze the molecular dynamics hydration data of a zwitterionic surface (polymer-brush). The surface
is made with a polymer chain of TMAO (pTMAO). The generated pG(r) for the pTMAO surface is published in the following JACS article.

Huang; Zhang et al. J. Am. Chem. Soc. 2021, 143, 40, 16786–16795 (See Figure 4a)
[!DOI]](https://doi.org/10.1021/jacs.1c08280)


## pG(r) 
pG(r) is a RDF for systems of irregular shape. It is calculated by taking the minimal distance between the solute-atoms 
and a refence atom (e.g., water-oxygen). This step is repeated for N solvent molecules. Since the solute is non-spherical, we
cannot obtain the local density by dividing the spherical shell volume as in standard RDF caculation. Here, we calculate the 
grid volume occupied by the solvent atoms in r and r+dr and divide the solvent atoms in that interval to get the local density. 
Finally, we normlize the local density by dividing the bulk density, averaged over distances between 10 and 15 angstrom.  

## How it works 
1. In the `test.xyz` file, there is a surface of polymer-TMAO ( or TMAO-polymer brush). The surface is attached to SiO2 substrate via -(CH2)10 connectors, 
   the atoms of which are not of our interest. We will only consider atoms of pTMAO, water, and ions.

2. This code will separate TMAO, water, ion-atoms from the SiO2 and connector atoms upon providing a reference index.
   The reference index is the begining of non-substrate/connector atoms assuming that the substrate/connector atoms 
   grouply indexed first, followed by TMAO/water/ion-atoms. For example, if the first 100 atoms in the trajectory file
   are substrate/connector atoms, the reference index would be 101. All atoms with index > 100 will be treated as
   surface/water/ion atoms.

3. We have to choose the distance (zdist = zmax - zmin)) in Z-direction (perpendicular distance to the surface)

4. The zval of the box dimesion will be initialized with 0. After determining zmin and zmax, zdist will set as zval.
   Alternatively, the user can also set customized zmin/zmax value to get user-defined zdist.

5. zmin = the lowest z of all TMAO/water atoms; zmax = the highest z of water atoms.

6. In the ion case, we may exclude the volumes occupied by the ion atoms.

7. The code will generate a plot as well as a 'surfPGR.dat' data file. 


## How to run

``` python surf_pRDF.py test.xyz ```


