#!/usr/bin/env python
# A python script to calculate surface hydration (pG(r)) written by Pranab Sarker
# Contact: srpranab@gmail.com
# See Figure 4a in Huang; Zhang et al. J. Am. Chem. Soc. 2021, 143, 40, 16786–16795
# DOI: https://doi.org/10.1021/jacs.1c08280

import sys
import numpy as np
from   numpy import arange, histogram, zeros, pi, where, sqrt
import matplotlib.pyplot as plt
import grid3D as gd3
from   itertools import chain
from   projectRDF import RDF

############### CALCULATION SETUP BEGIN ######################
#total atoms in the system 
tot  =  81138 
#vdW radii
rvdW     =  {
    'chemSP' : ['C', 'H', 'N', 'O', 'Na', 'Cl', 'Si'],
    'radii'  : [1.7, 1.2, 1.55, 1.52, 2.27, 1.75, 2.1]
}
#ions
ions     =  {
    'type' : ['Na', 'Cl'],
}
rdfType    =  'pgij'   # proximal rdf
typei      =  'solute' # solute molecule 
typej      =  'O'      # solvent atom, O -----> H2O
dr         =   0.2     # bin siz
#Assuming that substrate atoms are placed first in the trajectory file
#we will be ignoring those atoms, as those do not contribute to the hydration
#Reference atom must be provided to separate substrate+connector atoms from the brush/water/ion atoms
#refIndex = index of the reference atom (e.g., the last surface-atom); atoms with indices => refIndex are of our interest 
refIndex   =  33600 #tmao_5

#set box size and rmax here or set 'AUTO' to let the code automatically determine those
#set the maximum radial distance(rcut-off) to calculate RDF as a function of r
rMax_user  = 15   
#rMax_user = 'AUTO'   
#set the lower and upper cutoffs of zdist
zmin_user  =  36   #lower cutoff 
#zmin_user  = 'AUTO'   #lower cutoff 
zmax_user  = 68.4 #upper cutoff 
#zmax_user  = 'AUTO' #upper cutoff 
#_____________________________________________________________________________________
if rMax_user == 'AUTO':
    rMax = 0 # initialize now and update later after calculating zdist; rMax = 0.5*zdist
             # assuming that zdist => 30 angs.; if zdist < 30 angs., set rMax => 15 angs
             # this code will use bulk density of water within r = 10 and r = 15 angs
else:
    rMax = rMax_user

if zmin_user == 'AUTO' or zmax_user == 'AUTO':
    zdist_user = 0      #initialize now and update later after calculating zmin/zmax
else:
    zdist_user = zmax_user - zmin_user
#_____________________________________________________________________________________
#set cell size
xa, yb, zc =  (85.86, 79.38, zdist_user)
cell_size  =  {
    "a" : xa,
    "b" : yb,
    "c" : zc
}

#choose grid-size: dx/dy/dz ~vdW radii 
dx, dy, dz   =  (1.62, 1.62, 1.62)
#dx, dy, dz =  (2, 2, 2)
gridVol      = dx*dy*dz

#ion-case
#choose how the ions will be treated: visible or invisible?
#if visible,  the occupied volumes of ions will be considered
#if invisible,  the occupied volumes of ions will be ignored
#when ions are invisible
ionExcVol    = 'NO' #default
#when ions are visible
#ionExcVol    = 'YES'
############### CALCULATION SETUP END ########################

############### SIMULATION DATA PROCESSING ###################
#get the trajectory file from command line
try:
    f = sys.argv[1]
except:
    print('Error: no input file provided')
    sys.exit(1)
#load trajectory file
def divide_chunks(l, n):
     # looping till length l  
     for i in range(0, len(l), n):
          yield l[i:i + n]

trajFile = open(f).read().splitlines()
totAtom  = int(trajFile[0])
print('total atoms:', totAtom)
#check the input file
if tot  != totAtom:
    print('Wrong input file provided. The total atoms do not match')
    exit(1)
else:
    traj     = list(divide_chunks(trajFile, totAtom + 2))
    frames   = int((len(trajFile) / (totAtom + 2)))

#get ions coordinates
line  =  traj[0][2:totAtom+2]
data  =  [i.split() for i in line]
sp    =  [j[0] for j in data]
#get xyz for ions if there is any
crd             = RDF(refIndex, ions, rdfType, dr, rMax, cell_size)
for name in range(len(ions['type'])):
    if ions['type'][name] in sp:
        xyzions =  crd.get_coord_ions(traj, frames, totAtom) 
    else:
        xyzions = zeros(frames)
#if we want to ignore ions
xyzions_none = zeros(frames)

# separte atoms of polymer brush and water
sep             = RDF(refIndex, ions, rdfType, dr, rMax, cell_size)
poly, wat       = sep.get_coord_ij(traj, frames, totAtom)

############### SURFACE LOW/HIGH-CUTOFF FOR CALCULATION  ###################
#all_atom        = [0]*frames
zmax            = []
zwmax           = []
zmin            = []
for fm in range(frames):
    #all_atom[fm] = poly[fm] + wat[fm]
    zpval = []
    zwval = []
    for i in range(len(poly[fm])):
        zp   = poly[fm][i][3]       #poly-Z-coordinate
        zpval.append(zp)          
        sMAX = max(zpval)
    print(fm, 'surfMax', sMAX)     # max. Z
    for j in range(len(wat[fm])):
        zw   = wat[fm][j][3]       #wat-Z-coordinate
        zwval.append(zw)         
        wMIN = min(zwval)         
    print(fm, 'watMin', wMIN) 
    zmax.append(max(zpval))    
    zwmax.append(max(zwval)) 
    zmin.append(min(zwval)) 
#print('zmin:', zmin, 'zmax:', zmax[0])
surfLOW      = int(min(zmin))
surfHIGH     = int(max(zwmax))
polyMAX      = int(max(zmax))
print('zlow@', surfLOW, 'zhigh@', surfHIGH, 'interface@', polyMAX, 'all dimensions in angstrom')

############### UPDATE UPPER AND LOWER LIMIT OF Z  #############################
#update box-size and rmax if zmim/zmax are not user-defined
if zmin_user == 'AUTO' and  zmax_user != 'AUTO':
    zlow    = zmin_user
    zhigh   = surfHIGH
    zdist   = zhigh - zlow
    cell_size['c'] = zdist
    rMax    = 0.5*zdist
elif zmin_user != 'AUTO' and  zmax_user == 'AUTO':
    zlow    = surfLOW
    zhigh   = zmax_user
    zdist   = zhigh - zlow
    cell_size['c'] = zdist
    rMax    = 0.5*zdist
elif zmin_user == 'AUTO' and  zmax_user == 'AUTO':
    zlow    = surfLOW
    zhigh   = surfHIGH 
    zdist   = zhigh - zlow
    cell_size['c'] = zdist
    rMax    = 0.5*zdist
else:
    zlow    = zmin_user
    zhigh   = zmax_user

#ensure boxSize_a / dx = integer, boxSize_b / dy = integer, and zdist / dz = integer
rmdx1 = round(cell_size['a'] / dx, 2)
rmdx2 = round(cell_size['a'] / dx, 0)
diffx = rmdx1 - rmdx2
rmdy1 = round(cell_size['b'] / dy, 2)
rmdy2 = round(cell_size['b'] / dy, 0)
diffy  = rmdy1 - rmdy2
rmdz1 = round(cell_size['c'] / dz, 2)
rmdz2 = round(cell_size['c'] / dz, 0)
diffz = rmdz1 - rmdz2

#update cellsize to have integer grids in x, y, and z directions
cell_size_grid  =  {
    "a" : xa-diffx,
    "b" : yb-diffy,
    "c" : zc-diffz
}

############### ATOMS WITHIN THE SURFACE CUTOFF (Z-DISTANCE)  ###################
#selected atoms above zmin and below zmax
#set user-defined zmin/zmax or zmin/zmax from each snapshot
polys = []
wats  = []
IONs  = []
surf  = []
for fm in range(frames):
    zth = zlow 
    ps  = []
    ws  = []
    ins = []
    for i in range(len(poly[fm])):
        zp = poly[fm][i][3]
        if zp >= zth:
            ps.append(poly[fm][i])
    polys.append(ps)
    for j in range(len(wat[fm])):
        zw = wat[fm][j][3]
        sp = wat[fm][j][0]
        sp = sp.split('-')
        sp = sp[0]
        if zw <= zhigh  and sp == 'O':
            ws.append(wat[fm][j])
    wats.append(ws)

######### PG(R) CALCULATION  #################################
#grid set-up
grid         =  gd3.GridCenters(dx, dy, dz, rvdW, zlow, zhigh, cell_size_grid)
rdf          = RDF(refIndex, ions, rdfType, dr, rMax, cell_size)
if ionExcVol == 'NO':
    pd, pnt, pgd =  rdf.get_prox_rij(polys, wats, xyzions_none, grid)
#exclusion of ion-volume
elif ionExcVol == 'YES':
    pd, pnt, pgd    =  rdf.get_prox_rij(polys, wats, xyzions, grid)
pgr, uoccgrd =  rdf.get_timeAvg_pgij(pd, pgd, frames, len(wats[0]))
#print(pgr.tolist())
#print(uoccgrd.tolist())
pgr_norm, r  =  rdf.get_norm_pgij(pgr, uoccgrd, gridVol)
#print(pgr_norm)
#print(r)

######### WRITE PG(R) DATA  #################################
with open('surfPGR.dat', 'w') as f:
    f.write(str('#r') + "\t" + str('pG(r)'))
    f.write("\n")
    for i  in range(len(r)):
        f.write(str(r[i])+ "\t"+str(pgr_norm[i])) 
        f.write("\n")
f.close()

######### VISUALIZE  PG(R) ##################################
plt.figure()
plt.rcParams.update({'font.size': 22})
plt.plot(r, pgr_norm, color='black')
plt.xlabel('r ('r'$\AA$'')')
plt.ylabel('pg(r)')
plt.xlim((0, rMax))
plt.ylim((0, 1.05 * pgr_norm.max()))
plt.show()
