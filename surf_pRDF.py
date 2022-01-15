#!/usr/bin/env python
import sys
import numpy as np
from   numpy import arange, histogram, zeros, pi, where, sqrt
import matplotlib.pyplot as plt
import grid3D as gd3
from   itertools import chain
from   projectRDF import RDF

############### SYSTEM INFO  ###################
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

############### CALCULATION SETUP ###################
#rdfType    =  'gij'  #partial rdf
rdfType    =  'pgij'  #proximal rdf
typei      =  'O'    #modify this according to your need (only for rdfType = 'gij') 
typej      =  'O'
dr         =   0.2     #bin siz
#assuming that substrate atoms are placed first in the trajectory file
# we will be ignoring those atoms, as those do contribute to the hydration
#reference atom to separate substrate atoms from the brush/water/ion atoms
#refIndex = index to the reference atom; atoms with indices => refIndex are of our interest 
refIndex   =  33600 #tmao_5
#initialize box size and rmax; these will be updated later
rMax       =   0        #8.1      #maximum radial distance(rcut-off)
xa, yb, zc =  (85.86, 79.38, 0)
cell_size =  {
    "a" : xa,
    "b" : yb,
    "c" : zc
}

############### SIMULATION DATA PROCESSING ###################
#get the trajectory file from command line
try:
    f = sys.argv[1]
except:
    print('Error: wrong input file on the command line')
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
    print('Wrong input file. Check your system and input file.')
    exit(1)
else:
    traj     = list(divide_chunks(trajFile, totAtom + 2))
    frames   = int((len(trajFile) / (totAtom + 2)))

############### ION COORDINATES ###################
#get ions
line  =  traj[0][2:totAtom+2]
data  =  [i.split() for i in line]
sp    =  [j[0] for j in data]
#get xyz for ions if there is any
for name in range(len(ions['type'])):
    if ions['type'][name] in sp:
        xyzions =  rdf.get_coord_ions(traj, frames, totAtom) 
    else:
        xyzions = zeros(frames)
#print('ions:', xyzions[0])

############### POLYMER/WATER COORDINATES ###################
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
print('LOW', surfLOW, 'HIGH', surfHIGH, 'INTF', polyMAX)

############### SET UPPER AND LOWER LIMIT OF Z HERE #############################
############### Z-DISTANCE AND 3D-GRID SETUP ####################################
zlow    = 36   #or surfLOW
zhigh   = 68.4 #surfHIGH
zdist   = abs(zhigh-zlow) #surfHIGH - surfLOW
#update box-size and rmax
cell_size['c'] = zdist
rMax    = 0.5*zdist
#choose grid-size: dx/dy/dz ~vdW radii 
#and ensure boxSize_a / dx = integer, boxSize_b / dy = integer, and zdist / dz = integer
#dx, dy, dz = (0.81, 0.81, 0.81)
dx, dy, dz   =  (1.62, 1.62, 1.62)
#dx, dy, dz =  (2, 2, 2)
gridVol      = dx*dy*dz
grid         =  gd3.GridCenters(dx, dy, dz, rvdW, zlow, zhigh, cell_size)
#################################################################################
#################################################################################

############### ATOMS WITHIN THE SURFACE CUTOFF (Z-DISTANCE)  ###################
#selected atoms above zmin and below zmax
#set user defined zmin/zmax or zmin/zmax from each snapshot
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
#no ion consideration
#no exclusion of ion-volume
xyzions_none = zeros(frames)
rdf          = RDF(refIndex, ions, rdfType, dr, rMax, cell_size)
pd, pnt, pgd =  rdf.get_prox_rij(polys, wats, xyzions_none, grid)
# exclusion of ion-volume
#pd, pnt, pgd    =  rdf.get_prox_rij(polys, wats, xyzions, grid)
pgr, uoccgrd =  rdf.get_timeAvg_pgij(pd, pgd, frames, len(wats[0]))
#print(pgr.tolist())
#print(uoccgrd.tolist())
pgr_norm, r  =  rdf.get_norm_pgij(pgr, uoccgrd, gridVol)
#print(pgr_norm)
#print(r)

######### WRITE PG(R) DATA  #################################
with open('surfPGR.dat', 'w') as f:
    for i  in range(len(r)):
        f.write(str(r[i])+ "\t"+str(pgr_norm[i])) #+"\t"+str(n_r[i]))
        f.write("\n")
f.close()

######### VISUALIZE  #################################
plt.figure()
plt.rcParams.update({'font.size': 22})
plt.plot(r, pgr_norm, color='black')
plt.xlabel('r ('r'$\AA$'')')
plt.ylabel('pg(r)')
plt.xlim((0, rMax))
plt.ylim((0, 1.05 * pgr_norm.max()))
plt.show()
