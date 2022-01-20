import math
import numpy as np
from   numpy import arange, histogram, zeros, pi, where, sqrt
from   getDistWithPBC import pbcDistAB
from   getMinValIndex import minValIndex

class RDF(object):
    
    def __init__(self, rfIdx, ions, rdfType, dr, rMax, cellSize):
        self.rfIdx       =  rfIdx
        self.ions        =  ions
        self.rdfType     =  rdfType
        self.dr          =  dr
        self.rMax        =  rMax
        print('ZMAX', self.rMax)
        self.binEdges    =  arange(0., self.rMax + 0.5*self.dr, self.dr)
        self.numOf_dr    =  len(self.binEdges) - 1
        self.cellSize    =  cellSize
        self.dFrame      =  [[]]

    def get_coord_ij(self, traj, frames, totAtom):
       #list for ref and target atoms' coordinates
        atomi, xatomi, yatomi, zatomi, xyzi   = ([], [], [], [], [])
        atomj, xatomj, yatomj, zatomj, xyzj   = ([], [], [], [], [])
        for i in range(0, frames):
            spi, xi, yi, zi        = ([], [], [], [])
            spj, xj, yj, zj        = ([], [], [], [])
            for j in range(totAtom):
                xyz = traj[i][j+2].split()
                if j+1 <= self.rfIdx:
                    spi.append(xyz[0] + str('-') + str(j+1))
                    xi.append(float(xyz[1]))
                    yi.append(float(xyz[2]))
                    zi.append(float(xyz[3]))
                elif j+1 > self.rfIdx and xyz[0] not in self.ions['type']:
                    spj.append(xyz[0]+ str('-') + str(j+1))
                    xj.append(float(xyz[1]))
                    yj.append(float(xyz[2]))
                    zj.append(float(xyz[3]))

            atomi.append(spi)
            xatomi.append(xi)
            yatomi.append(yi)
            zatomi.append(zi)   
            zippedi = list(zip(atomi[i], xatomi[i], yatomi[i], zatomi[i]))
            xyzi.append(zippedi)

            atomj.append(spj)
            xatomj.append(xj)
            yatomj.append(yj)
            zatomj.append(zj)
            zippedj = list(zip(atomj[i], xatomj[i], yatomj[i], zatomj[i]))
            xyzj.append(zippedj)

        #print(xyzj[0]) 
        return(xyzi, xyzj) 

    def get_coord_ions(self, traj, frames, totAtom):
        ion, xion, yion, zion, xyzion   = ([], [], [], [], [])
        for i in range(0, frames):
            atom, x, y, z  = ([], [], [], [])
            for j in range(totAtom):
                xyz = traj[i][j+2].split()
                if xyz[0] in self.ions['type']:
                    atom.append(xyz[0])
                    x.append(float(xyz[1]))
                    y.append(float(xyz[2]))
                    z.append(float(xyz[3]))

            ion.append(atom)
            xion.append(x)
            yion.append(y)
            zion.append(z)   
            zippedion = list(zip(ion[i], xion[i], yion[i], zion[i]))
            xyzion.append(zippedion)

        #print(xyzion) 
        return(xyzion) 

    def get_rij(self, coordi, coordj):
        dist, pbcXYZ = ([], [])
        for frame in  range(len(coordi)):
            xi, yi, zi  = (0, 0, 0)
            dq, drq     = ([], [])
            spij, newxij, newyij, newzij, xyzij = ([], [], [], [], [])
            atomi = coordi[frame][11600][0]
            xi    = coordi[frame][11600][1]
            yi    = coordi[frame][11600][2]
            zi    = coordi[frame][11600][3]
            spj, newxj, newyj, newzj = ([], [], [], [])
            for q in range(len(coordj[frame])):
                atomj = coordj[frame][q][0]
                xj    = coordj[frame][q][1]
                yj    = coordj[frame][q][2]
                zj    = coordj[frame][q][3]
                #print('before', xj, yj, zj)
                xjPBC, yjPBC, zjPBC = (0, 0, 0) #, rij = pbcDistAB(xi, yi, zi, xj, yj, zj, self.cellSize)
                #print(atomj, rij)
                rij = sqrt((0-zj)**2)
                #print('after', xjPBC, yjPBC, zjPBC)
                #print(atomj, rij)
                spj.append(atomj)
                newxj.append(xjPBC)
                newyj.append(yjPBC)
                newzj.append(zjPBC)
                dq.append(rij)
            spij.append(spj)
            newxij.append(newxj)
            newyij.append(newyj)
            newzij.append(newzj)
            #pbcZip = list(zip(spij[frame], newxij[frame], newyij[frame], newzij[frame]))
            #pbcXYZ.append(pbcZip)
            #add all r for each frame    
            dist.extend(dq)
        #print(dist)
        #return(dist, pbcXYZ)
        return(dist)

    def get_prox_rij(self, coordi, coordj, coordIONS, grid):
        pdist, pgridDist, patoms, pmax, pgmax, pcount  = ([], [], [], [], [], [])
        for frame in  range(len(coordi)):
            xi, yi, zi           = (0, 0, 0)
            surfAtoms, dpq, pind = ([], [], [])
            #solvent
            print('WATER:', (len(coordj[frame])))
            print('TMAO:', (len(coordi[frame])))
            for p in range(len(coordj[frame])):
                atomj = coordj[frame][p][0]
                xj    = coordj[frame][p][1]
                yj    = coordj[frame][p][2]
                zj    = coordj[frame][p][3]
                r     = []
                #solute
                for q in range(len(coordi[frame])):
                    atomi = coordi[frame][q][0]
                    xi    = coordi[frame][q][1]
                    yi    = coordi[frame][q][2]
                    zi    = coordi[frame][q][3]
                    #impose PBCs on Target atoms
                    xjPBC, yjPBC, zjPBC, rij = pbcDistAB(xi, yi, zi, xj, yj, zj, self.cellSize)
                    r.append(float(rij))
                    #print(p, atomj, min(r))
                #get proximal distance for a solvent atom and the proximal solute atom
                dpq.append(min(r))
                #print(len(dpq))
                min_atomi = minValIndex(r)
                pind.append(min_atomi)
                surfAtoms.append(coordi[frame][min_atomi][0])
            pdist.append(dpq)
            #print(pdist[0])
            patoms.append(surfAtoms)
            pcount.append(len(np.unique(pind)))
            #get proximal grid distance (solute/ion-volume exclueded)
            pgridDist.append(grid.find_prox_dist(coordi[frame], coordIONS[frame]))
            pmax.append(max(pdist[frame]))
            pgmax.append(max(pgridDist[frame]))
            #update max. radial distance
        #if self.rdfType == 'pgij':
        #    self.rMax = min(min(pmax), min(pgmax))

        #check whethter grid is small enough
        for fm in range(len(pdist)):
           if min(pdist[fm]) < min(pgridDist[fm]):
               print('frame :', fm+1," ", 'prox_min = ', min(pdist[fm])," ", 'proxGrid_min = ', min(pgridDist[fm]))
               print('Larger grid size. Make it smaller. Otherwise, you might have an error in calculating pg(r)')
           #if max(pdist[fm]) > max(pgridDist[fm]):
           #    print(fm," ", max(pdist[fm])," ", max(pgridDist[fm]))

        return(pdist, pcount, pgridDist)

    def get_timeAvg_dist(self, dist, frames):
        g = zeros(self.numOf_dr)
        for frame in  range(0, frames):                
            (natom, bins) = histogram(dist[frame], bins=self.binEdges, normed=False)
            gTot   = np.add(g, natom)
            g      = gTot

        #time averqge    
        gAvg = gTot / frames
        return(gAvg)

    def get_timeAvg_pgij(self, pdist, pgridDist, frames, atomsj):
        self.binEdges = arange(0., self.rMax + 0.5 * self.dr, self.dr)
        self.numOf_dr = len(self.binEdges) - 1
        #print(self.binEdges)
        g = zeros(self.numOf_dr)
        n = zeros(self.numOf_dr)
        radii = zeros(self.numOf_dr)
        #constant density at r
        vol = round((self.cellSize["a"]* self.cellSize["b"]* self.cellSize["c"]),3)
        for frame in  range(0, frames):                
            (natom, bins) = histogram(pdist[frame], bins=self.binEdges, normed=False)
            (ngrid, bins) = histogram(pgridDist[frame], bins=self.binEdges, normed=False)
            #print('FRAME#',frame) 
            for i in range(self.numOf_dr):
                radii[i] = round((self.binEdges[i] + self.binEdges[i+1]) / 2.,3)
                #print('r:', radii[i], 'CN:', natom[i]) 
            #print(natom)
            #prho  = round((pcount[frame]*atomsj / vol),3)
            prho  = round((atomsj / vol),3)
            #gNorm = natom / prho
            gNorm = natom
            gTot  = np.add(g, gNorm)
            g     = gTot
            nTot  = np.add(n, ngrid)
            n     = nTot

        #time average    
        gAvg = gTot / frames
        nAvg = nTot / frames
        return(gAvg, nAvg)

    def get_norm_pgij(self, gAvg, nAvg, drVol):
        pgij     = zeros(self.numOf_dr)
        radii    = zeros(self.numOf_dr)
        bulkj    = []
        for i in range(self.numOf_dr):
            radii[i] = round((self.binEdges[i] + self.binEdges[i+1]) / 2.,3)
            voli     = nAvg[i] * drVol
            if gAvg[i] == 0:
               pgij[i] =  0
            else:
                pgij[i] = round((gAvg[i]) / voli, 3)
                print(radii[i], pgij[i])

            if radii[i] >= 0.7*self.rMax and radii[i] <= self.rMax:
                bulkj.append(pgij[i])

        tot = 0
        for val in bulkj:
            tot += val
        totAvg = tot / len(bulkj)
        pgij = pgij / totAvg

        return(pgij, radii)
