import numpy as np
from   numpy import arange, sqrt
from   getDistWithPBC import pbcDistAB

class GridCenters(object):

    def __init__(self, dx, dy, dz, rvdW, zlow, zhigh, cellSize):
        self.dx       = dx
        self.dy       = dy
        self.dz       = dz
        self.zlow     = zlow
        self.zhigh    = zhigh
        self.rvdW     = rvdW
        self.cellSize = cellSize

    def get_grid_centers(self):
        xCnt          = arange(0.5 * self.dx, self.cellSize["a"], self.dx)
        yCnt          = arange(0.5 * self.dy, self.cellSize["b"], self.dy)
        zCnt          = arange((0.5 * self.dz) + self.zlow, self.zhigh + 0.1, self.dz)
        xv, yv, zv    = np.meshgrid(xCnt, yCnt, zCnt, sparse=False, indexing='ij')
        return(xv, yv, zv)

    def find_prox_dist(self, coordA, coordIONS):
        xGrid, yGrid, zGrid = self.get_grid_centers()
        xMax = len(xGrid)
        yMax = len(yGrid)
        zMax = len(zGrid)

        dmin  = []
        for i in range(0, 53):
            for j in range(0, 49):
                for k in range(0, 20):
                     centX, centY, centZ = (xGrid[i][j][k], yGrid[i][j][k], zGrid[i][j][k])
                     d      = []
                     count1, count2  = (0, 0)
                     for atom in range(len(coordA)): 
                         xyz = coordA[atom] 
                         sp, x, y, z = (xyz[0], xyz[1], xyz[2], xyz[3])
                         #include PBC-image atoms, which are outside the box
                         xPBC, yPBC, zPBC, diffPBC = pbcDistAB(centX, centY, centZ, x, y, z, self.cellSize)
                         d.append(diffPBC)
                         sp = sp.split('-') 
                         sp = sp[0]
                         #exclude the grid vol. occupied a solute-atom
                         if sp in self.rvdW['chemSP']  and diffPBC <= self.rvdW['radii'][self.rvdW['chemSP'].index(sp)]:
                             count1 += 1
                         #exclude the grid vol. occupied an ion-atom
                     if coordIONS != 0:
                         for ion in range(len(coordIONS)): 
                             xyzion = coordIONS[ion] 
                             spion, xion, yion, zion = (xyzion[0], xyzion[1], xyzion[2], xyzion[3])
                             diff = sqrt((centX - xion)**2 + (centY - yion)**2 + (centZ -zion)**2)
                             #xiPBC, yiPBC, ziPBC, diffiPBC = pbcDistAB(centX, centY, centZ, xion, yion, zion, self.cellSize)
                             #print(ion, xyzion)
                             if spion in self.rvdW['chemSP']  and diff <= self.rvdW['radii'][self.rvdW['chemSP'].index(spion)]:
                         #if spion in self.rvdW['chemSP']  and diffiPBC <= self.rvdW['radii'][self.rvdW['chemSP'].index(spion)]:
                                 #print(ion, centX, centY, centZ, diff)
                                 #print(centX, centY, centZ, diffiPBC)
                                 count2 += 1
                     if count1 != 0 or count2 != 0: 
                         continue 
                     else: 
                         #print(min(d))
                         dmin.append(min(d))
                          
        #print('WatOccGrids:', len(dmin))
        return(dmin)

