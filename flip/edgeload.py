import math
import numpy as np
from flip.dofid import DofID 
from flip import load

class EdgeLoad(load.Load):
    def __init__(self, elems, edges, values):
        self.elementLabels = elems
        self.edges = edges
        self.values = values # dofID dict; global
    def _getLocationArray(self, elem, edge, domain):
        edgeNodes = domain.getElement(elem).getEdgeNodes(edge)
        dofs=domain.getElement(elem).getNodeDofs()
        loc=[]
        if edgeNodes:
            for n in edgeNodes:
                    for d in dofs:
                        loc.append(domain.equationNumbering[n][d])
        return loc
           
    def _getLocalIntensities(self, domain, element, edge, dofids):
        gl = [] # global values
        for i in dofids:
            gl.append(self.values.get(i, 0.0))
        # transform
        #print("gl:", gl)
        t = self.getT(dofids, domain.getElement(element).getEdgeLSC(edge))
        #print("t:",t)
        return np.matmul(t, gl)    

    def getT(self, dofids, lcs):
        # returns transformation matrix from global to local
        size = len(dofids)
        ans = np.zeros([size,size])
        i = j = 0
        #print("elcs:", lcs)
        for di in dofids:
            j=0
            for dj in dofids:
                if (di in (DofID.Dx, DofID.Dy, DofID.Dz)):
                    if (dj in (DofID.Dx, DofID.Dy, DofID.Dz)):
                        ans[i,j]=lcs[di-DofID.Dx][dj-DofID.Dx]
                elif (di in (DofID.Rx, DofID.Ry, DofID.Rz)):
                    if (dj in (DofID.Rx, DofID.Ry, DofID.Rz)):
                        ans[i,j]=lcs[di-DofID.Rx][dj-DofID.Rx]
                elif (di==dj):
                    ans[i,j]=1.0
                j = j+1
            i = i+1
        return ans

    def getLocationArrays(self, domain):
        loc = []
        for ielem, iedge in zip(self.elementLabels, self.edges):
            loc.append(self._getLocationArray(ielem, iedge, domain))
        return loc

    def getValue(self, domain, element, iedge):
        # return single element, edge contribution
        pass

    def getValues(self, domain):
        pass
        # aasemble load vector of the receiver
        ans = []
        for ielem, iedge in zip(self.elementLabels, self.edges):
            ival = self.getValue(domain, ielem, iedge)
            ans.append(ival)
        return ans
