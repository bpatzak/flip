from itertools import count
import math
import numpy as np
from flip.dofid import DofID 
from flip import edgeload

class Beam2dConstantEdgeLoad(edgeload.EdgeLoad):
    def __init__(self, elems, edges, values):
        edgeload.EdgeLoad.__init__(self, elems, edges, values)

    def getValue(self, domain, elem, edge):
        e = domain.getElement(elem)
        n1c = domain.getNode(e.nodes[0]).coordinates
        n2c = domain.getNode(e.nodes[1]).coordinates
        l = math.sqrt((n2c[0]-n1c[0])**2 + (n2c[2]-n1c[2])**2)
        lv = self._getLocalIntensities(domain, elem, edge, e.getNodeDofs())
        #print(lv)
        t = self.getT(domain.getElement(elem).getNodeDofs(), domain.getElement(elem).getEdgeLSC(edge))
       
        fl= np.array([0.5*l*lv[0], 0.5*l*lv[1], -lv[1]*l*l/12., 0.5*l*lv[0], 0.5*l*lv[1], lv[1]*l*l/12.])
        # transform to global cs
        return np.concatenate((np.matmul(t.T, fl[0:3]), np.matmul(t.T, fl[3:6])))


        

