from flip.element import Element
from flip.dofid import DofID
import numpy as np
import math
from flip import geometry

class PlaneStressTriLin2D(Element):
    def __init__(self, label, domain, nodes, props):
        Element.__init__(self,label, domain, nodes)
        self.E = props['E']
        self.nu = props['nu']
        self.t = props['t']
        c1 = domain.getNode(self.nodes[0]).coordinates
        c2 = domain.getNode(self.nodes[1]).coordinates
        c3 = domain.getNode(self.nodes[2]).coordinates
        self.A  = (1/2)*((c2[0]*c3[1]-c3[0]*c2[1])-(c1[0]*c3[1]-c3[0]*c1[1])+(c1[0]*c2[1]-c2[0]*c1[1]))
    def getNodeDofs (self):
        return (DofID.Dx, DofID.Dy)
    def computeT(self):
        pass
    def _computeB(self):
        c1 = self.domain.getNode(self.nodes[0]).coordinates
        c2 = self.domain.getNode(self.nodes[1]).coordinates
        c3 = self.domain.getNode(self.nodes[2]).coordinates

        y23 = c2[1]-c3[1]
        y31 = c3[1]-c1[1]
        y12 = c1[1]-c2[1]
        x32 = c3[0]-c2[0]
        x13 = c1[0]-c3[0]
        x21 = c2[0]-c1[0]

        be = (1/(2*self.A))*np.array([[ y23, 0, y31, 0, y12, 0],
                                [ 0, x32, 0, x13, 0, x21],
                                [x32,y23,x13,y31,x21,y12]])
        return be
    def _computeD(self):
        return (self.E/(1-self.nu*self.nu))*np.array([[ 1, self.nu, 0], [self.nu, 1, 0], [0, 0, (1-self.nu)/2 ]])
    def computeStiffness_local (self):
        b = self._computeB()
        d = self._computeD()
        return self.t*self.A*np.matmul(b.T, np.matmul(d, b))
    def computeStiffness(self):
        return self.computeStiffness_local()
    def computeEndForces(self, r):
        re = r[self.getLocationArray()]
        kl = self.computeStiffness()
        ans = np.matmul(kl, re)
        return ans
    def _computeIntVars(self, r):
        re = r[self.getLocationArray()]
        b = self._computeB()
        d = self._computeD()
        eps = np.matmul(b,re)
        sig = np.matmul(d, eps)
        return (eps, sig)
    def getEdgeNodes(self, iedge):
        if (iedge==0):
            return (self.nodes[0], self.nodes[1])
        elif (iedge==1):
            return (self.nodes[1], self.nodes[2])
        elif (iedge==2):
            return (self.nodes[2], self.nodes[0])
        else:
            raise IndexError("Element %d: edge %d out of range"%(self.label, iedge))

    def getEdgeLSC(self, iedge):
        return np.array([[self.cos, 0, self.sin],
                        [0        , 1, 0],
                        [-self.sin, 0, self.cos]])
    def getGeometry(self):
        return geometry.TriSurf(self.nodes, ((0,1,2),))
       
    def printState(self, context):
            (eps, sig)=self._computeIntVars(context.r)
            print (self.label, ":", "Eps=", eps, "Sig=", sig)

