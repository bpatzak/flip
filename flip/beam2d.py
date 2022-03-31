from flip import element
from flip.dofid import DofID
import numpy as np
import math
from flip import geometry

class Beam2d(element.Element):
    def __init__(self, label, domain, nodes, props):
        element.Element.__init__(self,label, domain, nodes)
        self.EA = props['E']*props['A']
        self.EI = props['E']*props['I']
        n1c = domain.getNode(self.nodes[0]).coordinates
        n2c = domain.getNode(self.nodes[1]).coordinates
        self.l = math.sqrt((n2c[0]-n1c[0])**2 + (n2c[2]-n1c[2])**2)
        self.cos = (n2c[0]-n1c[0])/self.l
        self.sin = (n2c[2]-n1c[2])/self.l
    def getNodeDofs (self):
        return (DofID.Dx, DofID.Dz, DofID.Ry)
    def computeT(self):
        # globalToLocal fl=T.fg
        c = self.cos
        s = self.sin
        return np.array([[c,s,0,0,0,0],
                    [-s,c,0,0,0,0],
                    [0,0,1,0,0,0],
                    [0,0,0,c,s,0],
                    [0,0,0,-s,c,0],
                    [0,0,0,0,0,1]])
    def computeStiffness_local (self):
        l = self.l
        l2 = l**2
        l3 = l2*l
        EA=self.EA
        EI=self.EI
        return np.array([[EA/l, 0, 0, -EA/l, 0, 0],
                    [0, 12*EI/l3, -6*EI/l2, 0, -12*EI/l3, -6*EI/l2],
                    [0, -6*EI/l2, 4*EI/l,   0,  6*EI/l2, 2*EI/l],
                    [-EA/l, 0, 0, EA/l, 0, 0],
                    [0, -12*EI/l3, 6*EI/l2, 0, 12*EI/l3, 6*EI/l2],
                    [0, -6*EI/l2, 2*EI/l,   0, 6*EI/l2, 4*EI/l]])

    def computeStiffness(self):
        kl = self.computeStiffness_local()
        t=self.computeT()
        return np.matmul(np.matmul(np.transpose(t),kl), t)
    def computeEndForces(self, r):
        re = r[self.getLocationArray()]
        kl = self.computeStiffness_local()
        t  = self.computeT()
        ans = np.matmul(kl, np.matmul(t, re))
        #  loop over element loads
        loads=self.domain.getLoadsOnElement(self.label)
        for l in loads:
            le = l.getValue(self.domain, self.label, 1)
            ans=np.subtract(ans,le)
        return ans
    def getEdgeNodes(self, iedge):
        return self.nodes
    def getEdgeLSC(self, iedge):
        return np.array([[self.cos, 0, self.sin],
                        [0        , 1, 0],
                        [-self.sin, 0, self.cos]])
    def getGeometry(self, r=None, defScale=0.0):
        n1c = list(self.domain.getNode(self.nodes[0]).coordinates)
        n2c = list(self.domain.getNode(self.nodes[1]).coordinates)
        if (r is None):
            return geometry.PolyLine((n1c[0], n2c[0]), (n1c[1], n2c[1]), (n1c[2], n2c[2]), "ko-")
        else:
            nseg=10
            T = self.computeT()
            rl = np.matmul(T, r[self.getLocationArray()])
            ug=[]
            wg=[]
            for iseg in range(nseg+1):
                xl=iseg/nseg; # [0,1]
                # components from end displacements
                wl = (1.0-3.0*xl*xl+2.0*xl*xl*xl)*rl[1]+self.l*(-xl+2.0*xl*xl-xl*xl*xl)*rl[2]+(3.0*xl*xl-2.0*xl*xl*xl)*rl[4]+self.l*(xl*xl-xl*xl*xl)*rl[5]
                ul = (1.-xl)*rl[0]+xl*rl[3]
                # add contributions of loads
                #for (let load of eloads) {
                #    let c = load.computeBeamDeflectionContrib(xl);
                #    wl += c.w;
                #    ul += c.u;
                #}`
                # position along centerline
                xg = n1c[0]+(n2c[0]-n1c[0])/(nseg)*iseg
                zg = n1c[2]+(n2c[2]-n1c[2])/(nseg)*iseg
                ug.append(xg+(ul*self.cos-wl*self.sin)*defScale)
                wg.append(zg+(wl*self.cos+ul*self.sin)*defScale)
            return geometry.PolyLine(ug, [0]*(nseg+1), wg, "r")       
    def printState(self, context):
            print (self.label, ":", *self.computeEndForces(context.r))