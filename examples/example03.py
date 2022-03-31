from flip import linearstatic
from flip import domain
from flip import node
from flip import planestresstrilin2d
from flip import nodalload
from flip.dofid import DofID
import numpy as np
from flip import graphics
from flip import util

E = 10
nu=0.1
t=0.1
d = domain.Domain()
# generate nodes
ny = 6
nx = 20
eprops = {'E':E, 'nu':nu, 't':t}
util.genRectangularTriMesh(d, 5, 1, nx, ny, planestresstrilin2d.PlaneStressTriLin2D, eprops)
# clamp vertical edge
for j in range(ny+1):
    d.getNode(j*(nx+1)+1).bcs=(DofID.Dx, DofID.Dy)
#add some load
d.addLoad(nodalload.NodalLoad(labels=(range(nx+1,(ny+1)*(nx+1)+1, nx+1)), values={DofID.Dy:-1.0}))
s = linearstatic.LinearStatic()
s.solve(d)

g = graphics.MatplotlibGraphis2D()
elemSet = list(range(1,len(d.elements)+1))
(nx, ny, cells) = util.makeTriSurf(d, elemSet, displacements=s.r, defScale=0.00005)
vals=[]
for e in elemSet:
    (eps,sig)=d.getElement(e)._computeIntVars(s.r)
    vals.append(sig[0])

g.addTriSurf(nx,ny,cells, cellVars = np.array(vals))
g.addTriSurf(nx,ny,cells)
g.plot()
