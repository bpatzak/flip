from flip import linearstatic
from flip import domain
from flip import node
from flip import beam2d
from flip import nodalload
from flip import beam2dedgeload
from flip.dofid import DofID
import numpy as np
from flip import graphics

E = 30e9
A=0.01
I=0.001
d = domain.Domain()
d.addNode(node.Node(label=1,coords=(0,0,0),bcs=(DofID.Dx, DofID.Dz, DofID.Ry)))
d.addNode(node.Node(label=2,coords=(3,0,0), bcs={}))
d.addElement(beam2d.Beam2d(label=1, domain=d, nodes=(1,2), props={'E':E, 'A':A, 'I':I}))
#d.addLoad(nodalload.NodalLoad(labels=(2,), values={DofID.Dz:10.0}))
d.addLoad(beam2dedgeload.Beam2dConstantEdgeLoad(elems=(1,), edges=(1,), values={DofID.Dz:1.0, DofID.Dx:3.0}))
s = linearstatic.LinearStatic()
s.solve(d)

g = graphics.MatplotlibGraphis2D()
g.add(d.getElement(1).getGeometry())
g.add(d.getElement(1).getGeometry(r=s.r, defScale=1.e6))
g.plot()

