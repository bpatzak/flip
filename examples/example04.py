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
eprops = {'E':E, 'A':A, 'I':I}
d = domain.Domain()
d.addNode(node.Node(label=1,coords=(0,0,0),bcs=(DofID.Dx, DofID.Dz, DofID.Ry)))
d.addNode(node.Node(label=2,coords=(0,0,4),bcs=()))
d.addNode(node.Node(label=3,coords=(4,0,4),bcs=()))
d.addNode(node.Node(label=4,coords=(4,0,0), bcs=(DofID.Dx, DofID.Dz, DofID.Ry)))
d.addElement(beam2d.Beam2d(label=1, domain=d, nodes=(1,2), props=eprops))
d.addElement(beam2d.Beam2d(label=2, domain=d, nodes=(2,3), props=eprops))
d.addElement(beam2d.Beam2d(label=3, domain=d, nodes=(3,4), props=eprops))

#d.addLoad(nodalload.NodalLoad(labels=(2,), values={DofID.Dz:10.0}))
d.addLoad(beam2dedgeload.Beam2dConstantEdgeLoad(elems=(2,), edges=(1,), values={DofID.Dz:-2.0, DofID.Dx:1.0}))
s = linearstatic.LinearStatic()
s.solve(d)

g = graphics.MatplotlibGraphis2D()
elemSet = list(range(1,len(d.elements)+1))
for e in elemSet:
    g.add(d.getElement(e).getGeometry())
    g.add(d.getElement(e).getGeometry(r=s.r, defScale=1.e6))
g.plot()

