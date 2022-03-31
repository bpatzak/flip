from flip import linearstatic
from flip import domain
from flip import node
from flip import planestresstrilin2d
from flip import nodalload
from flip.dofid import DofID
import numpy as np
from flip import graphics

E = 10
nu=0.1
t=0.1
d = domain.Domain()
d.addNode(node.Node(label=1,coords=(0,0,0),bcs=(DofID.Dx, DofID.Dy)))
d.addNode(node.Node(label=2,coords=(6,0,0), bcs=(DofID.Dy,)))
d.addNode(node.Node(label=3,coords=(6,3,0), bcs=(DofID.Dy,)))
d.addNode(node.Node(label=4,coords=(0,3,0), bcs=(DofID.Dx, DofID.Dy)))
d.addNode(node.Node(label=5,coords=(2,1,0), bcs={}))
eprops= {'E':E, 'nu':nu, 't':t}
d.addElement(planestresstrilin2d.PlaneStressTriLin2D(label=1, domain=d, nodes=(1,2,5), props=eprops))
d.addElement(planestresstrilin2d.PlaneStressTriLin2D(label=2, domain=d, nodes=(2,3,5), props=eprops))
d.addElement(planestresstrilin2d.PlaneStressTriLin2D(label=3, domain=d, nodes=(3,4,5), props=eprops))
d.addElement(planestresstrilin2d.PlaneStressTriLin2D(label=4, domain=d, nodes=(4,1,5), props=eprops))
d.addLoad(nodalload.NodalLoad(labels=(2,3), values={DofID.Dx:1.5}))
s = linearstatic.LinearStatic()
s.solve(d)

