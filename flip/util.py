from flip.dofid import DofID
from flip import domain
from flip import node
from flip import geometry

def genRectangularTriMesh(domain, lx, ly, nx, ny, elementClass, eprops):
    # generate nodes first
    label = 1
    dx = lx/nx
    dy = ly/ny
    for j in range(ny+1):
        for i in range(nx+1):
            domain.addNode(node.Node(label, (dx*i, dy*j,0), bcs={}))
            #print(domain.getNode(label))
            label=label+1
    #generate elements
    label = 1
    for j in range(ny):
        for i in range(nx):
            n1 = 1+ i + (nx+1)*j
            n2 = n1+1
            n3 = n1+nx+1
            n4=n3+1
            #print (label, n1,n2,n3,n4)
            domain.addElement(elementClass(label, domain, (n1, n2, n4), eprops ))
            label = label+1
            domain.addElement(elementClass(label, domain, (n1, n4, n3), eprops ))
            label=label+1

def makeTriSurf (domain, elementSet, displacements=None, defScale=0.0):
    #print(elementSet)
    # loop over eligble elements and add relavant nodes 
    nx = []
    ny = []
    nodeLabel2Indx={}
    index = 0
    cells = []
    for e in elementSet:
        elem = domain.getElement(e)
        eg = elem.getGeometry()
        eVertices = []
        if (isinstance(eg, geometry.TriSurf)):
            for node in eg.nodes: #loop over nodes of triangles
                if not node in nodeLabel2Indx:
                    # insert node
                    nodeLabel2Indx[node]=index
                    c=list(domain.getNode(node).coordinates)
                    if (not displacements is None):
                        if domain.hasNodeDof(node, DofID.Dx):
                            c[0]=c[0]+displacements[domain.getNodeDofEquation(node, DofID.Dx)]*defScale
                        if domain.hasNodeDof(node, DofID.Dy):
                            c[1]=c[1]+displacements[domain.getNodeDofEquation(node, DofID.Dy)]*defScale
                    nx.append(c[0])
                    ny.append(c[1])
                    index=index+1
                eVertices.append(nodeLabel2Indx[node])
            cells.append(eVertices)
    return (nx, ny, cells)

    # 