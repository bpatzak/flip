class Element:
    def __init__(self, label, domain, nodes, props={}):
        self.domain = domain
        self.label=label
        self.nodes=nodes
    def getNodeDofs (self):
        return ()
    def getLocationArray(self):
        loc = []
        for n in self.nodes:
            loc.extend(self.domain.getNodeEquations(n, self.getNodeDofs()))
        return loc
    def computeT(self):
        pass
    def computeStiffness(self):
        pass
    def computeEndForces(self, context):
        pass
    def printState(self, context):
        pass


    def getEdgeNodes(self, iedge):
        pass
    def getEdgeLSC(self, iedge):
        pass

    def getGeometry(self):
        pass
    