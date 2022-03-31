class Domain:
    def __init__(self):
        self.nodes = {}
        self.elements = {}
        self.loads = []
        self.equationNumbering = {} #key node label, value is dofID dictionary

    def getNode(self, label):
        return self.nodes.get(label)
    def getElement(self, label):
        return self.elements.get(label)
    def addNode(self, node):
        self.nodes[node.label]=node
    def addElement(self, element):
        element.domain = self
        self.elements[element.label]=element
    def addLoad(self, load):
        self.loads.append(load)
    def getNodeEquations(self, node, dofIDs):
        answer=[]
        for d in dofIDs:
            answer.append(self.equationNumbering[node][d])
        return answer
    def hasNodeDof(self, node, dof):
        return (dof in self.equationNumbering[node])
    def getNodeDofEquation(self, node, dof):
        if dof in self.equationNumbering[node]:
            return self.equationNumbering[node][dof]
        else:
            return None
    def getLoadsOnElement(self, elem):
        ans=[]
        for l in self.loads:
            if l.actsOnElement(elem):
                ans.append(l)
        return ans