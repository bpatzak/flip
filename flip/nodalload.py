from flip import load

class NodalLoad(load.Load):
    def __init__(self, labels, values):
        load.Load.__init__(self)
        self.nodalLabels = labels #list
        self.values = values #dict, keys dofIDs values intensities
    def getLocationArrays(self, domain):
        answer=[]
        for l in self.nodalLabels: # loop over nodes
            iloc = []
            for dof, val in self.values.items():
                if domain.hasNodeDof(l, dof):
                    iloc.append(domain.equationNumbering[l][dof])
            answer.append(iloc)
        return answer
    def getValues(self, domain):
        answer=[]
        for l in self.nodalLabels: # loop over nodes
            lval = []
            for dof, val in self.values.items():
                if domain.hasNodeDof(l, dof):
                    lval.append(val)
            answer.append(lval)
        return answer
