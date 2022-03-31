class Solver:
    def __init__(self):
        print("OOFEM.py (c) 2022 Borek Patzak")
        pass
    def numberEquations(self, domain):
        # fist loop over elements and collect dofID requirements in nodes
        domain.equationNumbering = {}
        nodalDofs = {}
        for n in domain.nodes:
            domain.equationNumbering[n]={} # empty dict
            nodalDofs[n] = set() # empty set
        for elabel, elem in domain.elements.items():
            dofs = elem.getNodeDofs()
            for n in elem.nodes:
                if (n in nodalDofs):
                    nodalDofs[n].update(dofs)
                    #for d in dofs:
                    #    nodalDofs[n].update({d:None})
                else:
                    raise KeyError ("Element node does not exist in domain")
        # count equations
        domain.neq=0
        domain.pneq = 0
        for nlabel,node in domain.nodes.items():
            for d in nodalDofs[nlabel]:
                if (d in node.bcs):
                    domain.pneq = domain.pneq+1
                else:
                    domain.neq = domain.neq+1
        # assign equtions to dofs
        un = 0
        pn = domain.neq
        for nlabel, node in domain.nodes.items():
            for d in nodalDofs[nlabel]:
                if (d in node.bcs):
                    domain.equationNumbering[nlabel][d]=pn
                    pn=pn+1
                else:
                    domain.equationNumbering[nlabel][d]=un
                    un=un+1   
        print ("Solver: neq=", domain.neq, ", pneq=", domain.pneq)             
    def assembleMatrix(self, k,loc,ke):
        for i, index in enumerate(loc):
            k[index, loc] += ke[i, :]
        return k
    def solve(self):
        pass
