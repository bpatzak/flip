from flip import solver
import numpy as np
from flip import dofid
import time

class LinearStatic(solver.Solver):
    def __init__(self):
        solver.Solver.__init__(self)
    def solve (self, domain):
        tstart = time.process_time()
        # number the equtions
        self.numberEquations(domain)
        #print(domain.equationNumbering)
        # assemble stifness matrix
        tneq=domain.neq+domain.pneq
        k = np.zeros([tneq, tneq])
        for e in domain.elements:
            elem = domain.getElement(e)
            ke = elem.computeStiffness()
            #print(ke)
            loc =  elem.getLocationArray()
            #print(loc)
            self.assembleMatrix(k, elem.getLocationArray(), ke)
        #print(k[0:domain.neq,0:domain.neq])

        # assemble load vector
        f = np.zeros(tneq)
        for l in domain.loads:
            locs = l.getLocationArrays(domain)
            vals = l.getValues(domain)
            #print(locs, vals)
            for iloc, ival in zip(locs, vals):
                # print (iloc, ival)
                f[iloc] += ival
        ru = np.linalg.solve(k[0:domain.neq,0:domain.neq], f[0:domain.neq])
        #print("ru:", ru)
        rp = np.zeros(domain.pneq)
        self.r = np.concatenate((ru, rp))
        self.R = np.matmul(k[domain.neq:tneq,:], self.r)-f[domain.neq:tneq]
        #print("R:", R)
        elapsedTime=time.process_time()-tstart
        print("LinearStatic: Solution done in %.2f [s]"%(elapsedTime))

        # print results
        print("Nodal displacements")
        for n in domain.nodes:
            print(n,end=' ')
            dofs = domain.equationNumbering[n]
            for d in dofs:
                print("{:s}:{:+10.6e}  ".format(str(d), self.r[dofs[d]]),end='')
            print()
        # print elements
        print("Element forces")
        for e in domain.elements:
            domain.getElement(e).printState(self)

        print("Nodal reactions")
        for n in domain.nodes:
            if (domain.getNode(n).bcs):
                print(n,end=' ')
                dofs = domain.equationNumbering[n]
                for d in domain.getNode(n).bcs:
                    print("{:s}:{:+10.6e}  ".format(str(d), self.R[dofs[d]-domain.neq]),end='')
                print()
        


