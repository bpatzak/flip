import numpy as np
from .core import DofID


class Domain:
    def __init__(self):
        self.nodes = {}
        self.elements = {}
        self.materials = {}
        self.cross_sections = {}

        self.nodal_loads = {}    # node_label → (fx, fz, My)
        self.element_loads = {}  # elem_label → [load objects]

        self.solver = None

    def add_node(self, node):
        self.nodes[node.label] = node

    def get_node(self, label):
        return self.nodes[label]

    def add_material(self, mat):
        self.materials[mat.label] = mat

    def get_material(self, label):
        return self.materials[label]

    def add_cs(self, cs):
        self.cross_sections[cs.label] = cs

    def get_cs(self, label):
        return self.cross_sections[label]

    def add_element(self, elem):
        self.elements[elem.label] = elem

    def get_element(self, label):
        return self.elements[label]

    def set_solver(self, solver):
        self.solver = solver
        solver.domain = self

    # loads
    def apply_nodal_load(self, node_label, fx=0.0, fz=0.0, My=0.0):
        self.nodal_loads[node_label] = (fx, fz, My)

    def apply_element_load(self, elem_label, load):
        self.element_loads.setdefault(elem_label, []).append(load)

    def get_element_loads(self, elem_label):
        return self.element_loads.get(elem_label, [])
    
import time


class Solver:
    def __init__(self):
        self.domain = Domain()
        self.domain.set_solver(self)

        self.neq = 0          # number of unknown DOFs
        self.pneq = 0         # number of prescribed DOFs
        self.code_numbers = {}  # node_label → {dof: code}
        self.code_generated = False

        self.K = None         # global stiffness matrix
        self.f = None         # global load vector
        self.r = None         # global displacement vector
        self.R = None         # reactions

    # ------------------------------------------------------------
    # Code number generation (faithful to TS version)
    # ------------------------------------------------------------
    def generate_code_numbers(self):
        nodal_dofs = {}

        # initialize
        for node_label, node in self.domain.nodes.items():
            self.code_numbers[node_label] = {}
            nodal_dofs[node_label] = set()

        # collect DOFs required by elements
        for elem in self.domain.elements.values():
            for nlabel in elem.nodes:
                dofs = elem.get_node_dofs(nlabel)
                nodal_dofs[nlabel].update(dofs)

        # count unknowns and prescribed
        self.neq = 0
        self.pneq = 0

        for node_label, node in self.domain.nodes.items():
            for d in nodal_dofs[node_label]:
                if d in node.bcs:
                    self.pneq += 1
                else:
                    self.neq += 1

        # assign code numbers
        eq = 0
        peq = self.neq

        for node_label, node in self.domain.nodes.items():
            for d in nodal_dofs[node_label]:
                if d in node.bcs:
                    self.code_numbers[node_label][d] = peq
                    peq += 1
                else:
                    self.code_numbers[node_label][d] = eq
                    eq += 1

        self.code_generated = True

    # ------------------------------------------------------------
    # Helper: get location array for a node
    # ------------------------------------------------------------
    def get_node_location_array(self, node_label, dofs):
        return [self.code_numbers[node_label][d] for d in dofs]

    # ------------------------------------------------------------
    # Assembly
    # ------------------------------------------------------------
    def assemble(self):
        ndof = self.neq + self.pneq
        self.K = np.zeros((ndof, ndof), dtype=float)
        self.f = np.zeros(ndof, dtype=float)

        # assemble stiffness
        for elem in self.domain.elements.values():
            ke = elem.compute_stiffness()
            loc = elem.get_location_array()

            for i in range(6):
                I = loc[i]
                for j in range(6):
                    J = loc[j]
                    self.K[I, J] += ke[i, j]

        # assemble nodal loads
        for node_label, (fx, fz, My) in self.domain.nodal_loads.items():
            loc = self.get_node_location_array(node_label, [DofID.Dx, DofID.Dz, DofID.Ry])
            vals = [fx, fz, My]
            for eq, val in zip(loc, vals):
                self.f[eq] += val

        # assemble element loads
        for elem in self.domain.elements.values():
            for load in self.domain.get_element_loads(elem.label):
                fe = load.get_load_vector(elem)
                loc = elem.get_location_array()
                for i in range(6):
                    I = loc[i]
                    self.f[I] += fe[i]

    # ------------------------------------------------------------
    # Solve system
    # ------------------------------------------------------------
    def solve(self):
        start = time.perf_counter()

        if not self.code_generated:
            self.generate_code_numbers()

        self.assemble()

        if self.neq == 0:
            self.r = np.zeros(self.neq + self.pneq)
            self.R = np.zeros(self.pneq)
            return self.r

        # unknown and prescribed index ranges
        unknowns = np.arange(0, self.neq)
        prescribed = np.arange(self.neq, self.neq + self.pneq)

        # prescribed displacements are zero (Dirichlet BC)
        rp = np.zeros(self.pneq)

        # compute RHS: fu - K_up * rp
        Kuu = self.K[np.ix_(unknowns, unknowns)]
        Kup = self.K[np.ix_(unknowns, prescribed)]
        fu = self.f[unknowns]

        b = fu - Kup @ rp

        # solve
        ru = np.linalg.solve(Kuu, b)

        # assemble full displacement vector
        self.r = np.zeros(self.neq + self.pneq)
        self.r[unknowns] = ru
        self.r[prescribed] = rp

        # reactions
        Kpu = self.K[np.ix_(prescribed, unknowns)]
        Kpp = self.K[np.ix_(prescribed, prescribed)]
        fp = self.f[prescribed]

        self.R = Kpu @ ru + Kpp @ rp - fp

        elapsed = time.perf_counter() - start
        print(f"Solution took {elapsed:.3f} sec")

        return self.r
