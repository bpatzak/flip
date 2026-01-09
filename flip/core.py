import numpy as np
from enum import IntEnum
from .internalforcecomputer import InternalForceComputer


class DofID(IntEnum):
    Dx = 0
    Dy = 1
    Dz = 2
    Rx = 3
    Ry = 4
    Rz = 5


class Material:
    def __init__(self, label, e=1.0, g=1.0, alpha=1.0, d=1.0):
        self.label = label
        self.e = e
        self.g = g
        self.alpha = alpha
        self.d = d


class CrossSection:
    def __init__(self, label, a=None, iy=None, iz=None, dyz=None,
                 h=None, k=1.e30, j=None, rho=None):
        self.label = label
        self.a = a
        self.iy = iy
        self.iz = iz
        self.dyz = dyz
        self.h = h
        self.k = k
        self.j = j
        self.rho = rho  # density for self-weight


class Node:
    def __init__(self, label, domain, coords=None, bcs=None, lcs=None):
        self.label = label
        self.domain = domain
        self.coords = np.array(coords if coords is not None else [0.0, 0.0, 0.0], dtype=float)
        self.bcs = set(bcs if bcs is not None else [])
        self.update_lcs(lcs)

    def update_lcs(self, lcs):
        if lcs is None:
            self.lcs = None
            return
        locx = np.array(lcs["lx"], dtype=float)
        locy = np.array(lcs["ly"], dtype=float)
        e1 = locx / np.linalg.norm(locx)
        e2 = locy / np.linalg.norm(locy)
        e3 = np.cross(e1, e2)
        self.lcs = np.vstack([e1, e2, e3])

    def has_lcs(self):
        return self.lcs is not None

    def get_transformation_matrix(self, dofs):
        size = len(dofs)
        if self.lcs is None:
            return np.eye(size)
        T = np.zeros((size, size), dtype=float)
        for i, dof_i in enumerate(dofs):
            if dof_i in (DofID.Dx, DofID.Dy, DofID.Dz):
                for j, dof_j in enumerate(dofs):
                    if dof_j in (DofID.Dx, DofID.Dy, DofID.Dz):
                        T[i, j] = self.lcs[dof_j][dof_i]
            elif dof_i in (DofID.Rx, DofID.Ry, DofID.Rz):
                for j, dof_j in enumerate(dofs):
                    if dof_j in (DofID.Rx, DofID.Ry, DofID.Rz):
                        T[i, j] = self.lcs[dof_j - DofID.Rx][dof_i - DofID.Rx]
        return T

    # convenience
    def apply_force(self, fx=0.0, fz=0.0, My=0.0):
        self.domain.apply_nodal_load(self.label, fx=fx, fz=fz, My=My)


class Element:
    def __init__(self, label, domain, nodes, mat, cs):
        self.label = label
        self.domain = domain
        self.nodes = nodes
        self.mat = mat
        self.cs = cs

    def get_material(self):
        return self.domain.get_material(self.mat)

    def get_cs(self):
        return self.domain.get_cs(self.cs)

    def get_node_dofs(self, node):
        raise NotImplementedError

    def compute_stiffness(self):
        raise NotImplementedError

    def get_location_array(self):
        raise NotImplementedError

    def compute_geo(self):
        raise NotImplementedError

    def compute_T(self):
        raise NotImplementedError


class Beam2D(Element):
    def __init__(self, label, domain, nodes, mat, cs, hinges=(False, False)):
        super().__init__(label, domain, nodes, mat, cs)
        self.hinges = hinges
        self.ifc = InternalForceComputer(self)

    def get_node_dofs(self, node):
        return [DofID.Dx, DofID.Dz, DofID.Ry]

    def get_location_array(self):
        loc = []
        for n in self.nodes:
            loc.extend(
                self.domain.solver.get_node_location_array(
                    n, [DofID.Dx, DofID.Dz, DofID.Ry]
                )
            )
        return np.array(loc, dtype=int)

    def compute_geo(self):
        c1 = self.domain.get_node(self.nodes[0]).coords
        c2 = self.domain.get_node(self.nodes[1]).coords
        dx = c2[0] - c1[0]
        dz = c2[2] - c1[2]
        L = np.sqrt(dx * dx + dz * dz)
        return {"l": L, "dx": dx, "dz": dz}

    def has_hinges(self):
        return self.hinges[0] or self.hinges[1]

    def compute_T(self):
        geo = self.compute_geo()
        L = geo["l"]
        c = geo["dx"] / L
        s = geo["dz"] / L
        T = np.array([
            [ c,  s, 0,  0, 0, 0],
            [-s,  c, 0,  0, 0, 0],
            [ 0,  0, 1,  0, 0, 0],
            [ 0,  0, 0,  c, s, 0],
            [ 0,  0, 0, -s, c, 0],
            [ 0,  0, 0,  0, 0, 1]
        ], dtype=float)
        n1 = self.domain.get_node(self.nodes[0])
        n2 = self.domain.get_node(self.nodes[1])
        if n1.has_lcs() or n2.has_lcs():
            Tn = np.zeros((6, 6), dtype=float)
            Tn[0:3, 0:3] = n1.get_transformation_matrix(self.get_node_dofs(self.nodes[0]))
            Tn[3:6, 3:6] = n2.get_transformation_matrix(self.get_node_dofs(self.nodes[1]))
            T = T @ Tn
        return T

    def compute_local_stiffness(self, ret_condensed=False):
        geo = self.compute_geo()
        mat = self.get_material()
        cs = self.get_cs()
        L = geo["l"]
        L2 = L * L
        L3 = L2 * L
        EA = mat.e * cs.a
        EI = mat.e * cs.iy
        fi = 12 * EI / (cs.k * mat.g * cs.a * L2)
        fi1 = 1 + fi
        k = np.array([
            [ EA/L,               0,                 0,  -EA/L,               0,                 0],
            [ 0,     12*EI/L3/fi1,   -6*EI/L2/fi1,     0,   -12*EI/L3/fi1,   -6*EI/L2/fi1],
            [ 0,     -6*EI/L2/fi1, (4+fi)*EI/L/fi1,    0,    6*EI/L2/fi1, (2-fi)*EI/L/fi1],
            [-EA/L,              0,                 0,   EA/L,               0,                 0],
            [ 0,    -12*EI/L3/fi1,    6*EI/L2/fi1,     0,   12*EI/L3/fi1,    6*EI/L2/fi1],
            [ 0,     -6*EI/L2/fi1, (2-fi)*EI/L/fi1,    0,    6*EI/L2/fi1, (4+fi)*EI/L/fi1]
        ], dtype=float)
        if not self.has_hinges():
            return {"answer": k}
        if self.hinges[0] and self.hinges[1]:
            a = [0, 1, 3, 4]
            b = [2, 5]
        elif self.hinges[0]:
            a = [0, 1, 3, 4, 5]
            b = [2]
        else:
            a = [0, 1, 2, 3, 4]
            b = [5]
        a = np.array(a, dtype=int)
        b = np.array(b, dtype=int)
        kaa = k[np.ix_(a, a)]
        kab = k[np.ix_(a, b)]
        kbb = k[np.ix_(b, b)]
        k_condensed = kaa - kab @ np.linalg.inv(kbb) @ kab.T
        k2 = np.zeros((6, 6), dtype=float)
        k2[np.ix_(a, a)] = k_condensed
        if ret_condensed:
            return {"answer": k2, "a": a, "b": b, "kaa": kaa, "kab": kab, "kbb": kbb}
        return {"answer": k2}

    def compute_stiffness(self):
        kl = self.compute_local_stiffness()
        T = self.compute_T()
        return T.T @ kl["answer"] @ T

    def _end_displacement_local(self):
        T = self.compute_T()
        loc = self.get_location_array()
        r_global = self.domain.solver.r[loc]
        r_local = T @ r_global
        if self.has_hinges():
            stiffrec = self.compute_local_stiffness(ret_condensed=True)
            a = stiffrec["a"]
            b = stiffrec["b"]
            kab = stiffrec["kab"]
            kbb = stiffrec["kbb"]
            bl = np.zeros(6, dtype=float)
            for load in self.domain.get_element_loads(self.label):
                bl += load.get_load_vector_for_clamped_beam(self)
            r_local[b] = np.linalg.inv(kbb) @ (-bl[b] - kab.T @ r_local[a])
        return r_local

    def _end_forces_local(self):
        T = self.compute_T()
        loc = self.get_location_array()
        r_global = self.domain.solver.r[loc]
        r_local = T @ r_global
        stiffrec = self.compute_local_stiffness(ret_condensed=True)
        k = stiffrec["answer"]
        
        fe = k @ r_local
        bl = np.zeros(6, dtype=float)
        for load in self.domain.get_element_loads(self.label):
            bl += load.get_load_vector_for_clamped_beam(self)
        if self.has_hinges():
            a = stiffrec["a"]
            b = stiffrec["b"]
            kab = stiffrec["kab"]
            kbb = stiffrec["kbb"]
            h1 = kab @ np.linalg.inv(kbb)
            if len(b) == 1:
                blb = bl[b][0]
                for i, ai in enumerate(a):
                    fe[ai] -= bl[ai] - h1[i, 0] * blb
            else:
                fe[a] -= bl[a] - h1 @ bl[b]
        else:
            fe -= bl
        return fe

    def compute_local_deflection(self, nseg):
        rl = self._end_displacement_local()
        geo = self.compute_geo()
        L = geo["l"]
        u_vals = []
        w_vals = []
        eloads = self.domain.get_element_loads(self.label)
        for iseg in range(nseg + 1):
            xl = iseg / nseg
            w = (
                (1 - 3*xl**2 + 2*xl**3) * rl[1]
                + L * (-xl + 2*xl**2 - xl**3) * rl[2]
                + (3*xl**2 - 2*xl**3) * rl[4]
                + L * (xl**2 - xl**3) * rl[5]
            )
            u = (1 - xl) * rl[0] + xl * rl[3]
            for load in eloads:
                c = load.compute_beam_deflection_contrib(self, xl)
                w += c["w"]
                u += c["u"]
            u_vals.append(u)
            w_vals.append(w)
        return {"u": np.array(u_vals), "w": np.array(w_vals)}

    def compute_global_deflection(self, nseg):
        ld = self.compute_local_deflection(nseg)
        geo = self.compute_geo()
        L = geo["l"]
        c = geo["dx"] / L
        s = geo["dz"] / L
        u_local = ld["u"]
        w_local = ld["w"]
        u_global = u_local * c - w_local * s
        w_global = w_local * c + u_local * s
        return {"u": u_global, "w": w_global}
    
    def get_element_internal_force_polynomials(self):
        """
            Returns the element's own internal force polynomials
            from FE nodal displacements.
        """

        F = self._end_forces_local()
    
        # Axial: N(x) = EA * u'(x)
        # u(x) = N1 u1 + N2 u2
        # u'(x) = (u2 - u1)/L
        A2 = 0.0
        A1 = 0.0
        A0 = -F[0]

        # Shear: V(x) = -EI w'''(x)
        # For Hermite cubic shape functions, w'''(x) is constant
        B2 = 0.0
        B1 = 0.0
        B0 = -F[1]

        # Moment: M(x) = EI w''(x)
        # w''(x) = linear function
        C2 = 0.0
        C1 = -F[1]
        C0 = -F[2]

        D3 = 0.0  # element itself never produces cubic M

        return (A2, A1, A0), (B2, B1, B0), (D3, C2, C1, C0)



    # convenience
    def apply_udl(self, w):
        from .loads import UniformDistributedLoad
        self.domain.apply_element_load(self.label, UniformDistributedLoad(w))
