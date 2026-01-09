import numpy as np

## Without point loads, internal forces are polynomials → extrema can be found analytically.
#With point loads, the internal force diagram becomes piecewise polynomial, with discontinuities in V and kinks in M.
#With multiple point loads, the number of pieces grows, and the extrema may occur:
#at segment boundaries (i.e., at point‑load locations),
#or inside a segment (roots of the derivative of the polynomial),
#or at the element ends.
#The good news: this is completely solvable in a clean, general way.

#Point loads do not break the polynomial nature — they only break continuity.
#So the internal force functions are, for example:
#N(x): piecewise linear
#V(x): piecewise linear
#M(x): piecewise quadratic
#The only complication is that each point load introduces a new polynomial segment.
#Thus the solution is:
#
#Split the element into segments at every point load location.
#On each segment, treat the loads as purely distributed → polynomial.
#Solve for extrema inside each segment.
#Also check the boundaries (point‑load locations).
#This gives exact extrema.




class NodalForce:
    def __init__(self, fx=0.0, fz=0.0, My=0.0):
        self.fx = fx
        self.fz = fz
        self.My = My


class UniformDistributedLoad:
    def __init__(self, fx=0.0, fz=0.0, local=False):
        """
        fx, fz = load components
        local = False → fx,fz are GLOBAL components
        local = True  → fx,fz are LOCAL components (u,w)
        """
        self.fx = fx
        self.fz = fz
        self.local = local

    # ------------------------------------------------------------
    # Convert global load vector → local (u,w)
    # ------------------------------------------------------------
    def _global_to_local(self, elem):
        """
        Convert global (fx,fz) to local (u,w) using element orientation.
        """
        geo = elem.compute_geo()
        L = geo["l"]
        cx = geo["dx"] / L
        cz = geo["dz"] / L

        # local axes:
        # u = along element axis
        # w = perpendicular to element axis
        u =  self.fx * cx + self.fz * cz
        w = -self.fx * cz + self.fz * cx
        return u, w

    # ------------------------------------------------------------
    # Local load vector for a clamped beam
    # ------------------------------------------------------------
    def get_load_vector_for_clamped_beam(self, elem):
        L = elem.compute_geo()["l"]

        # choose local components
        if self.local:
            u = self.fx
            w = self.fz
        else:
            u, w = self._global_to_local(elem)

        return np.array([
            u * L / 2.0,
            w * L / 2.0,
            -w * L**2 / 12.0,
            u * L / 2.0,
            w * L / 2.0,
            w * L**2 / 12.0
        ])

    # ------------------------------------------------------------
    # Global load vector (handles hinges + transformation)
    # ------------------------------------------------------------
    def get_load_vector(self, elem):
        T = elem.compute_T()
        f_local = self.get_load_vector_for_clamped_beam(elem)

        # hinge condensation
        if elem.has_hinges():
            stiffrec = elem.compute_local_stiffness(ret_condensed=True)
            a = stiffrec["a"]
            b = stiffrec["b"]
            kab = stiffrec["kab"]
            kbb = stiffrec["kbb"]

            h1 = kab @ np.linalg.inv(kbb)
            fb = f_local[b]
            fa = f_local[np.ix_(a)] - h1 @ fb

            f_condensed = np.zeros(6, dtype=float)
            f_condensed[np.ix_(a)] = fa

            return T.T @ f_condensed

        return T.T @ f_local

    # ------------------------------------------------------------
    # Internal force contributions (local)
    # ------------------------------------------------------------
    def compute_beam_deflection_contrib(self, elem, xl):
        return {"u": 0.0, "w": 0.0}

    def get_break_points(self, elem):
        return []

    def get_polynomial_contrib(self, elem):
        """
        Returns polynomial coefficients for N, V, M:
        N(x) = A1*x + A0
        V(x) = B1*x + B0
        M(x) = C2*x^2 + C1*x + C0
        """
        if self.local:
            u = self.fx
            w = self.fz
        else:
            u, w = self._global_to_local(elem)

        # N = -u x
        A1 = -u
        A0 = 0.0

        # V = -w x
        B1 = -w
        B0 = 0.0

        # M = -w x^2 / 2
        C2 = -w / 2.0
        C1 = 0.0
        C0 = 0.0

        return (0, A1, A0), (0, B1, B0), (0, C2, C1, C0)


class PointLoadOnElement:
    def __init__(self, fx=0.0, fz=0.0, a=0.0, local=False):
        """
        fx, fz = load components (global or local)
        a      = distance from node 1
        local  = False → fx,fz are GLOBAL components
        local  = True  → fx,fz are LOCAL components (u,w)
        """
        self.fx = fx
        self.fz = fz
        self.a = a
        self.local = local

    # ------------------------------------------------------------
    # Convert global load vector → local (u,w)
    # ------------------------------------------------------------
    def _global_to_local(self, elem):
        geo = elem.compute_geo()
        L = geo["l"]
        cx = geo["dx"] / L
        cz = geo["dz"] / L

        # local axes:
        # u = along element axis
        # w = perpendicular to element axis
        u =  self.fx * cx + self.fz * cz
        w = -self.fx * cz + self.fz * cx
        return u, w

    # ------------------------------------------------------------
    # Local load vector for a clamped beam
    # ------------------------------------------------------------
    def get_load_vector_for_clamped_beam(self, elem):
        L = elem.compute_geo()["l"]
        a = self.a
        b = L - a

        # choose local components
        if self.local:
            u = self.fx
            w = self.fz
        else:
            u, w = self._global_to_local(elem)

        # axial point load (u) → only axial DOFs
        # transverse point load (w) → bending DOFs
        return np.array([
            u * b / L,                         # Fx1
            w * b**2 * (3*a + b) / L**3,       # Fz1
            -w * a * b**2 / L**2,               # M1
            u * a / L,                         # Fx2
            w * a**2 * (a + 3*b) / L**3,       # Fz2
            w * a**2 * b / L**2               # M2
        ])

    # ------------------------------------------------------------
    # Global load vector (handles hinges + transformation)
    # ------------------------------------------------------------
    def get_load_vector(self, elem):
        T = elem.compute_T()
        f_local = self.get_load_vector_for_clamped_beam(elem)

        # hinge condensation
        if elem.has_hinges():
            stiffrec = elem.compute_local_stiffness(ret_condensed=True)
            a = stiffrec["a"]
            b = stiffrec["b"]
            kab = stiffrec["kab"]
            kbb = stiffrec["kbb"]

            h1 = kab @ np.linalg.inv(kbb)
            fb = f_local[b]
            fa = f_local[np.ix_(a)] - h1 @ fb

            f_condensed = np.zeros(6, dtype=float)
            f_condensed[np.ix_(a)] = fa

            return T.T @ f_condensed

        return T.T @ f_local

    # ------------------------------------------------------------
    # Internal force contributions (local)
    # ------------------------------------------------------------
    def compute_beam_deflection_contrib(self, elem, xl):
        return {"u": 0.0, "w": 0.0}

    def get_break_points(self, elem):
        return [self.a]

    def get_polynomial_contrib(self, elem):
        """
        Returns polynomial coefficients valid for x >= a.
        For x < a, contribution is zero.
        """
        if self.local:
            u = self.fx
            w = self.fz
        else:
            u, w = self._global_to_local(elem)

        # For x >= a:
        # N = -u
        A1 = 0.0
        A0 = -u

        # V = -w
        B1 = 0.0
        B0 = -w

        # M = -w (x - a) = -w x + w a
        C2 = 0.0
        C1 = -w
        C0 = w * self.a

        return (0, A1, A0), (0, B1, B0), (0, C2, C1, C0)


class SelfWeightLoad:
    def __init__(self, g=9.81, local=False):
        """
        g     = gravitational acceleration (positive value)
        local = False → weight acts in GLOBAL -Z direction
                True  → weight acts in LOCAL -w direction
        """
        self.g = g
        self.local = local

    # ------------------------------------------------------------
    # Compute distributed load components (u,w) in local system
    # ------------------------------------------------------------
    def _compute_local_components(self, elem):
        cs = elem.get_cs()
        mat = elem.get_material()

        if cs.rho is None or cs.a is None:
            return 0.0, 0.0  # no weight

        w_global = -cs.rho * cs.a * self.g  # global Z downward

        if self.local:
            # weight acts in local -w direction
            u = 0.0
            w = w_global
            return u, w

        # convert global (0, w_global) → local (u,w)
        geo = elem.compute_geo()
        L = geo["l"]
        cx = geo["dx"] / L
        cz = geo["dz"] / L

        # global load vector = (fx=0, fz=w_global)
        fx = 0.0
        fz = w_global

        # local components
        u =  fx * cx + fz * cz
        w = -fx * cz + fz * cx
        return u, w

    # ------------------------------------------------------------
    # Local equivalent nodal load vector
    # ------------------------------------------------------------
    def get_load_vector_for_clamped_beam(self, elem):
        L = elem.compute_geo()["l"]
        u, w = self._compute_local_components(elem)

        return np.array([
            u * L / 2.0,
            w * L / 2.0,
            -w * L**2 / 12.0,
            u * L / 2.0,
            w * L / 2.0,
            w * L**2 / 12.0
        ])

    # ------------------------------------------------------------
    # Global load vector (handles hinges + transformation)
    # ------------------------------------------------------------
    def get_load_vector(self, elem):
        T = elem.compute_T()
        f_local = self.get_load_vector_for_clamped_beam(elem)

        # hinge condensation
        if elem.has_hinges():
            stiffrec = elem.compute_local_stiffness(ret_condensed=True)
            a = stiffrec["a"]
            b = stiffrec["b"]
            kab = stiffrec["kab"]
            kbb = stiffrec["kbb"]

            h1 = kab @ np.linalg.inv(kbb)
            fb = f_local[b]
            fa = f_local[np.ix_(a)] - h1 @ fb

            f_condensed = np.zeros(6, dtype=float)
            f_condensed[np.ix_(a)] = fa

            return T.T @ f_condensed

        return T.T @ f_local

    # ------------------------------------------------------------
    # Internal force contributions (local)
    # ------------------------------------------------------------
    def compute_beam_deflection_contrib(self, elem, xl):
        return {"u": 0.0, "w": 0.0}

    def get_break_points(self, elem):
        return [self.a]

    def get_polynomial_contrib(self, elem):
        """
        Returns polynomial coefficients valid for x >= a.
        For x < a, contribution is zero.
        """
        if self.local:
            u = self.fx
            w = self.fz
        else:
            u, w = self._global_to_local(elem)

        # For x >= a:
        # N = -u
        A1 = 0.0
        A0 = -u

        # V = -w
        B1 = 0.0
        B0 = -w

        # M = -w (x - a) = -w x + w a
        C2 = 0.0
        C1 = -w
        C0 = w * self.a

        return (0, A1, A0), (0, B1, B0), (0, C2, C1, C0)


class TemperatureLoad:
    def __init__(self, dT=0.0, dT_top=0.0, dT_bottom=0.0, local=True):
        """
        dT        = uniform temperature change (°C)
        dT_top    = temperature at top fiber (°C)
        dT_bottom = temperature at bottom fiber (°C)
        local     = temperature loads are always local, but kept for consistency
        """
        self.dT = dT
        self.dT_top = dT_top
        self.dT_bottom = dT_bottom
        self.local = local

    # ------------------------------------------------------------
    # Local equivalent nodal load vector
    # ------------------------------------------------------------
    def get_load_vector_for_clamped_beam(self, elem):
        mat = elem.get_material()
        cs  = elem.get_cs()
        geo = elem.compute_geo()

        E = mat.e
        A = cs.a
        I = cs.iy
        alpha = mat.alpha
        L = geo["l"]

        # ----------------------------------------
        # 1) Uniform temperature → axial force
        # ----------------------------------------
        N_T = E * A * alpha * self.dT

        # axial equivalent nodal forces
        fN = np.array([
            +N_T/2,   # Fx1
            0.0,      # Fz1
            0.0,      # M1
            -N_T/2,   # Fx2
            0.0,      # Fz2
            0.0       # M2
        ])

        # ----------------------------------------
        # 2) Temperature gradient → bending moment
        # ----------------------------------------
        dTg = self.dT_top - self.dT_bottom

        if cs.h is not None and cs.h > 0:
            M_T = E * I * alpha * (dTg / cs.h)
        else:
            M_T = 0.0

        # bending equivalent nodal moments
        fM = np.array([
            0.0,
            0.0,
            +M_T,   # M1
            0.0,
            0.0,
            -M_T    # M2
        ])

        return fN + fM

    # ------------------------------------------------------------
    # Global load vector (handles hinges + transformation)
    # ------------------------------------------------------------
    def get_load_vector(self, elem):
        T = elem.compute_T()
        f_local = self.get_load_vector_for_clamped_beam(elem)

        # hinge condensation
        if elem.has_hinges():
            stiffrec = elem.compute_local_stiffness(ret_condensed=True)
            a = stiffrec["a"]
            b = stiffrec["b"]
            kab = stiffrec["kab"]
            kbb = stiffrec["kbb"]

            h1 = kab @ np.linalg.inv(kbb)
            fb = f_local[b]
            fa = f_local[np.ix_(a)] - h1 @ fb

            f_condensed = np.zeros(6, dtype=float)
            f_condensed[np.ix_(a)] = fa

            return T.T @ f_condensed

        return T.T @ f_local

    # ------------------------------------------------------------
    # Internal force contributions
    # ------------------------------------------------------------
    def compute_beam_deflection_contrib(self, elem, xl):
        return {"u": 0.0, "w": 0.0}

    def get_break_points(self, elem):
        return []

    def get_polynomial_contrib(self, elem):
        mat = elem.get_material()
        cs  = elem.get_cs()
        E = mat.e
        A = cs.a
        I = cs.iy
        alpha = mat.alpha

        # uniform temperature
        N_T = E * A * alpha * self.dT

        # gradient
        if cs.h:
            M_T = E * I * alpha * ((self.dT_top - self.dT_bottom) / cs.h)
        else:
            M_T = 0.0

        # N = -N_T
        A1 = 0.0
        A0 = -N_T

        # V = 0
        B1 = 0.0
        B0 = 0.0

        # M = M_T
        C2 = 0.0
        C1 = 0.0
        C0 = M_T

        return (0, A1, A0), (0, B1, B0), (0, C2, C1, C0)


class PrestressLoad:
    def __init__(self, N0=0.0, local=True):
        """
        N0    = axial prestress force (positive = tension, negative = compression)
        local = True  → N0 is along the element axis (recommended)
                False → N0 is given in global coordinates (fx,fz)
        """
        self.N0 = N0
        self.local = local

    # ------------------------------------------------------------
    # Convert global axial force → local (u)
    # ------------------------------------------------------------
    def _global_to_local(self, elem):
        geo = elem.compute_geo()
        L = geo["l"]
        cx = geo["dx"] / L
        cz = geo["dz"] / L

        # global axial force components
        fx = self.N0
        fz = 0.0

        # local axial component
        u = fx * cx + fz * cz
        return u

    # ------------------------------------------------------------
    # Local equivalent nodal load vector
    # ------------------------------------------------------------
    def get_load_vector_for_clamped_beam(self, elem):
        """
        Prestress produces only axial forces:
            +N0/2 at node 1
            -N0/2 at node 2
        """
        if self.local:
            u = self.N0
        else:
            u = self._global_to_local(elem)

        return np.array([
            +u/2,   # Fx1
            0.0,    # Fz1
            0.0,    # M1
            -u/2,   # Fx2
            0.0,    # Fz2
            0.0     # M2
        ])

    # ------------------------------------------------------------
    # Global load vector (handles hinges + transformation)
    # ------------------------------------------------------------
    def get_load_vector(self, elem):
        T = elem.compute_T()
        f_local = self.get_load_vector_for_clamped_beam(elem)

        # hinge condensation
        if elem.has_hinges():
            stiffrec = elem.compute_local_stiffness(ret_condensed=True)
            a = stiffrec["a"]
            b = stiffrec["b"]
            kab = stiffrec["kab"]
            kbb = stiffrec["kbb"]

            h1 = kab @ np.linalg.inv(kbb)
            fb = f_local[b]
            fa = f_local[np.ix_(a)] - h1 @ fb

            f_condensed = np.zeros(6, dtype=float)
            f_condensed[np.ix_(a)] = fa

            return T.T @ f_condensed

        return T.T @ f_local

    # ------------------------------------------------------------
    # Internal force contributions (local)
    # ------------------------------------------------------------
    def compute_beam_deflection_contrib(self, elem, xl):
        return {"u": 0.0, "w": 0.0}
    
    def get_break_points(self, elem):
        return []

    def get_polynomial_contrib(self, elem):
        if self.local:
            u = self.N0
        else:
            u = self._global_to_local(elem)

        # N = -u
        A1 = 0.0
        A0 = -u

        # V = 0
        B1 = 0.0
        B0 = 0.0

        # M = 0
        C2 = 0.0
        C1 = 0.0
        C0 = 0.0

        return (0, A1, A0), (0, B1, B0), (0, C2, C1, C0)

