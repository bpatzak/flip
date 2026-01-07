import numpy as np


class NodalForce:
    def __init__(self, fx=0.0, fz=0.0, My=0.0):
        self.fx = fx
        self.fz = fz
        self.My = My


class UniformDistributedLoad:
    def __init__(self, w):
        self.w = w

    def get_load_vector_for_clamped_beam(self, elem):
        L = elem.compute_geo()["l"]
        w = self.w
        return np.array([
            0.0,
            w * L / 2.0,
            -w * L**2 / 12.0,
            0.0,
            w * L / 2.0,
            w * L**2 / 12.0
        ])

    def get_load_vector(self, elem):
        T = elem.compute_T()
        f = self.get_load_vector_for_clamped_beam(elem)
        if (elem.has_hinges()):
            stiffrec = elem.compute_local_stiffness(ret_condensed=True)
            a = stiffrec["a"]
            b = stiffrec["b"]
            kab = stiffrec["kab"]
            kbb = stiffrec["kbb"]
            h1 = kab @ np.linalg.inv(kbb)
            if len(b) == 1:
                fb = f[b]
            else:
                fb = f[b]
            fa = f[np.ix_(a)] - h1 @ fb
            f_condensed = np.zeros(6, dtype=float)
            f_condensed[np.ix_(a)] = fa
            f_global = T.T @ f_condensed
            return f_global
        return self.get_load_vector_for_clamped_beam(elem)
    
    def compute_beam_deflection_contrib(self, elem, xl):
        # simple zero contribution placeholder
        return {"u": 0.0, "w": 0.0}

    def compute_beam_N_contrib(self, elem, x):
        return 0.0

    def compute_beam_V_contrib(self, elem, x):
        return -self.w * x

    def compute_beam_M_contrib(self, elem, x):
        return -self.w * x**2 / 2.0


class PointLoadOnElement:
    def __init__(self, P, a):
        self.P = P
        self.a = a  # distance from node 1

    def get_load_vector_for_clamped_beam(self, elem):
        L = elem.compute_geo()["l"]
        a = self.a
        b = L - a
        P = self.P
        return np.array([
            0.0,
            P * b**2 * (3*a + b) / L**3,
            P * a * b**2 / L**2,
            0.0,
            P * a**2 * (a + 3*b) / L**3,
            -P * a**2 * b / L**2
        ])

    def compute_beam_deflection_contrib(self, elem, xl):
        return {"u": 0.0, "w": 0.0}

    def compute_beam_N_contrib(self, elem, x):
        return 0.0

    def compute_beam_V_contrib(self, elem, x):
        # piecewise: left of load, right of load
        if x <= self.a:
            return 0.0
        return -self.P

    def compute_beam_M_contrib(self, elem, x):
        if x <= self.a:
            return 0.0
        return -self.P * (x - self.a)


class LinearDistributedLoad:
    def __init__(self, w1, w2):
        self.w1 = w1
        self.w2 = w2

    def get_load_vector_for_clamped_beam(self, elem):
        # approximate by equivalent uniform + linear part
        L = elem.compute_geo()["l"]
        w1, w2 = self.w1, self.w2
        weq = (w1 + w2) / 2.0
        # reuse UDL formula for equivalent load
        udl = UniformDistributedLoad(weq)
        return udl.get_load_vector_for_clamped_beam(elem)

    def compute_beam_deflection_contrib(self, elem, xl):
        return {"u": 0.0, "w": 0.0}

    def compute_beam_N_contrib(self, elem, x):
        return 0.0

    def compute_beam_V_contrib(self, elem, x):
        L = elem.compute_geo()["l"]
        w1, w2 = self.w1, self.w2
        # V(x) = -∫0^x w(s) ds
        # w(s) = w1 + (w2-w1)s/L
        return -(w1 * x + (w2 - w1) * x**2 / (2*L))

    def compute_beam_M_contrib(self, elem, x):
        L = elem.compute_geo()["l"]
        w1, w2 = self.w1, self.w2
        # M(x) = -∫0^x (x-s) w(s) ds (not exact here, simplified)
        return 0.0


class PrestressLoad:
    def __init__(self, N0):
        self.N0 = N0

    def get_load_vector_for_clamped_beam(self, elem):
        # handled via initial stress matrix usually; here zero
        return np.zeros(6)

    def compute_beam_deflection_contrib(self, elem, xl):
        return {"u": 0.0, "w": 0.0}

    def compute_beam_N_contrib(self, elem, x):
        return self.N0

    def compute_beam_V_contrib(self, elem, x):
        return 0.0

    def compute_beam_M_contrib(self, elem, x):
        return 0.0


class SelfWeightLoad:
    def __init__(self, g=9.81):
        self.g = g

    def get_load_vector_for_clamped_beam(self, elem):
        cs = elem.get_cs()
        if cs.rho is None or cs.a is None:
            return np.zeros(6)
        w = -cs.rho * cs.a * self.g
        return UniformDistributedLoad(w).get_load_vector_for_clamped_beam(elem)

    def compute_beam_deflection_contrib(self, elem, xl):
        return {"u": 0.0, "w": 0.0}

    def compute_beam_N_contrib(self, elem, x):
        return 0.0

    def compute_beam_V_contrib(self, elem, x):
        cs = elem.get_cs()
        if cs.rho is None or cs.a is None:
            return 0.0
        w = -cs.rho * cs.a * self.g
        return -w * x

    def compute_beam_M_contrib(self, elem, x):
        cs = elem.get_cs()
        if cs.rho is None or cs.a is None:
            return 0.0
        w = -cs.rho * cs.a * self.g
        return -w * x**2 / 2.0
