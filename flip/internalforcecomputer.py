import numpy as np


class InternalForceComputer:
    """
    Encapsulates:
      - breakpoint detection
      - polynomial construction
      - exact extrema
      - diagram sampling

    Conventions:
      N(x) = A2 x^2 + A1 x + A0
      V(x) = B2 x^2 + B1 x + B0
      M(x) = D3 x^3 + C2 x^2 + C1 x + C0
    """

    def __init__(self, elem):
        self.elem = elem

    # ------------------------------------------------------------
    # Breakpoint detection
    # ------------------------------------------------------------
    def get_breakpoints(self):
        L = self.elem.compute_geo()["l"]
        pts = [0.0, L]

        for load in getattr(self.elem, "loads", []):
            if hasattr(load, "get_break_points"):
                pts.extend(load.get_break_points(self.elem))

        pts = sorted(set(pts))
        return pts

    # ------------------------------------------------------------
    # Build piecewise polynomial representation
    # ------------------------------------------------------------
    def get_polynomials(self):
        L = self.elem.compute_geo()["l"]
        breaks = self.get_breakpoints()

        segments = []

        for i in range(len(breaks) - 1):
            x0 = breaks[i]
            x1 = breaks[i+1]

            # Start with base element polynomials
            (A2, A1, A0), (B2, B1, B0), (D3, C2, C1, C0) = self.elem.get_element_internal_force_polynomials()

            for load in self.elem.domain.get_element_loads(self.elem.label):
                if not hasattr(load, "get_polynomial_contrib"):
                    continue

                # Each load returns:
                # (A2,A1,A0), (B2,B1,B0), (D3,C2,C1,C0)
                (a2, a1, a0), (b2, b1, b0), (d3, c2, c1, c0) = load.get_polynomial_contrib(self.elem)

                # Point loads only apply for x >= a
                if hasattr(load, "a"):
                    if x1 <= load.a:
                        continue

                A2 += a2; A1 += a1; A0 += a0
                B2 += b2; B1 += b1; B0 += b0
                D3 += d3; C2 += c2; C1 += c1; C0 += c0

            segments.append({
                "x0": x0,
                "x1": x1,
                "N": (A2, A1, A0),
                "V": (B2, B1, B0),
                "M": (D3, C2, C1, C0)
            })

        return segments

    # ------------------------------------------------------------
    # Exact extrema
    # ------------------------------------------------------------
    def get_extrema(self):
        segs = self.get_polynomials()

        N_candidates = []
        V_candidates = []
        M_candidates = []

        for seg in segs:
            x0, x1 = seg["x0"], seg["x1"]
            A2, A1, A0 = seg["N"]
            B2, B1, B0 = seg["V"]
            D3, C2, C1, C0 = seg["M"]

            # -------- N: quadratic → derivative 2A2 x + A1 --------
            if abs(A2) > 1e-12:
                x_star = -A1 / (2.0 * A2)
                if x0 < x_star < x1:
                    N_candidates.append((x_star, A2*x_star**2 + A1*x_star + A0))
            for x in (x0, x1):
                N_candidates.append((x, A2*x**2 + A1*x + A0))

            # -------- V: quadratic → derivative 2B2 x + B1 --------
            if abs(B2) > 1e-12:
                x_star = -B1 / (2.0 * B2)
                if x0 < x_star < x1:
                    V_candidates.append((x_star, B2*x_star**2 + B1*x_star + B0))
            for x in (x0, x1):
                V_candidates.append((x, B2*x**2 + B1*x + B0))

            # -------- M: cubic → derivative 3D3 x^2 + 2C2 x + C1 ---
            a = 3.0 * D3
            b = 2.0 * C2
            c = C1

            if abs(a) > 1e-12:
                disc = b*b - 4.0*a*c
                if disc >= 0.0:
                    rdisc = disc**0.5
                    x_roots = [(-b + rdisc) / (2.0 * a),
                               (-b - rdisc) / (2.0 * a)]
                    for x_star in x_roots:
                        if x0 < x_star < x1:
                            M_candidates.append(
                                (x_star,
                                 D3*x_star**3 + C2*x_star**2 + C1*x_star + C0)
                            )
            elif abs(b) > 1e-12:
                # derivative linear: 2C2 x + C1 = 0
                x_star = -c / b
                if x0 < x_star < x1:
                    M_candidates.append(
                        (x_star,
                         D3*x_star**3 + C2*x_star**2 + C1*x_star + C0)
                    )

            for x in (x0, x1):
                M_candidates.append((x, D3*x**3 + C2*x**2 + C1*x + C0))

        # Guard against empty lists (no loads)
        if not N_candidates:
            N_candidates = [(0.0, 0.0)]
        if not V_candidates:
            V_candidates = [(0.0, 0.0)]
        if not M_candidates:
            M_candidates = [(0.0, 0.0)]

        return {
            "N_max": max(N_candidates, key=lambda t: t[1]),
            "N_min": min(N_candidates, key=lambda t: t[1]),
            "V_max": max(V_candidates, key=lambda t: t[1]),
            "V_min": min(V_candidates, key=lambda t: t[1]),
            "M_max": max(M_candidates, key=lambda t: t[1]),
            "M_min": min(M_candidates, key=lambda t: t[1]),
        }

    # ------------------------------------------------------------
    # Diagram sampling (for plotting)
    # ------------------------------------------------------------
    def get_diagram_points(self, n_per_segment=20):
        segs = self.get_polynomials()

        xs = []
        Ns = []
        Vs = []
        Ms = []

        for seg in segs:
            x0, x1 = seg["x0"], seg["x1"]
            A2, A1, A0 = seg["N"]
            B2, B1, B0 = seg["V"]
            D3, C2, C1, C0 = seg["M"]

            xs_loc = np.linspace(x0, x1, n_per_segment)
            xs.extend(xs_loc)

            Ns.extend(A2*xs_loc**2 + A1*xs_loc + A0)
            Vs.extend(B2*xs_loc**2 + B1*xs_loc + B0)
            Ms.extend(D3*xs_loc**3 + C2*xs_loc**2 + C1*xs_loc + C0)

        return np.array(xs), np.array(Ns), np.array(Vs), np.array(Ms)
