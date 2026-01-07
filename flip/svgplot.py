import math
import numpy as np
import drawsvg as draw

from flip import domain
from .core import DofID


# ------------------------------------------------------------
# Coordinate transform builder
# ------------------------------------------------------------
def _sx_sy_builder(xs, zs, width_px, height_px, margin):
    xs = np.array(xs)
    zs = np.array(zs)
    xmin, xmax = xs.min(), xs.max()
    zmin, zmax = zs.min(), zs.max()
    if xmax == xmin:
        xmax += 1.0
        xmin -= 1.0
    if zmax == zmin:
        zmax += 1.0
        zmin -= 1.0

    def sx(x):
        return margin + (x - xmin) / (xmax - xmin) * (width_px - 2*margin)

    def sy(z):
        #return height_px - (margin + (z - zmin) / (zmax - zmin) * (height_px - 2*margin))
        return (margin + (z - zmin) / (zmax - zmin) * (height_px - 2*margin))

    return sx, sy

# ------------------------------------------------------------
# Helper: scale X and Y to drawing coordinates
# ------------------------------------------------------------
def _build_scalers(xs, ys, width_px, height_px, margin):
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)

    dx = xmax - xmin if xmax > xmin else 1.0
    dy = ymax - ymin if ymax > ymin else 1.0

    sx = lambda x: margin + (x - xmin) / dx * (width_px - 2*margin)
    sy = lambda y: (margin + (y - ymin) / dy * (height_px - 2*margin))

    return sx, sy


# ------------------------------------------------------------
# Support symbols (Polygon-free)
# ------------------------------------------------------------

def draw_pinned_support(d, x, z, size=20):
    base = z + size

    # Base line
    d.append(draw.Line(x - size, base, x + size, base, stroke='black'))

    # Triangle using Path
    p = draw.Path(stroke='black', fill='lightgray')
    p.M(x - size, base)
    p.L(x + size, base)
    p.L(x, z)
    p.Z()
    d.append(p)

def draw_roller_support(d, x, z, size=20, lcs=None):
    """
    Roller support aligned with node LCS.
    Base line is offset outward from the node.
    Polygon-free version (Path, Line, Circle only).
    """

    # ------------------------------------------------------------
    # Determine roller free direction from node LCS
    # ------------------------------------------------------------
    if lcs is not None:
        # local z-axis projected to XZ plane
        nx = -lcs[2][0]
        nz = -lcs[2][2]
    else:
        nx, nz = 0.0, -1.0  # default: roller free in +X

    # Normalize
    L = (nx*nx + nz*nz)**0.5
    if L == 0:
        nx, nz = 1.0, 0.0
    else:
        nx /= L
        nz /= L

    # ------------------------------------------------------------
    # Support direction = opposite of roller free direction
    # ------------------------------------------------------------
    sx_dir = -nx
    sz_dir = -nz

    # Perpendicular direction (base line direction)
    tx, tz = -sz_dir, sx_dir

    # ------------------------------------------------------------
    # Offset the base line outward from the node
    # ------------------------------------------------------------
    offset = size * 0.6
    bx = x + sx_dir * offset
    bz = z + sz_dir * offset

    # Base endpoints
    bx1 = bx - tx * size
    bz1 = bz - tz * size
    bx2 = bx + tx * size
    bz2 = bz + tz * size

    d.append(draw.Line(bx1, bz1, bx2, bz2, stroke='black'))

    # ------------------------------------------------------------
    # Triangle (support body)
    # ------------------------------------------------------------
    p = draw.Path(stroke='black', fill='white')
    p.M(bx1, bz1)
    p.L(bx2, bz2)
    p.L(x, z)
    p.Z()
    d.append(p)

    # ------------------------------------------------------------
    # Rollers (two circles)
    # ------------------------------------------------------------
    r = size * 0.25
    roll_offset = size * 0.6
    bx = x + sx_dir * size * 0.85
    bz = z + sz_dir * size * 0.85

    # Roller 1
    rx1 = bx - tx * roll_offset
    rz1 = bz - tz * roll_offset
    d.append(draw.Circle(rx1, rz1, r, fill='white', stroke='black'))

    # Roller 2
    rx2 = bx + tx * roll_offset
    rz2 = bz + tz * roll_offset
    d.append(draw.Circle(rx2, rz2, r, fill='white', stroke='black'))




def draw_fixed_support(d, x, z, size=20):
    # Vertical wall
    d.append(draw.Line(x, z - size, x, z + size, stroke='black', stroke_width=3))

    # Hatch lines
    for i in range(-size, size, 6):
        d.append(draw.Line(x, z + i, x - size, z + i + 6, stroke='black'))


# ------------------------------------------------------------
# Loads (Polygon-free)
# ------------------------------------------------------------

def draw_force(d, x, z, fx, fz, scale=0.002):
    dx = fx * scale
    dz = -fz * scale
    d.append(draw.Line(x, z, x + dx, z + dz, stroke='blue', stroke_width=2))

    # Arrowhead (simple V shape)
    ah = 6
    d.append(draw.Line(x + dx, z + dz,
                       x + dx - ah*np.sign(dx+1e-9),
                       z + dz - ah*np.sign(dz+1e-9),
                       stroke='blue'))


def draw_moment(d, x, z, My, radius=20):
    if My == 0:
        return

    # Approximate circular arc using a Path
    p = draw.Path(stroke='green', fill='none', stroke_width=2)

    # Sweep direction
    sign = 1 if My > 0 else -1

    # Approximate 270° arc with 6 segments
    for i in range(7):
        ang = sign * (i / 6) * 1.5 * np.pi
        px = x + radius * np.cos(ang)
        pz = z + radius * np.sin(ang)
        if i == 0:
            p.M(px, pz)
        else:
            p.L(px, pz)

    d.append(p)


def draw_distributed_load(d, x1, z1, x2, z2, w, n_arrows=8, scale=0.002):
    dx = x2 - x1
    dz = z2 - z1
    L = np.sqrt(dx*dx + dz*dz)
    if L == 0:
        return

    tx, tz = dx/L, dz/L
    nx, nz = -tz, tx

    for i in range(n_arrows + 1):
        s = i / n_arrows
        px = x1 + s * dx
        pz = z1 + s * dz

        ax = px + nx * w * scale
        az = pz + nz * w * scale

        d.append(draw.Line(px, pz, ax, az, stroke='blue'))


# ------------------------------------------------------------
# Model plot
# ------------------------------------------------------------

def _draw_hinge(d, elem, hinge_indx, SX, SY, radius=6, color='black'): 
    n1 = elem.domain.get_node(elem.nodes[0])
    n2 = elem.domain.get_node(elem.nodes[1])
    x1, z1 = SX(n1.coords[0]), SY(n1.coords[2])
    x2, z2 = SX(n2.coords[0]), SY(n2.coords[2])
    l = math.hypot(x2 - x1, z2 - z1)
    if hinge_indx == 0:
        x = x1 + (radius / l) * (x2 - x1)
        y = z1 + (radius / l) * (z2 - z1)
    else:
        x = x2 - (radius / l) * (x2 - x1)
        y = z2 - (radius / l) * (z2 - z1) 
    d.append(draw.Circle(x, y, radius, fill='white', stroke=color, stroke_width=2))

def plot_model_drawsvg(domain, filename="model.svg",
                       width_px=800, height_px=600, margin=40,
                       show_node_labels=True,
                       show_element_labels=True,
                       show_deformed=False,
                       deform_scale=1.0,
                       nseg=40):
    """
    Draws the undeformed model and optionally the deformed shape.
    Deformed shape is drawn in real geometry using only Path/Line/Circle/Text.
    """

    # ------------------------------------------------------------
    # Collect coordinates for global scaling
    # ------------------------------------------------------------
    xs, zs = [], []
    for node in domain.nodes.values():
        xs.append(node.coords[0])
        zs.append(node.coords[2])

    sx, sy = _sx_sy_builder(xs, zs, width_px, height_px, margin)
    d = draw.Drawing(width_px, height_px, origin=(0, 0))

    # ------------------------------------------------------------
    # Draw elements (undeformed)
    # ------------------------------------------------------------
    for elem in domain.elements.values():
        n1 = domain.get_node(elem.nodes[0])
        n2 = domain.get_node(elem.nodes[1])

        x1, z1 = sx(n1.coords[0]), sy(n1.coords[2])
        x2, z2 = sx(n2.coords[0]), sy(n2.coords[2])

        d.append(draw.Line(x1, z1, x2, z2, stroke='black', stroke_width=2))

        if elem.hinges[0]:
            _draw_hinge(d, elem, 0, sx, sy)
        if elem.hinges[1]:
            _draw_hinge(d, elem, 1, sx, sy)

        if show_element_labels:
            xm, zm = (x1 + x2) / 2, (z1 + z2) / 2
            d.append(draw.Text(f"E{elem.label}", 14, xm, zm - 5,
                               center=True, fill='blue'))

        # Distributed loads
        for load in domain.get_element_loads(elem.label):
            if hasattr(load, "w"):
                draw_distributed_load(d, x1, z1, x2, z2, load.w)

    # ------------------------------------------------------------
    # Draw nodes, supports, nodal loads
    # ------------------------------------------------------------
    for node in domain.nodes.values():
        x, z = sx(node.coords[0]), sy(node.coords[2])
        d.append(draw.Circle(x, z, 4, fill='red'))

        if show_node_labels:
            d.append(draw.Text(f"N{node.label}", 12, x + 6, z - 6, fill='darkred'))

        # Supports
        bcs = node.bcs
        if DofID.Dx in bcs and DofID.Dz in bcs and DofID.Ry in bcs:
            draw_fixed_support(d, x, z)
        elif DofID.Dx in bcs and DofID.Dz in bcs:
            draw_pinned_support(d, x, z)
        elif DofID.Dz in bcs:
            draw_roller_support(d, x, z, lcs=node.lcs)

        # Nodal loads
        if node.label in domain.nodal_loads:
            fx, fz, My = domain.nodal_loads[node.label]
            if fx != 0 or fz != 0:
                draw_force(d, x, z, fx, fz)
            if My != 0:
                draw_moment(d, x, z, My)

    # ------------------------------------------------------------
    # Draw deformed shape (optional)
    # ------------------------------------------------------------
    if show_deformed:
        if domain.solver is None or domain.solver.r is None:
            raise RuntimeError("Solver must be run before drawing deformed shape.")

        # automatic scaling based on model size
        model_span = max(xs) - min(xs)
        auto_scale = deform_scale * 0.05 * model_span

        # collect all deformed coordinates for max/min detection
        all_dx = []
        all_dz = []
        all_u = []
        all_w = []

        # First pass: compute all deformed points
        elem_shapes = {}  # elem.label → list of (x_def, z_def)
        for elem in domain.elements.values():
            dsh = elem.compute_global_deflection(nseg)
            geo = elem.compute_geo()

            n1 = domain.get_node(elem.nodes[0])
            n2 = domain.get_node(elem.nodes[1])

            x1, z1 = n1.coords[0], n1.coords[2]
            dx = n2.coords[0] - x1
            dz = n2.coords[2] - z1
            L = geo["l"]

            pts = []
            for i in range(nseg + 1):
                s = i / nseg
                px = x1 + s * dx
                pz = z1 + s * dz

                # add deformation
                px_def = px + dsh["u"][i] * auto_scale
                pz_def = pz + dsh["w"][i] * auto_scale

                pts.append((px_def, pz_def))
                all_dx.append(px_def)
                all_dz.append(pz_def)
                all_u.append(dsh["u"][i])
                all_w.append(dsh["w"][i])

            elem_shapes[elem.label] = pts

        # --------------------------------------------------------
        # Draw deformed shape (semi-transparent)
        # --------------------------------------------------------
        for elem in domain.elements.values():
            pts = elem_shapes[elem.label]

            p = draw.Path(
                stroke='purple',
                fill='none',
                stroke_width=2,
                stroke_opacity=0.6   # transparency
            )

            for i, (px_def, pz_def) in enumerate(pts):
                if i == 0:
                    p.M(sx(px_def), sy(pz_def))
                else:
                    p.L(sx(px_def), sy(pz_def))

            d.append(p)

        # --------------------------------------------------------
        # Mark max/min displacement with numeric values
        # --------------------------------------------------------
        all_u = np.array(all_u)
        all_w = np.array(all_w)
        disp_mag = np.sqrt(all_u**2 + all_w**2)

        imax = np.argmax(disp_mag)
        imin = np.argmin(disp_mag)

        xmax = all_dx[imax]
        zmax = all_dz[imax]
        xmin = all_dx[imin]
        zmin = all_dz[imin]

        vmax_val = disp_mag[imax]
        vmin_val = disp_mag[imin]

        # Max marker
        d.append(draw.Circle(sx(xmax), sy(zmax), 5, fill='red'))
        d.append(draw.Text(f"{vmax_val:.4g}", 12,
                           sx(xmax), sy(zmax) - 12,
                           center=True, fill='red'))

        # Min marker
        d.append(draw.Circle(sx(xmin), sy(zmin), 5, fill='blue'))
        d.append(draw.Text(f"{vmin_val:.4g}", 12,
                           sx(xmin), sy(zmin) + 14,
                           center=True, fill='blue'))

    # ------------------------------------------------------------
    # Save SVG
    # ------------------------------------------------------------
    # d.save_svg(filename)
    return d



import drawsvg as draw
import numpy as np
import math


# ------------------------------------------------------------
# Helper: rotate + translate a point
# ------------------------------------------------------------
def _rot_trans(x, y, x0, y0, c, s):
    xr = x * c - y * s + x0
    yr = x * s + y * c + y0
    return xr, yr

def _draw_filled_diagram(d, xs, vals, sf, n1, c, s, SX, SY, color, alpha=0.25):
    """
    Draws a filled diagram polygon along the element.
    xs, vals: local coordinates and diagram values
    sf: scale factor
    n1: first node (global coords)
    c, s: cos/sin of element angle
    """
    x1 = n1.coords[0]
    z1 = n1.coords[2]
    # Build polygon points: baseline forward, diagram backward
    pts = []

    # 1) Baseline (element axis)
    for x in xs:
        Xg, Yg = _rot_trans(x, 0.0, x1, z1, c, s)
        pts.append((SX(Xg), SY(Yg)))

    # 2) Diagram curve (reverse order)
    for i in reversed(range(len(xs))):
        x = xs[i]
        v = vals[i] * sf
        Xg, Yg = _rot_trans(x, v, x1, z1, c, s)
        pts.append((SX(Xg), SY(Yg)))

    # Create polygon
    poly = draw.Path(fill=color, fill_opacity=alpha, stroke='none', stroke_opacity=alpha)
    x0, y0 = pts[0]
    poly.M(x0, y0)
    for (x, y) in pts[1:]:
        poly.L(x, y)
    poly.Z()

    d.append(poly)

# ------------------------------------------------------------
# Main function: plot N/V/M diagrams on the actual structure
# ------------------------------------------------------------
def plot_internal_forces_on_structure(domain,
                                      filename="structure_diagrams.svg",
                                      width_px=1200,
                                      height_px=900,
                                      margin=40,
                                      diagrams=("N", "V", "M"),
                                      scale=0.15,
                                      show_extrema=True,
                                      n_per_segment=40,
                                      flip={"N": -1, "V": -1, "M": 1}):
    """
    Draws internal force diagrams directly on the structure.

    diagrams = tuple/list of any of:
        "N" → normal force
        "V" → shear force
        "M" → bending moment

    scale = diagram height scaling factor (relative to element length)
    """

    diagrams = set(diagrams)

    # --------------------------------------------------------
    # Determine global bounding box of structure
    # --------------------------------------------------------
    xs_all = []
    ys_all = []
    for node in domain.nodes.values():
        xs_all.append(node.coords[0])
        ys_all.append(node.coords[2])

    Sx, Sy = _build_scalers(xs_all, ys_all, width_px, height_px, margin)
    

    # Drawing
    d = draw.Drawing(width_px, height_px, origin=(0, 0))


    # --------------------------------------------------------
    # Draw diagrams on each element
    # --------------------------------------------------------
    colors = {"N": "red", "V": "green", "M": "blue"}

    for elem in domain.elements.values():
        n1 = domain.get_node(elem.nodes[0])
        n2 = domain.get_node(elem.nodes[1])

        x1, z1 = n1.coords[0], n1.coords[2]
        x2, z2 = n2.coords[0], n2.coords[2]
        # Element geometry
        dx_e = x2 - x1
        dy_e = z2 - z1
        L = math.hypot(dx_e, dy_e)

        c = dx_e / L
        s = dy_e / L

        # Local → global offset direction (perpendicular)
        # For N: offset along normal direction
        # For V/M: offset perpendicular to element
        nx = -s
        ny =  c

        # Get polynomial diagram points
        xs, Ns, Vs, Ms = elem.ifc.get_diagram_points(n_per_segment)

        # Normalize diagrams for visual scaling
        # (use element length as reference)
        max_val = 1e-12
        if "N" in diagrams: max_val = max(max_val, max(abs(Ns)))
        if "V" in diagrams: max_val = max(max_val, max(abs(Vs)))
        if "M" in diagrams: max_val = max(max_val, max(abs(Ms)))

        if max_val < 1e-12:
            max_val = 1.0

        # scaling factor
        sf = scale * L / max_val

        # Filled diagrams 
        if "N" in diagrams: _draw_filled_diagram(d, xs, Ns, sf*flip["N"], n1, c, s, Sx, Sy, colors["N"])
        if "V" in diagrams: _draw_filled_diagram(d, xs, Vs, sf*flip["V"], n1, c, s, Sx, Sy, colors["V"])
        if "M" in diagrams: _draw_filled_diagram(d, xs, Ms, sf*flip["M"], n1, c, s, Sx, Sy, colors["M"])

        # ----------------------------------------------------
        # Draw each selected diagram
        # ----------------------------------------------------
        def draw_diagram(values, color, label):
            p = draw.Path(stroke=color, fill='none', stroke_width=2)

            for i in range(len(xs)):
                xloc = xs[i]
                vloc = values[i] * sf * flip[label]

                # Local coordinates:
                # point on element: (xloc, 0)
                # diagram offset:   (0, vloc)
                # Convert to global:
                Xg, Yg = _rot_trans(xloc, vloc, x1, z1, c, s)

                if i == 0:
                    p.M(Sx(Xg), Sy(Yg))
                else:
                    p.L(Sx(Xg), Sy(Yg))

            d.append(p)

        if "N" in diagrams:
            draw_diagram(Ns, colors["N"], "N")

        if "V" in diagrams:
            draw_diagram(Vs, colors["V"], "V")

        if "M" in diagrams:
            draw_diagram(Ms, colors["M"], "M")

        # ----------------------------------------------------
        # Draw extrema markers
        # ----------------------------------------------------
        if show_extrema:
            ext = elem.ifc.get_extrema()

            def draw_ext_point(x, val, color, label):
                Xg, Yg = _rot_trans(x, val * sf * flip[label], x1, z1, c, s)
                # Marker
                d.append(draw.Circle(Sx(Xg), Sy(Yg), 4, fill=color))
                # Value
                offset = 6 if val* flip[label] >= 0 else -6
                anchor = 'start' if val* flip[label] >= 0 else 'end'
                baseline = 'hanging' if val* flip[label] >= 0 else 'auto'
                if math.fabs(val) > 1000:
                    val_str = f"{val:.2e}"
                    offset *= 1.5
                else:
                    val_str = f"{val:.2f}"
                d.append(draw.Text(val_str, 12,
                                   Sx(Xg)+offset, Sy(Yg) + offset,
                                   text_anchor=anchor,
                                    dominant_baseline=baseline,
                                   fill=color))

            if "N" in diagrams:
                draw_ext_point(*ext["N_max"], colors["N"], "N")
                draw_ext_point(*ext["N_min"], colors["N"], "N")

            if "V" in diagrams:
                draw_ext_point(*ext["V_max"], colors["V"], "V")
                draw_ext_point(*ext["V_min"], colors["V"], "V")

            if "M" in diagrams:
                draw_ext_point(*ext["M_max"], colors["M"], "M")
                draw_ext_point(*ext["M_min"], colors["M"], "M")
                
    # --------------------------------------------------------
    # Draw structure (elements)
    # --------------------------------------------------------
    for elem in domain.elements.values():
        n1 = domain.get_node(elem.nodes[0])
        n2 = domain.get_node(elem.nodes[1])

        d.append(draw.Line(Sx(n1.coords[0]), Sy(n1.coords[2]),
                           Sx(n2.coords[0]), Sy(n2.coords[2]),
                           stroke='black', stroke_width=2))

    # Save
    # d.save_svg(filename)
    return d
