import numpy as np
import drawsvg as draw
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
    if zmax == zmin:
        zmax += 1.0

    def sx(x):
        return margin + (x - xmin) / (xmax - xmin) * (width_px - 2*margin)

    def sy(z):
        #return height_px - (margin + (z - zmin) / (zmax - zmin) * (height_px - 2*margin))
        return (margin + (z - zmin) / (zmax - zmin) * (height_px - 2*margin))

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


def draw_roller_support(d, x, z, size=20):
    base = z + size

    # Triangle
    p = draw.Path(stroke='black', fill='white')
    p.M(x - size, base)
    p.L(x + size, base)
    p.L(x, z)
    p.Z()
    d.append(p)

    # Rollers
    d.append(draw.Circle(x - size/2, base + size/3, size/4, fill='white', stroke='black'))
    d.append(draw.Circle(x + size/2, base + size/3, size/4, fill='white', stroke='black'))


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
def plot_model_drawsvg(domain, filename="model.svg",
                       width_px=800, height_px=600, margin=40,
                       show_node_labels=True, show_element_labels=True):

    xs, zs = [], []
    for node in domain.nodes.values():
        xs.append(node.coords[0])
        zs.append(node.coords[2])

    sx, sy = _sx_sy_builder(xs, zs, width_px, height_px, margin)
    d = draw.Drawing(width_px, height_px, origin=(0, 0))

    # Elements
    for elem in domain.elements.values():
        n1 = domain.get_node(elem.nodes[0])
        n2 = domain.get_node(elem.nodes[1])
        x1, z1 = sx(n1.coords[0]), sy(n1.coords[2])
        x2, z2 = sx(n2.coords[0]), sy(n2.coords[2])

        d.append(draw.Line(x1, z1, x2, z2, stroke='black', stroke_width=2))

        if show_element_labels:
            xm, zm = (x1 + x2)/2, (z1 + z2)/2
            d.append(draw.Text(f"E{elem.label}", 14, xm, zm - 5,
                               center=True, fill='blue'))

        # Distributed loads
        for load in domain.get_element_loads(elem.label):
            if hasattr(load, "w"):
                draw_distributed_load(d, x1, z1, x2, z2, load.w)

    # Nodes, supports, nodal loads
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
            draw_roller_support(d, x, z)

        # Nodal loads
        if node.label in domain.nodal_loads:
            fx, fz, My = domain.nodal_loads[node.label]
            if fx != 0 or fz != 0:
                draw_force(d, x, z, fx, fz)
            if My != 0:
                draw_moment(d, x, z, My)

    #d.save_svg(filename)
    #d.set_pixel_scale(2) # Set number of pixels per geometry unit
    return d


# ------------------------------------------------------------
# Diagram plot (Polygon-free)
# ------------------------------------------------------------
def plot_model_diagram(domain, diagram_type, filename,
                       width_px=900, height_px=900, margin=40,
                       scale=1.0,
                       fill_color_pos="rgba(0,120,255,0.35)",
                       fill_color_neg="rgba(255,80,80,0.35)",
                       line_color="black",
                       show_peaks=True,
                       nseg=40):
    # ------------------------------------------------------------
    # Collect node coordinates for global scaling
    # ------------------------------------------------------------
    xs = []
    zs = []
    for node in domain.nodes.values():
        xs.append(node.coords[0])
        zs.append(node.coords[2])

    sx, sy = _sx_sy_builder(xs, zs, width_px, height_px, margin)
    dsvg = draw.Drawing(width_px, height_px, origin=(0, 0))

    # ------------------------------------------------------------
    # Determine global diagram scale
    # ------------------------------------------------------------
    all_vals = []

    for elem in domain.elements.values():
        if diagram_type == "N":
            vals = elem.compute_normal_force(nseg)["N"]
        elif diagram_type == "V":
            vals = elem.compute_shear_force(nseg)["V"]
        elif diagram_type == "M":
            vals = elem.compute_bending_moment(nseg)["M"]
        all_vals.extend(vals)

    all_vals = np.array(all_vals)
    vmax = np.max(np.abs(all_vals))
    if vmax == 0:
        vmax = 1.0

    diag_scale = scale * 0.15 * (width_px + height_px) / vmax

    # ------------------------------------------------------------
    # Draw each element’s filled diagram
    # ------------------------------------------------------------
    for elem in domain.elements.values():
        n1 = domain.get_node(elem.nodes[0])
        n2 = domain.get_node(elem.nodes[1])

        x1, z1 = n1.coords[0], n1.coords[2]
        x2, z2 = n2.coords[0], n2.coords[2]

        dx = x2 - x1
        dz = z2 - z1
        L = np.sqrt(dx*dx + dz*dz)

        tx, tz = dx/L, dz/L
        nx, nz = -tz, tx

        # get diagram values
        if diagram_type == "N":
            data = elem.compute_normal_force(nseg)
            vals = data["N"]
        elif diagram_type == "V":
            data = elem.compute_shear_force(nseg)
            vals = data["V"]
        elif diagram_type == "M":
            data = elem.compute_bending_moment(nseg)
            vals = data["M"]

        xs_local = data["x"]

        # --------------------------------------------------------
        # Build filled area path
        # --------------------------------------------------------
        # First: upper curve (diagram)
        p = draw.Path(stroke=line_color, fill='none', stroke_width=2)

        # Positive and negative areas must be filled separately
        # so we split into segments of same sign
        def sign(x):
            return 1 if x >= 0 else -1

        # Identify sign-change indices
        segments = []
        start = 0
        for i in range(1, len(vals)):
            if sign(vals[i]) != sign(vals[i-1]):
                segments.append((start, i))
                start = i
        segments.append((start, len(vals)-1))

        # --------------------------------------------------------
        # Draw each same-sign segment as a filled polygon
        # --------------------------------------------------------
        for (i0, i1) in segments:
            seg_vals = vals[i0:i1+1]
            seg_xs = xs_local[i0:i1+1]

            if len(seg_vals) < 2:
                continue

            sgn = sign(seg_vals[0])
            fill_color = fill_color_pos if sgn > 0 else fill_color_neg

            pseg = draw.Path(stroke=line_color, fill=fill_color, stroke_width=1)

            # Upper curve (diagram)
            for k in range(len(seg_vals)):
                s = seg_xs[k] / L
                px = x1 + s * dx
                pz = z1 + s * dz
                offset = seg_vals[k] * diag_scale
                dxo = px + nx * offset
                dzo = pz + nz * offset

                if k == 0:
                    pseg.M(sx(dxo), sy(dzo))
                else:
                    pseg.L(sx(dxo), sy(dzo))

            # Lower curve (zero line) back to start
            for k in reversed(range(len(seg_vals))):
                s = seg_xs[k] / L
                px = x1 + s * dx
                pz = z1 + s * dz
                pseg.L(sx(px), sy(pz))

            pseg.Z()
            dsvg.append(pseg)

        # --------------------------------------------------------
        # Draw diagram outline (centerline)
        # --------------------------------------------------------
        p = draw.Path(stroke=line_color, fill='none', stroke_width=2)
        for i in range(len(xs_local)):
            s = xs_local[i] / L
            px = x1 + s * dx
            pz = z1 + s * dz
            offset = vals[i] * diag_scale
            dxo = px + nx * offset
            dzo = pz + nz * offset

            if i == 0:
                p.M(sx(dxo), sy(dzo))
            else:
                p.L(sx(dxo), sy(dzo))
        dsvg.append(p)

        # --------------------------------------------------------
        # Peaks
        # --------------------------------------------------------
        if show_peaks:
            imax = np.argmax(vals)
            imin = np.argmin(vals)

            for idx, color in [(imax, 'red'), (imin, 'blue')]:
                s = xs_local[idx] / L
                px = x1 + s * dx
                pz = z1 + s * dz
                offset = vals[idx] * diag_scale
                dxo = px + nx * offset
                dzo = pz + nz * offset

                dsvg.append(draw.Circle(sx(dxo), sy(dzo), 4, fill=color))
                dsvg.append(draw.Text(f"{vals[idx]:.2f}",
                                      12,
                                      sx(dxo),
                                      sy(dzo) - 10,
                                      center=True,
                                      fill=color))

    # ------------------------------------------------------------
    # Title
    # ------------------------------------------------------------
    title_map = {"N": "Normal Force", "V": "Shear Force", "M": "Bending Moment"}
    dsvg.append(draw.Text(title_map[diagram_type], 18,
                          width_px/2, margin/2, center=True))

    #dsvg.save_svg(filename)
    return dsvg
