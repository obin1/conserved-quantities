"""
Kinetic-invariant geometry for the RO2 + NO branching system (RONO2 / NO / NO2).

Coordinates:  x = [NO2],  y = [RONO2],  z = [NO] (vertical)
Invariants:
  (1) nitrogen conservation:  [NO2] + [RONO2] + [NO] = N_tot   -> green simplex (triangle)
  (2) kinetic invariant:      k2 [NO2] - k1 [RONO2] = L        -> magenta plane (parallel to NO axis)
Intersection = 1D chord of accessible states, drawn as a game-style laser.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
from matplotlib.patches import FancyBboxPatch
from scipy.ndimage import gaussian_filter

SCALE = 2          # 1 = preview, 2 = print resolution

# ---------------------------------------------------------------- chemistry
k1, k2, N, L = 0.7, 0.3, 1.0, 0.06

P_top = np.array([L/k2, 0.0, N - L/k2])                            # (0.20, 0.00, 0.80)
P_bot = np.array([(k1*N + L)/(k1+k2), (k2*N - L)/(k1+k2), 0.0])    # (0.76, 0.24, 0.00)
for P in (P_top, P_bot):
    assert abs(P.sum() - N) < 1e-12
    assert abs(k2*P[0] - k1*P[1] - L) < 1e-12
    assert (P >= -1e-12).all()

GREEN   = np.array([0.55, 0.85, 0.20])
MAGENTA = np.array([0.94, 0.12, 0.82])

# ------------------------------------------------- frosted / acrylic texture
def frosted_field(seed, sigma=13.0, res=128):     # sigma up -> translucent acrylic, not noise
    r = np.random.default_rng(seed).standard_normal((res, res))
    f = gaussian_filter(r, sigma=sigma, mode="reflect")
    f -= f.min(); f /= (f.max() + 1e-9)
    return f

GFIELD = frosted_field(11)
MFIELD = frosted_field(23)

def frost_color(base, t, lo=0.55, hi=1.28, alpha=0.50):
    b = lo + (hi - lo) * t
    rgb = np.clip(base * b, 0, 1)
    return (rgb[0], rgb[1], rgb[2], alpha)

# ------------------------------------------------------------- tessellation
A = np.array([N, 0, 0]); B = np.array([0, N, 0]); C = np.array([0, 0, N])

def tri_faces(n=11):
    faces, cols = [], []
    def bary(a, b):
        c = 1 - a - b
        return a*A + b*B + c*C
    for i in range(n):
        for j in range(n - i):
            p00 = bary(i/n, j/n); p10 = bary((i+1)/n, j/n); p01 = bary(i/n, (j+1)/n)
            faces.append([p00, p10, p01])
            cols.append(((i+0.33)/n, (j+0.33)/n))
            if j < n - i - 1:
                p11 = bary((i+1)/n, (j+1)/n)
                faces.append([p10, p11, p01])
                cols.append(((i+0.66)/n, (j+0.66)/n))
    return faces, np.array(cols)

anchor = np.array([L/k2, 0.0, 0.0])
dfloor = np.array([k1, k2, 0.0]); dfloor /= np.linalg.norm(dfloor)
zhat   = np.array([0, 0, 1.0])
A_MAX, H_MAX = 0.86, 1.06

def mag_faces(na=10, nh=8):
    faces, cols = [], []
    def pt(a, h):
        return anchor + a*A_MAX*dfloor + h*H_MAX*zhat
    for i in range(na):
        for j in range(nh):
            faces.append([pt(i/na, j/nh), pt((i+1)/na, j/nh),
                          pt((i+1)/na, (j+1)/nh), pt(i/na, (j+1)/nh)])
            cols.append(((i+0.5)/na, (j+0.5)/nh))
    return faces, np.array(cols)

gfaces, guv = tri_faces()
mfaces, muv = mag_faces()

gcolors = [frost_color(GREEN,   1, alpha=0.40) for (u, v) in guv]  # green denser
mcolors = [frost_color(MAGENTA, 1, alpha=0.30) for (u, v) in muv]

def edge_rgba(c, alpha, boost=1.0):
    rgb = tuple(np.clip(np.array(c[:3]) * boost, 0, 1))
    return (*rgb, alpha)

gedge = [edge_rgba(c, alpha=0.005,  boost=1.18) for c in gcolors]
medge = [edge_rgba(c, alpha=0.005, boost=1.22) for c in mcolors]



# ------------------------------------------------------------------- figure
fig = plt.figure(figsize=(9.6, 10.24), dpi=100*SCALE)
fig.patch.set_facecolor("black")
ax = fig.add_subplot(111, projection="3d")
ax.set_facecolor("black")
ax.set_proj_type("persp", focal_length=0.62)


# Split green simplex by the kinetic line
right_piece = [A, P_bot, P_top]              # k2*x - k1*y >= L
left_piece  = [B, C, P_top, P_bot]           # k2*x - k1*y <= L

# Split magenta plane by the same intersection line
mBL = anchor
mBR = anchor + A_MAX * dfloor
mTR = anchor + A_MAX * dfloor + H_MAX * zhat
mTL = anchor + H_MAX * zhat

mag_left_piece  = [mBL, P_bot, P_top]
mag_right_piece = [P_top, mTL, mTR, mBR, P_bot]

right_poly = Poly3DCollection(
    [right_piece],
    facecolors=[(*GREEN, 0.65)],
    edgecolors=[(*np.clip(GREEN * 1.35, 0, 1), 0.95)],
    linewidths=1.2,
    antialiaseds=True
)

left_poly = Poly3DCollection(
    [left_piece],
    facecolors=[(*GREEN, 0.65)],
    edgecolors=[(*np.clip(GREEN * 1.15, 0, 1), 0.65)],
    linewidths=1.0,
    antialiaseds=True
)

mag_left_poly = Poly3DCollection(
    [mag_left_piece],
    facecolors=[(*MAGENTA, 0.55)],
    edgecolors='none',
    linewidths=0,
    antialiaseds=True
)

mag_right_poly = Poly3DCollection(
    [mag_right_piece],
    facecolors=[(*MAGENTA, 0.55)],
    edgecolors='none',
    linewidths=0,
    antialiaseds=True
)

# Add fills
ax.add_collection3d(mag_left_poly)
ax.add_collection3d(right_poly)
ax.add_collection3d(mag_right_poly)
ax.add_collection3d(left_poly)

# zsort / zorder
mag_left_poly.set_zsort("min")
right_poly.set_zsort("min")
mag_right_poly.set_zsort("min")
left_poly.set_zsort("max")

mag_left_poly.set_zorder(80)
right_poly.set_zorder(90)
mag_right_poly.set_zorder(100)
left_poly.set_zorder(110)

# -------------------------------------------------------------
# Manual magenta outline: draw only visible edges
# Omit the edge you do NOT want inside the green overlap.
# -------------------------------------------------------------

def soft_edge_inward(points, color, *, inside_point=None, layers=16, shrink_step=0.008):
    """
    Draw a glowy outline by shrinking copies of the path toward an interior point.
    This makes the glow bleed inward instead of outward.

    points:
        Open or closed polyline, e.g. [P_top, mTL, mTR, mBR, P_bot]
        or [mBL, P_bot, P_top]
    inside_point:
        A point guaranteed to lie inside the plane piece.
        If None, the centroid of the points is used.
    layers:
        Number of glow layers.
    shrink_step:
        How much each outer layer is pulled inward.
    """
    pts = np.asarray(points, dtype=float)

    # If the path is explicitly closed, drop the duplicate end for scaling,
    # then close it again after shrinking.
    closed = np.allclose(pts[0], pts[-1])
    if closed:
        pts = pts[:-1]

    anchor = np.asarray(inside_point, dtype=float) if inside_point is not None else pts.mean(axis=0)

    glow = tuple(np.clip(np.array(color) * 1.35, 0, 1))

    # Draw widest / faintest first, then tighter / brighter layers on top
    for level in range(layers, 0, -1):
        scale = max(0.0, 1.0 - shrink_step * level)
        shrunk = anchor + (pts - anchor) * scale

        if closed:
            shrunk = np.vstack([shrunk, shrunk[0]])

        segs = np.stack([shrunk[:-1], shrunk[1:]], axis=1)

        # faint outer layers, stronger inner layers
        alpha = 0.007 * (layers - level + 1)

        lc = Line3DCollection(
            segs,
            colors=[(*glow, alpha)] * len(segs),
            linewidths=level,
            capstyle="round"
        )
        lc.set_zorder(18 + (layers - level))
        ax.add_collection3d(lc)

    # crisp core line on top
    core_pts = np.asarray(points, dtype=float)
    if np.allclose(core_pts[0], core_pts[-1]):
        core_pts = core_pts[:-1]
        core_pts = np.vstack([core_pts, core_pts[0]])
    core_segs = np.stack([core_pts[:-1], core_pts[1:]], axis=1)

    core = Line3DCollection(
        core_segs,
        colors=[(*color, 0.95)] * len(core_segs),
        linewidths=1.2,
        capstyle="round"
    )
    core.set_zorder(200)
    ax.add_collection3d(core)

# left piece: a V-shape inside the magenta plane
soft_edge_inward([mBL, P_bot, P_top], MAGENTA, inside_point=(mBL + P_bot + P_top) / 3)

# right piece: boundary path of the other magenta polygon
soft_edge_inward([P_top, mTL, mTR, mBR, P_bot], MAGENTA,
                 inside_point=(P_top + mTL + mTR + mBR + P_bot) / 5)

soft_edge_inward([mBL, mTL, P_top, P_bot, mBL], MAGENTA, shrink_step=0.018)


# --------------------------------------------------------- line geometry
elev, azim = 9, -22
ax.view_init(elev=elev, azim=azim)
el, az = np.radians(elev), np.radians(azim)
eye = np.array([np.cos(el)*np.cos(az), np.cos(el)*np.sin(az), np.sin(el)])

M = 160
ts   = np.linspace(0, 1, M)
pts  = np.outer(1-ts, P_top) + np.outer(ts, P_bot)
segs = np.stack([pts[:-1], pts[1:]], axis=1)
depth = (0.5*(pts[:-1] + pts[1:])) @ eye
dn = (depth - depth.min())/(np.ptp(depth) + 1e-9)     # 0 far .. 1 near
width_scale = 0.55 + 1.0*dn                           # taper toward viewer


def draw_line(target, base_w, color=(1,1,1), sparks=False, spark_s=60):
    rgba = (*color, 1.0)
    lc = Line3DCollection(
        segs,
        colors=[rgba]*len(segs),
        linewidths=base_w*width_scale,
        capstyle="round"
    )
    lc.set_zorder(50)
    target.add_collection3d(lc)
    if sparks:
        target.scatter(
            *np.stack([P_top, P_bot]).T,
            s=spark_s,
            c=[color],
            depthshade=False,
            zorder=60
        )
# NOTE: line is NOT drawn into the base scene; it is composited afterwards.

# ------------------------------------------------------------ soft edges

def soft_edge(loop, color, layers=16, shrink_step=0.014):
    loop = np.asarray(loop, dtype=float)

    # If the loop is already closed like [A, B, C, A], drop the repeated last point
    if np.allclose(loop[0], loop[-1]):
        loop = loop[:-1]

    centroid = loop.mean(axis=0)
    glow = tuple(np.clip(np.array(color) * 1.35, 0, 1))

    # Draw from the widest / faintest layer down to the narrowest / brightest layer
    for level in range(layers, 0, -1):
        # shrink inward toward the centroid
        scale = max(0.0, 1.0 - shrink_step * level)
        shrunk = centroid + (loop - centroid) * scale

        # close the polygon again
        closed = np.vstack([shrunk, shrunk[0]])
        e = np.stack([closed[:-1], closed[1:]], axis=1)

        # keep the wide layers faint and the inner layers stronger
        alpha = 0.007 * (layers - level + 1)

        lc = Line3DCollection(
            e,
            colors=[(*glow, alpha)] * len(e),
            linewidths=level,
            capstyle="round"
        )
        lc.set_zorder(18 + (layers - level))
        ax.add_collection3d(lc)

soft_edge([A, B, C, A], GREEN, shrink_step=0.008)


# ------------------------------------------------------------- floor grid
gl = []
ext, step = 1.5, 0.15
for x in np.arange(0, ext+1e-9, step):
    gl.append([[x, 0, 0], [x, ext, 0]]); gl.append([[0, x, 0], [ext, x, 0]])
grid = Line3DCollection(gl, colors=[(0.42, 0.5, 0.42, 0.15)]*len(gl), linewidths=0.7)
grid.set_zorder(5); ax.add_collection3d(grid)

# ------------------------------------------------------------------ axes
def draw_cone_head(ax, start, end, head_length=0.10, radius=0.03, n=24, color="black"):
    start = np.asarray(start, dtype=float)
    end   = np.asarray(end, dtype=float)

    d = end - start
    d /= np.linalg.norm(d)

    # choose a perpendicular basis
    ref = np.array([0.0, 0.0, 1.0])
    if abs(np.dot(ref, d)) > 0.95:
        ref = np.array([0.0, 1.0, 0.0])

    u = np.cross(d, ref)
    u /= np.linalg.norm(u)
    v = np.cross(d, u)

    # center of the cone's base
    base_center = end - head_length * d

    theta = np.linspace(0, 2*np.pi, n, endpoint=False)
    circle = np.array([
        base_center + radius * (np.cos(t) * u + np.sin(t) * v)
        for t in theta
    ])

    faces = []
    for i in range(n):
        faces.append([end, circle[i], circle[(i + 1) % n]])

    ax.add_collection3d(
        Poly3DCollection(
            faces,
            facecolor=color,
            edgecolor=color
        )
    )


def draw_axis(ax, start, end, color="black",
              lw=1.6, head_length=0.05, head_width=0.02):
    start = np.asarray(start, dtype=float)
    end   = np.asarray(end, dtype=float)

    # Shaft
    ax.plot(
        [start[0], end[0]],
        [start[1], end[1]],
        [start[2], end[2]],
        color=color,
        lw=lw
    )

    # Unit direction
    d = end - start
    d /= np.linalg.norm(d)

    # Choose a vector not parallel to d
    ref = np.array([0., 0., 1.])
    if abs(np.dot(ref, d)) > 0.95:
        ref = np.array([0., 1., 0.])

    # Perpendicular direction
    u = np.cross(d, ref)
    u /= np.linalg.norm(u)

    # Triangle vertices
    base = end - head_length * d

    draw_cone_head(ax, base, end, head_length=head_length, radius=head_width, n=28, color=color)


axis_len = {"x": 1.45, "y": 1.45, "z": 1.32}
draw_axis(ax, [0,0,0], [axis_len["x"],0,0], color="white", head_length=0.11, head_width=0.025)
draw_axis(ax, [0,0,0], [0,axis_len["y"],0], color="white", head_length=0.08, head_width=0.025)
draw_axis(ax, [0,0,0], [0,0,axis_len["z"]], color="white", head_length=0.09, head_width=0.025)

ax.text(axis_len["x"]+0.42, 0, 0.02, "[NO$_2$]", color="white", fontsize=28, ha="center")
ax.text(0, axis_len["y"]+0.30, 0.0, "[RONO$_2$]", color="white", fontsize=28, ha="center")
ax.text(0, 0, axis_len["z"]+0.06, "[NO]", color="white", fontsize=28, ha="center")

ax.set_xlim(0, 1.5); ax.set_ylim(0, 1.5); ax.set_zlim(0, 1.4)
ax.set_box_aspect((1.5, 1.5, 1.4)); ax.set_axis_off()

ax.text(
    0.50, 0.60, 0.45,
    "Nitrogen Conservation",
    color=GREEN,
    fontsize=15,
    fontweight="bold"
)
ax.text(
    0.50, 0.60, 0.38,
    "d$_t$[NO] + d$_t$[NO$_2$] + d$_t$[RONO$_2$] = 0",
    color=GREEN,
    fontsize=13
)

ax.text(
    0.90, 0.60, 0.08,
    "Kinetic Invariant",
    color=MAGENTA,
    fontsize=15,
    fontweight="bold"
)
ax.text(
    0.90, 0.60, 0,
    "k$_2$d$_t$[NO$_2$] - k$_1$d$_t$[RONO$_2$] = 0",
    color=MAGENTA,
    fontsize=13
)

ax.text(
    1.25, 0.10, 0,
    "Accessible states",
    color="white",
    fontsize=14,
    fontweight="bold"
)

# -------- curved arrow --------
t = np.linspace(0, 1, 40)

cx, cy, cz = 1.1, 0.18, 0
r = 0.08
theta = np.linspace(np.pi*0.3, np.pi*0.8, len(t))

x = cx + r*np.cos(theta)
y = cy - r*np.sin(theta)
z = cz + 0.03*t

def draw_curved_arrow_with_roll(ax, x, y, z, head_len=0.03, head_width=0.015, roll=0.0, color="black", lw=2):
    # shaft
    ax.plot(x, y, z, color=color, lw=lw)

    # tip direction from last segment
    p2 = np.array([x[-1], y[-1], z[-1]], dtype=float)
    p1 = np.array([x[-2], y[-2], z[-2]], dtype=float)
    d = p2 - p1
    d = d / (np.linalg.norm(d) + 1e-12)

    # make two perpendicular vectors to d
    a = np.array([0.0, 0.0, 1.0])
    if abs(np.dot(a, d)) > 0.9:
        a = np.array([0.0, 1.0, 0.0])

    u = np.cross(d, a)
    u = u / (np.linalg.norm(u) + 1e-12)
    v = np.cross(d, u)

    # roll the head around the shaft axis
    u2 = np.cos(roll) * u + np.sin(roll) * v

    # two wing endpoints
    left  = p2 - head_len * d + head_width * u2
    right = p2 - head_len * d - head_width * u2

    # draw head
    ax.plot([left[0], p2[0], right[0]], [left[1], p2[1], right[1]], [left[2], p2[2], right[2]],
            color=color, lw=lw)
    
draw_curved_arrow_with_roll(ax, x, y, z, roll=np.deg2rad(0))



fig.subplots_adjust(left=-0.02, right=1.02, top=1.04, bottom=-0.04)

# ============================================ game-laser line compositing
fig.canvas.draw()
buf = np.asarray(fig.canvas.buffer_rgba()).astype(np.float32)[..., :3] / 255.0

def line_buffer(base_w, color=(1,1,1), sparks=False, spark_s=60):
    """Render the line ALONE on black at identical camera -> RGB buffer."""
    lf = plt.figure(figsize=(9.6, 10.24), dpi=100*SCALE)
    lf.patch.set_facecolor("black")
    la = lf.add_subplot(111, projection="3d"); la.set_facecolor("black")
    la.set_proj_type("persp", focal_length=0.62); la.view_init(elev=elev, azim=azim)
    draw_line(la, base_w, color=color, sparks=sparks, spark_s=spark_s)
    la.set_xlim(0, 1.5); la.set_ylim(0, 1.5); la.set_zlim(0, 1.4)
    la.set_box_aspect((1.5, 1.5, 1.4)); la.set_axis_off()
    lf.subplots_adjust(left=-0.02, right=1.02, top=1.04, bottom=-0.04)
    lf.canvas.draw()
    b = np.asarray(lf.canvas.buffer_rgba()).astype(np.float32)[..., :3] / 255.0
    plt.close(lf)
    return b

def screen(a, b):
    return 1 - (1 - a) * (1 - np.clip(b, 0, 1))

WHITE   = (1.0, 1.0, 1.0)
MAGENTA = (0.94, 0.12, 0.82)


emissive = line_buffer(base_w=4.5, color=MAGENTA, sparks=False)   # glow source
crisp    = line_buffer(base_w=1.6, color=WHITE, sparks=False)    # sharp core


out = buf.copy()
for s, w in zip((2.5*SCALE, 8*SCALE, 20*SCALE), (1.1, 1.45, 1.5)):   # 2. blur  3. screen
    out = screen(out, np.stack([gaussian_filter(emissive[..., c], s) for c in range(3)], -1) * w)
out = screen(out, crisp)                                              # 4. crisp line on top
res = np.clip(out, 0, 1)
plt.close(fig)

GREEN_HEX = "#8fe23a"; MAG_HEX = "#f01fd0"; GREY = "#9a9a9a"
comp = plt.figure(figsize=(30.72, 20.48), dpi=1200)
comp.patch.set_facecolor("white")
ax3d = comp.add_axes([0.0, 0.0, 960/1536, 1.0]); ax3d.imshow(res); ax3d.axis("off")

ov = comp.add_axes([0, 0, 1, 1]); ov.set_xlim(0, 1536); ov.set_ylim(1024, 0)
ov.set_aspect("equal"); ov.axis("off")


comp.savefig("kinv_figure.png", facecolor="black", dpi=1200)
comp.savefig("kinv_figure.pdf", format="pdf", bbox_inches="tight", dpi=1200)
print("saved kinv_figure.png")
