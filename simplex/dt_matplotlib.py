import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
from matplotlib.patches import FancyBboxPatch
from scipy.ndimage import gaussian_filter

SCALE = 2          # 1 = preview, 2 = print resolution
fig_size = (32, 22)
dpi = 300
# ---------------------------------------------------------------- chemistry
k1, k2 = 0.7, 0.3
S = 1.0  

P_bot = np.array([
    S * k1 / (k1 + k2),
    S * k2 / (k1 + k2),
    -S
])

P_top = np.array([0.0, 0.0, 0.0])

GREEN   = np.array([0.2, 0.5, 0.1])
MAGENTA = np.array([0.85, 0, 0.67])

# -------------------------------------------------------------------
# Finite green patch for dt[NO] + dt[NO2] + dt[RONO2] = 0
# clipped to x > 0, y > 0, z < 0
#
# In the box 0 <= x <= S, 0 <= y <= S, -S <= z <= 0,
# the plane x + y + z = 0 becomes the triangle:
#   (S, 0, -S), (0, S, -S), (0, 0, 0)
# -------------------------------------------------------------------
G0 = np.array([S, 0.0, -S])
G1 = np.array([0.0, S, -S])
G2 = np.array([0.0, 0.0, 0.0])

# intersection line between:
#   x + y + z = 0
#   k2*x - k1*y = 0
# direction = (k1, k2, -(k1+k2))

# Split the green triangle into two pieces
right_piece = [G2, G0, P_bot]
left_piece  = [G2, P_bot, G1]

# -------------------------------------------------------------------
# Finite magenta patch for k2*x - k1*y = 0, extended vertically in z
# clipped to x > 0, y > 0, z < 0
#
# Rectangle in the box:
#   (0,0,-S) -> (S,(k2/k1)S,-S) -> (S,(k2/k1)S,0) -> (0,0,0)
#
# One edge is exactly on the dt[NO] axis (x = y = 0).
# -------------------------------------------------------------------
mBL = np.array([0.0, 0.0, -S])                 # bottom-left
mTL = np.array([0.0, 0.0, 0.0])                 # top-left, on z-axis
mTR = np.array([S, (k2 / k1) * S, 0.0])         # top-right
mBR = np.array([S, (k2 / k1) * S, -S])          # bottom-right

# Split the magenta rectangle by the same intersection line
mag_left_piece = [mBL, P_bot, P_top]
mag_right_piece = [P_top, mTL, mTR, mBR, P_bot]

# ------------------------------------------------------------------- figure
fig = plt.figure(figsize=fig_size, dpi=dpi)
fig.patch.set_facecolor("white")
ax = fig.add_subplot(111, projection="3d")
ax.set_facecolor("white")
ax.set_proj_type("persp", focal_length=0.9)

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
    facecolors=[(*MAGENTA, 0.65)],
    edgecolors='none',
    linewidths=0,
    antialiaseds=True
)

mag_right_poly = Poly3DCollection(
    [mag_right_piece],
    facecolors=[(*MAGENTA, 0.65)],
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

# only the magenta boundary that does NOT lie on the z-axis
soft_edge_inward(
    [P_top, mTR, mBR, P_bot, P_top],
    MAGENTA,
    inside_point=(P_top + mTR + mBR + P_bot) / 4,
    shrink_step=0.022
)


# --------------------------------------------------------- line geometry

elev, azim = 18, 45
ax.view_init(elev=elev, azim=azim)
el, az = np.radians(elev), np.radians(azim)
eye = np.array([np.cos(el)*np.cos(az), np.cos(el)*np.sin(az), np.sin(el)])

line_top = np.array([0.0, 0.0, -0.66])
line_bot = np.array([
    (S * k1 / (k1 + k2))-0.14,
    (S * k2 / (k1 + k2))+0.05,
    -S-0.56
])

factor = 0.9

line_bot_short = line_top + factor*(line_bot-line_top)

M = 160
ts   = np.linspace(0, 1, M)
pts  = np.outer(1-ts, line_top) + np.outer(ts, line_bot)
shift = np.array([0.0, 0.0, 0.35])
pts = pts + shift


segs = np.stack([pts[:-1], pts[1:]], axis=1)
depth = (0.5*(pts[:-1] + pts[1:])) @ eye
dn = (depth - depth.min())/(np.ptp(depth) + 1e-9)     # 0 far .. 1 near
width_scale = 0.55 + 1.0*dn                           # taper toward viewer


def draw_line(
    target,
    p0,
    p1,
    base_w,
    color=(1,1,1),
    sparks=False,
    spark_s=60
):
    p0 = np.asarray(p0, dtype=float)
    p1 = np.asarray(p1, dtype=float)

    M = 160
    ts = np.linspace(0, 1, M)

    pts = np.outer(1-ts, p0) + np.outer(ts, p1)
    segs = np.stack([pts[:-1], pts[1:]], axis=1)

    depth = (0.5*(pts[:-1] + pts[1:])) @ eye
    dn = (depth - depth.min()) / (np.ptp(depth) + 1e-9)
    width_scale = 0.55 + 1.0 * dn

    rgba = (*color, 1.0)

    lc = Line3DCollection(
        segs,
        colors=[rgba] * len(segs),
        linewidths=base_w * width_scale,
        capstyle="round"
    )

    lc.set_zorder(50)
    target.add_collection3d(lc)

    if sparks:
        target.scatter(
            *np.stack([p0, p1]).T,
            s=spark_s,
            c=[color],
            depthshade=False,
            zorder=60
        )

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

soft_edge([G0, G1, G2, G0], GREEN, shrink_step=0.008)


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


    draw_cone_head(ax, base, end, head_length=head_length, radius=head_width * 1.6, n=28, color=color)



axis_len = {"x": 1.35, "y": 1.45, "z": 1.5}
draw_axis(ax, [0,0,0], [axis_len["x"],0,0], color="black", head_length=0.09, head_width=0.013)
draw_axis(ax, [0,0,0], [0,axis_len["y"],0], color="black", head_length=0.09, head_width=0.013)
draw_axis(ax, [0,0,0], [0,0,-axis_len["z"]], color="black", head_length=0.09, head_width=0.013)


ax.text(axis_len["x"]+0.32, 0, 0.02, "d$_t$[NO$_2$]", color="black", fontsize=34, ha="center")
ax.text(0, axis_len["y"]+0.42, 0.0, "d$_t$[RONO$_2$]", color="black", fontsize=34, ha="center")
ax.text(0, 0, -axis_len["z"]-0.16, "d$_t$[NO]", color="black", fontsize=34, ha="center")

ax.set_xlim(-1.5, 1.5); ax.set_ylim(-1.5, 1.5); ax.set_zlim(-1.5, 1.5)
ax.set_box_aspect((1.5, 1.5, 1.5)); ax.set_axis_off()

ax.text2D(
    0.17, 0.38,
    "Kinetic Invariant",
    transform=ax.transAxes,
    fontsize=24,
    color=MAGENTA
)

ax.text2D(
    0.17, 0.35,
    "k$_2$d$_t$[NO$_2$] - k$_1$d$_t$[RONO$_2$] = 0",
    transform=ax.transAxes,
    fontsize=24,
    color=MAGENTA
)

ax.text2D(
    0.64, 0.38,
    "Nitrogen Conservation",
    transform=ax.transAxes,
    fontsize=24,
    color=GREEN
)

ax.text2D(
    0.64, 0.35,
    "d$_t$[NO] + d$_t$[NO$_2$] + d$_t$[RONO$_2$] = 0",
    transform=ax.transAxes,
    fontsize=24,
    color=GREEN
)

ax.text2D(
    0.3, 0.24,
    "Accessible states",
    transform=ax.transAxes,
    fontsize=24,
    color="black"
)

# -------- curved arrow --------
t = np.linspace(0, 1, 40)

cx, cy, cz = 1.2, 0.78, -0.88
r = 0.15
phi = np.deg2rad(110)
theta = np.linspace(np.pi*1.1, np.pi*0.3, len(t)) + phi

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



fig.subplots_adjust(left=0.02, right=0.98, top=0.98, bottom=0.02)

# ============================================ game-laser line compositing
fig.canvas.draw()
buf = np.asarray(fig.canvas.buffer_rgba()).astype(np.float32)[..., :3] / 255.0

leg_left_a = np.array([0.9, 0.7, -1.22])
leg_left_b = np.array([0.91, 0.7, -1.12])
leg_right_b = np.array([0.86, 0.7, -1.12])


def line_buffer(base_w, color=(1,1,1), sparks=False, spark_s=60):
    """Render the line ALONE on black at identical camera -> RGB buffer."""
    lf = plt.figure(figsize=fig_size, dpi=300)
    lf.patch.set_facecolor("black")
    la = lf.add_subplot(111, projection="3d"); la.set_facecolor("black")
    la.set_proj_type("persp", focal_length=0.62); la.view_init(elev=elev, azim=azim)
    draw_line(la, line_top, line_bot_short, base_w, color=color, sparks=sparks, spark_s=spark_s)
    la.set_xlim(0, 1.5); la.set_ylim(0, 1.5); la.set_zlim(-1.4, 1.4)
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
out = screen(out, crisp)          

res = np.clip(out, 0, 1)
plt.close(fig)

GREEN_HEX = "#8fe23a"; MAG_HEX = "#f01fd0"; GREY = "#9a9a9a"
comp = plt.figure(figsize=fig_size, dpi=dpi)
comp.patch.set_facecolor("white")
ax3d = comp.add_axes([0.0, 0.0, 1.0, 1.0])
ax3d.imshow(res)
ax3d.axis("off")

cone_ax = comp.add_axes([0.0, 0.0, 1.0, 1.0], projection="3d")
cone_ax.patch.set_alpha(0.0)
cone_ax.set_proj_type("persp", focal_length=0.62)
cone_ax.view_init(elev=elev, azim=azim)
cone_ax.set_xlim(-1.5, 1.5)
cone_ax.set_ylim(-1.5, 1.5)
cone_ax.set_zlim(-1.5, 1.5)
cone_ax.set_box_aspect((1.5, 1.5, 1.5))
cone_ax.set_axis_off()
cone_ax.set_zorder(500)


cone_tip = np.array([1.01, 0.64, -0.77])
cone_back = np.array([0.98, 0.64, -0.67])


draw_cone_head(
    cone_ax,
    cone_back,
    cone_tip,
    head_length=0.08,
    radius=0.02,
    n=32,
    color="white"
)



comp.savefig("dt.png", facecolor="white", bbox_inches="tight", dpi=dpi, pad_inches=0)
comp.savefig("dt.pdf", format="pdf", bbox_inches="tight", dpi=dpi)
print("saved dt.pdf")
