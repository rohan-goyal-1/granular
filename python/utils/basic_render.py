import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.patches import Circle, Rectangle, Polygon
from matplotlib.collections import LineCollection
import seaborn as sns

# === CONFIG ===
FILE = "run.h5"
GROUP = "/frames"
INTERVAL = 5
SAVE = False
OUTFILE = "pbc_animation.mp4"

# === LOAD FRAMES AND SIGMA ===
with h5py.File(FILE, 'r') as f:
    frame_keys = sorted(f[GROUP].keys(), key=lambda x: int(x))
    sigmas = [f[GROUP][k].attrs['sigma'] for k in frame_keys]
    frames = [f[GROUP][k][()] for k in frame_keys]

# === BASIC INFO ===
num_particles, num_verts, _ = frames[0].shape
box_length = 9.13391852902415557

# === FIGURE SETUP ===
fig, ax = plt.subplots(figsize=(8, 8))
ax.set_aspect('equal')
xlims = (-box_length * 0.5, box_length * 1.5)
ylims = (-box_length * 0.5, box_length * 1.5)
ax.set_xlim(*xlims)
ax.set_ylim(*ylims)

# === COMPUTE NECESSARY PBC SHIFTS TO COVER THE SCREEN ===
n_shift_x = int(np.ceil((xlims[1] - xlims[0]) / box_length))
n_shift_y = int(np.ceil((ylims[1] - ylims[0]) / box_length))
range_x = range(-n_shift_x, n_shift_x + 1)
range_y = range(-n_shift_y, n_shift_y + 1)
shifts = list((dx, dy) for dx in range_x for dy in range_y)
num_shifts = len(shifts)

# === DRAW SIMULATION BOX ===
box = Rectangle((0, 0), box_length, box_length, fill=False, edgecolor='black', linewidth=1, alpha = 0.8)
ax.add_patch(box)

# === COLORS ===
colors = sns.color_palette('muted', n_colors=num_particles)

# === CIRCLES ===
circle_groups = []
for pi in range(num_particles):
    group = []
    for vi in range(num_verts):
        for _ in shifts:
            circle = Circle((0, 0), radius=sigmas[0] / 2, facecolor=colors[pi], alpha=1.0)
            ax.add_patch(circle)
            group.append(circle)
    circle_groups.append(group)

# === POLYGONS ===
polygon_groups = []
for pi in range(num_particles):
    group = []
    for _ in shifts:
        poly = Polygon(np.zeros((num_verts, 2)), closed=True,
                       facecolor=colors[pi], alpha=1.0)
        ax.add_patch(poly)
        group.append(poly)
    polygon_groups.append(group)

# === CONTACT LINE COLLECTION ===
contact_lines = LineCollection([], colors='red', linewidths=1.0)
ax.add_collection(contact_lines)

def init():
    for group in circle_groups:
        for circle in group:
            circle.center = (0, 0)
    for group in polygon_groups:
        for poly in group:
            poly.set_xy(np.zeros((num_verts, 2)))
    contact_lines.set_segments([])
    return sum(circle_groups, []) + sum(polygon_groups, []) + [contact_lines]

def animate(frame_idx):
    print(f"{frame_idx} / {len(frames)}")
    verts = frames[frame_idx]  # shape: (P, V, 2)
    sigma_frame = sigmas[frame_idx]
    segments = []

    all_image_positions = {}

    # Step 1: Circles and polygons
    for pi in range(num_particles):
        radius = sigma_frame / 2
        for vi in range(num_verts):
            base = verts[pi, vi]
            for shift_idx, (dx, dy) in enumerate(shifts):
                shifted = base + np.array([dx * box_length, dy * box_length])
                circle_idx = vi * num_shifts + shift_idx
                circle = circle_groups[pi][circle_idx]
                circle.center = shifted
                circle.set_radius(radius)
                all_image_positions[(pi, vi, shift_idx)] = shifted

        # Polygon per shift
        for shift_idx, (dx, dy) in enumerate(shifts):
            shift = np.array([dx * box_length, dy * box_length])
            shifted_poly = verts[pi] + shift
            polygon_groups[pi][shift_idx].set_xy(shifted_poly)

    # Step 2: Contacts between all particle images (excluding same particle)
    keys = list(all_image_positions.keys())
    contact_count = {pi: 0 for pi in range(num_particles)}  # NEW

    for i in range(len(keys)):
        pi1, vi1, s1 = keys[i]
        p1 = all_image_positions[keys[i]]
        for j in range(i + 1, len(keys)):
            pi2, vi2, s2 = keys[j]
            if pi1 == pi2:
                continue  # skip same particle, even across periodic images
            p2 = all_image_positions[keys[j]]
            if np.linalg.norm(p1 - p2) < sigma_frame:
                segments.append([p1, p2])
                contact_count[pi1] += 1
                contact_count[pi2] += 1

    # NEW: Print warning for undercoordinated particles
    for pi, count in contact_count.items():
        if count < 3:
            print(f"WARNING: Particle {pi} has only {count} contacts at frame {frame_idx}")

    contact_lines.set_segments(segments)
    return sum(circle_groups, []) + sum(polygon_groups, []) + [contact_lines]

anim = animation.FuncAnimation(fig, animate, init_func=init, frames=len(frames),
                               interval=INTERVAL, blit=True, repeat=False)

if SAVE:
    anim.save(OUTFILE, fps=1000 // INTERVAL, dpi=200)
    print(f"Saved to {OUTFILE}")
else:
    plt.show()

