import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.patches import Circle, Rectangle
from matplotlib.collections import LineCollection
import matplotlib
import itertools
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
box_length = 5.96727
shifts = list(itertools.product([-1, 0, 1], [-1, 0, 1]))  # 9 periodic shifts

# === FIGURE SETUP ===
fig, ax = plt.subplots(figsize=(8, 8))
ax.set_aspect('equal')
# ax.set_xlim(-0.5, 1.5)
# ax.set_ylim(-0.5, 1.5)
ax.set_xlim(0 - box_length * 0.5, box_length + box_length * 0.5)
ax.set_ylim(0 - box_length * 0.5, box_length + box_length * 0.5)
ax.set_title("PBC Particle System with Correct Overlap Lines")

# === DRAW SIMULATION BOX ===
box = Rectangle((0, 0), box_length, box_length, fill=False, edgecolor='black', linewidth=2)
ax.add_patch(box)

# === COLORS ===
colors = sns.color_palette('muted', n_colors=num_particles)

# === CIRCLES ===
# Store all 9 image circles per vertex
circle_groups = []
for pi in range(num_particles):
    group = []
    for vi in range(num_verts):
        for dx, dy in shifts:
            circle = Circle((0, 0), radius=sigmas[0] / 2, facecolor=colors[pi], edgecolor='black', alpha=0.8)
            ax.add_patch(circle)
            group.append(circle)
    circle_groups.append(group)  # shape: (num_particles, 3 * 9)

# === CONTACT LINE COLLECTION ===
contact_lines = LineCollection([], colors='red', linestyles='dotted', linewidths=1.0)
ax.add_collection(contact_lines)

def init():
    for group in circle_groups:
        for circle in group:
            circle.center = (0, 0)
    contact_lines.set_segments([])
    return sum(circle_groups, []) + [contact_lines]

def animate(frame_idx):
    print(f"{frame_idx} / {len(frames)}")
    verts = frames[frame_idx]  # shape: (P, 3, 2)
    sigma_frame = sigmas[frame_idx]  # assume shape: (P,)
    segments = []

    # Step 1: Position and resize all periodic image circles
    all_image_positions = {}
    for pi in range(num_particles):
        radius = sigma_frame / 2
        for vi in range(num_verts):
            base = verts[pi, vi]
            for shift_idx, (dx, dy) in enumerate(shifts):
                shifted = base + np.array([dx * box_length, dy * box_length])
                circle_idx = vi * len(shifts) + shift_idx
                circle = circle_groups[pi][circle_idx]
                circle.center = shifted
                circle.set_radius(radius)
                all_image_positions[(pi, vi, shift_idx)] = shifted

    # Step 2: Find overlaps — compare central vertex (shift = (0, 0)) to all periodic images of *other* particles
    center_shift_index = shifts.index((0, 0))
    for pi1 in range(num_particles):
        for vi1 in range(num_verts):
            p1 = all_image_positions[(pi1, vi1, center_shift_index)]
            for pi2 in range(num_particles):
                if pi1 == pi2:
                    continue  # skip same particle
                for vi2 in range(num_verts):
                    for shift_idx, _ in enumerate(shifts):
                        p2 = all_image_positions[(pi2, vi2, shift_idx)]
                        if np.linalg.norm(p1 - p2) < sigma_frame:
                            segments.append([p1, p2])

    contact_lines.set_segments(segments)
    ax.set_title(f"Frame {frame_idx:04d}")
    return sum(circle_groups, []) + [contact_lines]

anim = animation.FuncAnimation(fig, animate, init_func=init, frames=len(frames),
                               interval=INTERVAL, blit=True, repeat=False)

if SAVE:
    anim.save(OUTFILE, fps=1000 // INTERVAL, dpi=200)
    print(f"Saved to {OUTFILE}")
else:
    plt.show()

