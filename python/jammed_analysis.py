import numpy as np
import argparse
import os
import sys
import math
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Polygon, Rectangle
from system_reader import SystemReader

EPS = 1e-6
PHI_EPS = 1e-15

show_states = False
filename = ""

def compute_eigenvalues(matrix):
    eigenvalues = np.linalg.eigvalsh(matrix)
    return np.sort(eigenvalues)

def are_eigenvalues_close(a, b, eps=EPS):
    return a.shape == b.shape and np.allclose(a, b, atol=eps)

def is_phi_close(phi1, phi2, eps=PHI_EPS):
    return abs(phi1 - phi2) <= eps

def draw_state(ax, vertex_positions, phi, sigma, count):
    ax.set_xticks([])
    ax.set_yticks([])
    L = 1.0  # Assume unit box as L isn't explicitly stored in the new format
    ax.set_xlim(-0.25 * L, 1.25 * L)
    ax.set_ylim(-0.25 * L, 1.25 * L)
    ax.set_aspect('equal')

    box = Rectangle((0, 0), L, L, edgecolor='lightgray', facecolor='none', linewidth=1)
    ax.add_patch(box)

    translations = [(dx * L, dy * L) for dx in [-1, 0, 1] for dy in [-1, 0, 1]]
    npart, num_vertices, _ = vertex_positions.shape

    all_centers = []
    all_sigmas = []

    for p in range(npart):
        verts = vertex_positions[p, :, :2]
        for dx, dy in translations:
            shifted = verts + np.array([dx, dy])
            center = np.mean(shifted, axis=0)
            all_centers.append(shifted)
            all_sigmas.append(sigma)

            for v in shifted:
                circle = Circle(v, radius=sigma / 2, fc='#dddddd', alpha=0.3)
                ax.add_patch(circle)

            polygon = Polygon(shifted, closed=True, edgecolor='lightgray', facecolor='lightgray', alpha=1)
            ax.add_patch(polygon)

    # Draw contact lines
    for i in range(len(all_centers)):
        for j in range(i + 1, len(all_centers)):
            for c1 in all_centers[i]:
                for c2 in all_centers[j]:
                    dist = np.linalg.norm(c1 - c2)
                    if dist <= (all_sigmas[i] + all_sigmas[j]) / 2:
                        ax.plot([c1[0], c2[0]], [c1[1], c2[1]], color='gray', linewidth=0.75, alpha=1.0)

    ax.set_title(f"ϕ={phi:.3f}, p={count}", fontsize=8)

def main():
    reader = SystemReader(filename)
    unique_states = []
    total_states = 0
    states_info = []

    for system in reader.systems.values():
        for frame_idx, frame in enumerate(system.frames):
            phi = frame.phi
            sigma = frame.sigma
            vertex_positions = frame.vertices
            adj_m = frame.adj_contacts

            eigenvalues = compute_eigenvalues(adj_m)
            is_unique = True

            for state in unique_states:
                existing_eigen, _, existing_phi, _ = state
                if are_eigenvalues_close(eigenvalues, existing_eigen) and is_phi_close(phi, existing_phi):
                    state[3] += 1
                    is_unique = False
                    break

            if is_unique:
                unique_states.append([eigenvalues, (system.system_id, frame_idx), phi, 1])

            states_info.append((system.system_id, frame_idx))
            total_states += 1

    print("\nSummary:")
    print(f"Total states processed: {total_states}")
    print(f"Unique states found: {len(unique_states)}\n")

    for i, (eigen, (sys_id, frame_idx), phi, count) in enumerate(unique_states):
        print(f"[{i}] STATE system_{sys_id}/frame_{frame_idx} :\tϕ={phi:.3f}\tp={count / total_states:.4f}")

    # Distribution plot
    phi_vals = dict()
    for _, _, phi, freq in unique_states:
        phi_vals[phi] = phi_vals.get(phi, 0) + freq

    plt.bar(list(phi_vals.keys()), [v / total_states for v in phi_vals.values()],
            width=0.0, align='center', edgecolor='black')
    plt.xlabel('ϕ')
    plt.ylabel('Probability')
    plt.title('Distribution of jammed states across ϕ')
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()

    if show_states:
        n = len(unique_states)
        cols = min(n, 6)
        rows = (n + cols - 1) // cols
        fig, axs = plt.subplots(rows, cols, figsize=(2.5 * cols, 2.5 * rows))

        # Ensure axs is always a 1D array
        if isinstance(axs, np.ndarray):
            axs = axs.flatten()
        else:
            axs = np.array([axs])  # Single subplot case

        for ax in axs[n:]:
            ax.axis('off')

        for i, (_, (sys_id, frame_idx), phi, count) in enumerate(unique_states):
            frame = reader.systems[sys_id].frames[frame_idx]
            draw_state(axs[i], frame.vertices, frame.phi, frame.sigma, count / total_states)

        plt.tight_layout(rect=[0, 0, 1, 0.97])

    plt.show()

if __name__ == "__main__":
    try:
        parser = argparse.ArgumentParser()
        parser.add_argument("filename", help="HDF5 file to process")
        parser.add_argument("--show-states", action="store_true", help="Whether to display individual states")

        args = parser.parse_args()

        filename = args.filename
        show_states = args.show_states

        main()
    except Exception as e:
        print("Error:", e)
        sys.exit(1)
