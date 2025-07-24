import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Polygon, Rectangle
from matplotlib.widgets import Button, TextBox
import argparse
import sys
import math

from system_reader import SystemReader

class SystemDrawer:
    def __init__(self, file_path):
        self.decoder = SystemReader(file_path)
        self.systems = self.decoder.systems
        self.state_idx = 0
        self.frame_idx = 0

        self.fig, self.ax = plt.subplots(figsize=(8, 8))
        self.meta_ax = self.fig.add_axes([0.75, 0.5, 0.2, 0.4])
        self.meta_ax.axis('off')

        self.next_ax = self.fig.add_axes([0.8, 0.05, 0.1, 0.05])
        self.prev_ax = self.fig.add_axes([0.65, 0.05, 0.1, 0.05])
        self.next_button = Button(self.next_ax, 'Next')
        self.prev_button = Button(self.prev_ax, 'Previous')
        self.next_button.on_clicked(self.next_state)
        self.prev_button.on_clicked(self.prev_state)

        # Create input box to jump to a specific state
        self.state_input_ax = self.fig.add_axes([0.8, 0.15, 0.1, 0.05])
        self.state_input = TextBox(self.state_input_ax, 'Jump to State:', initial="0")
        self.state_input.on_submit(self.jump_to_state)

        self.load_state()

    def load_state(self):
        self.frames = self.systems[self.state_idx].frames
        self.frame_idx = 0

    def draw_frame(self):
        self.ax.clear()
        self.meta_ax.clear()
        self.meta_ax.axis('off')

        system = self.systems[self.state_idx]
        frame = self.frames[self.frame_idx]

        num_p = system.num_p
        num_v = system.num_v
        ndim = system.ndim

        sigma = frame.sigma
        vertices = frame.vertices
        adj_contacts = frame.adj_contacts
        phi = frame.phi
        N_c = int(np.sum(adj_contacts) // 2)  # Each contact counted twice

        self.ax.set_xlim(-0.5, 1.5)
        self.ax.set_ylim(-0.5, 1.5)
        self.ax.set_aspect('equal')
        self.ax.set_xticks([])
        self.ax.set_yticks([])

        boundary = Rectangle((0, 0), 1, 1, linewidth=1, edgecolor='gray', facecolor='none')
        self.ax.add_patch(boundary)

        translations = [(dx, dy) for dx in [-1, 0, 1] for dy in [-1, 0, 1]]
        all_shifted_verts = []

        for p in range(num_p):
            verts = vertices[p]

            for dx, dy in translations:
                shifted_verts = verts + np.array([dx, dy])
                all_shifted_verts.append(shifted_verts)

                for v in shifted_verts:
                    circle = Circle(v, radius=sigma / 2, fc='lightblue', alpha=0.5)
                    self.ax.add_patch(circle)

                polygon = Polygon(shifted_verts, closed=True, ec='gray', fc='gray')
                self.ax.add_patch(polygon)

        # Draw contacts (red dashed lines)
        n_total = len(all_shifted_verts)
        for i in range(n_total):
            for j in range(i + 1, n_total):
                for c1 in all_shifted_verts[i]:
                    for c2 in all_shifted_verts[j]:
                        dist = np.linalg.norm(c1 - c2)
                        if dist <= sigma:
                            self.ax.plot([c1[0], c2[0]], [c1[1], c2[1]],
                                         color='red', linewidth=1, alpha=0.5, linestyle='dashed')

        # Display metadata
        meta_text = (
            f"$N_p$: {num_p}\n"
            f"$N_v$: {num_v}\n"
            f"$N_c$: {N_c}\n"
            f"$\\phi$: {phi:.16f}\n"
        )
        self.meta_ax.text(0.5, 1, meta_text, fontsize=10, verticalalignment='top', horizontalalignment='center')

        self.ax.set_title(f"STATE: {self.state_idx}     FRAME: {self.frame_idx}")

    def draw(self):
        self.draw_frame()
        self.timer = self.fig.canvas.new_timer(interval=500)
        self.timer.add_callback(self.update)
        self.timer.start()
        plt.show()

    def update(self):
        self.frame_idx += 1
        if self.frame_idx >= len(self.frames):
            self.frame_idx -= 1
        self.draw_frame()
        self.fig.canvas.draw_idle()

    def next_state(self, event):
        self.state_idx = (self.state_idx + 1) % len(self.states)
        self.load_state()
        self.draw_frame()

    def prev_state(self, event):
        self.state_idx = (self.state_idx - 1) % len(self.states)
        self.load_state()
        self.draw_frame()

    def jump_to_state(self, text):
        try:
            target_state = int(text)
            if 0 <= target_state < len(self.states):
                self.state_idx = target_state
                self.load_state()
                self.draw_frame()
            else:
                print(f"Invalid state index. Please enter a number between 0 and {len(self.states) - 1}.")
        except ValueError:
            print("Invalid input. Please enter a valid state index.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Visualize particle simulation log from HDF5 file.")
    parser.add_argument("file_path", type=str, help="Path to the HDF5 file")
    args = parser.parse_args()

    try:
        drawer = SystemDrawer(args.file_path)
        drawer.draw()
    except Exception as e:
        print("Error:", e)
        sys.exit(1)
