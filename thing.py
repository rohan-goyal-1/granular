import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

# Parameters
num_vertices = 11
side_length = 1.0
circle_radius = 0.5
spacing = 4.42  # Horizontal distance between centers of shapes

# Calculate radius of circumcircle of the polygon
circumradius = side_length / (2 * np.sin(np.pi / num_vertices))

def get_circle_centers(center_x, center_y):
    vertices = []
    for i in range(num_vertices):
        theta = 2 * np.pi * i / num_vertices
        x = center_x + circumradius * np.cos(theta)
        y = center_y + circumradius * np.sin(theta)
        vertices.append((x, y))
    return vertices

# Plot setup
fig, ax = plt.subplots()
ax.set_aspect('equal')
ax.set_facecolor('white')

# Draw three shapes side-by-side
for i in range(3):
    cx = i * spacing
    cy = 0.0
    centers = get_circle_centers(cx, cy)
    for x, y in centers:
        circle = Circle((x, y), radius=circle_radius, edgecolor='black', facecolor='white', linewidth=1.5)
        ax.add_patch(circle)

# Set limits and remove axes
total_width = 3 * spacing
padding = circumradius + circle_radius + 0.5
ax.set_xlim(-padding, total_width - spacing + padding)
ax.set_ylim(-padding, padding)
ax.axis('off')

plt.show()

