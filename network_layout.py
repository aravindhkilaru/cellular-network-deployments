import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection

def make_grid(nx, ny, min_diam, n, crop_circ, rotate_deg, align_to_origin):
    """
    Computes the coordinates of hexagon centers for a hexagonal grid.
    """
    ratio = np.sqrt(3) / 2
    if n > 0:  # Overwrite nx, ny if n is provided
        ny = int(np.sqrt(n / ratio))
        nx = n // ny

    coord_x, coord_y = np.meshgrid(np.arange(nx), np.arange(ny), sparse=False, indexing='xy')
    coord_y = coord_y * ratio
    coord_x = coord_x.astype('float')
    coord_x[1::2, :] += 0.5  # Shift every other row
    coord_x = coord_x.reshape(-1, 1)
    coord_y = coord_y.reshape(-1, 1)

    coord_x *= min_diam  # Scale to ISD
    coord_y = coord_y.astype('float') * min_diam

    mid_x = (np.ceil(nx / 2) - 1) + 0.5 * (np.ceil(ny/2) % 2 == 0)
    mid_y = (np.ceil(ny / 2) - 1) * ratio
    mid_x *= min_diam
    mid_y *= min_diam

    if crop_circ > 0:
        rad = ((coord_x - mid_x)**2 + (coord_y - mid_y)**2)**0.5
        coord_x = coord_x[rad.flatten() <= crop_circ, :]
        coord_y = coord_y[rad.flatten() <= crop_circ, :]

    if not np.isclose(rotate_deg, 0):
        RotMatrix = np.array([[np.cos(np.deg2rad(rotate_deg)), np.sin(np.deg2rad(rotate_deg))],
                              [-np.sin(np.deg2rad(rotate_deg)), np.cos(np.deg2rad(rotate_deg))]])
        rot_locs = np.hstack((coord_x - mid_x, coord_y - mid_y)) @ RotMatrix.T
        coord_x, coord_y = np.hsplit(rot_locs + np.array([mid_x, mid_y]), 2)

    if align_to_origin:
        coord_x -= mid_x
        coord_y -= mid_y

    return coord_x, coord_y

def plot_single_lattice(coord_x, coord_y, face_color, edge_color, min_diam, plotting_gap, rotate_deg, h_ax=None):
    """
    Plots a hexagonal lattice on the given axes.
    """
    if face_color is None:
        face_color = (1, 1, 1, 0)  # Transparent face
    if edge_color is None:
        edge_color = 'k'
    if h_ax is None:
        h_fig = plt.figure(figsize=(5, 5))
        h_ax = h_fig.add_axes([0.1, 0.1, 0.8, 0.8])

    patches = []
    for curr_x, curr_y in zip(coord_x, coord_y):
        polygon = mpatches.RegularPolygon((curr_x, curr_y), numVertices=6,
                                          radius=min_diam / np.sqrt(3) * (1 - plotting_gap),
                                          orientation=np.deg2rad(-rotate_deg))
        patches.append(polygon)
    collection = PatchCollection(patches, edgecolor=edge_color, facecolor=face_color)
    h_ax.add_collection(collection)
    h_ax.set_aspect('equal')
    h_ax.axis([coord_x.min() - 1 * min_diam, coord_x.max() + 1 * min_diam,
               coord_y.min() - 1 * min_diam, coord_y.max() + 1 * min_diam])
    return h_ax

def generate_perfect_honeycomb(ISD, n_cells=19):
    """
    Generate a honeycomb hexagonal grid with given ISD, cropped to n_cells.
    """
    nx, ny = 5, 5  # Grid size sufficient for 19 cells
    coord_x, coord_y = make_grid(nx, ny, min_diam=ISD, n=0, crop_circ=0, rotate_deg=0, align_to_origin=True)

    # Compute distances from origin
    dist = np.sqrt(coord_x**2 + coord_y**2)
    sort_idx = np.argsort(dist.flatten())
    coord_x = coord_x[sort_idx][:n_cells]
    coord_y = coord_y[sort_idx][:n_cells]

    return coord_x, coord_y

def plot_honeycomb(coord_x, coord_y, ISD):
    """
    Plot the hexagonal grid with BS positions and sector arrows.
    """
    fig, ax = plt.subplots(figsize=(8, 8))
    # Plot hexagons
    plot_single_lattice(coord_x, coord_y, face_color=None, edge_color='k', min_diam=ISD,
                        plotting_gap=0.01, rotate_deg=0, h_ax=ax)
    # Add BS positions and sector arrows
    for idx, (x, y) in enumerate(zip(coord_x.flatten(), coord_y.flatten())):
        ax.plot(x, y, 'bo')  # BS position
        '''
        # ha - Horizontal Alignment - 'left', 'center', 'right'
        # va - Vertical Alignment - 'top', 'center', 'bottom', 'baseline'
        '''
        ax.text(x, y, str(idx), ha='left', va='top', color='blue', fontsize=10) 
        for sector in range(3):
            angle_deg = sector * 120 #120
            angle_rad = np.deg2rad(angle_deg)
            r = 0.25 * ISD
            text_x = x + r * np.cos(angle_rad)
            text_y = y + r * np.sin(angle_rad)
            sector_number = 3 * idx + sector
            ax.text(text_x, text_y, str(sector_number), ha='right', va='baseline', color='green', fontsize=8)
        for sector in range(3):
            angle_deg = sector * 120
            angle_rad = np.deg2rad(angle_deg)
            dx_arrow = 0.4 * ISD * np.cos(angle_rad)
            dy_arrow = 0.4 * ISD * np.sin(angle_rad)
            ax.arrow(x, y, dx_arrow, dy_arrow, head_width=ISD/12, color='r', linestyle='dashed')

    ax.set_xlim(coord_x.min() - ISD, coord_x.max() + ISD)
    ax.set_ylim(coord_y.min() - ISD, coord_y.max() + ISD)
    ax.set_xlabel('X Position [m]')
    ax.set_ylabel('Y Position [m]')
    ax.set_title('19-Cell Hexagonal Honeycomb Grid')
    plt.grid(False)
    plt.show()

# Example usage
ISD = 1200  # Inter-Site Distance
coord_x, coord_y = generate_perfect_honeycomb(ISD, n_cells=19)
plot_honeycomb(coord_x, coord_y, ISD)
