import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
import itertools

# === Hexagonal Grid Utilities ===
def make_grid(nx, ny, min_diam, n, crop_circ, rotate_deg, align_to_origin):
    ratio = np.sqrt(3) / 2
    if n > 0:
        ny = int(np.sqrt(n / ratio))
        nx = n // ny
    coord_x, coord_y = np.meshgrid(np.arange(nx), np.arange(ny), sparse=False, indexing='xy')
    coord_y = coord_y * ratio
    coord_x = coord_x.astype('float')
    coord_x[1::2, :] += 0.5
    coord_x = coord_x.reshape(-1, 1)
    coord_y = coord_y.reshape(-1, 1)
    coord_x *= min_diam
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
    if face_color is None:
        face_color = (1, 1, 1, 0)
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
    nx, ny = 5, 5
    coord_x, coord_y = make_grid(nx, ny, min_diam=ISD, n=0, crop_circ=0, rotate_deg=0, align_to_origin=True)
    dist = np.sqrt(coord_x**2 + coord_y**2)
    sort_idx = np.argsort(dist.flatten())
    coord_x = coord_x[sort_idx][:n_cells]
    coord_y = coord_y[sort_idx][:n_cells]
    return coord_x, coord_y

def plot_honeycomb(coord_x, coord_y, ISD, ues=None):
    fig, ax = plt.subplots(figsize=(8, 8))
    plot_single_lattice(coord_x, coord_y, face_color=None, edge_color='k', min_diam=ISD,
                        plotting_gap=0.01, rotate_deg=0, h_ax=ax)
    for idx, (x, y) in enumerate(zip(coord_x.flatten(), coord_y.flatten())):
        ax.plot(x, y, 'bo')
        ax.text(x, y, str(idx), ha='left', va='top', color='blue', fontsize=10)
        for sector in range(3):
            angle_deg = sector * 120
            angle_rad = np.deg2rad(angle_deg)
            r = 0.25 * ISD
            text_x = x + r * np.cos(angle_rad)
            text_y = y + r * np.sin(angle_rad)
            sector_number = 3 * idx + sector
            ax.text(text_x, text_y, str(sector_number), ha='right', va='baseline', color='green', fontsize=8)
            dx_arrow = 0.4 * ISD * np.cos(angle_rad)
            dy_arrow = 0.4 * ISD * np.sin(angle_rad)
            ax.arrow(x, y, dx_arrow, dy_arrow, head_width=ISD/12, color='r', linestyle='dashed')
    if ues is not None:
        for u in ues:
            ax.plot(u.x, u.y, 'g.', markersize=3)
    ax.set_xlim(coord_x.min() - ISD, coord_x.max() + ISD)
    ax.set_ylim(coord_y.min() - ISD, coord_y.max() + ISD)
    ax.set_xlabel('X Position [m]')
    ax.set_ylabel('Y Position [m]')
    ax.set_title('4G Network Layout: Base Stations, Sectors, and UEs for Urban Micro (UMi)')
    plt.grid(False)
    plt.show()

# === LTE RSRP and RSRQ Calculation ===
def calculate_lte_rsrp_rsrq(RSSI, BW):
    BW_to_N = {1.4e6: 6, 3e6: 15, 5e6: 25, 10e6: 50, 15e6: 75, 20e6: 100}
    N = BW_to_N.get(BW, None)
    if N is None:
        raise ValueError("Unsupported Bandwidth")
    rsrp_lte = RSSI - 10 * np.log10(12 * N)
    rsrq_lte = 10 * np.log10(N) + rsrp_lte - RSSI
    return rsrp_lte, rsrq_lte

# === Entities ===
class BaseStation:
    bs_id = itertools.count()
    def __init__(self, x, y, height=25, tx_power=23):
        self.ID = next(BaseStation.bs_id)
        self.x, self.y, self.height, self.tx_power = x, y, height, tx_power
        self.sectors = [Sector(self.ID, tx_power, i * 120) for i in range(3)]

class Sector:
    sector_id = itertools.count()
    def __init__(self, bs_id, tx_power, orientation):
        self.ID = next(Sector.sector_id)
        self.bs_id, self.tx_power, self.orientation = bs_id, tx_power, orientation

class UE:
    id_iter = itertools.count()
    def __init__(self, x, y, sector_id):
        self.ID = next(UE.id_iter)
        self.x, self.y, self.sector_id = x, y, sector_id

# === Network ===
class Network:
    def __init__(self, ISD=500, seed=42):
        np.random.seed(seed)
        self.ISD = ISD
        self.BSs, self.UEs = [], []
        self.results = []

    def add_bs(self, x, y):
        self.BSs.append(BaseStation(x, y))

    def add_ue(self, x, y, sector_id):
        self.UEs.append(UE(x, y, sector_id))

    def compute_metrics(self):
        c = 3e8
        fc = 3.5e9
        h_bs, h_ue = 25, 1.5
        d_bp = 4 * (h_bs - 1.5) * (h_ue - 1.5) * fc / c
        calibration_offset = 30
        BW = 10e6
        self.results = []
        for ue in self.UEs:
            bs_index = ue.sector_id // 3
            serving_bs = self.BSs[bs_index]
            dx = ue.x - serving_bs.x
            dy = ue.y - serving_bs.y
            d2D = np.hypot(dx, dy)
            d3D = np.sqrt(dx**2 + dy**2 + (h_bs - h_ue)**2)
            if d2D <= d_bp:
                pl = 28 + 22 * np.log10(d3D) + 20 * np.log10(fc / 1e9)
            else:
                pl = 28 + 40 * np.log10(d3D) + 20 * np.log10(fc / 1e9) - 9 * np.log10(d_bp**2 + (h_bs-h_ue)**2)
            raw_rsrp = serving_bs.tx_power - pl
            RSSI = raw_rsrp + calibration_offset
            lte_rsrp, lte_rsrq = calculate_lte_rsrp_rsrq(RSSI, BW)
            self.results.append((d2D, pl, lte_rsrp, lte_rsrq))

    def compute_SINR_matrix(self):
        nBS = len(self.BSs)
        nSectors = nBS * 3
        nUE = len(self.UEs)
        self.RSRP_Matrix = np.zeros((nUE, nSectors))
        c = 3e8
        fc = 3.5e9
        h_bs, h_ue = 25, 1.5
        d_bp = 4 * (h_bs - 1.5) * (h_ue - 1.5) * fc / c
        calibration_offset = 30
        BW = 10e6
        for ue in self.UEs:
            ue_idx = ue.ID
            for sector_idx in range(nSectors):
                bs_idx = sector_idx // 3
                bs = self.BSs[bs_idx]
                dx = ue.x - bs.x
                dy = ue.y - bs.y
                d2D = np.hypot(dx, dy)
                d3D = np.sqrt(dx**2 + dy**2 + (h_bs - h_ue)**2)
                if d2D <= d_bp:
                    pl = 28 + 22 * np.log10(d3D) + 20 * np.log10(fc / 1e9)
                else:
                    pl = 28 + 40 * np.log10(d3D) + 20 * np.log10(fc / 1e9) - 9 * np.log10(d_bp**2 + (h_bs-h_ue)**2)
                raw_rsrp = bs.tx_power - pl
                RSSI_sec = raw_rsrp + calibration_offset
                lte_rsrp_sec, _ = calculate_lte_rsrp_rsrq(RSSI_sec, BW)
                self.RSRP_Matrix[ue_idx, sector_idx] = lte_rsrp_sec
        self.SINR_Matrix = np.zeros(nUE)
        self.Des_RSRP_Matrix = np.zeros(nUE)
        SEC_list = np.arange(nSectors)
        N_dBm = -95  # Noise power in dBm for 10 MHz bandwidth
        for ue in self.UEs:
            ue_idx = ue.ID
            serving_sector = ue.sector_id
            S_dBm = self.RSRP_Matrix[ue_idx, serving_sector]
            S_lin = 10**(S_dBm / 10)
            other_sectors = [sec for sec in SEC_list if sec != serving_sector]
            I_lin = sum(10**(self.RSRP_Matrix[ue_idx, sec] / 10) for sec in other_sectors)
            N_lin = 10**(N_dBm / 10)
            SINR_lin = S_lin / (I_lin + N_lin + 1e-12)
            self.SINR_Matrix[ue_idx] = 10 * np.log10(SINR_lin)
            self.Des_RSRP_Matrix[ue_idx] = S_dBm

# === Simulation Execution ===
ISD = 500
net = Network(ISD)
coord_x, coord_y = generate_perfect_honeycomb(ISD)
for x, y in zip(coord_x.flatten(), coord_y.flatten()):
    net.add_bs(x, y)

for sid in range(57):
    bx, by = coord_x[sid // 3], coord_y[sid // 3]
    theta = np.deg2rad((sid % 3) * 120)
    for _ in range(10):
        r = np.random.uniform(10, ISD / 2)
        ang = np.random.uniform(theta - 1, theta + 1)
        net.add_ue(bx + r * np.cos(ang), by + r * np.sin(ang), sid)

net.compute_metrics()
net.compute_SINR_matrix()

SINR_values = net.SINR_Matrix
des_rsrp_values = net.Des_RSRP_Matrix
res = np.array(net.results)

# === Plots ===
# Combined plot of base stations, sectors, and UEs
plot_honeycomb(coord_x, coord_y, ISD, ues=net.UEs)

plt.scatter(res[:, 0], SINR_values, s=10, c='blue', label="SINR")
plt.xlabel('Distance to BS (m)')
plt.ylabel('SINR (dB)')
plt.title('SINR vs Distance')
plt.grid(True)
plt.legend()
plt.show()

plt.scatter(res[:, 0], des_rsrp_values, s=10, c='green', label="RSRP")
plt.xlabel('Distance to BS (m)')
plt.ylabel('RSRP (dBm)')
plt.title('RSRP vs Distance')
plt.grid(True)
plt.legend()
plt.show()

plt.scatter(res[:, 0], res[:, 1], s=10, c='orange', label="Path Loss")
plt.xlabel('Distance to BS (m)')
plt.ylabel('Path Loss (dB)')
plt.title('Path Loss vs Distance')
plt.grid(True)
plt.legend()
plt.show()

spec_eff = np.log2(1 + 10**(SINR_values/10))
cdf_spec = np.sort(spec_eff)
p_spec = np.arange(1, len(cdf_spec)+1) / len(cdf_spec)
plt.plot(cdf_spec, p_spec, label="Spectral Efficiency CDF")
plt.xlabel('Spectral Efficiency (bps/Hz)')
plt.ylabel('CDF')
plt.title('Spectral Efficiency CDF')
plt.grid(True)
plt.legend()
plt.show()

countSINR, bins_countSINR = np.histogram(SINR_values, bins=100)
pdfSINR = countSINR / countSINR.sum()
cdfSINR = np.cumsum(pdfSINR)
plt.plot(bins_countSINR[1:], cdfSINR, color='blue', label="SINR CDF")
plt.xlabel("SINR (dB)")
plt.ylabel("CDF")
plt.title("CDF of SINR")
plt.grid(True)
plt.legend()
plt.show()

def sinr_to_cqi(sinr):
    if sinr < 0: return 1
    elif sinr < 5: return 2
    elif sinr < 10: return 4
    elif sinr < 15: return 7
    elif sinr < 20: return 10
    else: return 15

spectral_efficiency_table = {1: 0.15, 2: 0.3, 4: 0.6, 7: 1.0, 10: 2.0, 15: 3.5}
bandwidth_hz = 10e6
throughput_list = [spectral_efficiency_table.get(sinr_to_cqi(sinr), 0) * bandwidth_hz / 1e6 for sinr in SINR_values]
throughput_sorted = sorted(throughput_list)
cdf_throughput = [i / len(throughput_sorted) for i in range(1, len(throughput_sorted) + 1)]
plt.plot(throughput_sorted, cdf_throughput, color='red', label="Throughput CDF")
plt.xlabel("Throughput (Mbps)")
plt.ylabel("CDF")
plt.title("CDF of Throughput")
plt.grid(True)
plt.legend()
plt.show()

# === Final Metrics ===
area_km2 = ((coord_x.max() - coord_x.min() + 2 * ISD) * (coord_y.max() - coord_y.min() + 2 * ISD)) / 1e6
print(f"Connection Density: {len(net.UEs)/area_km2:.2f} UEs/km²")
print(f"Average Spectral Efficiency: {spec_eff.mean():.2f} bps/Hz")
print(f"Area Capacity: {spec_eff.mean() * len(net.UEs)/area_km2:.2f} bps/Hz/km²")
