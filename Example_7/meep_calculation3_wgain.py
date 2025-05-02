import meep as mp
import numpy as np
import matplotlib.pyplot as plt
import gdstk
from shapely.geometry import Polygon, Point
from shapely.ops import unary_union

# --- Parameter Setup ---
pixel_size = 0.04
resolution = 10
wavelength = 1.55
fcen = 1 / wavelength
df = 0.0
dpml = 1.0
z_thickness = 0.22
padding_above = 2.0
padding_below = 2.0

sio2 = mp.Medium(index=1.44)

# === Load gain spectrum from TiberCAD ===
file_path = "/home/tbaba/ドキュメント/Tibercad/examples/Example_7/output/spectrum_spectrum_Vb_4.dat"
data = np.loadtxt(file_path, comments='#')
energy_eV = data[:, 0]
gain_py = data[:, 8]

# === Constants and target λ ===
h = 4.135667696e-15     # eV·s
c = 299792458           # m/s
e = 1.602176634e-19     # C
c_um_fs = 299.792458    # um/fs

wavelength_um = (h * c) / (energy_eV * e) * 1e6
idx = np.argmin(np.abs(wavelength_um - wavelength))
lambda_sel = wavelength_um[idx]
gain_sel = gain_py[idx]  # [1/cm]

f0 = c_um_fs / lambda_sel
omega_0 = 2 * np.pi * f0
gamma = 0.05 * f0
n_real = 3.4
eps_inf = n_real ** 2
Im_eps = gain_sel * lambda_sel * n_real / (2 * np.pi)
sigma = Im_eps * gamma * omega_0

gain_medium = mp.Medium(
    epsilon=eps_inf,
    E_susceptibilities=[
        mp.LorentzianSusceptibility(
            frequency=f0,
            gamma=gamma,
            sigma=sigma
        )
    ]
)

# --- Read and merge GDS polygons ---
gds_file = "straight_waveguide2.gds"
lib = gdstk.read_gds(gds_file)
cell = lib.top_level()[0]
polygons = []
target_layer = (1, 0)
for poly in cell.polygons:
    if (poly.layer, poly.datatype) == target_layer:
        polygons.append(Polygon(poly.points))
merged = unary_union(polygons)

bbox = merged.bounds
xrange = bbox[0], bbox[2]
yrange = bbox[1], bbox[3]
grid_x = int((xrange[1] - xrange[0]) / pixel_size)
grid_y = int((yrange[1] - yrange[0]) / pixel_size)

xx, yy = np.meshgrid(
    np.linspace(xrange[0], xrange[1], grid_x),
    np.linspace(yrange[0], yrange[1], grid_y)
)
bitmap = np.array([
    merged.contains(Point(x, y))
    for x, y in zip(xx.ravel(), yy.ravel())
]).astype(float).reshape((grid_y, grid_x))

# --- MaterialGrid with gain medium ---
grid_size = mp.Vector3(bitmap.shape[1], bitmap.shape[0])
material_grid = mp.MaterialGrid(
    grid_size,
    sio2,
    gain_medium,   # ← 利得媒質に差し替え
    weights=bitmap.T
)

# --- Domain and geometry ---
sx = (xrange[1] - xrange[0]) + 2 * dpml
sy = (yrange[1] - yrange[0]) + 2 * dpml
sz = z_thickness + padding_above + padding_below + 2 * dpml
cell_size = mp.Vector3(sx, sy, sz)

z_center = -0.5 * sz + dpml + padding_below + 0.5 * z_thickness
x_center = -0.5 * sx + dpml + 0.5 * (xrange[1] - xrange[0])
y_center = -0.5 * sy + dpml + 0.5 * (yrange[1] - yrange[0])

geometry = [mp.Block(
    size=mp.Vector3(xrange[1] - xrange[0], yrange[1] - yrange[0], z_thickness),
    center=mp.Vector3(x_center, y_center, z_center),
    material=material_grid
)]

# --- Source ---
src_x = -0.5 * sx + dpml + 0.1
sources = [mp.EigenModeSource(
    src=mp.ContinuousSource(frequency=fcen),
    center=mp.Vector3(src_x, 0, z_center),
    size=mp.Vector3(0, sy - 2 * dpml, z_thickness),
    direction=mp.X,
    eig_band=1,
    eig_parity=mp.NO_PARITY,
    eig_match_freq=True
)]

# --- Monitor ---
monitor_x = 0.5 * sx - dpml - 0.1
flux_region = mp.FluxRegion(
    center=mp.Vector3(monitor_x, 0, z_center),
    size=mp.Vector3(0, sy - 2 * dpml, z_thickness),
    direction=mp.X
)

input_flux_region = mp.FluxRegion(
    center=mp.Vector3(src_x + 0.05, 0, z_center),
    size=mp.Vector3(0, sy - 2 * dpml, z_thickness),
    direction=mp.X
)

pml_layers = [mp.PML(dpml)]

# --- Simulation ---
sim = mp.Simulation(
    cell_size=cell_size,
    boundary_layers=pml_layers,
    geometry=geometry,
    sources=sources,
    resolution=resolution,
    default_material=sio2
)

trans = sim.add_flux(fcen, 0, 1, flux_region)
refl = sim.add_flux(fcen, 0, 1, input_flux_region)

sim.run(until=500)

trans_flux = mp.get_fluxes(trans)[0]
input_flux = mp.get_fluxes(refl)[0]
T = trans_flux / input_flux if input_flux != 0 else 0
print(f"λ = {wavelength:.3f} μm → Transmission = {T:.4f}")

# --- Dielectric distribution ---
eps_data = sim.get_array(center=mp.Vector3(), size=cell_size, component=mp.Dielectric)
z_center_index = eps_data.shape[2] // 2
eps_xy = eps_data[:, :, z_center_index]
x = np.linspace(-0.5 * sx, 0.5 * sx, eps_data.shape[0])
y = np.linspace(-0.5 * sy, 0.5 * sy, eps_data.shape[1])
plt.figure(figsize=(10, 8))
plt.imshow(eps_xy.T, extent=[x[0], x[-1], y[0], y[-1]], origin='lower', cmap='viridis')
plt.colorbar(label='Dielectric Constant (ε)')
plt.title('ε Distribution in xy-plane (z center)')
plt.xlabel('x (μm)')
plt.ylabel('y (μm)')
plt.tight_layout()
plt.show()

# --- Ey field ---
Ey_data = sim.get_array(center=mp.Vector3(), size=cell_size, component=mp.Ey)
Ey_xy = Ey_data[:, :, z_center_index]
plt.figure(figsize=(10, 8))
plt.imshow(Ey_xy.T, extent=[x[0], x[-1], y[0], y[-1]], origin='lower', cmap='RdBu')
plt.colorbar(label='Ey Field')
plt.title('Ey Field Distribution in xy-plane (z center)')
plt.xlabel('x (μm)')
plt.ylabel('y (μm)')
plt.tight_layout()
plt.show()
