import meep as mp
import numpy as np
import matplotlib.pyplot as plt

# === Load gain spectrum file ===
file_path = "/home/tbaba/ドキュメント/Tibercad/examples/Example_7/output/spectrum_spectrum_Vb_4.dat"
data = np.loadtxt(file_path, comments='#')
energy_eV = data[:, 0]
gain_py = data[:, 8]

# === Physical constants ===
h = 4.135667696e-15
c = 299792458
e = 1.602176634e-19
c_um_fs = 299.792458

# === Spectrum to wavelength [μm]
wavelength_um = (h * c) / (energy_eV * e) * 1e6

# === Pick gain at λ=1.55μm
target_lambda = 1.55
idx = np.argmin(np.abs(wavelength_um - target_lambda))
lambda_sel = wavelength_um[idx]
gain_sel = gain_py[idx]  # [1/cm]

# === Frequency setup
f0 = c_um_fs / lambda_sel
omega_0 = 2 * np.pi * f0
gamma = 0.05 * f0
n_real = 3.4
eps_inf = n_real ** 2
Im_eps = gain_sel * lambda_sel * n_real / (2 * np.pi)
sigma = Im_eps * gamma * omega_0

# === Gain medium (Lorentz model)
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

# === Simulation setup ===
dpml = 2.0
length = 50  # µm
cell = mp.Vector3(length)
geometry = [
    mp.Block(
        size=mp.Vector3(length - 2 * dpml),  # avoid PML overlap
        center=mp.Vector3(),
        material=gain_medium
    )
]

source_pos = -0.4 * (length - 2 * dpml)

sources = [
    mp.Source(
        mp.ContinuousSource(frequency=f0),
        component=mp.Ex,
        center=mp.Vector3(source_pos)
    )
]

sim = mp.Simulation(
    cell_size=cell,
    resolution=50,
    geometry=geometry,
    sources=sources,
    boundary_layers=[mp.PML(dpml)],
    dimensions=1
)

# === Run simulation ===
sim.add_dft_fields(
    fields=[mp.Ex],       # ← キーワード引数で指定！
    center=mp.Vector3(),  # ← これも必要
    size=mp.Vector3(46, 1, 1),
    fcen=f0,
    df=0,
    nfreq=1
)

sim.run(until=200)
ex_data = sim.get_dft_array(mp.Ex, 0)


#
x_vals = np.linspace(-23, 23, len(ex_data))

#
plt.figure(figsize=(10, 4))
plt.plot(x_vals, ex_data)
plt.xlabel("Position x [µm]")
plt.ylabel("Ex field")
plt.title("Electric Field Ex(x) at Final Time Step")
plt.grid(True)
plt.tight_layout()
plt.show()
