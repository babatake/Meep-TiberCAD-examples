# Meep-TiberCAD-examples
example 1

🧪 Gain-assisted Waveguide Simulation with Meep + TiberCAD

This repository demonstrates how to perform optical FDTD simulations in Meep by integrating semiconductor optical gain data calculated from TiberCAD via a Lorentzian susceptibility model.
📌 Highlights

    ✅ GDS-based waveguide layout (loaded via gdstk)

    ✅ MaterialGrid based structure definition

    ✅ Gain modeling using TiberCAD spectra:

        Gain spectrum extracted from TiberCAD .dat file

        Converted to complex dielectric function using:
        Im(ε)=g(λ)⋅λ⋅n2π
        Im(ε)=2πg(λ)⋅λ⋅n​

        Mapped into Meep's Lorentzian susceptibility model
    
    ✅ Field visualization (ε-distribution and Ey field)

    ✅ Gain effect observed: Transmission > 1.0

📂 Project Structure

.
├── straight_waveguide2.gds     # GDS layout of the waveguide
├── meep_calculation3_wgain.py  # Main FDTD simulation with gain
├── spectrum_spectrum_Vb_4.dat  # Gain spectrum output from TiberCAD
└── README.md                   # ← You're reading this!

📈 Example Result

At λ = 1.55 µm, the transmission ratio was:

λ = 1.550 μm → Transmission = 1.0473

This confirms that amplification due to optical gain is working correctly.
📌 Requirements

    Python 3.9+

    Meep (Python bindings)

    gdstk, shapely, matplotlib, numpy

🧠 Notes

This workflow is useful for simulating active photonic devices like:

    SOAs (semiconductor optical amplifiers)

    Gain-guided waveguides

    Laser cavity sections

🔬 Lorentzian Susceptibility Model for Optical Gain

In order to simulate optical gain in semiconductors, we model the complex permittivity ε(ω) using the Lorentz oscillator model, which captures frequency-dependent dispersion and amplification.
📘 Mathematical Formulation

The permittivity ε(ω) is defined as:
ε(ω)=ε∞+σω02−ω2−iγω
ε(ω)=ε∞​+ω02​−ω2−iγωσ​

Where:

    ε∞: High-frequency background permittivity (e.g., n² of passive Si = 3.4²)

    ω₀: Resonance (gain peak) angular frequency

    γ: Damping term (related to gain bandwidth)

    σ: Oscillator strength (controls amplitude of gain or absorption)

This expression yields both real and imaginary parts of ε(ω). The imaginary part, Im(ε), represents either absorption (if negative) or gain (if positive).
🧪 From TiberCAD Gain Spectrum to Lorentz Parameters

TiberCAD provides a gain spectrum as a function of photon energy:

    Gain is extracted at a target wavelength (e.g. λ = 1.55 µm)

    Gain unit is [1/cm]; we convert it into dielectric function via:

Im(ε)=g(λ)⋅λ⋅n2π
Im(ε)=2πg(λ)⋅λ⋅n​

This value is then used to compute sigma:
σ=Im(ε)⋅γ⋅ω0
σ=Im(ε)⋅γ⋅ω0​

Which is passed into Meep’s Lorentzian model:

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

This setup causes Meep to simulate gain-enhanced wave propagation, where the electromagnetic energy is amplified as it travels through the active medium.
🔍 Remarks

    A positive Im(ε) leads to energy amplification (gain)

    The Lorentzian profile ensures frequency-selective behavior

    You can use multiple LorentzianSusceptibility terms to model broader or more complex gain spectra
