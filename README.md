# pulseprop — Nonlinear optical pulse propagation

 ultrafast pulse propagation

Time‑domain simulations for second‑order (χ²) nonlinear optics: second‑harmonic generation (SHG), sum‑frequency generation (SFG), cascaded tripling, general three‑wave mixing, and optical parametric amplification (OPA). Uses Sellmeier dispersion and integrates the coupled envelopes along the crystal.

---

## Install

Requires Python ≥ 3.9, NumPy, and SciPy, plus two external libraries:

```bash
pip install numpy scipy
pip install git+https://github.com/decoherer/wavedata
pip install git+https://github.com/decoherer/sellmeier
```

- **wavedata** supplies the `Wave` container and plotting helpers.
- **sellmeier** supplies refractive index models, group velocity, and QPM utilities.

---

## Quick start

```python
from pulseprop import triplerpulseprop

# 0.2 ps pump at 1555 nm, 1 mm crystal, SHG η1=30 %/W/cm², SFG η2=250 %/W/cm²
_ = triplerpulseprop(τ=.2, L=1, P1=1000,
                     η1=30, η2=250, λ1=1555,
                     sell='mglnridgewg', plot=1)
```

Run `python pulseprop.py` for more demos.

---

## Example

![Pump, SHG, and SFG powers vs time](sfgpulsepropagation.png)

Plot of pump, SHG, and SFG powers vs time at the output of a 1mm length MgLN ridge waveguide with interleaved QPM poling for both 1555nm SHG doubling and 1555nm+777.5nm SFG tripling given 0.2ps Gaussian FWHM input pulse with 1000W peak power at 1555nm. This figure is produced by:

```python
triplerpulseprop(τ=.2, L=1, P1=1000, η1=30, η2=250,
                 λ1=1555, sell='mglnridgewg', plot=1)
```

Dashed lines mark group‑delay times \(t=L/v_g(\lambda_j)\) for \(\lambda_1\), \(\lambda_2=\lambda_1/2\), and \(\lambda_3=\lambda_1/3\).

---

## Units and conventions

- Time \(t\): picoseconds.
- Length \(L\): millimeters.
- Power \(P\): watts. Peak power of a Gaussian envelope.
- Pulse width \(τ\): full width at half maximum (FWHM) of **power**.
- Wavelengths \(\lambda\): nanometers.
- Material model: `sell` selects a dataset key from `sellmeier`.
- Quasi‑phase matching: `Λ` in micrometers (µm). If provided, the code applies a QPM mismatch term \(k_g\).
- Nonlinear “efficiency” parameters \(η\) are in %/W/cm². Internally converted to a coupling constant that scales with effective \(d_\mathrm{eff}\), modal overlap, and normalization.

---

## API

### `singlepulseprop(τ, L, P, λ0=1550, sell='ktp', dt=0.005, plot=False)`

Linear propagation sanity check. Propagates a Gaussian at \(\lambda_0\) and a copy at \(\lambda_0/2\) through a dispersive medium. No nonlinear coupling. Returns the output field envelopes as `Wave` objects. Optional plot overlays group‑delay markers.

---

### `shgpulseprop(τ, L, P, η=0, λ0=1550, sell='ktp', Type='zzz', Λ=None, dt=0.005, d=None, rtol=1e-3, atol=1e-6, plot=False)`

Second‑harmonic generation of a single Gaussian pump. Two coupled envelopes: fundamental \(\lambda_1=\lambda_0\) and second harmonic \(\lambda_3=\lambda_0/2\).

- `η` (%/W/cm²): SHG coupling strength.
- `Λ` (µm): optional QPM period override.
- `d(z)`: optional complex longitudinal poling profile (defaults to 1).

Returns `A1, A3` as `Wave` objects at \(z=L\).

---

### `sfgpulseprop(τ, L, P1, P2, η=0, λ1=1550, λ2=775, sell='ktp', Type='zzz', Λ=None, dt=0.005, rtol=1e-3, atol=1e-6, plot=False)`

Sum‑frequency generation of two Gaussian inputs at \(\lambda_1\) and \(\lambda_2\), producing \(\lambda_3 = (1/\lambda_1 + 1/\lambda_2)^{-1}\). Three coupled envelopes with dispersion and optional QPM. Returns `A1, A2, A3` at \(z=L\). Plots include dashed lines at \(L/v_g(\lambda_j)\).

---

### `triplerpulseprop(τ, L, P1, η1=0, η2=0, λ1=1550, sell='ktp', Type='zzz', Λ1=None, Λ2=None, d1=None, d2=None, nres=40, rtol=1e-3, atol=1e-6, plot=False)`

Cascaded tripling: SHG \(\lambda_1	o \lambda_2=\lambda_1/2\) then SFG \(\lambda_1+\lambda_2	o \lambda_3=\lambda_1/3\). Independent SHG and SFG couplings `η1`, `η2`. Optional distinct QPM periods `Λ1`, `Λ2` and poling profiles `d1(z)`, `d2(z)`. Returns `A0, A1, A2, A3` where `A0` is the input pump and `A1..A3` are outputs at \(z=L\). The demo that creates the figure uses this routine.

---

### `threewavepulseprop(P1=0, P2=0.01, P3=5000, τ1=50, τ2=50, τ3=0.5, L=1, η=600, λ1=920, λ2=1170, sell='mglnridgewg', Type='zzz', Λ=None, d=None, dt=0.01, rtol=1e-3, atol=1e-6, plot=False)`

General three‑wave mixing engine. Covers SFG or DFG by sign conventions. You set peak powers and durations for each band \((\lambda_1, \lambda_2, \lambda_3)\) with \(\lambda_3=(1/\lambda_1 + 1/\lambda_2)^{-1}\). Includes dispersion, optional QPM, and a longitudinal `d(z)` profile. Returns `A1, A2, A3, e1, e2, e3` (outputs at \(z=L\) and the corresponding inputs).

---

### `opapulseprop(P3=5e3, P2=0.01, τ3=0.5, τ2=50, L=1, η=600, λ3=515, λ2=1170, sell='mglnridgewg', Type='zzz', Λ=None, d=None, dt=0.01, rtol=1e-3, atol=1e-6, plot=False)`

OPA convenience wrapper built on `threewavepulseprop`. You pass the pump (`P3, τ3, λ3`) and seed (`P2, τ2, λ2`); the idler wavelength is computed from \(1/\lambda_1 = 1/\lambda_3 - 1/\lambda_2\). Returns `A1, A2, A3, e2, e3`.

---

### `temporal2spectral(A, λ0, vsfreq=False, correctphaseramp=False, plot=False, logplot=True)`

Fourier transform utility mapping a time‑domain envelope `A(t)` to a spectrum `B(λ)` or `B(f)` with consistent energy scaling. Useful for verifying time–bandwidth products and visualizing spectra. If `vsfreq=True`, returns `B(f)`; otherwise returns `B(λ)`.

---

## Model sketch

For each band \(j \in \{1,2,3\}\):
\[\partial_z E_j(t) = \mathcal{D}_j[E_j](t) + S_j(E_1,E_2,E_3; z)\]
where \(\mathcal{D}_j\) is a frequency‑domain linear operator from the Sellmeier \(n(\lambda)\), and \(S_j\) are the χ² coupling terms with optional QPM factor \(e^{\pm i k_g z}\). Integration along \(z\) uses `scipy.integrate.solve_ivp` with FFT‑based application of \(\mathcal{D}_j\).

---

## Repository layout

- `pulseprop.py` — core functions above.
- `README.md` — this file.
- Example PNG in the repo root as shown above.
