# Thin Airfoil Solver

A Python implementation of classical **linearized thin-airfoil theory** for NACA
4-digit cambered airfoils. Given a camber distribution and an angle of attack, the
solver computes the bound vortex sheet strength via a Glauert (Fourier) expansion,
integrates the induced velocity field with Biot–Savart, and visualizes the
resulting pressure coefficient and streamline pattern.


## Features

- Exact (symbolic) Fourier coefficients `A0..AN` for any NACA 4-digit camber line,
  derived analytically with SymPy rather than numerical quadrature.
- Circulation `Γ` and lift coefficient `c_l` from the classical Glauert results.
- Full induced velocity field `(u, v)` on an arbitrary grid via Biot–Savart
  integration of the vortex sheet.
- Pressure coefficient field `C_p`.
- Stream function `ψ` reconstructed from the velocity field for streamline
  visualization.

## Theory

The airfoil's mean camber line is linearized onto the chord (`y = 0`), and the
flow-tangency boundary condition is expanded as a Fourier series in the
transformed coordinate `x = (c/2)(1 − cos θ)`:

```
γ(θ) = 2 V∞ [ A0 (1 + cos θ)/sin θ + Σ An sin(nθ) ]
```

with

```
Γ = c V∞ [ π A0 + (π/2) A1 ]
c_l = π (2 A0 + A1)
```

This is the standard first-order approximation used throughout classical
aerodynamics (Glauert 1926). Because the boundary condition is enforced on the
flat chord line rather than on the true cambered surface, the linear theory is
most accurate for small camber and small angles of attack; see
[Limitations](#limitations) below.

## Project structure

```
config.py                 # All user-defined parameters (geometry, flow, grid)
fourier_coefficients.py   # Symbolic derivation of A0..AN (SymPy)
vortex.py                 # Vortex sheet strength, circulation, lift coefficient
velocity.py                # Induced/total velocity field, pressure coefficient
streamlines.py             # ODE-based streamline tracing (alternative to ψ contours)
visualize.py                # Plotting: Cp field, stream function heatmap
main.py                     # Driver script — runs the full pipeline end to end
```

## Installation

```bash
pip install numpy scipy sympy matplotlib
```

## Usage

1. Edit the parameters in `config.py`:
   - `M`, `P`, `C` — NACA camber, camber position, and chord length
   - `ALPHA_DEG`, `V_INF` — angle of attack and freestream speed
   - `N_FOURIER` — number of Fourier coefficients
   - `X_MIN/X_MAX/Y_MIN/Y_MAX`, `NX/NY` — plotting grid extent and resolution

2. Run the solver:

   ```bash
   python main.py
   ```

   This prints the Fourier coefficients, lift coefficient, and circulation,
   computes the velocity/pressure fields, and saves `flow_field.png` with two
   panels: the pressure coefficient field and the streamline pattern.

## Example output

```
Computing Fourier coefficients...
A0 = 0.082774, A1 = 0.081495, A2 = 0.013861
Lift coefficient: c_l = 0.776106
Circulation: Γ = 0.388053
```

## Limitations

- **The camber line is not an exact streamline.** The vortex sheet's boundary
  condition is linearized onto the chord line `y = 0`, not the true cambered
  surface. As a result, the plotted camber curve will not exactly coincide
  with a streamline of the computed field — the mismatch is small but most
  visible near the leading edge. This is an inherent property of first-order
  thin-airfoil theory, not a numerical error (see below for how it was
  verified). Enforcing tangency on the true camber curve would require a
  vortex panel method instead.
- **Leading-edge singularity.** `γ(θ)` has a `1/sin θ` singularity at `θ = 0`
  (the leading edge), which is why a masking radius is used when plotting the
  fields near `x/c = 0`.
- The Biot–Savart velocity integral (`scipy.integrate.quad`) is evaluated
  point-by-point on the grid, so increasing `NX`/`NY` significantly will slow
  down `main.py`.

