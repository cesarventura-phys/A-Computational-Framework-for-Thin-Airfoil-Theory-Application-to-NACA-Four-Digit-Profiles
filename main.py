"""
main.py
Driver script for the thin airfoil solver.
"""

import numpy as np
import matplotlib.pyplot as plt
from config import *
from fourier_coefficients import fourier_coefficients
from vortex import circulation, lift_coefficient
from velocity import compute_field_grid
from visualize import (
    plot_pressure_field,
    plot_streamfunction_heatmap,
    y_camber
)

# ------------------------------------------------------------
# 1. Compute Fourier coefficients
# ------------------------------------------------------------
print("Computing Fourier coefficients...")
A = fourier_coefficients(M, P, ALPHA, N_FOURIER)
print(f"A0 = {A[0]:.6f}, A1 = {A[1]:.6f}, A2 = {A[2]:.6f}")

# ------------------------------------------------------------
# 2. Compute lift and circulation
# ------------------------------------------------------------
cl = lift_coefficient(A)
Gamma = circulation(V_INF, C, A)
print(f"Lift coefficient: c_l = {cl:.6f}")
print(f"Circulation: Γ = {Gamma:.6f}")

# ------------------------------------------------------------
# 3. Generate grid
# ------------------------------------------------------------
print("Generating grid...")
X = np.linspace(X_MIN, X_MAX, NX)
Y = np.linspace(Y_MIN, Y_MAX, NY)
X_grid, Y_grid = np.meshgrid(X, Y)

# ------------------------------------------------------------
# 4. Compute velocity and pressure fields
# ------------------------------------------------------------
print("Computing velocity and pressure fields (this may take a moment)...")
field = compute_field_grid(X_grid, Y_grid, V_INF, ALPHA, C, A,
                           include_pressure=True, include_velocity=True,
                           epsrel=1e-6)

U = field['U']
V = field['V']
Cp = field['Cp']

print("Field computation complete.")

# ------------------------------------------------------------
# 5. Visualize results (1x2 layout with camber and chord lines)
# ------------------------------------------------------------
print("Creating visualizations...")

# Generate points for the camber line (only if cambered)
x_camber = np.linspace(0, C, 200)
if M > 0:
    y_camber_vals = y_camber(x_camber, C, M, P)
else:
    y_camber_vals = np.zeros_like(x_camber)   # symmetric airfoil: no camber

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# --- Left: Pressure field ---
plot_pressure_field(X_grid, Y_grid, Cp, ax=ax1)
ax1.set_title(r'Pressure Coefficient $C_p$')   # raw string to avoid escape warning
ax1.plot(x_camber, y_camber_vals, 'k-', linewidth=1.5)
ax1.plot([0, C], [0, 0], 'k--', alpha=0.4, linewidth=1)
ax1.set_xlim(X_MIN, X_MAX)
ax1.set_ylim(Y_MIN, Y_MAX)

# --- Right: Streamline heatmap ---
plot_streamfunction_heatmap(X_grid, Y_grid, U, V, ax=ax2, mask_radius=0.05)
ax2.set_title(r'Streamlines (Stream Function $\psi$)')   # raw string
ax2.plot(x_camber, y_camber_vals, 'k-', linewidth=1.5)
ax2.plot([0, C], [0, 0], 'k--', alpha=0.4, linewidth=1)
ax2.set_xlim(X_MIN, X_MAX)
ax2.set_ylim(Y_MIN, Y_MAX)

plt.tight_layout()
plt.savefig('flow_field.png', dpi=300, bbox_inches='tight')
plt.show()

print("Done!")
print(f"Results saved to flow_field.png")