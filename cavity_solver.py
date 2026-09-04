import numpy as np
import matplotlib.pyplot as plt

# ==========================================================
# 1. GRID & SIMULATION PARAMETERS
# ==========================================================
Nx = 41            # Number of grid points along the X axis
Ny = 41            # Number of grid points along the Y axis
Re = 100.0         # Reynolds Number (Controls fluid thickness)
dt = 0.001         # Time step size (How far we jump forward per frame)
max_steps = 2000   # Total number of frames/loops to compute
U_lid = 1.0        # Velocity of the sliding top wall

# Calculate the physical distance between grid rows and columns
dx = 1.0 / (Nx - 1)
dy = 1.0 / (Ny - 1)
# ==========================================================
# 2. INITIALIZE ARRAYS WITH ZEROS
# ==========================================================
omega = np.zeros((Ny, Nx)) # Vorticity (How fast parts of the fluid spin)
psi = np.zeros((Ny, Nx))   # Stream Function (Tracks fluid lines)
u = np.zeros((Ny, Nx))     # Real velocity going Left-to-Right (X)
v = np.zeros((Ny, Nx))     # Real velocity going Up-and-Down (Y)

# --- Apply Boundary Condition ---
# Tell the very top row of our grid (index -1) to move at U_lid speed
u[-1, :] = U_lid 

print("🚀 Grid and variables initialized successfully!")
# ==========================================================
# 3. MAIN SOLVER LOOP (Time-Stepping)
# ==========================================================
for step in range(max_steps):
    # Keep a snapshot of the previous step's vorticity to calculate changes
    omega_old = omega.copy()
    
    # --- A. Wall Boundary Conditions (Thom's Formula) ---
    # These lines calculate fluid spin against the solid boundaries
    omega[0, :] = -2.0 * psi[1, :] / (dy**2)                             # Bottom wall
    omega[-1, :] = -2.0 * (psi[-2, :] - psi[-1, :]) / (dy**2) - 2.0 * U_lid / dy # Top moving wall
    omega[:, 0] = -2.0 * psi[:, 1] / (dx**2)                             # Left wall
    omega[:, -1] = -2.0 * psi[:, -2] / (dx**2)                            # Right wall

    # --- B. Solve Vorticity Transport (Interior Grid Points) ---
    for i in range(1, Ny - 1):
        for j in range(1, Nx - 1):
            # Convective terms: How the moving fluid carries rotation around
            u_grad_w = u[i, j] * (omega_old[i, j+1] - omega_old[i, j-1]) / (2.0 * dx)
            v_grad_w = v[i, j] * (omega_old[i+1, j] - omega_old[i-1, j]) / (2.0 * dy)
            
            # Diffusive terms: How the fluid friction smooths out rotation
            laplacian_w = (
                (omega_old[i, j+1] - 2.0 * omega_old[i, j] + omega_old[i, j-1]) / (dx**2) +
                (omega_old[i+1, j] - 2.0 * omega_old[i, j] + omega_old[i-1, j]) / (dy**2)
            )
            
            # March forward in time (Forward Euler integration)
            omega[i, j] = omega_old[i, j] + dt * (-u_grad_w - v_grad_w + (1.0 / Re) * laplacian_w)

    # --- C. Solve Stream Function (Enforce Continuity/Mass Conservation) ---
    # We run relaxation loops to smooth out fluid path lines instantly
    for relaxation in range(15):
        for i in range(1, Ny - 1):
            for j in range(1, Nx - 1):
                psi[i, j] = 0.25 * (psi[i, j+1] + psi[i, j-1] + psi[i+1, j] + psi[i-1, j] + (dx**2) * omega[i, j])

    # --- D. Update Velocity Fields (Calculate u and v from psi) ---
    for i in range(1, Ny - 1):
        for j in range(1, Nx - 1):
            u[i, j] = (psi[i+1, j] - psi[i-1, j]) / (2.0 * dy)
            v[i, j] = -(psi[i, j+1] - psi[i, j-1]) / (2.0 * dx)

    # Print out progress tracking updates every 500 steps
    if step % 500 == 0:
        print(f"   Iteration {step:4d} / {max_steps} computed successfully...")

print("\n🎉 Simulation complete! Processing data into visual plots...")
# ==========================================================
# 4. POST-PROCESSING & VISUALIZATION
# ==========================================================
# Create coordinate points from 0 to 1 for our X and Y axes
x = np.linspace(0, 1, Nx)
y = np.linspace(0, 1, Ny)
X, Y = np.meshgrid(x, y)

# Initialize a clean canvas for our plot
plt.figure(figsize=(6, 5))

# Generate a filled contour plot of the fluid pathways (Stream Function)
contour_plot = plt.contourf(X, Y, psi, levels=20, cmap='jet')
plt.colorbar(contour_plot, label='Stream Function ($\psi$)')

# Add descriptive titles and physical labels
plt.title(f'2D Lid-Driven Cavity Fluid Vortex (Re = {Re})', fontsize=12, fontweight='bold')
plt.xlabel('Cavity Width (X)', fontsize=10)
plt.ylabel('Cavity Height (Y)', fontsize=10)
plt.tight_layout()

# Save the plot directly into our folder so we can upload it to GitHub later!
plt.savefig('cavity_vortex_result.png', dpi=300)
print("💾 Visual plot cleanly saved as 'cavity_vortex_result.png'!")

# Display the interactive plot on your screen
plt.show()
