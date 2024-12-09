---
title: '3D - Vortical Flow'
parent: Sample Programs
nav_order: 4
---

# 3D - Vortical Flow
[**Click to view the program's source code**](https://github.com/UCBCFD/MLegS/blob/modern/src/apps/vortical_flow_3d.f90){:target="_blank"}

## Compilation and Run
```bash
#! bash
np=$(nproc)
make vortical_flow_3d
mpirun.openmpi -n $np --oversubscribe ./build/bin/vortical_flow_3d
# If using IntelMPI
# mpiexec -n $np ./build/bin/vortical_flow_3d
```

## Description
The sample program in `[root_dir]/src/apps/vortical_flow_3d.f90` simulates the temporal evolution of a vortex with an axial jet. The base vortex model is the \\( q \\)-vortex, with a velocity field defined as:

$$
\mathbf{V}_q (r',\phi',z, t=0) = \frac{1 - \exp(-r'^2)}{r'} \hat{\mathbf{e}}_{\phi'} + \frac{\exp(-r'^2)}{q}\hat{\mathbf{e}}_{z},
$$

where \\( q \\) is the swirling parameter, determining the relative strength of the swirl. \\( r' \\) and \\( \phi' \\) represent *translated* polar variables centered at the \\( q \\)-vortex's center. In this example, \\( q = 1 \\). Initially, two \\( q \\)-vortices form a co-rotating vortex pair, with a distance of 4 between their centers \\( (x,y) = (\pm 2, 0) \\). The program assumes a uniform background flow \\( \mathbf{U} = U_z \hat{\mathbf{e}}_z \\), resulting in a total velocity field \\( \mathbf{U} + \mathbf{V} \\). Here, \\( U_z \\) is set to \\( -0.5 \\).

The evolution of the vortex field \\( \mathbf{V}(r,\phi,z,t) \\) is governed by:

$$
\frac{\partial \mathbf{V}}{\partial t} = - \nabla \varphi - \left( \nabla \times \mathbf{V} \right) \times (\mathbf{U} + \mathbf{V}) + {\rm{VISC}} \cdot \nabla^2 \mathbf{V},
$$

$$
\nabla \cdot \mathbf{V} = 0,
$$

where \\( \varphi \\) is the specific energy term. The viscosity constant \\( {\rm{VISC}} (\equiv {\rm{Re}}^{-1}) \\) represents the inverse Reynolds number of the flow. Using the Toroidal-Poloidal (TP) decomposition:

$$
\mathbf{V} = \nabla \times \left\{\psi \hat{\mathbf{e}}_z \right\} + \nabla \times \left[ \nabla \times \left\{ \chi \hat{\mathbf{e}}_z \right\} \right],
$$

the equations reduce to:

$$
\frac{\partial}{\partial t} \begin{bmatrix} \psi \\ \chi \end{bmatrix} = - \mathbb{P} \left[ \left( \nabla \times \mathbf{V} \right) \times (\mathbf{U} + \mathbf{V}) \right] + {\rm{VISC}} \cdot \nabla^2 \begin{bmatrix} \psi \\ \chi \end{bmatrix},
$$

where \\( \mathbb{P} \\) is the toroidal-poloidal projection operator. This operator extracts the toroidal and poloidal scalars from the vector field \\( \left( \nabla \times \mathbf{V} \right) \times (\mathbf{U} + \mathbf{V}) \\).

The default parameters are specified in `[root_dir]/input.params`:


```
!!! COMPUTATIONAL DOMAIN INFO !!!
# ---------- NR ----------- NP ----------- NZ ----------------------------------
            200            128            128    
# ------ NRCHOP ------- NPCHOP ------- NZCHOP ----------------------------------
            200             65             65    
# --------- ELL --------- ZLEN ------ ZLENxPI ---(IF ZLENxPI==T, ZLEN=ZLEN*PI)-- 
           4.D0           2.D0              T
#
!!! TIME STEPPING INFO !!!
# ---------- DT ----------- TI --------- TOTT ----------- NI --------- TOTN ----
          1.D-2           0.D0           5.D1              0           5000
#
!!! FIELD PROPERTY INFO !!!
# -------- VISC ----- HYPERPOW ---- HYPERVISC ----------------------------------
          1.D-4              8          5.D-7
#
!!! FIELD SAVE INFO !!!
# --------------------- FLDDIR -------------------------------------------------
                  ./output/fld
# ---- ISFLDSAV -- FLDSAVINTVL ---(IF ISFLDSAV!=T, FIELDS ARE NOT SAVED)--------
              T            100
#
!!! DATA LOGGING INFO !!!
# --------------------- LOGDIR -------------------------------------------------
                  ./output/dat
# ---- ISLOGSAV -- LOGSAVINTVL ---(IF ISLOGSAV!=T, LOGS ARE NOT GENERATED)------
              F              1
/* ------------------------------ END OF INPUT ------------------------------ */
```

From these parameters, the simulation is set as follows:
- Collocation points in physical space: \\( NR = 200 \\), \\( NP = 128 \\), \\( NZ = 128 \\).
- Spectral elements in spectral space: \\( NRCHOP = 200 \\), \\( NPCHOP = 65 \\), \\(NZCHOP = 65\\).
- Mapping parameter: \\( L = 4.0 \\) (one-half of \\(NR \\) are in \\( 0 \le r \le 4.0 \\)).
- Domain Length in \\( z \\): \\( 2 \pi \\).
- Time integration parameters:
  - Time step: \\( dt = 0.01 \\).
  - Initial time: \\( t_0 = 0 \\) (0th step).
  - Final time: \\( t_f = 50 \\) (5,000th step).
- Scalar field properties: \\( \rm{VISC} = 10^{-4} \\) (i.e.,  \\( {\rm{Re}} = 10^4 \\)) and \\( \rm{HYPERVISC} = 5 \times 10^{-7} \\) with \\( \rm{HYPERPOW} = 8 \\).
- Field information is stored at every 100 steps in `[root_dir]/output/fld/`.
- No log information is stored.

Note that a non-zero *hyperviscous* dissipation term \\( - {\rm{HYPERVISC}} \cdot (-\nabla^2)^{({\rm{HYPERPOW}}/2)} \\), along with the *physical* dissipation term \\( {\rm{VISC}} \cdot \nabla^2 \\), is additionally taken into account. The hyperviscosity term does not represent physical dissipation but prevents numerical instability caused by high-frequency oscillations at small scales. It should remain minimal to avoid affecting the physical accuracy of the results.

## Animated Results

The program outputs 51 field sets of vorticity magnitude (i.e., \\( \lVert \nabla \times \mathbf{V} \rVert \\)) from \\( t = t_0 = 0 \\) to \\( t = t_f = 50 \\) at intervals of \\( 1 \\). The following Python `matplotlib` code helps visualize the vortex merging process, using a 3-dimensional iso-contour (i.e., iso-surface) of vorticity magnitude.

```python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.patches as mpatches
from skimage.measure import marching_cubes as mc
from scipy.interpolate import RegularGridInterpolator

# Load static data outside the animate function
NR = 200; NP = 128; NZ = 128
dt = 1.E-2; dn = 1E2
coords_r = np.loadtxt('../output/dat/coords_r.dat', skiprows=1, max_rows=NR)
coords_p = np.loadtxt('../output/dat/coords_p.dat', skiprows=1, max_rows=NP+1)
coords_z = np.loadtxt('../output/dat/coords_z.dat', skiprows=1, max_rows=NZ+1)

# Initialize the figure and 3D axis once
fig = plt.figure(figsize=(10, 8), facecolor='none')
ax = fig.add_subplot(projection='3d', facecolor='none')
xlim = 5; ylim = 5; zlim = 2*np.pi
ax.set_xlabel("Horizontal (x)"); ax.set_ylabel("Vertical (y)"); ax.set_zlabel("Axial (z)")
ax.set_xlim(-xlim, xlim); ax.set_ylim(-ylim, ylim); ax.set_zlim(0, zlim)

# Animation function
def animate(frame):
    nn = frame * dn
    f_val = []
    for z in range(NZ):
        f_val.append(np.loadtxt('../output/fld/vormagfld' + f"{nn:06n}" + '_PPP', skiprows=z+1+(NR+8)*z, max_rows=NR))
    f_val.append(f_val[0])
    f_val = np.swapaxes(np.swapaxes(np.array(f_val), 0, 1), 1, 2)
    f_val = f_val[0:NR, 0:NP+1, 0:NZ+1]
    f_val[:,NP,:] = f_val[:,0,:] # phi of 2pi == phi of 0
    f_val[:,:,NZ] = f_val[:,:,0] # z of ZLEN == z of 0

    # Interpolation grid in Cartesian coordinates
    cartesian_x = np.linspace(-xlim, xlim, 100)  # Number of points can be adjusted
    cartesian_y = np.linspace(-ylim, ylim, 100)
    cartesian_z = np.linspace(0, zlim, 10)
    X, Y, Z = np.meshgrid(cartesian_x, cartesian_y, cartesian_z, indexing="ij")

    # Interpolation to create a smooth isosurface
    R = np.sqrt(X**2 + Y**2)
    P = np.arctan2(Y, X) % (2 * np.pi)  # Wrap angles to [0, 2*pi]
    interpolator = RegularGridInterpolator( (coords_r, coords_p, coords_z), 
        f_val, bounds_error=False, fill_value=0 )
    cylindrical_points = np.stack((R.ravel(), P.ravel(), Z.ravel()), axis=-1)
    f_cartesian = interpolator(cylindrical_points).reshape(X.shape)

    # Isosurface extraction
    iso_value = 1.0
    verts, faces, normals, values = mc(f_cartesian, level=iso_value)

    # Transform vertices directly in Cartesian coordinates
    verts_cartesian = verts * np.array(
        [2 * xlim / (len(cartesian_x) - 1), 2 * ylim / (len(cartesian_y) - 1), zlim / (len(cartesian_z) - 1)]
    ) + np.array([-xlim, -ylim, 0])

    # Add labels and adjust view
    ax.clear()
    ax.set_xlabel("Horizontal (x)")
    ax.set_ylabel("Vertical (y)")
    ax.set_zlabel("Axial (z)")
    ax.set_xlim(-xlim, xlim)
    ax.set_ylim(-ylim, ylim)
    ax.set_zlim(0, zlim)

    ax.set_title(f"Isosurface (Vorticity Magnitude={iso_value}) at t = {dt*nn}")

    # Plot the surface with shading and transparency
    ax.plot_trisurf( verts_cartesian[:, 0], verts_cartesian[:, 1],
        faces, verts_cartesian[:, 2],
        lw=.5, alpha=.9, shade=True )
    plt.tight_layout()

# Run the animation
anim = animation.FuncAnimation(fig, animate, frames=range(0, 51), interval=50)
anim.save('animation.mp4', writer='ffmpeg', fps=15, dpi=150)
```

<video controls loop class="d-block mx-auto" style="width:100%; max-width:480px">
  <source src="{{ '/assets/videos/vortical_flow_result.mp4' | relative_url }}" type="video/mp4">
</video>

### Bonus

Adding a slight perturbation in the axial direction can significantly alter the dynamics of vortex evolution. Consider superposing a sinusoidally perturbed axial jet initially located at the center, with a velocity field defined as:

$$
\mathbf{V} (r, \phi, z) = \exp(-r^2) \left( \frac{1 + \cos z}{2} \right) \hat{\mathbf{e}}_{z}.
$$

As a result, the vortex system undergoes more drastic evolution, with the emergence of pronounced axially varying flow structures. This scenario demonstrates that slight perturbations can significantly affect the flow dynamics, leading to complex vortical flow patterns.

<video controls loop class="d-block mx-auto" style="width:100%; max-width:480px">
  <source src="{{ '/assets/videos/vortical_flow_result_bonus.mp4' | relative_url }}" type="video/mp4">
</video>