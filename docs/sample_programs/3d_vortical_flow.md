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
The sample program in `[root_dir]/src/apps/vortical_flow_3d.f90` simulates the temporal evolution of a vortex with an axial jet. The base vortex model is the \\( q \\)-vortex, with an initial velocity field defined as:

$$
\mathbf{V}(r,\phi,z, t=0) = \frac{1 - \exp(-r^2)}{r} \hat{\mathbf{e}}_{\phi} + \frac{\exp(-r^2)}{q}\hat{\mathbf{e}}_{z},
$$

where \\( q \\) is the swirling parameter, determining the relative strength of the swirl. In this example, \\( q = 1 \\). The vortex is embedded in a uniform background flow \\( \mathbf{U} = U_z \hat{\mathbf{e}}_z \\), resulting in a total velocity field \\( \mathbf{U} + \mathbf{V} \\). Here, \\( U_z \\) is set to \\( -0.5 \\), which induces vortex breakdown based on the analysis of linear absolute instability (Olendraru *et al.*, 1999)[^1].

The evolution of the vortex field \\( \mathbf{V}(r,\phi,z,t) \\) is governed by:

$$
\frac{\partial \mathbf{V}}{\partial t} = - \nabla \varphi - \left( \nabla \times \mathbf{V} \right) \times (\mathbf{U} + \mathbf{V}) + {\rm{VISC}} \cdot \nabla^2 \mathbf{V},
$$

$$
\nabla \cdot \mathbf{V} = 0,
$$

where \\( \varphi \\) is the specific energy. The viscosity constant \\( {\rm{VISC}} (\equiv {\rm{Re}}^{-1}) \\), in the current context, indicates the inverse Reynolds number of the flow. Using the Toroidal-Poloidal (TP) decomposition:

$$
\mathbf{V} = \nabla \times \left\{\psi \hat{\mathbf{e}}_z \right\} + \nabla \times \left[ \nabla \times \left\{ \chi \hat{\mathbf{e}}_z \right\} \right],
$$

the equations reduce to:

$$
\frac{\partial}{\partial t} \begin{bmatrix} \psi \\ \chi \end{bmatrix} = - \mathbb{P} \left[ \left( \nabla \times \mathbf{V} \right) \times (\mathbf{U} + \mathbf{V}) \right] + {\rm{VISC}} \cdot \nabla^2 \begin{bmatrix} \psi \\ \chi \end{bmatrix},
$$

where \\( \mathbb{P} \\) is the toroidal-poloidal projection operator. In other words, \\( \mathbb{P} \left[ \left( \nabla \times \mathbf{V} \right) \times (\mathbf{U} + \mathbf{V}) \right] \\) yields the toroidal and poloidal scalars of the vector field \\( \left( \nabla \times \mathbf{V} \right) \times (\mathbf{U} + \mathbf{V})\\).

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
          1.D-2           0.D0           1.D2              0          10000
#
!!! FIELD PROPERTY INFO !!!
# -------- VISC ----- HYPERPOW ---- HYPERVISC ----------------------------------
          1.D-4              4          1.D-8
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
  - Final time: \\( t_f = 100 \\) (10,000th step).
- Scalar field properties: \\( \rm{VISC} = 10^{-4} \\) (i.e.,  \\( {\rm{Re}} = 10^4 \\)) and \\( \rm{HYPERVISC} = 10^{-8} \\) with \\( \rm{HYPERPOW} = 4 \\).
- Field information is stored at every 100 steps in `[root_dir]/output/fld/`.
- No log information is stored.

Note that a non-zero *hyperviscous* dissipation term \\( - {\rm{HYPERVISC}} \cdot (-\nabla^2)^{({\rm{HYPERPOW}}/2)} \\), along with the *physical* dissipation term \\( {\rm{VISC}} \cdot \nabla^2 \\), is additionally taken into account. The hyperviscosity term does not represent physical dissipation but prevents numerical instability caused by high-frequency oscillations at small scales. It should remain minimal to avoid affecting the physical accuracy of the results.

## Animated Results

The program outputs 101 scalar field data sets from \\( t = t_0 = 0 \\) to \\( t = t_f = 100 \\) at intervals of \\( 1 \\). The following Python `matplotlib` code helps visualize the vortex breakdown process, using 3-dimensional iso-contours of vorticity magnitude (\\( \left\lVert\nabla \times \mathbf{V}\right\rVert \\)).

```python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.patches as mpatches
from matplotlib import cm
from matplotlib.colors import Normalize

# Load static data outside the animate function
NR = 32; NP = 48
dt = 1.E-3; dn = 1E2
coords_r = np.loadtxt('../output/dat/coords_r.dat', skiprows=1, max_rows=NR)
coords_p = np.loadtxt('../output/dat/coords_p.dat', skiprows=1, max_rows=NP)

# Initialize the figure and 3D axis once
fig = plt.figure(figsize=(6, 8), facecolor='none')
ax = fig.add_subplot(facecolor='none')
xlim = 5; ylim = 5
ax.set_xlim(-xlim, xlim)
ax.set_xticks(np.arange(-xlim,xlim,5))
ax.set_ylim(-ylim, ylim)
ax.set_yticks(np.arange(-ylim,ylim,5))
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_title('Scalar Transport (2D)')
ax.set_aspect('equal')

# Set up a color norm and color map for the colorbar
norm = Normalize(vmin=0, vmax=1)
cmap = cm.viridis
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])

# Add the static colorbar
cbar = fig.colorbar(sm, ax=ax, orientation='horizontal')
cbar.set_label('Scalar Value')

# Animation function
def animate(frame):
    nn = frame * dn
    f_val = np.loadtxt('../output/fld/sfld' + f"{nn:06n}" + '_PPP', skiprows=1, max_rows=NR)
    
    X_d, Y_d, Z_d = [], [], []
    for index, value in np.ndenumerate(f_val[:NR-1, :NP]):
        x = coords_r[index[0]] * np.cos(coords_p[index[1]])
        y = coords_r[index[0]] * np.sin(coords_p[index[1]])
        z = value
        
        if -xlim <= x <= xlim and -ylim <= y <= ylim:
            X_d.append(x)
            Y_d.append(y)
            Z_d.append(z)
    
    ax.clear()
    ax.set_xlim(-xlim, xlim)
    ax.set_xticks(np.linspace(-xlim,xlim,11))
    ax.set_ylim(-ylim, ylim)
    ax.set_yticks(np.linspace(-ylim,ylim,11))
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_title('2-Dimensional Scalar Transport (Contour View)')
    ax.set_aspect('equal')
   
    surface_patch = mpatches.Patch(color=plt.cm.viridis(0.))
    ax.tricontourf(X_d, Y_d, Z_d, cmap=plt.cm.viridis, vmin=0, vmax=1, levels=np.linspace(0,1,51), extend='both')
    ax.scatter([0.5], [0.5], c='r', s=20)
    
    ax.legend(handles=[surface_patch], labels=['t='+f"{nn*dt:10.3f}"])

    # Use tight layout to reduce whitespace around the plot
    fig.tight_layout()
    
# Run the animation
anim = animation.FuncAnimation(fig, animate, frames=range(0, 101), interval=50)
anim.save('animation.mp4', writer='ffmpeg', fps=30, dpi=150)
```

<video controls loop class="d-block mx-auto" style="width:100%; max-width:480px">
  <source src="{{ '/assets/videos/scalar_transport_result.mp4' | relative_url }}" type="video/mp4">
</video>

---
[^1]: Olendraru, C., Sellier, A., Rossi, M., & Huerre, P. (1999). Inviscid instability of the Batchelor vortex: Absolute-convective transition and spatial branches. *Physics of Fluids*, 11(7), 1805–1820. [https://doi.org/10.1063/1.870045](https://doi.org/10.1063/1.870045)