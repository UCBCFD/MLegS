---
title: '2D - Scalar Diffusion'
parent: Sample Programs
nav_order: 2
---

# 2D - Scalar Diffusion
[**Click to view the program's source code**](https://github.com/UCBCFD/MLegS/blob/modern/src/apps/scalar_diffusion_2d.f90){:target="_blank"}

## Compilation and Run
```bash
#! bash
np=$(nproc)
make scalar_diffusion_2d
mpirun.openmpi -n $np --oversubscribe ./build/bin/barebone_template
# if using IntelMPI
# mpiexec -n $np ./build/bin/barebone_template
```

## Description
The sample program in `[root_dir]/src/apps/scalar_diffusion_2d.f90` demonstrates how to simulate a scalar diffusion process in a transient manner. It solves for a scalar field \\( s(r,\phi,t) \\), with the initial condition \\( s(r,\phi, t=0) = s_0 \\), governed by the following equation:

$$
\frac{\partial s}{\partial t} = \rm{VISC} \nabla^2 s
$$

In the sample program, the initial scalar distribution \\( s_0(r, \phi) \\) is defined as:

$$
\displaylines{ 
s_0 (r, \phi) = \max \left(1 - \left[ (r \cos \phi - 0.5)^2 + (r \sin \phi - 0.5)^2 \right], 0\right) \\ 
+ \max \left(0.5 - 0.5\left[ (r \cos \phi + 0.5)^2 + (r \sin \phi + 0.5)^2 \right], 0\right),
}
$$

representing two scalar clusters:
- A cluster centered at \\( (x, y) = (0.5, 0.5) \\) with a peak value of 1.
- A cluster centered at \\( (x, y) = (-0.5, -0.5) \\) with a peak value of 0.5.

Each cluster is approximately circular with a diameter of 1.

The default input parameters are defined in `[root_dir]/input_2d.params`:

```
!!! COMPUTATIONAL DOMAIN INFO !!!
# ---------- NR ----------- NP ----------- NZ ----------------------------------
             32             48              1    
# ------ NRCHOP ------- NPCHOP ------- NZCHOP ---------------------------------- 
             32             25              1    
# --------- ELL --------- ZLEN ------ ZLENxPI ---(IF ZLENxPI==T, ZLEN=ZLEN*PI)--
           1.D0           1.D0              F
#
!!! TIME STEPPING INFO !!!
# ---------- DT ----------- TI --------- TOTT ----------- NI --------- TOTN ----
          1.D-3           0.D0          1.D+1              0          10000
#
!!! FIELD PROPERTY INFO !!!
# -------- VISC ----- HYPERPOW ---- HYPERVISC ----------------------------------
          5.D-3              0           0.D0
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
              T            100
/* ------------------------------ END OF INPUT ------------------------------ */
```

From these parameters, the simulation is set as follows:
- The diffusion constant is \\( \rm{VISC} = 0.005 \\).
- Collocation points in physical space: \\( NR = 32 \\), \\( NP = 48 \\).
- Spectral elements in spectral space: \\( NRCHOP = 32 \\), \\( NPCHOP = 25 \\).
- Mapping parameter: \\( L = 1.0 \\).
- Time integration parameters:
  - Time step: \\( dt = 0.001 \\).
  - Initial time: \\( t_0 = 0 \\) (0th step).
  - Final time: \\( t_f = 10 \\) (10,000th step).
- Scalar field properties: \\( \rm{VISC} = 0.005 \\) and no hyperdiffusion[^1].
- Field information is stored at every 100 steps in `[root_dir]/output/fld/`.
- Log information is stored at every 100 steps in `[root_dir]/output/dat/`.

## Animated Results

The program outputs 101 scalar field data sets from \\( t = t_0 = 0 \\) to \\( t = t_f = 10 \\) at intervals of \\( 0.1 \\). The following Python `matplotlib` code helps visualize the diffusion process of the scalar field over time:

```python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.patches as mpatches
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.colors import Normalize

# Ensure the animation will display correctly in Jupyter
plt.rcParams["animation.html"] = "jshtml"
plt.rcParams["animation.embed_limit"] = 30  # Limit in MB
plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'

# Load static data outside the animate function
NR = 32; NP = 48
dt = 1.E-3; dn = 1E2
coords_r = np.loadtxt('../output/dat/coords_r.dat', skiprows=1, max_rows=NR)
coords_p = np.loadtxt('../output/dat/coords_p.dat', skiprows=1, max_rows=NP)

# Initialize the figure and 3D axis once
fig = plt.figure(figsize=(10, 6), facecolor='none')
ax = fig.add_subplot(projection='3d', facecolor='none')
xlim = 5; ylim = 5
ax.set_xlim(-xlim, xlim)
ax.set_ylim(-ylim, ylim)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('Scalar Diffusion (2D)')

# Set view to be perpendicular to the x-axis
ax.view_init(elev=30, azim=-45)

# Set up a color norm and color map for the colorbar
norm = Normalize(vmin=0, vmax=1)
cmap = cm.viridis
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])

# Add the static colorbar
cbar = fig.colorbar(sm, ax=ax, orientation='vertical')
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
    ax.set_ylim(-ylim, ylim)
    ax.set_zlim(0, 1)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('2-Dimensional Scalar Diffusion (3D View)')

    surface_patch = mpatches.Patch(color=plt.cm.viridis(0.))
    ax.plot_trisurf(X_d, Y_d, Z_d, cmap=plt.cm.viridis, antialiased=True, vmin=0, vmax=1)

    # Add the reference lines indicating the center
    ax.plot([0.5, 0.5], [0.5, 0.5], [0, 1], color='red', linewidth=2, label="Reference Line")
    ax.plot([-0.5, -0.5], [-0.5, -0.5], [0, 1], color='red', linewidth=2, label="Reference Line")
    
    ax.legend(handles=[surface_patch], labels=['t='+f"{nn*dt:10.3f}"])

    # Use tight layout to reduce whitespace around the plot
    fig.tight_layout()
    
# Run the animation
anim = animation.FuncAnimation(fig, animate, frames=range(0, 101), interval=50)
anim.save('animation.mp4', writer='ffmpeg', fps=30, dpi=150)
```

---
[^1]: In the absence of a physical diffusion term, numerical solutions often become unstable because high-wavenumber motions do not dissipate and instead accumulate at the highest wavenumbers due to the discretization limits of computation. Hyperdiffusion (or hyperviscosity) mitigates these numerical artifacts by selectively dissipating energy at high wavenumbers. Mathematically, it introduces an additional term: \\( \nu_h \nabla^{p} \\), where \\( \nu_h \\) and \\( p \\) are determined by the parameters `HYPERVISC` and `HYPERPOW`, respectively.
