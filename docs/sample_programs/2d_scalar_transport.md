---
title: '2D - Scalar Transport'
parent: Sample Programs
nav_order: 3
---

# 2D - Scalar Transport
[**Click to view the program's source code**](https://github.com/UCBCFD/MLegS/blob/modern/src/apps/scalar_transport_2d.f90){:target="_blank"}

## Compilation and Run
```bash
#! bash
np=$(nproc)
make scalar_transport_2d
mpirun.openmpi -n $np --oversubscribe ./build/bin/scalar_transport_2d
# if using IntelMPI
# mpiexec -n $np ./build/bin/scalar_transport_2d
```

## Description
The sample program in `[root_dir]/src/apps/scalar_transport_2d.f90` demonstrates how to simulate a scalar transport process where a scalar quantity is passively transported under a given background convection. It solves for a scalar field \\( s(r,\phi,t) \\), with the initial condition \\( s(r,\phi, t=0) = s_0 \\) and a time-varying source term \\( f(r,\phi, t) \\), governed by the following convection-diffusion equation:

$$
\frac{\partial s}{\partial t} + \mathbf{u} \nabla s = {\rm{VISC}} \nabla^2 s + f
$$

In this equation:
- \\( \mathbf{u} \\) denotes the background *swirling* velocity field (Lamb-Oseen vortex), defined as: 

$$ \mathbf{u} = \hat{\mathbf{e}}_\phi \left( 1 - \exp(-r^2)\right) / r  $$

- \\( {\rm{VISC}} \\) is the diffusion constant.
- \\( f(r, \phi, t) \\) is a source term introducing scalar quantities into the domain.

The initial scalar distribution \\( s_0(r, \phi) \\) is defined as:

$$
s_0 (r, \phi) = 0,
$$

indicating that the domain is empty at the beginning. The time-varying source term \\( f(r, \phi, t) \\) is defined as:

$$
\displaylines{ 
f (r, \phi, t ) = \max \left(1 - \left[ \left\{4(r \cos \phi - 0.5)\right\}^8 + \left\{4(r \sin \phi - 0.5)\right\}^8 \right], 0\right) \cdot \frac{1 - \cos (\pi t)}{2},
}
$$

representing a periodic injection of scalar quantities in a patch centered at \\( (x, y) = (0.5, 0.5) \\) with a diameter of 0.25. The source term oscillates with a temporal period of 2, smoothly transitioning between injection and no injection over time.

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

The program outputs 101 scalar field data sets from \\( t = t_0 = 0 \\) to \\( t = t_f = 10 \\) at intervals of \\( 0.1 \\). The following Python `matplotlib` code helps visualize the transport process of the scalar field along the frozen Lamb-Oseen vortex over time:

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
[^1]: When a physical diffusion term is weak or absent, unlike in the current sample case, numerical solutions can become unstable. This is especially true in systems with advection (e.g., convection), as high-wavenumber motions fail to dissipate and instead accumulate at the highest wavenumbers due to computational discretization limits. Hyperdiffusion (or hyperviscosity) helps resolve these issues by selectively dissipating energy at high wavenumbers. Mathematically, it adds a term involving the following operator: \\( (-1)^{p/2} \nu_h \nabla^{p} \\), where \\( \nu_h \\) and \\( p \\) are determined by the parameters `HYPERVISC` and `HYPERPOW`, respectively.