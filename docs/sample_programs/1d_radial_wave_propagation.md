---
title: '1D - Radial Wave Propagation'
parent: Sample Programs
nav_order: 1
---
# 1D - Radial Wave Propagation
[**Click to view the program's source code**](https://github.com/UCBCFD/MLegS/blob/modern/src/apps/wave_propagation_1d.f90){:target="_blank"}

## Compilation and Run
For 1-dimensional problems, MLegS operates in serial mode only. This is because all scalar information (e.g., radial collocation points or spectral elements) is necessarily loaded on a single processor during the entire computation.

```bash
#! bash
make wave_propagation_1d
mpirun.openmpi -n 1 ./build/bin/wave_propagation_1d
# if using IntelMPI
# mpiexec -n 1 ./build/bin/wave_propagation_1d
```

## Description
The sample program in `[root_dir]/src/apps/wave_propagation_1d.f90` demonstrates axisymmetric wave propagation in the radial direction. It solves for a wave scalar function \\( s(r,t) \\), with initial conditions \\( s(r, t=0) = s_0 \\) and \\( \dot{s} (r, t=0) = \dot{s}_0 \\), governed by the radial wave propagation equation of the unit wave propagation speed:

$$
\frac{\partial^2 s}{\partial t^2} = \nabla^2 s = \frac{1}{r} \frac{\partial}{\partial r} \left( r \frac{\partial s}{\partial r} \right),
$$

where the initial conditions are:

$$
\displaylines{ 
s_0 (r) = \exp \left[ - \left(\frac{r}{4}\right)^2 \right], \\
\dot{s}_0 (r) = 0.
}
$$

This setup represents a Gaussian pulse propagation, where the wave expands outward radially.

The default input parameters are defined in `[root_dir]/input_1d.params`:

```
!!! COMPUTATIONAL DOMAIN INFO !!!
# ---------- NR ----------- NP ----------- NZ ----------------------------------
            192              1              1    
# ------ NRCHOP ------- NPCHOP ------- NZCHOP ----------------------------------
            192              1              1    
# --------- ELL --------- ZLEN ------ ZLENxPI ---(IF ZLENxPI==T, ZLEN=ZLEN*PI)--
           1.D1           1.D0              F
#
!!! TIME STEPPING INFO !!!
# ---------- DT ----------- TI --------- TOTT ----------- NI --------- TOTN ----
          1.D-4           0.D0          5.D+1              0         500000
#
!!! FIELD PROPERTY INFO !!!
# -------- VISC ----- HYPERPOW ---- HYPERVISC ----------------------------------
           0.D0              0           0.D0
#
!!! FIELD SAVE INFO !!!
# --------------------- FLDDIR -------------------------------------------------
                  ./output/fld
# ---- ISFLDSAV -- FLDSAVINTVL ---(IF ISFLDSAV!=T, FIELDS ARE NOT SAVED)--------
              T           1000
#
!!! DATA LOGGING INFO !!!
# --------------------- LOGDIR -------------------------------------------------
                  ./output/dat
# ---- ISLOGSAV -- LOGSAVINTVL ---(IF ISLOGSAV!=T, LOGS ARE NOT GENERATED)------
              T            100
/* ------------------------------ END OF INPUT ------------------------------ */
```

From these parameters, the simulation is set as follows:
- Collocation points in physical space: \\( NR = 192 \\).
- Spectral elements in spectral space: \\( NRCHOP = 192 \\).
- Mapping parameter: \\( L = 10.0 \\) (one-half of \\(NR \\) are in \\( 0 \le r \le 10.0 \\)).
- Time integration parameters:
  - Time step: \\( dt = 0.0001 \\).
  - Initial time: \\( t_0 = 0 \\) (0th step).
  - Final time: \\( t_f = 50 \\) (500,000th step).
- Scalar field properties: no diffusion (this program does not simulate a diffusion process).
- Field information is stored at every 1000 steps in `[root_dir]/output/fld/`.
- Log information is stored at every 100 steps in `[root_dir]/output/dat/`.

## Animated Results

The program outputs 501 scalar field data sets from \\( t = t_0 = 0 \\) to \\( t = t_f = 50 \\) at intervals of \\( 0.1 \\). The following Python `matplotlib` code helps visualize the diffusion process of the scalar field over time:

```python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Load static data outside the animate function
NR = 192; NP = 1
dt = 1.E-4; dn = 1E3
coords_r = np.loadtxt('../output/dat/coords_r.dat', skiprows=1, max_rows=NR)

# Initialize the figure
fig = plt.figure(figsize=(6, 6), facecolor='none')
ax = fig.add_subplot(facecolor='none')
xlim = 50
ax.set_xlim(0, xlim)
ax.set_xticks(np.linspace(0,xlim,11))
ax.set_ylim(-1.1, 1.1)
ax.set_yticks(np.linspace(-1,1,11))
ax.set_xlabel('X')
ax.set_ylabel('Scalar')
ax.set_title('Radial Wave Propagation')

# Animation function
def animate(frame):
    nn = frame * dn
    f_val = np.loadtxt('../output/fld/sfld' + f"{nn:06n}" + '_PPP', skiprows=1, max_rows=NR)
    
    R_d, Z_d = [], []
    for index, value in np.ndenumerate(f_val[:NR-1, :NP]):
        R_d.append(coords_r[index[0]])
        Z_d.append(value)
    
    ax.clear()
    ax.set_xlim(0, xlim)
    ax.set_xticks(np.linspace(0,xlim,11))
    ax.set_ylim(-1.1, 1.1)
    ax.set_yticks(np.linspace(-1,1,11))
    ax.set_xlabel('r')
    ax.set_ylabel('Scalar')
    ax.set_title('Radial Wave Propagation')

    ax.plot(R_d, Z_d, '-k')
    
    ax.legend(['t='+f"{nn*dt:10.3f}"])

    # Use tight layout to reduce whitespace around the plot
    fig.tight_layout()
    
# Run the animation
anim = animation.FuncAnimation(fig, animate, frames=range(0, 501), interval=50)
anim.save('animation.mp4', writer='ffmpeg', fps=30, dpi=150)
```

<video controls loop class="d-block mx-auto" style="width:100%; max-width:480px">
  <source src="{{ '/assets/videos/radial_wave_propagation_result.mp4' | relative_url }}" type="video/mp4">
</video>
