---
title: '#3. Spectral Operation'
parent: Tutorials
nav_order: 4
---

# MLegS Tutorial 03: Spectral Operation
*Disclaimer: This MLegS tutorial assumes a Linux or other Unix-based environment that supports bash terminal commands. If you are using Windows, consider installing the [Windows Subsystem for Linux (WSL)](https://learn.microsoft.com/en-us/windows/wsl/install).*

In this tutorial, you will apply a differentiation operation to a scalar field in spectral space. Spectral space operations that convert one set of spectral coefficients to another generally provide more accurate results than physical space operations, such as finite differences or other collocation point-based methods. Many differentiation operations can be performed directly within the spectral space, especially for linear operations like the Laplacian \\( \nabla^2 \equiv \frac{1}{r} \frac{\partial}{\partial r} \left( r \frac{\partial}{\partial r} \right) + \frac{1}{r^2} \frac{\partial^2}{\partial \phi^2} + \frac{\partial^2}{\partial z^2} \\).

MLegS, in its most recent version, provides essential spectral operations commonly used in physical simulations, like the Laplacian \\( \nabla^2 \\), the radial derivative operator \\( r \frac{\partial}{\partial r} \\), and the Helmholtz operator \\( \nabla^2 + \alpha I \\), along with their inverse operations when inversion is feasible. This tutorial explains how these spectral space operations are implemented, using the Laplacian \\( \nabla^2 \\) and its inverse as examples, so that you can learn the following:

1. **Obtain Laplacian of a Scalar**
   - Given a scalar field \\(s \\), its Laplacian, \\( \nabla^2 s \\), is obtainable via the MLegS built-in subroutine `del2()`.
   - The operation does not require the scalar's physical form and is directly converted in the spectral space.
2. **Get Acquainted with Other Built-In MLegS Spectral Differentiation Operators**
   - MLegS also provides the 2-dimensional Laplacian operator \\( \nabla_\perp^2 \equiv \frac{1}{r} \frac{\partial}{\partial r} \left( r \frac{\partial}{\partial r} \right) \\), the Helmholtz operator \\( \nabla^2 + \alpha I \\), and the powered Helmholtz operator \\( \nabla^{h} + \beta \nabla^2 + \alpha I \\) for \\( h = 4, 6, 8 \\).
   - Their mathematical definitions and typical usage for simulations are presented.

Completing this tutorial will help you apply spectral differentiation operators in MLegS with confidence.

## Under Construction