# Ebsolve

[![Fortran](https://img.shields.io/badge/Fortran-734F96?logo=fortran&logoColor=fff)](https://fortran-lang.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-red.svg)](https://opensource.org/licenses/MIT)

* [PRESENTATION](#presentation)
* [USAGE AND OPTIONS](#usage-and-options)
* [REFERENCES](#references)

## PRESENTATION

Ebsolve is a Fortran 2008 program to solve the matrix transport equation of *radiant field theory* [[1](#1), [2](#2)] in a rectilinear waveguide geometry.
The solution of this equation directly provides the distribution of transmission eigenvalues through a disordered waveguide.

The matrix transport equation reads [[1](#1), [2](#2)]

<p>$$ \mathbf{\Omega}\cdot\nabla_{\mathbf{r}}\mathsf{g} = 
 -\frac{1}{2\ell_{\rm s}} [\tilde{\mathsf{Q}}(\mathbf{r}), \mathsf{g}] - \frac{\varepsilon}{2k} [\Lambda_3, \mathsf{g}] 
 + \frac{\mathrm{i}}{T}\Omega_{x}\Theta(m_{\rm a}-|\mathbf{\Omega}_\perp|)\delta(x-x_{\rm a})[\Lambda_+, \mathsf{g}] 
 + \mathrm{i}          \Omega_{x}\Theta(m_{\rm b}-|\mathbf{\Omega}_\perp|)\delta(x-x_{\rm b})[\Lambda_-, \mathsf{g}] $$</p>

and is similar to the *Eilenberger equation* for type-II superconductivity in presence of impurities [[3](#3), [4](#4)], hence the name of the program, a contraction of *Eilenberger solver*.

List of notations:
* $\mathsf{g}(\mathbf{\Omega},\mathbf{r})$: *Matrix radiance*, a 2-by-2 complex matrix which is the central variable and the main unknown of the transport equation.
* $\mathbf{\Omega}$: Unit vector pointing in the direction of propagation, such that $\mathbf{\Omega}=(\Omega_x,\Omega_{\perp})$.
* $\Omega_x$: Longitudinal component of the direction vector.
* $\mathbf{\Omega}_{\perp}$: Transverse component of the direction vector.
* $\mathbf{r}$: Position vector in usual space, such that $\mathbf{r}=(x,\mathbf{y})$.
* $x$: Longitudinal component of the position $\mathbf{r}$ in the waveguide direction.
* $\nabla_{\mathbf{x}}$: Gradient with respect to variable $\mathbf{x}$.
* $\ell_{\rm s}$: Scattering mean free path.
* $\varepsilon$: Imaginary shift of the Green function wavenumber. This parameter can be used to simulate absorption with $$. 
* $k$: The central wavenumber.
* $\mathrm{i}$: The imaginary unit.
* $T$: Transmission eigenvalue. Note that it must be shifted in the imaginary direction to get positive transmission eigenvalue distribution: $T=T_{\rm r}-\mathrm{i}\epsilon$ for $\epsilon>0$.
* $[\cdot,\cdot]$: Matrix commutator defined by $[\mathsf{A},\mathsf{B}]=\mathsf{A}\mathsf{B}-\mathsf{B}\mathsf{A}$.
* $\delta(x)$: Dirac delta distribution.
* $\Theta(x)$: Heaviside unit-step function.
* $x_{\rm a}$, $x_{\rm b}$: Positions of the contact interactions serving to measure the observables. $x_{\rm a}$ is the start point, and $x_{\rm b}$ the end point.
* $m_{\rm a}$, $m_{\rm b}$: Numerical aperatures of the contact interactions 'a' and 'b', respectively.
* $\tilde{\mathsf{Q}}(\mathbf{r})$: The *matrix field*, a 2-by-2 complex matrix given by the directional integral:
$$\tilde{\mathsf{Q}}(\mathbf{r}) = \oint \frac{d\mathbf{\Omega}}{S_d} \mathsf{g}(\mathbf{\Omega},\mathbf{r})$$
where $S_d$ is the surface of the unit $d$-ball.
* $\Lambda_1,\Lambda_2,\Lambda_3$: Standard 2-by-2 Pauli matrices.
* $\Lambda_+,\Lambda_-$: Raising/lowering Pauli matrices defined by $\Lambda_\pm=(\Lambda_1\pm\mathrm{i}\Lambda_2)/2$.


## USAGE AND OPTIONS

### `ebsolve` executable

The syntax of the executable `ebsolve` reads:
```
ebsolve settings.nml
```
This program reads and executes the instructions given by the configuration file `settings.nml`.

This file contains several [Fortran namelists](https://www.intel.com/content/www/us/en/docs/fortran-compiler/developer-guide-reference/2024-2/namelist.html) :

* `main_settings`: The only entry of this namelist is `task` which determines the quantity to compute. See [Tasks](#tasks) for the list of available tasks.

## REFERENCES

<a id="1">[1]</a>
David Gaspard and Arthur Goetschy,
*Radiant Field Theory: A Transport Approach to Coherent Control of Transmission through Disordered Media*,
[https://arxiv.org/abs/2411.10360](https://arxiv.org/abs/2411.10360).

<a id="2">[2]</a>
David Gaspard and Arthur Goetschy,
*Transmission eigenvalue distribution in disordered media from radiant field theory*,
[https://arxiv.org/abs/2411.10355](https://arxiv.org/abs/2411.10355).

<a id="3">[3]</a>
Gert Eilenberger,
*Transformation of Gorkov's equation for type II superconductors into transport-like equations*,
[Z. Phys. A **214**, 195-213 (1968)](https://doi.org/10.1007/BF01379803)

<a id="4">[4]</a>
Klaus D. Usadel,
*Generalized Diffusion Equation for Superconducting Alloys*,
[Phys. Rev. Lett. **25**, 507-509 (1970)](https://doi.org/10.1103/PhysRevLett.25.507)
