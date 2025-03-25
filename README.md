# Ebsolve

[![Fortran](https://img.shields.io/badge/Fortran-734F96?logo=fortran&logoColor=fff)](https://fortran-lang.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-red.svg)](https://opensource.org/licenses/MIT)

* [PRESENTATION](#presentation)
    - [List of symbols](#list-of-symbols)
    - [Transmission eigenvalue distribution](#transmission-eigenvalue-distribution)
* [USAGE AND OPTIONS](#usage-and-options)
* [REFERENCES](#references)

## PRESENTATION

Ebsolve is a Fortran 2008 program to solve the matrix transport equation of *radiant field theory* [[1](#1), [2](#2)] in a rectilinear [waveguide](https://en.wikipedia.org/wiki/Waveguide) geometry.
The solution of this equation provides the distribution of [transmission eigenvalues](https://en.wikipedia.org/wiki/Landauer_formula) through a disordered waveguide.

The matrix transport equation reads [[1](#1), [2](#2)]

<p>$$ \mathbf{\Omega}\cdot\nabla_{\mathbf{r}}\mathsf{g} = 
 -\frac{1}{2\ell_{\rm s}} [\tilde{\mathsf{Q}}(\mathbf{r}), \mathsf{g}] - \frac{\varepsilon}{2k} [\Lambda_3, \mathsf{g}] 
 + \mathrm{i} \Omega_{x} \gamma_{\rm a} \Theta_{\rm a} \delta(x-x_{\rm a}) [\Lambda_+, \mathsf{g}] 
 + \mathrm{i} \Omega_{x} \gamma_{\rm b} \Theta_{\rm b} \delta(x-x_{\rm b}) [\Lambda_-, \mathsf{g}] $$</p>

where $\tilde{\mathsf{Q}}(\mathbf{r})$ is the *matrix field*, a 2-by-2 complex matrix given by the directional integral

<p>$$ \tilde{\mathsf{Q}}(\mathbf{r}) = \oint \frac{d\mathbf{\Omega}}{S_d} \mathsf{g}(\mathbf{\Omega},\mathbf{r}) $$</p>

$S_d$ being the surface of the unit $d$-ball.

The matrix transport equation is similar to the [Boltzmann equation](https://en.wikipedia.org/wiki/Boltzmann_equation) in [radiative transport theory](https://en.wikipedia.org/wiki/Radiative_transfer) but, in contrast to the latter, it is able to capture [coherent effects](https://en.wikipedia.org/wiki/Coherence_(physics)) described by the [transmission matrix](https://en.wikipedia.org/wiki/Landauer_formula).
This equation is also very similar to the *Eilenberger equation* for [type-II superconductivity](https://en.wikipedia.org/wiki/Type-II_superconductor) in presence of impurities [[3](#3), [4](#4)], hence the name of the program, a contraction of *Eilenberger solver*.

The matrix transport equation is supplemented by the boundary conditions at infinity ($\mathbf{r}\rightarrow\infty$) [[1](#1), [2](#2)]

<p>$$ \mathsf{g}_{\rm out}(\mathbf{\Omega},\mathbf{r}) = \begin{pmatrix}1 & g_{12}\\ 0 & -1\end{pmatrix} , \qquad
\mathsf{g}_{\rm in}(\mathbf{\Omega},\mathbf{r}) = \begin{pmatrix}1 & 0\\ g_{21} & -1\end{pmatrix} $$</p>

where the indices 'in' and 'out' denote the incoming and outgoing directions $\mathbf{\Omega}$ on the disordered region, respectively.
Note that the elements $g_{12}$ and $g_{21}$ are not fixed by the boundary conditions.

### List of symbols

Here is a comprehensive list of symbols encountered in the matrix transport equation above:

* $\mathsf{g}(\mathbf{\Omega},\mathbf{r})$: The *matrix radiance*, a 2-by-2 complex matrix which is the central quantity and the main unknown of the transport equation.
This quantity also implicitly depends on the [transmission eigenvalue](https://en.wikipedia.org/wiki/Landauer_formula) $T$.
It has the properties of being traceless, $\mathrm{Tr}\mathsf{g}=0$, and normalized to 1 according to $\mathsf{g}^2=\mathsf{1}$, due to the boundary conditions.
* $\mathbf{\Omega}$: Unit vector pointing in the direction of propagation, such that $\mathbf{\Omega}=(\Omega_x,\mathbf{\Omega}_{\perp})$.
* $\Omega_x$: Longitudinal component of the direction vector in the direction of the waveguide.
* $\mathbf{\Omega}_{\perp}$: Transverse component of the direction vector with respect to the waveguide.
* $\mathbf{r}$: Position vector in usual space, such that $\mathbf{r}=(x,\mathbf{y})$.
* $x$: Longitudinal component of the position $\mathbf{r}$ in the waveguide direction.
* $\nabla_{\mathbf{x}}$: [Gradient](https://en.wikipedia.org/wiki/Del) with respect to the generic variable $\mathbf{x}$.
* $\mathrm{i}$: The [imaginary unit](https://en.wikipedia.org/wiki/Imaginary_unit).
* $k$: The central [wavenumber](https://en.wikipedia.org/wiki/Wavenumber).
* $\ell_{\rm s}$: The scattering [mean free path](https://en.wikipedia.org/wiki/Mean_free_path) of the disordered region of the waveguide.
* $\varepsilon$: Imaginary shift of the Green function wavenumber. This parameter can be used to simulate absorption with the relation $\varepsilon=\frac{k}{\ell_{\rm a}}$ where $\ell_{\rm a}$ is the [absorption length](https://en.wikipedia.org/wiki/Attenuation_length).
* $[\cdot,\cdot]$: [Matrix commutator](https://en.wikipedia.org/wiki/Commutator) defined by $[\mathsf{A},\mathsf{B}]=\mathsf{A}\mathsf{B}-\mathsf{B}\mathsf{A}$.
* $\delta(x)$: [Dirac delta distribution](https://en.wikipedia.org/wiki/Dirac_delta_function).
* $x_{\rm a}, x_{\rm b}$: Positions of the contact interactions serving to measure the observables related to the transmission matrix. $x_{\rm a}$ is the start point, and $x_{\rm b}$ the end point.
* $\gamma_{\rm a}, \gamma_{\rm b}$: Parameters of the contact interactions. They are related to the [transmission eigenvalues](https://en.wikipedia.org/wiki/Landauer_formula) by $\gamma_{\rm a}\gamma_{\rm b}=\gamma=\frac{1}{T}+\mathrm{i}\epsilon$ for arbitrarily small $\epsilon>0$.
In the numerical simulations, they are chosen as equal ($\gamma_{\rm a}=\gamma_{\rm b}$), but as explained in the paper [[2](#2)], it is somehow arbitrary as long as they are related to $T$ by the previous expression.
* $\Theta_{\rm a}, \Theta_{\rm b}$: Window functions limiting the [numerical aperture](https://en.wikipedia.org/wiki/Numerical_aperture) of the contact interactions 'a' and 'b'.
They are defined by $\Theta_{\rm a}=\Theta(m_{\rm a}-|\mathbf{\Omega}_\perp|)$, where $\Theta(x)$ stands for the [Heaviside unit-step function](https://en.wikipedia.org/wiki/Heaviside_step_function).
* $m_{\rm a}, m_{\rm b}$: [Numerical apertures](https://en.wikipedia.org/wiki/Numerical_aperture) of the contact interactions 'a' and 'b', respectively.
* $\Lambda_1,\Lambda_2,\Lambda_3$: The standard 2-by-2 [Pauli matrices](https://en.wikipedia.org/wiki/Pauli_matrices).
* $\Lambda_+,\Lambda_-$: The raising/lowering Pauli matrices defined by $\Lambda_\pm=(\Lambda_1\pm\mathrm{i}\Lambda_2)/2$.

### Transmission eigenvalue distribution

Once the matrix transport equation is solved for $\mathsf{g}(\mathbf{\Omega},\mathbf{r})$, the transmission eigenvalue distribution can be obtained from the generating function

<p>$$ F(\gamma) = \mathrm{i} \left( \gamma_{\rm a}' \tilde{J}_x^{21}(x_{\rm a}) + \gamma_{\rm b}' \tilde{J}_x^{12}(x_{\rm b}) \right) $$</p>

where $\tilde{\mathsf{J}}_x$ is the longitudinal component of the matrix current given by

<p>$$ \tilde{\mathsf{J}}_x(x) = \oint \frac{d\mathbf{\Omega}}{2V_{d-1}} \Omega_x \mathsf{g}(\mathbf{\Omega},\mathbf{r}) $$</p>

The prime refers to derivative with respect to $\gamma$, and $V_d$ is the volume of the unit $d$-ball.
The sought transmission eigenvalue distribution can be obtained by

<p>$$ \rho(T) = \frac{1}{\pi T^2} \mathrm{Im} F(\tfrac{1}{T} + \mathrm{i}\epsilon) $$</p>

for arbitrarily small $\epsilon>0$.

## INSTALLATION

Compiler requirements...

Use `git clone` and compile with `make all`...

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
