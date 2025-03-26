# Ebsolve

[![Fortran](https://img.shields.io/badge/Fortran-734F96?logo=fortran&logoColor=fff)](https://fortran-lang.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-red.svg)](https://opensource.org/licenses/MIT)

* [PRESENTATION](#presentation)
    - [List of symbols](#list-of-symbols)
    - [Transmission eigenvalue distribution](#transmission-eigenvalue-distribution)
* [INSTALLATION](#installation)
    - [Dependencies](#dependencies)
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

The source files can be downloaded using the `git clone` command:
```
git clone https://github.com/<name_of_repository>.git
```
To compile the program, call the `make` utility in the root directory:
```
make all
```

### Dependencies

The program has the following dependencies:

* The `gfortran` compiler, or any other compiler compliant with the [Fortran 2008](https://en.wikipedia.org/wiki/Fortran#Fortran_2008) standard and providing support for [OpenMP](https://en.wikipedia.org/wiki/OpenMP).
* The [LAPACK](https://en.wikipedia.org/wiki/LAPACK) Library.
* The `mkdir` command is called to create subdirectories to store the output data.
* [Python 3](https://en.wikipedia.org/wiki/Python_(programming_language)) scripts are called to generate the [TikZ](https://en.wikipedia.org/wiki/PGF/TikZ) codes for the plots.
* The LaTeX compiler `pdflatex` with the [PGFPlots package](https://www.ctan.org/pkg/pgfplots) is called to compile the plots.

## USAGE AND OPTIONS

To call the program, enter the following command:
```
ebsolve settings.nml
```
The `ebsolve` executable reads and executes the instructions given by the configuration file `settings.nml`.
This file contains several [Fortran namelists](https://www.intel.com/content/www/us/en/docs/fortran-compiler/developer-guide-reference/2024-2/namelist.html) containing the simulations parameters.
These namelists are described in details below.
See also the `settings.nml` file provided with the program.

### `main_settings`

The only entry in this namelist is `task` which specifies the task to be accomplished or the quantity to compute.
This entry points to another namelist of the same title in the `settings.nml` file.

### `eilenberger_system`

This namelist specifies the physical parameters of the disordered waveguide. The disordered aveguide is assumed to have length $L$.

* `modetype`: Defines the transverse profile of the waveguide, the shape and the boundary conditions. Either [`periodic`](#periodic) or [`infinite`](#infinite). See the namelists below.
* `nxdiso`: Number of 'x' points on the mesh in the bulk. Should be at least 2. Recommended value is `nxdiso=100`.
* `nxfree`: Number of 'x' points on the mesh at one edge. Should be at least 1. Recommended value is `nxfree=1` for `distrib` simulations.
* `dscat`: Scattering thickness of the disordered region, $L/\ell_{\rm s}$, where $\ell_{\rm s}$ is the scattering mean free path.
* `dabso`: Absorption thickness of the disordered region, $L/\ell_{\rm a}$, where $\ell_{\rm a}$ is the absorption length.
* `xa`: Position of the contact point $x_{\rm a}$ on the x mesh. The disordered region is for $x\in[0, 1]$.
* `xb`: Position of the contact point $x_{\rm b}$ on the x mesh. The disordered region is for $x\in[0, 1]$.
* `naper_a`: Numerical aperture of the input lead, $m_{\rm a}$. In 2D, fraction of excited modes. It must be in $[0, 1]$.
* `naper_b`: Numerical aperture of the output lead, $m_{\rm b}$. In 2D, fraction of observed modes. It must be in $[0, 1]$.

### `periodic`

This namelist specifies a squared-cross-section waveguide with *periodic* boundary conditions in the transverse direction.

* `d`: Total number of dimensions of the waveguide. Typically, d=2 for a 2D system.
* `wol`: Width-to-wavelength ratio, $W/\lambda$, defining the number of modes. In 2D, thresholds of opening of modes occur at *integer* values of `wol`.
Therefore, integer values of `wol` are forbidden because of division by zero. The standard value of the short paper [[1]](#1) is `wol=50.5`.

### `infinite`

This namelist specifies an infinite slab (infinitely wide waveguide). Integration path over the directions cosines $\mu=\cos\theta$ is deformed in the complex plane of $\mu$.

* `d`: Total number of dimensions of the waveguide. Typically, d=2 for a 2D system.
* `nmu`: Number of $\mu$ points, ideally very large. In the interval `nmu=50..200` is generally enough.
* `ash`: Shift factor of the integration path in the complex $\mu$ plane. The integration path is: $\mu(t) = t + \mathrm{i} a t(1-t)(1+t)$ for $t\in[0, 1]$.

### `distrib`

This namelist specifies the `distrib` task ordered by the `main_settings` namelist.
Compute the probability density function, $\rho(T)$, for the transmission eigenvalues over an interval of transmission eigenvalue.
The points of $T$ are taken as the [Chebyshev nodes](https://en.wikipedia.org/wiki/Chebyshev_nodes) of the first kind in order to refine the mesh at the two edges of the distribution (where the density usually increases).

* `ntm`: Number of desired samples for the distribution $\rho(T)$. Typically, a multiple of the thread number. Recommended value is `ntm=256`.
* `tmin`: Minimum transmission eigenvalue. This bound is never reached exactly by the Chebyshev nodes.
* `tmax`: Maximum transmission eigenvalue. This bound is never reached exactly by the Chebyshev nodes.
* `geps`: Shift of the gamma values to avoid being exactly on the real gamma axis. In principle, this value can be very small but nevertheless positive. Recommended value is `geps=1.0e-15`.
* `nthreads`: Number of threads used by OpenMP to parallelize the computations of the distribution points. Recommended value is the number of CPU cores.

### `fields`



### `solver`


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
