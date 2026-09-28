# Calculation methods and interpretation

## Fixed-condition evolution

For a selected network and constant conditions, the application solves

$$
\frac{d\mathbf N}{dt}=\mathbf A\mathbf N,\qquad
\mathbf N(t)=\exp(\mathbf A t)\mathbf N_0.
$$

A matrix column corresponds to a parent isotope. Diagonal losses and off-diagonal production are assembled from the enabled transitions. A decay branch contributes $\lambda_j b_{j\to i}$; a neutron reaction contributes $10^{-24}\Phi\sigma_{j\to i}$ when cross sections are in barns and flux in neutrons per square centimetre per second. Registry reaction links already include the isomeric branch yield in their `sigma_barn` value.

Decay-type and half-life filters act on channels. Manual decay channels are explicit overrides. Disabling a channel is not equivalent to removing its parent isotope.

## Solver and time grid

The main method is adaptive sparse CRAM-16. Step doubling checks internal accuracy, with additional concentration and applicable inventory checks. The solver refines internal steps independently of the saved output grid. It does not fix a trajectory by clipping significant negative values or rescaling its total inventory.

Padé uses SciPy's matrix exponential. The Python API also offers Taylor and BDF paths. The `cram48` method name is accepted for existing callers but routes to adaptive CRAM-16; comparing these two names is not independent solver validation. There is no MMPA implementation in the current code.

Evolution normally saves zero plus 399 geometrically spaced positive times. A time slider selects one of these saved points; it does not recompute at an arbitrary time. More output points improve temporal quadrature resolution, but do not substitute for solver error control.

## Flux and neutron density

The density option converts neutron number density to a flux using the entered energy:

$$
\Phi=n_n v,\qquad v=\sqrt{2E/m_n}.
$$

The implementation converts the speed to centimetres per second. This is a single-energy relation; it is not a complete thermal spectrum treatment. Keep the meaning of the energy and the selected MACS/pointwise data consistent. Equilibrium and sensitivity have their own flux settings.

## Inventories and mass-number sums

The input determines the units of $N_i$. For example, Fe-56 = 1 represents one unit of initial seed inventory, not automatically one nucleus per cubic centimetre.

$$
N(A,t)=\sum_{i:A_i=A}N_i(t),\qquad
f(A,t)=\frac{N(A,t)}{\sum_i N_i(t)}.
$$

A capture-weighted mass curve is

$$
S(A)=\sum_{i:A_i=A}\sigma_i N_i.
$$

It is not an arbitrary single $\sigma(A)$ multiplied by the mass-number sum. A sum over all isobars is also not the same selection as a classical stable s-only isotope comparison.

For averaged results, capture cross sections are the saved MT=102 snapshot, summed over isomeric branches once. Missing values are unavailable, not zero, and can produce gaps in the graph. The ordinary time-slice plotting path uses its stored isotope links and currently takes the first matching capture channel; for branched captures its weighting can differ from the averaged snapshot. Do not treat that display difference as a change in the physical inventory.

## Uniform physical-time average

For an interval inside the stored trajectory,

$$
\langle N_i\rangle_t=\frac{1}{t_b-t_a}\int_{t_a}^{t_b}N_i(t)\,dt.
$$

The implementation inserts linearly interpolated endpoint values if needed, then integrates the saved piecewise-linear trajectory with the trapezoid rule. Equal durations receive equal weight. It is not an arithmetic mean of output rows on a logarithmic grid.

The output retains the input inventory units. This is a mean over the history of one irradiation, not automatically a model for mixing stellar zones.

## Exponential irradiation-duration mixture

$$
p(t)=\frac{1}{t_0}e^{-t/t_0},\qquad
\overline N_i=\int_0^{t_{\max}}N_i(t)p(t)\,dt.
$$

Here $t_0$ is the mean duration of the complete exponential distribution. The implementation integrates the product at saved points by trapezoids. It combines endpoints of different irradiation durations from the same initial inventory; it does not successively irradiate a mixture.

The covered distribution weight is

$$
W=1-e^{-t_{\max}/t_0}.
$$

There is no automatic division by $W$. The missing tail is reported. At $t_{\max}=7t_0$, about 99.91% of the distribution weight is included. Small missing weight is not, by itself, a bound on the error of every rare isotope. Refine the saved grid and extend the calculation to check convergence.

## Exponential neutron-exposure mixture

For a positive constant flux,

$$
\tau[\mathrm{mbarn}^{-1}]=10^{-27}\Phi t,\qquad
\rho(\tau)=\frac{1}{\tau_0}e^{-\tau/\tau_0}.
$$

The concentrations are not rescaled: $N_i(\tau_k)=N_i(t_k)$. The engine computes the finite-range integral of $N_i(\tau)\rho(\tau)$ with the exposure interval widths. Its temporal equivalent has

$$
t_0=\frac{\tau_0}{10^{-27}\Phi}.
$$

Temporal averaging supports zero-flux decay trajectories; exposure averaging requires positive flux and fixed-condition provenance. The GUI does not offer a time-varying matrix history.

After any averaging, the relative mode divides by the total averaged inventory. That normalization is a composition display, not a correction for lost material. The three averaging modes share the same result tables, graph renderer and session storage.

## Heat

Heat is a dot product of the parent heat coefficients and the inventory. Depending on the configured network, terms may include decay and reaction energy. The GUI's volumetric conversion uses

$$
Q[\mathrm{W/cm^3}]=Q[\mathrm{MeV/s/initial\ nucleus}]
\,1.602\times10^{-13}\,\frac{\rho_m N_A}{M}.
$$

This conversion assumes a normalized initial inventory and uses a representative initial mass. It is not valid to apply the density multiplier a second time to an inventory already expressed as nuclei per cubic centimetre. Missing channel energies limit the meaning of the heat prediction; the model does not solve radiation transport or energy deposition.

## Equilibrium and sensitivity

Equilibrium mode searches for a concentration plateau or evaluates a user-specified time for each independent flux. A stopped search is not proof of equilibrium: inspect the returned times, convergence and boundary losses.

For a target response $F$ and fractional parameter perturbation $\delta$, central OAT sensitivity is

$$
S_p\approx\frac{F(p(1+\delta))-F(p(1-\delta))}{2\delta F(p)}.
$$

The scanned parameters are reaction cross sections, decay constants, flux/energy and data-library choice. Library spread is a discrete comparison, not a derivative. A zero baseline makes relative sensitivity undefined; a displayed zero in that case should not be interpreted as insensitivity. OAT is local and does not capture joint correlated uncertainties. Verify the sensitivity configuration separately; it does not inherit every Evolution override/filter.

## Network boundaries and scientific validation

An existing daughter object outside the matrix is a sink in the current assembly. For reactions with no linked daughter, the open-channel settings determine boundary loss. Decays outside the selected matrix can also remove material. Closed one-heavy-daughter networks can conserve nucleus count; arbitrary reaction networks need not.

Check the initial inventory, enabled channels, data provenance, endpoints, signed concentrations and the total inventory before interpreting a curve. Averaging cannot recover missing nuclei or make an incomplete network reproduce a stellar abundance pattern. Observational comparisons require matching isotope selection, normalization, metallicity, initial seeds and irradiation assumptions.
