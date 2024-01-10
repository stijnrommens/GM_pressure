# Gibbs-Marangoni pressure
This model is made to numerically solve the **modified Poisson-Boltzmann equation** from Duignan (2021) in Julia. This equation is modified as it not only includes a electrostatic potential, but also an adsorption energy function of the present ions [1].

<img src="https://github.com/stijnrommens/GM_pressure/assets/113170925/cc816c5d-426a-4f06-9b4b-0b75e1747da0" width="400">

\
The model needs the following input:
- ionic strength;
- hydrated radius and charge of each ion;
- adsorption energy function as a function of distance from the gas-liquid interface (also for each ion).

\
As output, the electrostatic potential of salts at the interface is obtained. This can be used to calculate the surface excess of ions, which results in the **Gibbs-Marangoni pressure**.

\
Sources:\
[1] Duignan, T. T. (2021). The surface potential explains ion specific bubble coalescence inhibition. _Journal of Colloid and Interface Science, 600_, 338–343. https://doi.org/10.1016/j.jcis.2021.04.144

# Prerequisites
Julia needs to be installed. The attached _.tomls_ in this repo should automatically install the required packages.

# How it works
The input can be given in the file _input.jl_:
- At the top, the ionic strenght (**ionS**) needs to be given. This value, together with the predefined constants, is then used to calculate the required parameters.
- The next part of the input file asks for the adsorption energy function (**U**), hydrated radius (**r**) and charge (**q**) of each ion. The U-function for alpha ions is already given, together with a couple of input examples. The _misc.jl_ file also contains U-functions for the beta ions H<sup>+</sup> and ClO<sub>4</sub><sup>-</sup>.
- Lastly, the total ionic composition (**ion_tot**) and the number of salts (**n_salts**) need to be filled in. As an example, a mixture conposed of Na<sub>2</sub>SO<sub>4</sub> and NH<sub>4</sub>Cl has ion_tot = (Na, Na, SO<sub>4</sub>, NH<sub>4</sub>, Cl) and n_salts = 2. The ion_tot tupple can be extended as much as possible.

# Debugging
It may happen that certain scenarios (high ionic strenght) result in warnings/errors from the BVP-solver. This will most likely return values which go to infinity. The following solver parameters in _solver.jl_ can be tweaked:
- absolute tolerance (**abstol**)
- relative tolerance (**reltol**)
- time span (**tspan**)
  
