# Gibbs-Marangoni pressure
Numerical solver of the **modified Poisson-Boltzmann equation** from Duignan (2021) in Julia. This equation is modified as it not only includes an electrostatic potential but also an adsorption energy function of the present ions [1].

<img src="https://github.com/stijnrommens/GM_pressure/assets/113170925/cc816c5d-426a-4f06-9b4b-0b75e1747da0" width="400">

The model needs the following input:
- ionic strength;
- hydrated radius and charge of each ion;
- adsorption energy function as a function of distance from the gas-liquid interface, for each ion.

As output, the electrostatic potential of salts at the interface is obtained. This can be used to calculate the surface excess of ions, which results in the **Gibbs-Marangoni pressure** (at 10 nm separation).

Sources:\
[1] Duignan, T. T. (2021). The surface potential explains ion-specific bubble coalescence inhibition. _Journal of Colloid and Interface Science, 600_, 338–343. https://doi.org/10.1016/j.jcis.2021.04.144

# Prerequisites
Julia (v1.9) needs to be installed. The attached _.tomls_ in this repo should automatically install the required packages.

# How it works
The input can be given in the file _input.jl_:
- At the top, the ionic strength (**ionS**) needs to be given. This value and the predefined constants are then used to calculate the required parameters. \
  <img src="https://github.com/stijnrommens/GM_pressure/blob/main/ionS_fig.PNG" width="400">
- The next part of the input file asks for the adsorption energy function (**U**), hydrated radius (**r**) and charge (**q**) of each ion. The U-function for α-ions is already given, with a couple of input examples. The _misc.jl_ file also contains U-functions for the β-ions H<sup>+</sup> and ClO<sub>4</sub><sup>-</sup>. \
  <img src="https://github.com/stijnrommens/GM_pressure/blob/main/Uqr_fig.PNG" width="400">
- Lastly, the total ionic composition (**ion_tot**) and the number of salts (**n_salts**) need to be filled in. As an example, a mixture conposed of Na<sub>2</sub>SO<sub>4</sub> and NH<sub>4</sub>Cl has ion_tot = (Na, Na, SO<sub>4</sub>, NH<sub>4</sub>, Cl) and n_salts = 2. The ion_tot tuple can be extended as much as possible. \
  <img src="https://github.com/stijnrommens/GM_pressure/blob/main/ntot_fig.PNG" width="500">
- The file also prints the results and allows the creation of plots. To make plots, the code at the bottom of the script can be uncommented (also don't forget the import at the top of the file!).


# Debugging & Tweaking
Certain scenarios (e.g. high ionic strengths) may result in warnings/errors from the BVP-solver. This will most likely return values which go to infinity. The following solver arguments in _solver.jl_ can be changed:
- absolute tolerance (**abstol**); currently at 1e-12 (default is 1e-6)
- relative tolerance (**reltol**); currently at 1e-6 (default is 1e-3)
- time span (**tspan**); currently set at 10x Debye length

Changing these arguments also increases speed, as long as the results are accurate enough. Accuracy/validity can be tested by checking if the surface excess of all ions fulfils the **electroneutrality condition** ($\sum$ surface excesses = 0).
