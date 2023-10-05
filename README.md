# GM_pressure
This model is made to numerically solve the **modified Poisson-Boltzmann equation** from Duignan (2021) in Julia. This equation is modified as it not only includes a electrostatic potential, but also an adsorption energy function of the present ions [1].

The model needs the following input:
- concentration of salt;
- hydrated radius of each ion;
- adsorption energy function as a function of distance from the gas-liquid interface (for each ion).

As output, the electrostatic potential of salts at the interface is obtained. This can be used to calculate the surface excess of ions, which results in the **Gibbs-Marangoni pressure**.

Sources:\
[1] Duignan, T. T. (2021). The surface potential explains ion specific bubble coalescence inhibition. _Journal of Colloid and Interface Science, 600_, 338–343. https://doi.org/10.1016/j.jcis.2021.04.144
