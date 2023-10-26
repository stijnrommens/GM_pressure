import numpy as np
from pHcalc import Acid, System

MW_NH3 = 17.031 # g/mol
H_nh3  = 5.9e-4 # mol/L.Pa
p_tot  = 101325 # Pa (1 atm)
R = 8.31446261815324 # m3.Pa/K.mol
T = 298.15 # K

# Airflow through column
F_air = 40 /1000/60 # L/min -> m3/s

# Concentrations in column
c_salt   = 0.5 # M
c_carbon = 0   # M Can be adjusted to have real values, will not have significant impact

# Dimensions column & workplace
h_column = 0.8   # m
d_column = 0.192 # m

h_kast = 3.0 # m
w_kast = 1.0 # m
l_kast = 1.5 # m
n_vent = 8   # 1/h, ventilatievoud


""" NH3 concentration in column """
sul = Acid(pKa=[-2.80, 1.99], charge=0, conc=c_salt  ) # H2SO4
nh4 = Acid(pKa=9.25         , charge=1, conc=2*c_salt) # 1 mol of salt contains 2 mol NH4+
co2 = Acid(pKa=[6.35, 10.33], charge=0, conc=c_carbon) # Open system, CO2 can dissolve in water

system = System(sul, nh4, co2)
system.pHsolve()
pH = system.pH
c_nh3 = 2*c_salt*nh4.alpha(system.pH)[1] # mol/L, 1 mol salt contains 2 mol NH4+

print('In column:')
print(f'   pH = {pH:.2f}')
print(f'   NH3 concentration = {c_nh3:.3e} mol/L \n')


""" NH3 flux to gas phase """
### Non-equilibrium
Ac = np.pi*(d_column/2)**2 # m3
V_water = np.pi*h_column*(d_column/2)**2 # m3

V_kast = h_kast*w_kast*l_kast # m3

v_Gs  = F_air/Ac   # m/s
kLa   = 2*v_Gs     # 1/s, homogeneous

c_nh3_air = 0 # mol/L, initial NH3 concentration in air
J_water = -kLa*(c_nh3_air - c_nh3) * MW_NH3  # g/m3.s, with respect to water
J_air = J_water*V_water/V_kast * 1000 * 3600 # g/m3.s -> mg/m3.h, with respect to gas

### Equilibrium
p_nh3 = c_nh3/H_nh3 # Pa

print('In air:')
print(f'   NH3 (kLa) = {J_air/n_vent:.3f} mg/m3')
print(f'   NH3 (eq.) = {p_nh3/(R*T)*MW_NH3*1000:.3f} mg/m3')