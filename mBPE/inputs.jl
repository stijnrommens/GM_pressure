cd(@__DIR__)
using Pkg
Pkg.activate(".")
using Printf, QuadGK #,Plots


### ---------- Concentration input ----------
ionS = 0.1; # Ionic strength [M]



### ---------- Constants ----------
kB = 1.3806503e-23;         # Boltzmann constant [J/K]k
elc = 1.6021765e-19;        # Elementary charge  [C]
Avog = 6.0221415e23;        # Avogadro constant  [1/mol]
epsilon_o = 8.85418782e-12; # Vacuum's permittivity [F/m]
epsilon_w = 78.3;           # Water's permittivity [F/m]
T = 297.15;                 # Temperature [K]

rho_ion = ionS * Avog * 1000;   # Ion density [particles/L] -> [particles/m3]
beta = 1 / (kB * T);            # Thermodynamic beta [1/J] = [s2/kg.m2] = [1/N.m]
STconst = -1e-6*1e10*beta/Avog; # Surfce tension constant [mol/N.m] -> [M.m.Å/mN]
kappa = sqrt(2 * elc^2 * Avog * ionS*1000 * beta / (epsilon_w*epsilon_o)); # Inverse Debye lenght [1/m]
boundary = 10e10/kappa;         # 10x Debye lenght [Å]


f(k, a) = (k*(sqrt(kappa^2 + k^2)*cosh(k*a) - k*sinh(k*a))) / (sqrt(kappa^2 + k^2)*(sqrt(kappa^2 + k^2)*cosh(k*a) + k*sinh(k*a)));
function W(q, ah)
    factor = beta * (q*elc)^2 / (2epsilon_w * 4pi * epsilon_o)
    W = quadgk(k -> f(k,ah), 0, 10.0e10)[1] * factor
end;



### ---------- Ion specific functions ----------
function U_alpha(t, ah, kappa, W)
    if t < 1e10ah # while inside hydrated ion radius
        U = 1000
    else          # outside of radius
        U = W*1e10ah/t * exp(-2kappa * (1e-10t - ah))
    end
end;
q_Na, r_Na = +1, 2.50e-10; # Levin (2009)
W_Na   = W(q_Na, r_Na);
U_Na(t) = U_alpha(t, r_Na, kappa, W_Na);
Na   = (q_Na, U_Na);

q_Cl, r_Cl = -1, 2.00e-10; # Levin (2010)
W_Cl   = W(q_Cl, r_Cl);
U_Cl(t) = U_alpha(t, r_Cl, kappa, W_Cl);
Cl   = (q_Cl, U_Cl);

q_NH4, r_NH4 = +1, 2.50e-10; # Kielland (1937)
W_NH4 = W(q_NH4, r_NH4);
U_NH4(t) = U_alpha(t, r_NH4, kappa, W_NH4);
NH4 = (q_NH4, U_NH4);

q_SO4, r_SO4 = -2, 3.79e-10; # Levin (2010)
W_SO4 = W(q_SO4, r_SO4);
U_SO4(t) = U_alpha(t, r_SO4, kappa, W_SO4);
SO4 = (q_SO4, U_SO4);



### ---------- Ions/salts input ----------
ion_tot = (Na, Na, SO4, NH4, Cl); # The ions composing the salts
n_salts = 2;                      # Number of salts



### ---------- Calculate ----------
include("solver.jl")

print("\n\nResults:")
print("\n   • Surface excess × charge  = ", round.(sol2[1]; digits=3), " = ", round.(sum(sol2[1]); digits=3), " Å")
@printf("\n   • Surface tension          = %.3f mN/m.M", sol2[2])
@printf("\n   • Electrostatic potential  = %.3e mV", 1000*last(sol.u)[1])
@printf("\n   • Gibbs-Marangoni pressure = %.3f Pa", sol2[3])



### ---------- Plot ----------
# plot(sol.t, reduce(hcat, sol.u)[1,:],
#     ylims=(-0.45,0.45),
#     # xlims=(0,60),
#     ylabel="Potential [V]",
#     xlabel="z [Å]"
# )

# plot(time_range, sol2[1],
#     labels=["H" "Cl"],
#     ylims=(0,0.8),
#     xlims=(0,12),
#     ylabel="Concentration [M]",
#     xlabel="z [Å]")
# scatter!(time_range, sol2[1],
#     ylims=(-2,8),
#     xlims=(0,12))
