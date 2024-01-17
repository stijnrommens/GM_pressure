cd(@__DIR__)
using Pkg
Pkg.activate(".")
using Printf, QuadGK ,Plots


### ---------- Constants ----------
kB = 1.3806503e-23;         # Boltzmann constant [J/K]
elc = 1.6021765e-19;        # Elementary charge  [C]
Avog = 6.0221415e23;        # Avogadro constant  [1/mol]
epsilon_o = 8.85418782e-12; # Vacuum's permittivity [F/m]
epsilon_w = 78.3;           # Water's permittivity [F/m]
T = 297.15;                 # Temperature [K]
h = 10e-9;                  # Separation distance [nm] where pGM is calculated at in Duignan (2021)
beta = 1 / (kB * T);            # Thermodynamic beta [1/J] = [s2/kg.m2] = [1/N.m]
STconst = -1e-6*1e10*beta/Avog; # Surfce tension constant [mol/N.m] -> [M.m.Å/mN]



### ---------- Input ----------
# function give_ions(ion_list)
#     # ionS = 0.5 *
#     ion_tot = zeros(1, length(ion_list))
#     for (index, ion) in enumerate(ion_list)
#         W = W(ion[2], ion[2])
#         U(t) = U_alpha(t, ion[3], kappa, W)
#         ion_tot[index] = [ion[1], ion[2], U]

#     return ion_tot
# end;
# ion_tot = give_ions(ion_list)


# ionC = 0.1; # Ionic concentration [M]
# ionS = 0.1; # Ionic strength [M] CALCULATE FROM ionC FOR EACH SALT!
n_salts = 2; # Number of salts [-]
C = 0.1;
print_flag = false

c_Na,   q_Na,   r_Na   = C*2, +1, 2.50e-10; # Ionic concentration [M], Charge [-], Hydrated radius [m] from Levin (2009)
c_Cl,   q_Cl,   r_Cl   = 0, -1, 2.00e-10; # Levin (2010)
c_NH4,  q_NH4,  r_NH4  = 0, +1, 2.50e-10; # Kielland (1937)
c_SO4,  q_SO4,  r_SO4  = C, -2, 3.79e-10; # Levin (2010)
c_H,    q_H,    r_H    = 0.02/2, +1, 1.97e-10; # Levin (2009)
c_ClO4, q_ClO4, r_ClO4 = 0.02/2, -1, 2.83e-10; # Levin (2009)

ionS = 0.5 * (c_Na*q_Na^2 + c_SO4*q_SO4^2) # Ionic strength [M]
kappa = sqrt(2 * elc^2 * Avog * ionS*1000 * beta / (epsilon_w*epsilon_o)); # Inverse Debye-Hückel lenght [1/m]
boundary = 10e10/kappa;         # 10x Debye-Hückel lenght [Å]


f(k, a) = (k*(sqrt(kappa^2 + k^2)*cosh(k*a) - k*sinh(k*a))) / (sqrt(kappa^2 + k^2)*(sqrt(kappa^2 + k^2)*cosh(k*a) + k*sinh(k*a))); # [1/m2 / 1/m2] -> [-]?
function W(q, ah)
    """
    Energy to bring ion out of Gibbs dividing surface.
    """

    factor = beta * (q*elc)^2 / (2epsilon_w * 4pi * epsilon_o) # Squared radius of sphere [m2]?
    W = quadgk(k -> f(k,ah), 0, 10.0e10)[1] * factor # [m]?
end;



### ---------- Ion specific functions ----------
function U_alpha(t, ah, kappa, W)
    """ 
    Gibbs adsorption energy of an α-ion. 
    """

    if t < 1e10ah # While inside hydrated ion radius
        U = 1000 # Gibbs adsorption energy [J]
    else          # Outside of radius
        U = W*1e10ah/t * exp(-2kappa * (1e-10t - ah)) # Gibbs adsorption energy [J]
    end
end;
W_Na   = W(q_Na, r_Na);
U_Na(t) = U_alpha(t, r_Na, kappa, W_Na);
Na   = (c_Na, q_Na, U_Na);

W_Cl   = W(q_Cl, r_Cl);
U_Cl(t) = U_alpha(t, r_Cl, kappa, W_Cl);
Cl   = (c_Cl, q_Cl, U_Cl);

W_NH4 = W(q_NH4, r_NH4);
U_NH4(t) = U_alpha(t, r_NH4, kappa, W_NH4);
NH4 = (c_NH4, q_NH4, U_NH4);

W_SO4 = W(q_SO4, r_SO4);
U_SO4(t) = U_alpha(t, r_SO4, kappa, W_SO4);
SO4 = (c_SO4, q_SO4, U_SO4);


function U_beta1(z, ah, kappa, elc, beta, epsilon_o, epsilon_w)
    if z < 1e10ah # while inside hydrated ion radius
        U = 1/(4pi*epsilon_o) * elc^2*beta / (1e-10z*4epsilon_w) * exp(-2kappa*1e-10z) - 3.05
    else # outside of radius
        U = 1/(4pi*epsilon_o) * elc^2*beta / (1e-10z*4epsilon_w) * exp(-2kappa*1e-10z)
    end
end;
U_H(t) = U_beta1(t, r_H, kappa, elc, beta, epsilon_o, epsilon_w)
W_H    = W(q_H, r_H);
H      = (c_H, q_H, U_H);


function U_beta2(z, ah, kappa, W)
    if z < 1e10ah # while inside hydrated ion radius
        U = W * 1e10ah/z * exp(-2kappa * (1e-10z - ah)) - 2.1
    else # outside of radius
        U = W * 1e10ah/z * exp(-2kappa * (1e-10z - ah))
    end
end;
W_ClO4 = W(q_ClO4, r_ClO4);
U_ClO4(t) = U_beta2(t, r_ClO4, kappa, W_ClO4)
ClO4 = (c_ClO4, q_ClO4, U_ClO4);


### ---------- Ions/salts input ----------
ion_tot = (Na, SO4)#, NH4, Cl); # The ions composing the salts



### ---------- Calculate ----------
include("solver.jl")

if print_flag == true
    print("\n\nResults:")
    print("\n   • Surface excess × charge  = ", round.(sol2[2]; digits=3), " = ", round.(sum(sol2[2]); digits=3), " Å")
    @printf("\n   • Surface tension          = %.3f mN/m.M", sol2[3])
    @printf("\n   • Electrostatic potential  = %.3e mV", 1000*last(sol.u)[1])
    @printf("\n   • Gibbs-Marangoni pressure = %.3f Pa", sol2[4])
    @printf("\n \U1FAE7")
end



### ---------- Plot ----------
# plot(sol.t, reduce(hcat, sol.u)[1,:],
#     ylims=(-0.45,0.45),
#     # xlims=(0,60),
#     ylabel="Potential [V]",
#     xlabel="z [Å]"
# )

# plot(time_range, sol2[1],
#     labels=["Na" "Na" "SO4"],
#     # ylims=(-0.005,0.005),
#     xlims=(0,boundary),
#     ylabel="Concentration [M]",
#     xlabel="z [Å]")
# scatter!(time_range, sol2[1],
    # ylims=(-2,8),
    # xlims=(0,boundary))
