cd(@__DIR__)
using Pkg
Pkg.activate(".")
Pkg.resolve()
Pkg.instantiate()
using Printf, QuadGK, BenchmarkTools, Plots


### ---------- Initialize for solver ----------
function ionic_strenght(ion_list, ionS)
    """Calculate the ionic strength

    Args:
    ion_list (array): 
    1. c (float): Concentration of all ions [M]
    2. q (float): Valency of all ions [-]
    ionS (float): Initialized ionic strength [M]
    
    Returns:
    ionS (float): Ionic strength for the given concentrations and valencies [M]
    """
    for ion in ion_list
        ionS += ion[1] * ion[2]^2
    end
    return 0.5*ionS
end;


f(k, a) = (k*(sqrt(kappa^2 + k^2)*cosh(k*a) - k*sinh(k*a))) / (sqrt(kappa^2 + k^2)*(sqrt(kappa^2 + k^2)*cosh(k*a) + k*sinh(k*a))); # [1/m2 / 1/m2] -> [-]?

function W(q, ah, beta, elc, epsilon_w, epsilon_o, factor)
    """Calculate the energy to bring ion out of Gibbs dividing surface
    
    Args:
    q (float): Valency of the ion [-]
    ah (float): Hydrated radius of the ion [m]
    
    Returns:
    W (float): Energy to bring ion out of Gibbs dividing surface [m?]
    """
    factor = beta * (q*elc)^2 / (2epsilon_w * 4pi * epsilon_o) # Squared radius of sphere [m2]?
    W = quadgk(k -> f(k,ah), 0, 10.0e10)[1] * factor # [m]?
    return W
end;


function U_alpha(t, ah, kappa, W, U)
    """Calculate the Gibbs adsorption energy of an α-ion
    
    Args:
    t (float): Distance from G/L-interface [m]
    ah (float): Hydrated radius of the ion [m]
    kappa (float): Inverse Debye-Hückel lenght [1/m]
    W (float): Energy to bring ion out of Gibbs dividing surface [m?]
    
    Returns:
    U (float): Gibbs adsorption energy [J]
    """
    if t < 1e10ah # While inside hydrated ion radius
        U = 1000 # Gibbs adsorption energy [J]
    else          # Outside of radius
        U = W*1e10ah/t * exp(-2kappa * (1e-10t - ah)) # Gibbs adsorption energy [J]
    end
    return U
end;


function give_ions(ion_list, kappa, beta, elc, epsilon_w, epsilon_o, factor, U0, ion_tot)
    """Combine all information needed for solver.jl
    
    Args:
    ion_list (array): 
    1. c (float): Concentration of the ion [M]
    2. q (float): Valency of the ion [-]
    3. ah (float): Hydrated radius of the ion [m]
    kappa (float): Inverse Debye-Hückel lenght [1/m]
    
    Returns:
    ion_tot (tuple):
    1. q (float): Valency of the ion [-]
    2. ah (float): Hydrated radius of the ion [m]
    3. U(t) (function): Gibbs adsorption energy as function over distance [J/m]
    """
    for ion in ion_list
        W_cal = W(ion[2], ion[3], beta, elc, epsilon_w, epsilon_o, factor)
        U_func(t) = U_alpha(t, ion[3], kappa, W_cal, U0)
        ion_tot = tuple(ion_tot..., (ion[1], ion[2], U_func))
    end
    return ion_tot
end;


function main(ion_list=false, print_flag::Bool = true)
    ### ---------- Constants & Input ----------
    kB = 1.3806503e-23;         # Boltzmann constant [J/K]
    elc = 1.6021765e-19;        # Elementary charge  [C]
    Avog = 6.0221415e23;        # Avogadro constant  [1/mol]
    epsilon_o = 8.85418782e-12; # Vacuum's permittivity [F/m]
    epsilon_w = 78.3;           # Water's permittivity [F/m]
    T = 297.15;                 # Temperature [K]
    h = 10e-9;                  # Separation distance [nm] where pGM is calculated at in Duignan (2021)

    ionS = 0.0;                 # Initialize ionic strength [M]
    factor = 0.0;
    U0 = 0.0;
    ion_tot = ();

    C = 1.00142; # M

    if ion_list == false
        ion_list = [[C, +1, 2.50e-10],  # Na, Ionic concentration [M], Charge [-], Hydrated radius [m] from Levin (2009)
                    [C, -1, 2.00e-10],  # Cl, Levin (2010)
                    [2C, +1, 2.50e-10], # NH4, Kielland (1937)
                    [C, -2, 3.79e-10]]; # SO4, Levin (2010)
    end

    n_salts = length(ion_list)/2; # Number of salts [-]

    ionS = ionic_strenght(ion_list, ionS);

    beta = 1 / (kB * T);            # Thermodynamic beta [1/J] = [s2/kg.m2] = [1/N.m]
    STconst = -1e-6*1e10*beta/Avog; # Surfce tension constant [mol/N.m] -> [M.m.Å/mN]
    kappa = sqrt(2 * elc^2 * Avog * ionS*1000 * beta / (epsilon_w*epsilon_o)); # Inverse Debye-Hückel lenght [1/m]
    boundary = 10e10/kappa;         # 10x Debye-Hückel lenght [Å]

    ion_tot = give_ions(ion_list, kappa, beta, elc, epsilon_w, epsilon_o, factor, U0, ion_tot);


    ### ---------- Calculate ----------
    include("solver.jl")


    if print_flag
        print("\n\nResults:")
        print("\n   • Surface excess × q × c   = ", round.(sol_pGM[2]; digits=3), " = ", round.(sum(sol_pGM[2]); digits=3), " Å")
        @printf("\n   • Surface tension          = %.3e mN/m", sol_pGM[3]*C)
        @printf("\n   • Surface tension          = %.3e mN/m.M", sol_pGM[3])
        @printf("\n   • Electrostatic potential  = %.3e mV", 1000*last(sol_PBE.u)[1])
        @printf("\n   • Gibbs-Marangoni pressure = %.3f Pa", sol_pGM[4])
        # @printf("\n \U1FAE7")
    end
end

main()

### ---------- Plot ----------
# plot(sol.t, reduce(hcat, sol.u)[1,:],
#     ylims=(-0.45,0.45),
#     # xlims=(0,60),
#     ylabel="Potential [V]",
#     xlabel="z [Å]"
# )

# plot(time_range, sol2[1],
#     labels=["Na" "Cl"],
#     # ylims=(-0.005,0.005),
#     xlims=(0,30),
#     ylabel="Concentration [M]",
#     xlabel="z [Å]")
# scatter!(time_range, sol2[1],
#     ylims=(-2,8),
#     xlims=(0,boundary))
