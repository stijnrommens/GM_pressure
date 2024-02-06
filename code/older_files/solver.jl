using BoundaryValueDiffEq, OrdinaryDiffEq, Trapz


### ---------- Check input ----------
if print_flag == true
    print("Solver messages:")
end

function check_input(ion_tot, n_salts, print_flag)
    """Check if sum of charges = 0  and if number of salts is correct.

    Args:
        ion_tot (tuple):
            1. q (float): Valency of the ion [-]
            2. ah (float): Hydrated radius of the ion [m]
            3. U(t) (function): Gibbs adsorption energy as function over distance [J/m]
        n_salts (integer): Number of salts present in the system [-]
        print_flag (boolean): Statement to (not) print the checks [-]

    Returns:
        nothing
    """
    if print_flag == true
        charge_tot = 0
        for ion in ion_tot
            charge_tot += ion[1]*ion[2] # Sum of charges [-]
        end
        if charge_tot != 0
            print("\n   ✖ Warning: Sum of charges is NOT 0!")
            # error("Inputs NOT correct!")
        else
            print("\n   ✔ Sum of charges is 0.")
        end
        if length(ion_tot)/n_salts > 2
            print("\n   ✖ Warning: Number of salts may NOT be correct!")
        else
            print("\n   ✔ Number of salts is correct.")
        end
    end
    nothing
end;
check_input(ion_tot, n_salts, print_flag);
# @btime check_input(ion_tot, n_salts, print_flag)


### ---------- modified Poisson-Boltzmann equation ----------
function mPBE!(du, u, p, t)
    """Initialize the modified Poisson-Boltzmann equation.

        Args:
            du (array):
                1. Second derivative [?]
                2. First derivative [?]
            u (array):
                1. Electrostatic potential [V]
                2. Derivative of electrostatic potential [V/m]
            p (tuple):
                1. constants (tuple):
                    1.1 beta (float): Thermodynamic beta [1/J]
                    1.2 elc (float): Elementary charge [C]
                    1.3 epsilon_o (float): Vacuum's permittivity [F/m]
                    1.4 epsilon_w (float): Water's permittivity [F/m]
                    1.5 Avog (float): Avogadro constant [1/mol]
                2. ion_tot (tuple):
                    2.1. q (float): Valency of the ion [-]
                    2.2. ah (float): Hydrated radius of the ion [m]
                    2.3. U(t) (function): Gibbs adsorption energy as function over distance [J/m]
                3. n_salts (integer): Number of salts present in the system [-]
            t (float): Distance from G/L-interface [m]
    
        Returns:
            nothing
        """
    constants, ion_tot, n_salts = p
    beta, elc, epsilon_o, epsilon_w, Avog, kappa = constants
    Gads_epot = 0
    # print(ion_tot)
    for ion in ion_tot
        rho_ion = ion[1] * Avog * 1000 # Ion density [particles/L] -> [particles/m3]
        if t < 1e10*ion[3]
            U = 1000
        else
            U = ion[4]*1e10*ion[3]/t * exp(-2kappa * (1e-10t - ion[3]))
        end
        Gads_epot += rho_ion * ion[2] * exp(-U - ion[2] * u[1]) # Sum of exponentials [-]
        # Gads_epot += rho_ion * ion[2] * exp(-(ion[3])(t) - ion[2] * u[1]) # Sum of exponentials [-]
    end
    du[1] = -u[2] # Second derivative [?]
    du[2] = 1e-20beta * elc^2 / (epsilon_o*epsilon_w) * Gads_epot #/ n_salts # First derivative [?] IS n_salts actually needed, or did Tim used it to use comparable numbers for mixtures?
    nothing
end;

function bc!(residual, u, p, t)
    """Initialize boundary conditions

        Args:
            residual (array):
                1. Boundary condition at G/L-interface [V/m]
                2. Boundary condition in bulk liquid [V]
            u (array):
                1. Electrostatic potential [V]
                2. Derivative of electrostatic potential [V/m]
            p (tuple):
                1. constants (tuple):
                    1.1 beta (float): Thermodynamic beta [1/J]
                    1.2 elc (float): Elementary charge [C]
                    1.3 epsilon_o (float): Vacuum's permittivity [F/m]
                    1.4 epsilon_w (float): Water's permittivity [F/m]
                    1.5 Avog (float): Avogadro constant [1/mol]
                2. ion_tot (tuple):
                    2.1. q (float): Valency of the ion [-]
                    2.2. ah (float): Hydrated radius of the ion [m]
                    2.3. U(t) (function): Gibbs adsorption energy as function over distance [J/m]
                3. n_salts (integer): Number of salts present in the system [-]
            t (float): Distance from G/L-interface [m]  
    
        Returns:
            nothing
        """
    residual[1] = u[end][2] -0.0 # Zero electrostatic field in vacuum
    residual[2] = u[1][1] -0.0 # Zero electrostatic potential in bulk liquid
    nothing
end;
tspan = (boundary, 0.01); # Distance [Å]
u0  = [0.0, 0.0];
mBPE_constants = (beta, elc, epsilon_o, epsilon_w, Avog, kappa);
mPBE_param = (mBPE_constants, ion_tot, n_salts);
bvp = BVProblem(mPBE!, bc!, u0, tspan, mPBE_param);
# @btime BVProblem(mPBE!, bc!, u0, tspan, mPBE_param)
sol = solve(bvp, Shooting(RadauIIA5(autodiff=false)), abstol=1e-15, reltol= 1e-15);
print(sol.retcode)
# @btime solve(bvp, Shooting(RadauIIA5(autodiff=false)), abstol=1e-12, reltol= 1e-6)



### ---------- Gibbs-Marangoni pressure ----------
function pGM!(time, potential, p, print_flag)
    """
    Calculate the surface excesses, total surface tension & total Gibbs-Marangoni pressure at separation 'h'

    Args:
        time (array):
            1. t (float): Distance from G/L-interface [m]
        potential (array):
            1. Electrostatic potential at distance t [V]
        p (tuple):
                1. constants (tuple):
                    1.1 STconst (float): Surface tension constant [M.m.Å/mN]
                    1.2 Avog (float): Avogadro constant [1/mol]
                    1.3 ionS (float): Ionic strength [M]
                    1.4 beta (float): Thermodynamic beta [1/J]
                    1.5 h (float): Separation distance [nm]
                2. ion_tot (tuple):
                    2.1. q (float): Valency of the ion [-]
                    2.2. ah (float): Hydrated radius of the ion [m]
                    2.3. U(t) (function): Gibbs adsorption energy as function over distance [J/m]
                3. n_salts (integer): Number of salts present in the system [-]
        print_flag (boolean): Statement to (not) print the checks [-]

    Returns:
        conc_matrix (array):
            1. conc (float): Ion's local concentration fraction at distance t from G/L-interface [-]
        gammacon_qc (array):
            1. gammacon_qc (float): Ion's surface excess × charge × concentration [Å.M]
        tension (float): Total surface tension [mN/m.M]
        pGM (float): Total Gibbs-Marangoni pressure [Pa]
    """
    constants, ion_tot, n_salts = p
    STconst, Avog, ionS, beta, h, C = constants
    conc_matrix = zeros(Float64, size(time)[1], length(ion_tot))
    for (index, t) in enumerate(time)
        for (index2, ion) in enumerate(ion_tot)
            conc_matrix[index, index2] = exp(-(ion[3])(t)-ion[2]*potential(t)[1]) - 1 # Fractionate ionic concentration profile [M/M] = [-], Eq.6
            # conc_matrix[index, index2] = C*exp(-(ion[3])(t)-ion[2]*potential(t)[1]) - C 
        end
    end
    gammacon = zeros(1, length(ion_tot))
    gammacon_qc = zeros(1, length(ion_tot))
    for (index, ion) in enumerate(ion_tot)
        gammacon[index] = trapz(time, conc_matrix[:, index]) # Surface excess over cocentration [M.Å/M] = [Å], Eq.6
        gammacon_qc[index] = gammacon[index]*ion[1]*ion[2] # Surface excess × charge × concentration [Å.M]
    end
    if print_flag == true
        if sum(gammacon_qc) > 0.05
            print("\n   ✖ Warning: Electroneutrality condition is NOT satisfied!")
        else
            print("\n   ✔ Electroneutrality condition is satisfied.")
        end
    end
    # print(gammacon)
    # tension = sum(gammacon)/(STconst*n_salts)# Surface tension [Å] / [M.m.Å/mN] -> [mN/m.M]
    tension = sum(gammacon)/STconst # Surface tension [Å] / [M.m.Å/mN] -> [mN/m.M]
    pGM = 0
    for (index, ion) in enumerate(ion_tot)
        # gammacon_squared = sum(map(x->x^2, gammacon)) # Surface excess squared [Å^2]
        gammacon_squared = gammacon[index]^2 # Surface excess squared [Å^2]
        pGM += 4Avog * ion[1]*1000 * gammacon_squared*1e-20 / (beta*h^2) # Gibbs-Marangoni pressure [m4.mol.kg/m5.mol.s2] -> [Pa], Eq.5
    end
    return (conc_matrix, gammacon, tension, pGM)
end;
time_range = range(0.01, boundary, 100001);
pGM_constants = (STconst, Avog, ionS, beta, h, C);
pGM_param = (pGM_constants, ion_tot, n_salts);
# sol2 = pGM!(time_range, sol, pGM_param, print_flag);
# @btime pGM!(time_range, sol, pGM_param, print_flag)
