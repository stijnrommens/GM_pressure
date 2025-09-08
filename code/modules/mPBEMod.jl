using Parameters

function mPBE!(du, u, params, t)
    """Initialize the modified Poisson-Boltzmann equation.

        Args:
            du (array):
                1. First derivative [?]
                2. Second derivative [?]
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
    @unpack kB, T, elc, epsilon_o, epsilon_w, Avog, kappa, boundary = params
    @unpack ion_conc, ion_charges, ion_hyd_radii, ion_Wcal, ion_types = params
    # constants, ion_tot, n = p
    # beta, elc, epsilon_o, epsilon_w, Avog, kappa = constants
    Gads_epot = 0.0

    for (c_i, ch_i, hr_i, W_i, typ_i) in zip(ion_conc, ion_charges, ion_hyd_radii, ion_Wcal, ion_types)
        # TODO calculate kappa for the individual ion
        # TODO calculate W_i for the individual ion
        if typ_i == "beta"
            if t < hr_i
                U = W_i * hr_i / t * exp(-2kappa * (t - hr_i)) - 2.1 * kB * T
            else
                U = W_i * hr_i / t * exp(-2kappa * (t - hr_i))
            end
        elseif typ_i == "proton"
            if t < hr_i
                U = 1 / (4pi * t) * (elc^2) / (4epsilon_o * epsilon_w) * exp(-2kappa * t) - 3.05 * kB * T
            else
                U = 1 / (4pi * t) * (elc^2) / (4epsilon_o * epsilon_w) * exp(-2kappa * t)
            end
        else
            if t < hr_i
                U = 1000 * kB * T
            else
                U = W_i * hr_i / t * exp(-2kappa * (t - hr_i))
            end
            # Gads_epot += rho_ion * ch_i * exp(-U - ch_i * u[1])        
        end
        c_i_m = c_i * Avog # Ion density [mol/m3] -> [m-3]
        Gads_epot += ch_i * c_i_m * exp(-1 / (kB * T) * (U + ch_i * elc * u[1]))
    end
    
    du[1] = -u[2] # First derivative [V/Å]
    du[2] = elc / (epsilon_o * epsilon_w) * Gads_epot # Second derivative [V/Å2]
    nothing
end

function bc!(residual, sol, params, t)
    """Initialize boundary conditions

        Args:
            residual (array):
                1. Boundary condition at G/L-interface [V/m]
                2. Boundary condition in bulk liquid [V]
            sol (OdeSol):
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
    @unpack boundary = params
    residual[1] = sol(0.0)[2] - 0.0 # Zero electrostatic field in vacuum
    residual[2] = sol(boundary)[1] - 0.0 # Zero electrostatic potential in bulk liquid
    nothing
end
