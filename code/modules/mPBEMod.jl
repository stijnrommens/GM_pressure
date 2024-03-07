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
    constants, ion_tot, n = p
    beta, elc, epsilon_o, epsilon_w, Avog, kappa = constants
    Gads_epot = 0.0

    for ion in ion_tot
        rho_ion = ion[1] * Avog * 1000 # Ion density [particles/L] -> [particles/m3]
        if ion[4] == "beta"
            if t < 1e10*ion[3]
                U = ion[5]*1e10*ion[3]/t * exp(-2kappa * (1e-10t - ion[3])) - 2.1
            else
                U = ion[5]*1e10*ion[3]/t * exp(-2kappa * (1e-10t - ion[3]))
            end
        elseif ion[4] == "proton"
            if t < 1e10*ion[3]
                U = beta/(4pi*1e-10t) * (elc^2)/(4epsilon_o*epsilon_w) * exp(-2kappa*1e-10t) - 3.05
            else
                U = beta/(4pi*1e-10t) * (elc^2)/(4epsilon_o*epsilon_w) * exp(-2kappa*1e-10t)
            end
        else
            if t < 1e10*ion[3]
                U = 1000
            else
                U = ion[5]*1e10*ion[3]/t * exp(-2kappa * (1e-10t - ion[3]))
            end
            # Gads_epot += rho_ion * ion[2] * exp(-U - ion[2] * u[1])        
        end
        Gads_epot += rho_ion * ion[2] * exp(-U - ion[2] * u[1])
    end
    du[1] = -u[2] # Second derivative [?]
    du[2] = 1e-20beta * elc^2 / (epsilon_o*epsilon_w) * Gads_epot #/ n # First derivative [?] IS n_salts actually needed, or did Tim used it to use comparable numbers for mixtures?
    nothing
end

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
end
