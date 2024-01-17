using BoundaryValueDiffEq, OrdinaryDiffEq, Trapz #,BenchmarkTools


### ---------- Check input ----------
if print_flag == true
    print("Solver messages:")
end

function check_input(ion_tot, n_salts, print_flag)
    """
    Check if sum of charges = 0  and if number of salts is correct.
    """
    if print_flag == true
        charge_tot = 0

        for (index, ion) in enumerate(ion_tot)
            charge_tot += ion[1]*ion[2] # Sum of charges [-]
        end

        if charge_tot != 0
            print("\n   ✖ Warning: Sum of charges is NOT 0!")
            # error("Inputs NOT correct!")
        else
            print("\n   ✔ Sum of charges is 0.")
        end

        if length(ion_tot)/n_salts > 3
            print("\n   ✖ Warning: Number of salts may NOT correct!")
        else
            print("\n   ✔ Number of salts is correct.")
        end
    end
    nothing
end;
check_input(ion_tot, n_salts, print_flag)


### ---------- modified Poisson-Boltzmann equation ----------
function mPBE!(du, u, p, t)
    """ 
    Modified Poisson-Boltzmann equation. 
    This equation is solved as a BVP to obtain the electrostatic potential. 
    """

    constants, ion_tot, n_salts = p
    beta, elc, epsilon_o, epsilon_w, Avog = constants

    Gads_epot = 0
    for ion in ion_tot
        rho_ion = ion[1] * Avog * 1000 # Ion density [particles/L] -> [particles/m3]
        Gads_epot += rho_ion * ion[2] * exp(-(ion[3])(t) - ion[2] * u[1]) # Sum of exponentials [-]
    end

    du[1] = -u[2] # Second derivative [?]
    du[2] = 1e-20beta * elc^2 / (epsilon_o*epsilon_w) * Gads_epot #/ n_salts # First derivative [?] IS n_salts actually needed, or did Tim used it to use comparable numbers for mixtures?
    nothing
end;

function bc!(residual, u, p, t)
    """ 
    Boundary conditions of BVP.
    """

    residual[1] = u[end][2] -0.0 # Zero electrostatic field in vacuum
    residual[2] = u[1][1] -0.0 # Zero electrostatic potential in bulk liquid
    nothing
end;
tspan = (boundary, 0.01); # Distance [Å]
u0  = [0.0, 0.0];
mBPE_constants = (beta, elc, epsilon_o, epsilon_w, Avog);
mPBE_param = (mBPE_constants, ion_tot, n_salts);
bvp = BVProblem(mPBE!, bc!, u0, tspan, mPBE_param);
sol = solve(bvp, Shooting(RadauIIA5(autodiff=false)), abstol=1e-12, reltol= 1e-6);
# @btime sol = solve(bvp, Shooting(RadauIIA5(autodiff=false)))



### ---------- Gibbs-Marangoni pressure ----------
function pGM!(time, potential, p, print_flag)
    """
    Use electrostatic potential to obtain the Gibbs-Marangoni pressure (at 10 nm separation).
    """

    constants, ion_tot, n_salts = p
    STconst, Avog, ionS, beta, h = constants

    conc_matrix = zeros(Float64, size(time)[1], length(ion_tot))
    for (index, t) in enumerate(time)
        for (index2, ion) in enumerate(ion_tot)
            conc_matrix[index, index2] = exp(-(ion[3])(t)-ion[2]*potential(t)[1]) - 1 # Fractionate ionic concentration profile [M/M] = [-], Eq.6
        end
    end

    gammacon = zeros(1, length(ion_tot))
    gammacon_q = zeros(1, length(ion_tot))
    for (index, ion) in enumerate(ion_tot)
        gammacon[index] = trapz(time, conc_matrix[:, index]) # Surface excess over cocentration [M.Å/M] = [Å], Eq.6
        gammacon_q[index] = gammacon[index]*ion[2] # Surface excess × charge [Å]
    end

    if print_flag == true
        if sum(gammacon_q) > 0.05
            print("\n   ✖ Warning: Electroneutrality condition is NOT satisfied!")
        else
            print("\n   ✔ Electroneutrality condition is satisfied.")
        end
    end

    tension = sum(gammacon)/(STconst*n_salts) # Surface tension [Å] / [M.m.Å/mN] -> [mN/m.M]
    # tension = sum(gammacon)/STconst # Surface tension [Å] / [M.m.Å/mN] -> [mN/m.M]
    
    pGM = 0
    for (index, ion) in enumerate(ion_tot)
        # gammacon_squared = sum(map(x->x^2, gammacon)) # Surface excess squared [Å^2]
        gammacon_squared = gammacon[index]^2 # Surface excess squared [Å^2]
        pGM += 4Avog * ion[1]*1000 * gammacon_squared*1e-20 / (beta*h^2) # Gibbs-Marangoni pressure [m4.mol.kg/m5.mol.s2] -> [Pa], Eq.5
    end

    return (conc_matrix, gammacon_q, tension, pGM)
end;
time_range = range(0.01, boundary, 100001);
pGM_constants = (STconst, Avog, ionS, beta, h);
pGM_param = (pGM_constants, ion_tot, n_salts);
sol2 = pGM!(time_range, sol, pGM_param, print_flag);
# @btime sol2 = pGM!(time_range, sol, pGM_param)
