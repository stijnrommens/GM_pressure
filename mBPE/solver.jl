cd(@__DIR__)
using Pkg
Pkg.activate(".")
using DifferentialEquations, Trapz #,BenchmarkTools


### ---------- modified Poisson-Boltzmann equation ----------
function mPBE!(du, u, p, t)
    constants, ion_tot, n_salts = p
    beta, elc, epsilon_o, epsilon_w, rho_ion = constants

    Gads_epot = 0
    for ion in ion_tot
        Gads_epot += ion[1]*exp(-(ion[2])(t) - ion[1] * u[1]) 
    end

    du[1] = -u[2]
    du[2] = 1e-20beta * elc^2 / (epsilon_o*epsilon_w) * (rho_ion*Gads_epot) / n_salts
    nothing
end;

function bc!(residual, u, p, t)
    residual[1] = u[end][2] -0.0
    residual[2] = u[1][1] -0.0
    nothing
end;
tspan = (boundary, 0.01);
u0  = [0.0, 0.0];
mBPE_constants = (beta, elc, epsilon_o, epsilon_w, rho_ion);
mPBE_param = (mBPE_constants, ion_tot, n_salts);
bvp = BVProblem(mPBE!, bc!, u0, tspan, mPBE_param);
sol = solve(bvp, Shooting(RadauIIA5(autodiff=false)));
# @btime sol = solve(bvp, Shooting(RadauIIA5(autodiff=false)))



### ---------- Gibbs-Marangoni pressure ----------
function pGM!(time, potential, p)
    constants, ion_tot, n_salts = p
    STconst, Avog, ionS, beta = constants

    conc_matrix = zeros(Float64, size(time)[1], length(ion_tot))
    for (index, t) in enumerate(time)
        for (index2, ion) in enumerate(ion_tot)
            conc_matrix[index, index2] = exp(-(ion[2])(t)-ion[1]*potential(t)[1]) - 1 # Concentration profile [M/M] = [-]

        end
    end

    gammacon = zeros(1, length(ion_tot))
    for (index, ion) in enumerate(ion_tot)
        gammacon[index] = trapz(time, conc_matrix[:, index]) # Surface excess over cocentration [M.Å/M] = [Å]
    end

    tension = sum(gammacon)/(STconst*n_salts) # [Å] / [M.m.Å/mN] -> [mN/m.M]
    gammacon_squared = sum(map(x->x^2, gammacon))
    pGM = 4Avog * ionS/n_salts*1000 * gammacon_squared*1e-20 / (beta*(10e-9)^2)

    return (gammacon, tension, pGM);
end;
time_range = range(0.01, boundary, 100001);
pGM_constants = (STconst, Avog, ionS, beta);
pGM_param = (pGM_constants, ion_tot, n_salts);
sol2 = pGM!(time_range, sol, pGM_param);
# @btime sol2 = pGM!(time_range, sol, param, ionS)
