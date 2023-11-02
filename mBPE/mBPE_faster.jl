cd(@__DIR__)
using Pkg
Pkg.activate(".")
# Pkg.instantiate()
using QuadGK, DifferentialEquations, BenchmarkTools, Trapz, Plots#, NumericalIntegration, FastGaussQuadrature#, StaticArrays, Profile, ProfileView, SnoopCompile

kB = 1.3806503e-23;  # Boltzmann constant [J/K]
elc = 1.6021765e-19; # Elementary charge  [C]
Avog = 6.0221415e23; # Avogadro constant  [1/mol]
epsilon_o = 8.85418782e-12; # Vacuum's permittivity [F/m]
epsilon_w = 78.3;           # Water's permittivity [F/m]
ionS = 0.1; # M
T = 297.15; # K

rho_ion = ionS * Avog * 1000;
beta = 1 / (kB * T);
STconst = -1e-6*1e10*beta/Avog;
kappa = sqrt(2000*elc^2*Avog*ionS*beta/(epsilon_w*epsilon_o));
boundary = 10e10/kappa;
tspan = (boundary, 0.01);
u0  = [0.0, 0.0];
constants = (kappa, beta, elc, epsilon_o, epsilon_w, rho_ion);

f(k, a) = (k*(sqrt(kappa^2 + k^2)*cosh(k*a) - k*sinh(k*a))) / (sqrt(kappa^2 + k^2)*(sqrt(kappa^2 + k^2)*cosh(k*a) + k*sinh(k*a)));
function W(q, ah)
    factor = beta * (q*elc)^2 / (2epsilon_w * 4pi * epsilon_o)
    W = quadgk(k -> f(k,ah), 0, 10.0e10)[1] * factor
end;

function U_Na(z, ah, kappa, W)
    if z < 1e10ah # while inside hydrated ion radius
        U = 1000
    else # outside of radius
        U = W*1e10ah/z * exp(-2kappa * (1e-10z - ah))
    end
end;

function U_Cl(z, ah, kappa, W)
    if z < 1e10ah # while inside hydrated ion radius
        U = 1000
    else # outside of radius
        U = W*1e10ah/z * exp(-2kappa * (1e-10z - ah))
    end
end;

function U_H(z, ah, kappa, elc, beta, epsilon_o, epsilon_w)
    if z < 1e10ah # while inside hydrated ion radius
        U= 1/(4pi*epsilon_o) * elc^2*beta / (1e-10z*4epsilon_w) * exp(-2kappa*1e-10z) - 3.05
    else # outside of radius
        U = 1/(4pi*epsilon_o) * elc^2*beta / (1e-10z*4epsilon_w) * exp(-2kappa*1e-10z)
    end
end;

function U_ClO4(z, ah, kappa, W)
    if z < 1e10ah # while inside hydrated ion radius
        U = W * 1e10ah/z * exp(-2kappa * (1e-10z - ah)) - 2.1
    else # outside of radius
        U = W * 1e10ah/z * exp(-2kappa * (1e-10z - ah))
    end
end;


function mPBE!(du, u, p, t)
    constants, ion_tot, n_salts = p
    kappa, beta, elc, epsilon_o, epsilon_w, rho_ion = constants
    ioni, ionj, ionk, ionl = ion_tot
    qi, ahi, Wi = ioni
    qj, ahj, Wj = ionj
    qk, ahk, Wk = ionk
    ql, ahl, Wl = ionl

    # Ui = U_Na(t, ahi, kappa, Wi)
    Uj = U_Cl(t, ahj, kappa, Wj)
    Ui = U_H(t, ahk, kappa, elc, beta, epsilon_o, epsilon_w)
    # Uj = U_ClO4(t, ahl, kappa, Wl)

    du[1] = -u[2]
    if n_salts == 1
        du[2] = 1e-20beta * elc^2 / (epsilon_o*epsilon_w) * rho_ion/n_salts*( qi*exp(-Ui - qi * u[1]) + qj*exp(-Uj - qj * u[1]) )
    else
        du[2] = 1e-20beta * elc^2 / (epsilon_o*epsilon_w) * rho_ion/n_salts*( qi*exp(-Ui - qi * u[1]) + qj*exp(-Uj - qj * u[1]) + qk*exp(-Uk - qk * u[1]) + ql*exp(-Ul - ql * u[1]) )
    end
    nothing
end;

function bc!(residual, u, p, t)
    residual[1] = u[end][2] -0.0
    residual[2] = u[1][1] -0.0
    nothing
end;

Na   = (+1, 2.50e-10, W(+1, 2.50e-10));
Cl   = (-1, 2.00e-10, W(-1, 2.00e-10));
H    = (+1, 1.97e-10, W(+1, 1.97e-10));
ClO4 = (-1, 2.83e-10, W(-1, 2.83e-10));

ion_tot = (Na, Cl, H, ClO4); # Only change Gads functions in mBPE! and pGM!
n_salts = 1;

param = (constants, ion_tot, n_salts);
bvp = BVProblem(mPBE!, bc!, u0, tspan, param);
sol = solve(bvp, Shooting(RadauIIA5(autodiff=false)))#, dtmax=0.1)
# @btime sol = solve(bvp, Shooting(RadauIIA5(autodiff=false)))

# plot(sol.t, reduce(hcat, sol.u)[1,:],
#     ylims=(-0.45,0.45),
#     xlims=(0,60),
#     ylabel="Potential [V]",
#     xlabel="z [Å]"
# )


function pGM!(time, potential, p, c_bulk)
    constants, ion_tot, n_salts = p
    kappa, beta, elc, epsilon_o, epsilon_w, rho_ion = constants
    ioni, ionj, ionk, ionl = ion_tot
    qi, ahi, Wi = ioni
    qj, ahj, Wj = ionj
    qk, ahk, Wk = ionk
    ql, ahl, Wl = ionl

    conc_matrix = zeros(Float64, size(time)[1], 2n_salts)
    for (index, t) in enumerate(time)
        # Ui = U_Na(t, ahi, kappa, Wi)
        Uj = U_Cl(t, ahj, kappa, Wj)
        Ui = U_H(t, ahk, kappa, elc, beta, epsilon_o, epsilon_w)
        # Uj = U_ClO4(t, ahl, kappa, Wl)
        
        if n_salts == 1
            conc_matrix[index,1] = exp(-Ui-qi*potential(t)[1]) - 1
            conc_matrix[index,2] = exp(-Uj-qj*potential(t)[1]) - 1
        else
            conc_matrix[index,1] = exp(-Ui-qi*potential[index][1])-1
            conc_matrix[index,2] = exp(-Uj-qj*potential[index][1])-1
            conc_matrix[index,3] = exp(-Uk-qk*potential[index][1])-1
            conc_matrix[index,4] = exp(-Ul-ql*potential[index][1])-1
        end
    end

    if n_salts == 1
        gammaconi = trapz(time, conc_matrix[:,1])
        gammaconj = trapz(time, conc_matrix[:,2])
        gammacon = gammaconi, gammaconj
        tension = sum(gammacon)/(STconst*n_salts)
        pGM = 4Avog * ionS/n_salts*1000 * (gammaconi^2 + gammaconj^2)*1e-20 / (beta*(10e-9)^2)
    
    else
        gammaconi = trapz(time, conc_matrix[:,1])
        gammaconj = trapz(time, conc_matrix[:,2])
        gammaconk = trapz(time, conc_matrix[:,3])
        gammaconl = trapz(time, conc_matrix[:,4])
        gammacon = gammaconi, gammaconj, gammaconk, gammaconl
        tension = sum(gammacon)/(STconst*n_salts)
        pGM = 4Avog * ionS/n_salts*1000 * (gammaconi^2 + gammaconj^2 + gammaconk^2+ gammaconl^2)*1e-20 / (beta*(10e-9)^2)
    end
    return gammacon, tension, pGM
end;
time_range = range(0.01, boundary, 1000001);
sol2 = pGM!(time_range, sol, param, ionS)
# @btime sol2 = pGM!(time_range, sol, param, ionS)

# plot(time_range, sol2[1],
#     labels=["H" "Cl"],
#     ylims=(0,0.8),
#     xlims=(0,12),
#     ylabel="Concentration [M]",
#     xlabel="z [Å]")
# scatter!(time_range, sol2[1],
#     ylims=(-2,8),
#     xlims=(0,12))
