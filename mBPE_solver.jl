cd(@__DIR__)
using Pkg
Pkg.activate(".")
using QuadGK, DifferentialEquations, Plots

kB = 1.3806503e-23;  # Boltzmann constant [J/K]
elc = 1.6021765e-19; # Elementary charge  [C]
Avog = 6.0221415e23; # Avogadro constant  [1/mol]
epsilon_o = 8.85418782e-12; # Vacuum's permittivity [F/m]
epsilon_w = 78.3;           # Water's permittivity [F/m]

ionS = 0.1; # M

T = 297.15; # K

     #     Na+,      Cl-,       H+,    ClO4-
q  = [      +1,       -1,       +1,       -1] * elc; # Charge [C]
ah = [2.50e-10, 2.00e-10, 1.97e-10, 2.83e-10];       # Hydrated radius

kappa = sqrt(ionS) / 0.304e-9;
rho_ion = ionS * Avog * 1000;
boundary = 10e10/kappa;
beta = 1 / (kB * T);
STconst = -1e-6*1e10*beta/Avog;


f(k, a) = (k*(sqrt(kappa^2 + k^2)*cosh(k*a) - k*sinh(k*a))) / (sqrt(kappa^2 + k^2)*(sqrt(kappa^2 + k^2)*cosh(k*a) + k*sinh(k*a)));
function W(i)
    factor = beta * q[i]^2 / (2epsilon_w * 4pi * epsilon_o)
    W = quadgk(k -> f(k,ah[i]), 0, 10.0e10)[1] * factor
    return W
end;

function U_Na(z)
    if z < 1e10ah[1] # while inside hydrated ion radius
        U = 1000
    else # outside of radius
        U = W(1)*1e10ah[1]/z * exp(-2kappa * (1e-10z - ah[1]))
    end
    return U
end;

function U_Cl(z)
    if z < 1e10ah[2] # while inside hydrated ion radius
        U = 1000
    else # outside of radius
        U = W(2)*1e10ah[2]/z * exp(-2kappa * (1e-10z - ah[2]))
    end
    return U
end;

function U_H(z)
    if z < 1e10ah[3] # while inside hydrated ion radius
        U= 1/(4pi*epsilon_o) * elc^2*beta / (1e-10z*4epsilon_w) * exp(-2kappa*1e-10z) - 3.05
    else # outside of radius
        U = 1/(4pi*epsilon_o) * elc^2*beta / (1e-10z*4epsilon_w) * exp(-2kappa*1e-10z)
    end
    return U
end;

function U_ClO4(z)
    if z < 1e10ah[4] # while inside hydrated ion radius
        U = W(4) * 1e10ah[4]/z * exp(-2kappa * (1e-10z - ah[4])) - 2.1
    else # outside of radius
        U = W(4) * 1e10ah[4]/z * exp(-2kappa * (1e-10z - ah[4]))
    end
    return U
end;

function mPBE(du, u, p, t)
    i, j = 3, 4
    Ui = U_H(t)
    Uj = U_ClO4(t)
   
    du[1] = -1*u[2]
    du[2] = -1*-1e-20beta * elc^2 / (epsilon_o*epsilon_w) * ( q[i]/elc*rho_ion*exp(-Ui - q[i]/elc * u[1]) + q[j]/elc*rho_ion*exp(-Uj - q[j]/elc * u[1]))
end;

function bc(residual, u, p, t)
    residual[1] = u[end][2] -0.0
    residual[2] = u[1][1] -0.0
end;

tspan = (boundary, 0.1);
bvp = BVProblem(mPBE, bc, [0.0, 0.0], tspan)
sol = solve(bvp, Shooting(Vern7()))

HCl_list    = zeros(Float64, size(sol.u)[1], 2);
function pre_plot!(u_list, sol)
    for (i, solu) in enumerate(sol.u)
        u_list[i,1] = sol.t[i]
        u_list[i,2] = solu[1]
    end
    return u_list
end;
HCl_list    = pre_plot!(HCl_list, sol);

plot(HCl_list[:,1], HCl_list[:,2],
    label=["HCl"],
    ylims=(-0.45,0.45),
    xlims=(0,60),
    ylabel="Potential [V]",
    xlabel="z [Å]"
)