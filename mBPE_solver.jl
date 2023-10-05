cd(@__DIR__)
using Pkg
Pkg.activate(".")
Pkg.instantiate()
using QuadGK, DifferentialEquations, Plots

kB = 1.3806503e-23;  # Boltzmann constant [J/K]
elc = 1.6021765e-19; # Elementary charge  [C]
Avog = 6.0221415e23; # Avogadro constant  [1/mol]
epsilon_o = 8.85418782e-12; # Vacuum's permittivity [F/m]
epsilon_w = 78.3;           # Water's permittivity [F/m]

ionS = 0.1; # M

T = 297.15; # K

     #     Na+,      Cl-,       H+,    ClO4-
q  = [      +1,       -1,       +1,       -1]; # Charge [C]
ah = [2.50e-10, 2.00e-10, 1.97e-10, 2.83e-10];       # Hydrated radius

# kappa = sqrt(ionS) / 0.304e-9
rho_ion = ionS * Avog * 1000;
beta = 1 / (kB * T);
STconst = -1e-6*1e10*beta/Avog;
kappa = sqrt(2000*elc^2*Avog*ionS*beta/(epsilon_w*epsilon_o));
boundary = 10e10/kappa;
tspan = (boundary, 0.1);


f(k, a) = (k*(sqrt(kappa^2 + k^2)*cosh(k*a) - k*sinh(k*a))) / (sqrt(kappa^2 + k^2)*(sqrt(kappa^2 + k^2)*cosh(k*a) + k*sinh(k*a)));
function W(i)
    factor = beta * (q[i]*elc)^2 / (2epsilon_w * 4pi * epsilon_o)
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
    i, j, k, l = p

    if i == 1
        Ui = U_Na(t)
    end
    if k == 1
        Uk = U_Na(t)
    end

    if i == 3
        Ui = U_H(t)
    end
    if k == 3
        Uk = U_H(t)
    end

    if j == 2
        Uj = U_Cl(t)
    end
    if l == 2
        Ul = U_Cl(t)
    end

    if j == 4
        Uj = U_ClO4(t)
    end
    if l == 4
        Ul = U_ClO4(t)
    end

    du[1] = -u[2]
    if k==0 && l==0
        du[2] = 1e-20beta * elc^2 / (epsilon_o*epsilon_w) * rho_ion/1*( q[i]*exp(-Ui - q[i] * u[1]) + q[j]*exp(-Uj - q[j] * u[1]))
    else
        du[2] = 1e-20beta * elc^2 / (epsilon_o*epsilon_w) * rho_ion/2*( q[i]*exp(-Ui - q[i] * u[1]) + q[j]*exp(-Uj - q[j] * u[1]) + q[k]*exp(-Uk - q[k] * u[1]) + q[l]*exp(-Ul - q[l] * u[1]))
    end
end;

function bc(residual, u, p, t)
    residual[1] = u[end][2] -0.0
    residual[2] = u[1][1] -0.0
end;

function solver(param)
    bvp = BVProblem(mPBE, bc, [0.0, 0.0], tspan, param);
    sol = solve(bvp, Shooting(RadauIIA5(autodiff=false)))#, reltol=1e-5, abstol=1e-6)
    return sol
end

NaCl   = 1,2,0,0;
HCl    = 3,2,0,0;
NaClO4 = 1,4,0,0;
HClO4  = 3,4,0,0;
HCl_NaClO4 = 3,2,1,4;

NaCl_sol   = solver(NaCl)
HCl_sol    = solver(HCl)
NaClO4_sol = solver(NaClO4)
HClO4_sol  = solver(HClO4)
HCl_NaClO4_sol = solver(HCl_NaClO4)

NaCl_list   = zeros(Float64, size(NaCl_sol.u)[1], 2);
HCl_list    = zeros(Float64, size(HCl_sol.u)[1], 2);
NaClO4_list = zeros(Float64, size(NaClO4_sol.u)[1], 2);
HClO4_list  = zeros(Float64, size(HClO4_sol.u)[1], 2);
HCl_NaClO4_list = zeros(Float64, size(HCl_NaClO4_sol.u)[1], 2);
function pre_plot!(u_list, sol)
    for (i, solu) in enumerate(sol.u)
        u_list[i,1] = sol.t[i]
        u_list[i,2] = solu[1]
    end
    return u_list
end;
NaCl_list   = pre_plot!(NaCl_list, NaCl_sol);
HCl_list    = pre_plot!(HCl_list, HCl_sol);
NaClO4_list = pre_plot!(NaClO4_list, NaClO4_sol);
HClO4_list  = pre_plot!(HClO4_list, HClO4_sol);
HCl_NaClO4_list = pre_plot!(HCl_NaClO4_list, HCl_NaClO4_sol);


plot([NaCl_list[:,1], HCl_list[:,1], NaClO4_list[:,1], HClO4_list[:,1], HCl_NaClO4_list[:,1]], 
     [NaCl_list[:,2], HCl_list[:,2], NaClO4_list[:,2], HClO4_list[:,2], HCl_NaClO4_list[:,2]],
    label=["NaCl" "HCl" "NaClO4" "HClO4" "HCl + NaClO4"],
    ylims=(-0.45,0.45),
    xlims=(0,60),
    ylabel="Potential [V]",
    xlabel="z [Å]"
)