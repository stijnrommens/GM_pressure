cd(@__DIR__)
using Pkg
Pkg.activate(".")
Pkg.instantiate()
using QuadGK, DifferentialEquations, Plots, BenchmarkTools, StaticArrays, Profile

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

factor_i, Wi, Ui = 0.0, 0.0, 0.0
factor_j, Wj, Uj = 0.0, 0.0, 0.0

function mPBE!(du, u, p, t)
    i, j = 1,2
    Ui = U_Na(t)
    Uj = U_Cl(t)

    # factor_i = beta * (q[i]*elc)^2 / (2epsilon_w * 4pi * epsilon_o)
    # Wi = quadgk(k -> f(k,ah[i]), 0, 10.0e10)[1] * factor_i
    # if t < 1e10ah[i] # while inside hydrated ion radius
    #     Ui = 1000
    # else # outside of radius
    #     Ui = Wi*1e10ah[i]/t * exp(-2kappa * (1e-10t - ah[i]))
    # end

    # factor_j = beta * (q[j]*elc)^2 / (2epsilon_w * 4pi * epsilon_o)
    # Wj = quadgk(k -> f(k,ah[j]), 0, 10.0e10)[1] * factor_j
    # if t < 1e10ah[j] # while inside hydrated ion radius
    #     Uj = 1000
    # else # outside of radius
    #     Uj = Wj*1e10ah[j]/t * exp(-2kappa * (1e-10t - ah[j]))
    # end
    Ui = 0.0
    Uj = 0.0

    du[1] = -u[2]
    du[2] = 1e-20beta * elc^2 / (epsilon_o*epsilon_w) * rho_ion/2*( q[i]*exp(-Ui - q[i] * u[1]) + q[j]*exp(-Uj - q[j] * u[1]))
    nothing
end;

function bc(residual, u, p, t)
    residual[1] = u[end][2] -0.0
    residual[2] = u[1][1] -0.0
end;

param = 1,2;
bvp = BVProblem(mPBE!, bc, [0.0, 0.0], tspan);
sol = solve(bvp, Shooting(RadauIIA5(autodiff=false)))#, save_everystep = false)
@btime solve(bvp, Shooting(RadauIIA5(autodiff=false)), save_everystep = false)

function test()
    bvp = BVProblem(mPBE, bc, [0.0, 0.0], tspan, param);
    sol = solve(bvp, Shooting(RadauIIA5(autodiff=false)))
end

@allocated test()


###
NaCl_list   = zeros(Float64, size(NaCl_sol.u)[1], 2);
function pre_plot!(u_list, sol)
    for (i, solu) in enumerate(sol.u)
        u_list[i,1] = sol.t[i]
        u_list[i,2] = solu[1]
    end
    return u_list
end;
NaCl_list   = pre_plot!(NaCl_list, NaCl_sol);

plot(NaCl_list[:,1], NaCl_list[:,2],
    label=["NaCl" "HCl" "NaClO4" "HClO4" "HCl + NaClO4"],
    ylims=(-0.45,0.45),
    xlims=(0,60),
    ylabel="Potential [V]",
    xlabel="z [Å]"
)