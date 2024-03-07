cd(@__DIR__)
using Pkg
Pkg.activate(".")
using BoundaryValueDiffEq, OrdinaryDiffEq, Printf#, BenchmarkTools

include("welcome.jl")
include("constants.jl")
include("modules/CheckMod.jl")
include("modules/WMod.jl")
include("modules/mPBEMod.jl")
include("modules/pGMMod.jl")

using .Check_Mod

function main(list=false; 
    print_flag::Bool=true, T::Float64=297.15, h::Float64=10e-9, ionS::Float64=0.0, ion_tot::Tuple=(),
    elc::Float64=elc, Avog::Float64=Avog, epsilon_w::Float64=epsilon_w, epsilon_o::Float64=epsilon_o,
    abstol::Float64=1e-15, reltol::Float64=1e-15, fun=mPBE!, bc=bc!)
    
    if list == false
        list = [[0.1, +1, 2.5e-10], [0.1, -1, 2.0e-10]]
    end
    n = length(list)/2 # Number of salts [-]

    beta = 1 / (kB * T)            # Thermodynamic beta [1/J] = [s2/kg.m2] = [1/N.m]
    STconst = -1e-6*1e10*beta/Avog # Surfce tension constant [mol/N.m] -> [M.m.Å/mN]

    for ion in list
        ionS += ion[1] * ion[2]^2
    end
    ionS = 0.5*ionS # Ionic strength [M]

    kappa = sqrt(2 * elc^2 * Avog * ionS*1000 * beta / (epsilon_w*epsilon_o)) # Inverse Debye-Hückel lenght [1/m]
    boundary = 10e10/kappa         # 10x Debye-Hückel lenght [Å]

    for ion in list
        Wcal = W(ion[2], ion[3]; kappa=kappa, beta=beta, elc=elc, epsilon_w=epsilon_w, epsilon_o=epsilon_o)
        if last(ion) isa String
            ion_tot = tuple(ion_tot..., (ion[1], ion[2], ion[3], ion[4], Wcal))
        else
            ion_tot = tuple(ion_tot..., (ion[1], ion[2], ion[3], "alpha", Wcal))
        end
    end
    
    # --- Electrostatic potential ---
    mPBE_constants = (beta, elc, epsilon_o, epsilon_w, Avog, kappa)
    mPBE_param     = (mPBE_constants, ion_tot, n)
    tspan = (boundary, 0.0)
    u0  = [0.0, 0.0]
    bvp_problem = BVProblem(fun, bc, u0, tspan, mPBE_param)
    bvp_sol     = solve(bvp_problem, Shooting(RadauIIA5(autodiff=false)), abstol=abstol, reltol=reltol)

    # --- Gibbs-Marangoni pressure ---
    time_range = range(0.0, boundary, 100001)
    pGM_constants = (STconst, Avog, ionS, beta, h, kappa)
    pGM_param = (pGM_constants, ion_tot, n)
    pGM_sol = pGM!(time_range, bvp_sol, pGM_param)

    if print_flag == true
        println("\nChecks:")
        Check_Mod.check_input(list)  
        Check_Mod.check_mPBE(bvp_sol)
        Check_Mod.check_pGM(pGM_sol[2])

        println("\nResults:")
        println("   • Surface excess × q × c   = ", round.(pGM_sol[2]; digits=3), " = ", round.(sum(pGM_sol[2]); digits=3), " Å")
        @printf("   • Surface tension          = %.3e mN/m.M", pGM_sol[3])
        @printf("\n   • Electrostatic potential  = %.3f mV", 1000*last(bvp_sol.u)[1])
        @printf("\n   • Gibbs-Marangoni pressure = %.3f Pa", pGM_sol[4])
    end

    return pGM_sol[4]
end
