cd(@__DIR__)
using Pkg
Pkg.activate(".")
using BoundaryValueDiffEq, OrdinaryDiffEq, Printf#, BenchmarkTools
using Parameters
import BoundaryValueDiffEq: BVProblem

include("welcome.jl")
include("constants.jl")
include("modules/CheckMod.jl")
include("modules/WMod.jl")
include("modules/mPBEMod.jl")
include("modules/pGMMod.jl")

using .Check_Mod

@with_kw struct mPBE_Para
    beta::Float64
    elc::Float64 = 1.6021765e-19    # Elementary charge [C]
    epsilon_o::Float64 = 8.85418782e-12 # Vacuum's permitivity [F/m]
    epsilon_w::Float64 = 78.3       # Water's permitivity [F/m]
    Avog::Float64 = 6.0221415e23    # Avogadro's number [mol-1]
    h::Float64 = 10e-9  
    kappa::Float64
    n::Int64
    STconst::Float64
    ionS::Float64
    ion_conc::Vector{Float64}
    ion_charges::Vector{Int64}
    ion_hyd_radii::Vector{Float64}
    ion_Wcal::Vector{Float64}
    ion_types::Array{String}
end


function main(list=false; 
    print_flag::Bool=true, T::Float64=297.15, h::Float64=10e-9, ionS::Float64=0.0, ion_tot::Tuple=(),
    elc::Float64=elc, Avog::Float64=Avog, epsilon_w::Float64=epsilon_w, epsilon_o::Float64=epsilon_o,
    abstol::Float64=1e-12, reltol::Float64=1e-12, fun=mPBE!, bc=bc!)
    
    if list == false
        # contains molar concentration, charge, hydrated radius, type {'alpha', 'beta', 'proton'}
        list = [[0.1, +1, 2.5e-10], [0.1, -1, 2.0e-10]]
    end
    n = length(list)/2 # Number of salts [-]

    beta::Float64 = 1 / (kB * T)            # Thermodynamic beta [1/J] = [s2/kg.m2] = [1/N.m]
    STconst = -1e-6 * 1e10 * beta / Avog    # Surface tension constant [mol/N.m] -> [M.m.Å/mN]

    for ion in list
        ionS += ion[1] * ion[2]^2
    end
    ionS = 0.5*ionS # Ionic strength [M]

    kappa::Float64 = sqrt(2 * elc^2 * Avog * ionS*1000 * beta / (epsilon_w*epsilon_o)) # Inverse Debye-Hückel lenght [1/m]
    if kappa == 0.0
        return [0, 0, 0, 0]
    end
    boundary = 10e10/kappa         # 10x Debye-Hückel lenght [Å]

    ion_conc = zeros(Float64, size(list, 1))
    ion_charges = zeros(Int64, size(list, 1))
    ion_hyd_radii = zeros(Float64, size(list, 1))
    ion_Wcal = zeros(Float64, size(list, 1))
    ion_types = Array{String}(undef, size(list, 1))

    for (i, ion) in enumerate(list)
        Wcal = W(ion[2], ion[3]; kappa=kappa, beta=beta, elc=elc, epsilon_w=epsilon_w, epsilon_o=epsilon_o)
        ion_conc[i] = ion[1]
        ion_charges[i] = ion[2]
        ion_hyd_radii[i] = ion[3]
        ion_Wcal[i] = Wcal
        if last(ion) isa String
            ion_types[i] = ion[4]
        else
            ion_types[i] = "alpha"
        end
    end
    
    # --- Electrostatic potential ---
    # mPBE_constants = (beta, elc, epsilon_o, epsilon_w, Avog, kappa)
    # mPBE_param = (mPBE_constants, ion_tot, n)
    mPBE_param = mPBE_Para(
        beta = beta,
        kappa = kappa,
        n = n,
        STconst = STconst,
        ionS = ionS,
        ion_conc = ion_conc,
        ion_charges = ion_charges,
        ion_hyd_radii = ion_hyd_radii,
        ion_Wcal = ion_Wcal,
        ion_types = ion_types)
    # println(eltype(mPBE_param))

    paramized_fun(du, u, p, t) = fun(du, u, mPBE_param, t)
    tspan = (boundary, 0.0)
    u0  = [0.0, 0.0]
    # TODO rewrite for use of MultipleShooting - Need mPBE_param to have concrete 
    # type in eltype(). That means rewriting mPBE and maybe bc to have a less
    # convoluted input in parameters
    bvp_problem = BVProblem(paramized_fun, bc, u0, tspan)
    bvp_sol     = solve(
        bvp_problem,
        MultipleShooting(;ode_alg = Rodas4P(autodiff = false), nshoots = 10),
        # Shooting(;ode_alg = AutoVern7(Rodas4(autodiff = false))),
        abstol = abstol,
        reltol = reltol)

    # println(bvp_sol[1], bvp_sol[-1])
    # --- Gibbs-Marangoni pressure ---
    time_range = range(0.0, boundary, 100001)
    pGM_constants = (beta, elc, epsilon_o, epsilon_w, Avog, kappa, STconst, ionS, h)
    pGM_param = (pGM_constants, ion_tot, n)
    pGM_sol = pGM!(time_range, bvp_sol, mPBE_param)

    # println("\nChecks:")
    Check_Mod.check_input(list, print_flag=print_flag)  
    Check_Mod.check_mPBE(bvp_sol, print_flag=print_flag)
    Check_Mod.check_pGM(pGM_sol[2], print_flag=print_flag)
    
    if print_flag == true
        println("\nResults:")
        println("   • Surface excess × q × c   = ", round.(pGM_sol[2]; digits=3), " = ", round.(sum(pGM_sol[2]); digits=3), " Å")
        @printf("   • Surface tension          = %.3e mN/m.M", pGM_sol[3])
        @printf("\n   • Electrostatic potential  = %.3f mV", 1000*last(bvp_sol.u)[1])
        @printf("\n   • Gibbs-Marangoni pressure = %.3f Pa\n", pGM_sol[4])
    end

    return pGM_sol[4], pGM_sol[3], SciMLBase.successful_retcode(bvp_sol)
end
