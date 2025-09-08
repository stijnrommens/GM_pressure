# cd(@__DIR__)
# using Pkg
# Pkg.activate(".")
using BoundaryValueDiffEq, OrdinaryDiffEq, Printf#, BenchmarkTools
using Parameters
# import BoundaryValueDiffEq: BVProblem

# TODO Rewrite to use angstrom length scale in all parameters.
# Currently, the system is defined in meter scale, but solved in angstrom scale.
# This means we're always operating close to the floating point limit, resulting
# in unreliable outcomes. Operating directly with angstrom scale units should
# make this easier

include("welcome.jl")
include("constants.jl")
include("modules/CheckMod.jl")
include("modules/WMod.jl")
include("modules/mPBEMod.jl")
include("modules/pGMMod.jl")

using .Check_Mod

@with_kw struct Para
    # beta::Float64
    kB::Float64 = 1.3806503e-23     # Boltzmann constant [kg.m2/s2.K]
    T::Float64
    elc::Float64 = 1.6021765e-19    # Elementary charge [s.A]
    epsilon_o::Float64 = 8.85418782e-12 # Vacuum's permitivity [s4.A2/kg.m3]
    epsilon_w::Float64 = 78.3       # Water's permitivity [-]
    Avog::Float64 = 6.0221415e23    # Avogadro's number [mol-1]
    h::Float64 = 10e-9              # Film thickness between two bubbles [m]
    kappa::Float64                  # Inverse Debye-Hueckel length [1/m]
    # n::Int64
    STconst::Float64                # Constant for surface tension calculations [mol.m.Å/mN.dm3] = 
    ionS::Float64                   # Ionic strength [mol/m3]
    ion_conc::Vector{Float64} = Vector{Float64}[]   # Ion concentration [mol/m3]
    ion_charges::Vector{Int64}  = Vector{Float64}[]
    ion_hyd_radii::Vector{Float64}  = Vector{Float64}[]
    ion_Wcal::Vector{Float64}  = Vector{Float64}[]
    ion_types::Array{String}  = Vector{String}[]
    boundary::Float64
end

function Para_to_angstrom(params)
    @unpack kB, epsilon_o, h, kappa, STconst, boundary = params
    @unpack ionS, ion_conc, ion_hyd_radii, ion_Wcal = params
    angstrom_per_meter = 1e10       # [Å/m]
    meter_per_angstrom = 1e-10      # [m/Å]
    new_params = Para(params;
        boundary = boundary,    # boundary already in Å
        kB = kB * angstrom_per_meter^2, # kB originally in [kg.m2/s.K]
        h = h * angstrom_per_meter,     # h originally in [m]
        kappa = kappa * meter_per_angstrom,     # κ originally in [1/m]
        ionS = ionS * meter_per_angstrom^3,     # ionS originally in [mol/m3]
        STconst = STconst * meter_per_angstrom^2,   # STconst originally in [mol.s2/kg.m2]
        ion_conc = ion_conc * meter_per_angstrom^3,  # ion_conc originally in [mol/m3]
        epsilon_o = epsilon_o * meter_per_angstrom^3,   # epsilon originally in [s4.A2/kg.m3]
        ion_Wcal = ion_Wcal * angstrom_per_meter^2,     # Wcal originally in [kg.m2/s2]
        ion_hyd_radii = ion_hyd_radii * angstrom_per_meter, # ion_hyd_radii originally in [m]
    )
    return new_params
end

const kB = 1.3806503e-23         # Boltzmann constant [J/K]
const elc = 1.6021765e-19        # Elementary charge  [C]
const Avog = 6.0221415e23        # Avogadro constant  [1/mol]
const epsilon_o = 8.85418782e-12 # Vacuum's permittivity [F/m] = [s4.A2/kg.m3]
const epsilon_w = 78.3           # Water's relative permittivity [-] 

function main(list=false; 
    print_flag::Bool=true, T::Float64=297.15, abstol::Float64=1e-3, reltol::Float64=1e-3)
    
    # TODO consider using concentration as [mol/m3]. Makes everything more SI-consistent
    if list == false
        # contains molar concentration, charge, hydrated radius, type {'alpha', 'beta', 'proton'}
        # list = [[0.1, +1, 2.5e-10], [0.1, -1, 2.0e-10]]
        
        # Extra tough concentration of NH4Cl from other work
        # list = [[0.0053, +1, 2.5e-10], [0.0053, -1, 2.0e-10]]
        
        # Extra tough concentration of NaAc from other work
        list = [[0.95e3, +1, 2.5e-10, "alpha"], [0.95e3, -1, 3.22e-10, "beta"]]
    end
    # n = length(list)/2 # Number of salts [-]

    # STbase in [mol/N.m] = [mol.s2/kg.m2]
    STbase = - 1 / (Avog * kB * T)
    # STbase * [N/mN] * [m2/m2] * [Å/m] * [dm3/m3] -> [M.m.Å/mN] = [mol.Å.s2/dm3.kg]
    # STconst = STbase * 1e-3 * 1 * 1e10 * 1e-3
    STconst = STbase

    ionS = 0.0
    for ion in list
        ionS += ion[1] * ion[2]^2
    end
    ionS = 0.5 * ionS   # Ionic strength [mol/m3]
    
    # Inverse Debye-Hückel lenght [1/m]
    # κ is based on ionic strength as indicated in the Wikipedia article
    kappa::Float64 = sqrt(
        elc^2 * Avog * 2 * ionS
        / (kB * T * epsilon_w * epsilon_o)
    )
    
    if kappa == 0.0
        return [0, 0, 0, 0]
    end

    # Boundary defined as 10x Debye-Hückel length [Å]
    boundary = 10 * 1e10 / kappa

    params = Para(
        T = T,
        kappa = kappa,
        STconst = STconst,
        ionS = ionS,
        boundary = boundary
    )

    n_present = 0
    for (i, ion) in enumerate(list)
        if ion[1] > 0.0
            n_present += 1
        end
    end
    ion_conc = zeros(Float64, n_present)
    ion_charges = zeros(Int64, n_present)
    ion_hyd_radii = zeros(Float64, n_present)
    ion_Wcal = zeros(Float64, n_present)
    ion_types = Array{String}(undef, n_present)

    i = 1
    for (_, ion) in enumerate(list)
        if ion[1] == 0.0
            continue
        end
        Wcal = W(ion[2], ion[3], params)
        ion_conc[i] = ion[1]
        ion_charges[i] = ion[2]
        ion_hyd_radii[i] = ion[3]
        ion_Wcal[i] = Wcal
        if last(ion) isa String
            ion_types[i] = ion[4]
        else
            ion_types[i] = "alpha"
        end
        i += 1
    end
    # println(ion_types)
    
    # --- Electrostatic potential ---
    # mPBE_constants = (beta, elc, epsilon_o, epsilon_w, Avog, kappa)
    # mPBE_param = (mPBE_constants, ion_tot, n)
    full_param = Para(params,
        ion_conc = ion_conc,
        ion_charges = ion_charges,
        ion_hyd_radii = ion_hyd_radii,
        ion_Wcal = ion_Wcal,
        ion_types = ion_types)
    # println(eltype(mPBE_param))
    
    full_param = Para_to_angstrom(full_param)

    paramized_fun!(du, u, p, t) = mPBE!(du, u, full_param, t)
    paramized_bc!(residual, u, p, t) = bc!(residual, u, full_param, t)
    tspan = (boundary, 0.0)     # [Å]
    u0  = [0, 0]
    
    bvp_problem = BVProblem(paramized_fun!, paramized_bc!, u0, tspan)
    bvp_sol     = solve(
        bvp_problem,
        # Rodas4P worked for the 'tough' NH4Cl concentration at tolerances of 1e-9
        MultipleShooting(;ode_alg = Rodas4P(autodiff = false), nshoots = 20),   
        # MultipleShooting(;ode_alg = TRBDF2(autodiff = false), nshoots = 20),
        # MultipleShooting(;ode_alg = Rodas5P(autodiff = false), nshoots = 20),
        # MultipleShooting(;ode_alg = ROS3(autodiff = false), nshoots = 20),
        # MultipleShooting(;ode_alg = DFBDF(autodiff = false), nshoots = 20),
        # MultipleShooting(;ode_alg = FBDF(autodiff = false), nshoots = 20),
        # MultipleShooting(;ode_alg = QNDF(autodiff = false), nshoots = 20),
        # Shooting(;ode_alg = AutoVern7(Rodas4(autodiff = false))),
        abstol = abstol,
        reltol = reltol,
        maxiters = 1e7)

    # println(bvp_sol[1], bvp_sol[-1])
    # --- Gibbs-Marangoni pressure ---
    # TODO rewrite time to be distance. There is no time things happening here.
    time_range = range(0.0, boundary, 100001)
    # pGM_constants = (beta, elc, epsilon_o, epsilon_w, Avog, kappa, STconst, ionS, h)
    # pGM_param = (pGM_constants, ion_tot, n)
    pGM_sol = pGM!(time_range, bvp_sol, full_param)

    # println("\nChecks:")
    Check_Mod.check_input(list, print_flag=print_flag)  
    Check_Mod.check_mPBE(bvp_sol, print_flag=print_flag)
    Check_Mod.check_pGM(pGM_sol[2], print_flag=print_flag)
    
    if print_flag == true
        println("\nResults:")
        println("   • Surface excess × q × c   = ", round.(pGM_sol[2]; digits=3), " = ", round.(sum(pGM_sol[2]); digits=3), " Å.mol/m3")
        @printf("   • Surface tension          = %.3e mN/m.M", pGM_sol[3])
        @printf("\n   • Electrostatic potential  = %.3f mV", 1000*last(bvp_sol.u)[1])
        @printf("\n   • Gibbs-Marangoni pressure = %.3f Pa\n", pGM_sol[4])
    end

    return pGM_sol[4], pGM_sol[3], SciMLBase.successful_retcode(bvp_sol)
end
