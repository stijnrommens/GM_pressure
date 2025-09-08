using Parameters
using Plots
using Optimization
using Zygote
using Logging

include("modules/WMod.jl")

function U_p_base(x, p)
    z, ah, q = p
    x = x[1]
    # h = max(0, ah - z)
    # Avoid acos of value < -1 to avoid errors
    acos_arg = max(-1, 1 - z / ah)
    theta = acos(acos_arg)
    g = 10

    # Bjerrum length in water
    l_B = (q * elc)^2 / (epsilon_w * epsilon_o)
    # Electrostatic self energy costs
    essec_w = pi * x^2 / theta
    if isapprox(theta, pi)
        essec_a = 0
    else
        essec_a = (pi * (1 - x)^2 * epsilon_w * epsilon_o) / (epsilon_o * (pi - theta))
    end
    essec = essec_w + essec_a
    # induced surface inhomogeneity energy cost
    isiec = g * (x - (1 - cos(theta)) / 2)^2
    # if isapprox(essec, isiec, rtol=0.5)
    #     println("Contributions close to equal (z=$z):")
    #     # println("ESSEC:\t$essec")
    #     # println("ISIEC:\t$isiec")
    # end
    U_p = l_B / (2 * ah) * essec + isiec
    return U_p
end

function Utot(distance, kappa, q, ah)
    v = 0.3 * kB * T * 1e30
    U_tot = zeros(Float64, size(distance))
    for (i, z) in enumerate(distance)
        x0 = [0.5]
        p = [z, ah, q]
        optf = OptimizationFunction(U_p_base, AutoZygote())
        prob = OptimizationProblem(optf, x0, p, lb=[0.0], ub=[1.0])
        sol = solve(prob, Optimization.LBFGS())
        if !SciMLBase.successful_retcode(sol)
            @warn "Unsuccesful x-optimization for z = $z"
            println(sol.retcode)
        end
        xmin = sol[1]
        W_factor =  (q * elc)^2 / (2 * epsilon_o * epsilon_w)
        if z >= ah
            # Part for z >= ah
            W_int, W_err = quadgk(k -> Wintegrator(k, kappa, ah, z), 0, 40e10)
            W_z = W_factor * W_int
            # U_p = U_p_base(xmin, p)
            U_s = ((q * elc)^2 / (2 * epsilon_o * epsilon_w * ah))
            U_cav = v * ah^3
        elseif z < ah
            W_z = 1
            U_p = 0
            U_s = U_p
            U_cav = 0
        elseif 0 < z < ah
            # Part for 0 < z < a0
            W_int, W_err = quadgk(k -> Wintegrator(k, kappa, ah, ah), 0, 40e10)
            W_z = W_factor * W_int * z / ah
            # U_p = U_p_base(xmin, p)
            # U_s = U_p
            U_s = (((q * elc)^2 / (2 * epsilon_o * epsilon_w * ah)) * (1 ./ ((acos.(z ./ ah) ./ pi) .+ 1 / 2epsilon_w)))
            U_cav = 1/4 * v * ah^3 * (z ./ ah .+ 1).^2 .* (2 .- z ./ ah)
        end
        U_tot[i] = W_z + U_s + U_cav
    end

    return U_tot
end

const kB = 1.3806503e-23         # Boltzmann constant [J/K]
const elc = 1.6021765e-19        # Elementary charge  [C] = [s.A]
const Avog = 6.0221415e23        # Avogadro constant  [1/mol]
const epsilon_o = 8.85418782e-12 # Vacuum's permittivity [F/m] = [s4.A2/kg.m3]
const epsilon_w = 78.3           # Water's relative permittivity [-] 

# NaAc at high concentration
list = [
    [0.95e3, +1, 2.5e-10, "alpha"],
    [0.95e3, -1, 3.22e-10, "beta"]
]

T = 293.15

ionS = 0.0
for ion in list
    global ionS += ion[1] * ion[2]^2
end
ionS = 0.5 * ionS   # Ionic strength [mol/m3]

kappa::Float64 = sqrt(
    elc^2 * Avog * 2 * ionS
    / (kB * T * epsilon_w * epsilon_o)
)

boundary = 10 / kappa
boundary = 10e-10

# println(boundary)

z = range(1e-20, boundary, length=1000)

UtotNa = Utot(z, kappa, list[1][2], list[1][3]) / (kB * T)
UtotAc = Utot(z, kappa, list[2][2], list[2][3]) / (kB * T)
# println(maximum(UtotNa))
p2 = plot(z*1e10, [UtotNa UtotAc], label=["Na" "Ac"])
xlabel!(p2, "Distance [Å]")
ylabel!(p2, "U(r) / kB T")
ylims!(p2, 0, 70)
