# module W_Mod
using QuadGK
using Logging
using Printf

function Wintegrator(k::Union{Vector{Float64},Float64}, kappa::Float64,
                     ah::Float64, z::Union{Vector{Float64}, Float64, Nothing})
    s = sqrt.(kappa^2 .+ k.^2)
    if z === nothing
        z = 1 ./ k
    end
    # [1/m2 / 1/m2] -> [-]
    return (
        exp.(-2s .* (z .- ah)) .*
        (k .* (s .* cosh.(k .* ah) .- k .* sinh.(k .* ah)))
        ./ (s .* (s .* cosh.(k .* ah) .+ k .* sinh.(k .* ah))
    ))
end

function W(q::Real, ah::Float64, c_b::Float64, params)
    """Calculate the energy to bring ion out of Gibbs dividing surface

    Args:
        q (float): Valency of the ion [-]
        ah (float): Hydrated radius of the ion [m]

    Returns:
        W (float): Energy to bring ion out of Gibbs dividing surface [m?]
    """
    @unpack kappa, elc, epsilon_w, epsilon_o = params
    
    # factor in [kg.m3/s2]
    factor = (q * elc)^2 / (2 * epsilon_w * epsilon_o)
    # TODO should this 10e10 be based on some boundary value?
    Wint, Werr = quadgk(k -> Wintegrator(k, kappa, ah), 0, 20.0e10; rtol=1e-10, atol=1e-10)   # [1/m]
    if Werr > 1e-10
        @warn "W integration inaccurate: err= $Werr > 1e-10"
        println("Wint value: $(@sprintf("%.3f", Wint)) J")
        println("Wint value: $(@sprintf("%.3f", Wint*1e-20)) kg.Å2/s2")
    end
    W = factor * Wint   # [kg.m2/s2]
    return W
end

# end # end of module