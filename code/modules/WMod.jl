# module W_Mod
using QuadGK
using Logging

function Wintegrator(k::Float64, kappa::Float64, ah::Float64)
    s = sqrt(kappa^2 + k^2)
    # [1/m2 / 1/m2] -> [-]
    return (
        (k * (s * cosh(k * ah) - k * sinh(k * ah)))
        / (s * (s * cosh(k * ah) + k * sinh(k * ah))
    ))
end

function W(q::Real, ah::Float64, params)
    """Calculate the energy to bring ion out of Gibbs dividing surface

    Args:
        q (float): Valency of the ion [-]
        ah (float): Hydrated radius of the ion [m]

    Returns:
        W (float): Energy to bring ion out of Gibbs dividing surface [m?]
    """
    @unpack kappa, elc, epsilon_w, epsilon_o = params
    
    # kappa = sqrt(
    #     8pi * q * c_ion * 1000 / 
    #     (epsilon_w * epsilon_o * kB * T)
    # )
    # factor in [kg.m3/s2]
    factor = (q * elc)^2 / (8pi * epsilon_w * epsilon_o)
    # TODO should this 10e10 be based on some boundary value?
    Wint, Werr = quadgk(k -> Wintegrator(k, kappa, ah), 0, 10.0e10; rtol=1e-10)   # [1/m]
    if Werr > 1e-10
        @warn "W integration inaccurate: err > 1e-10"
        # println("Warning: W integration inaccurate")
    end
    W = factor * Wint   # [kg.m2/s2]
    return W
end

# end # end of module