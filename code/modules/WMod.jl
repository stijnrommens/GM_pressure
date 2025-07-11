# module W_Mod
using QuadGK

function Wintegrator(k::Float64, kappa::Float64, ah::Float64)
    # [1/m2 / 1/m2] -> [-]?
    s = sqrt(kappa^2 + k^2)
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
    
    factor = (q * elc)^2 / (8pi * epsilon_w * epsilon_o) # Squared radius of sphere [m2]?
    Wint = quadgk(k -> Wintegrator(k, kappa, ah), 0, 10.0e10)[1]
    W = factor * Wint   # [m]?
    return W
end

# end # end of module