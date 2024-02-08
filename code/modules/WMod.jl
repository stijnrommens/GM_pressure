# module W_Mod
using QuadGK

function W(q::Real, ah::Float64; kappa::Float64, beta::Float64, elc::Float64, epsilon_w::Float64, epsilon_o::Float64)
    """Calculate the energy to bring ion out of Gibbs dividing surface

    Args:
        q (float): Valency of the ion [-]
        ah (float): Hydrated radius of the ion [m]

    Returns:
        W (float): Energy to bring ion out of Gibbs dividing surface [m?]
    """
    factor = beta * (q*elc)^2 / (2epsilon_w * 4pi * epsilon_o) # Squared radius of sphere [m2]?
    f(k) = (k*(sqrt(kappa^2 + k^2)*cosh(k*ah) - k*sinh(k*ah))) / (sqrt(kappa^2 + k^2)*(sqrt(kappa^2 + k^2)*cosh(k*ah) + k*sinh(k*ah))) # [1/m2 / 1/m2] -> [-]?
    # W = quadgk(k -> f(k), 0, 10.0e10)[1] * factor # [m]?
    W = quadgk(k -> f(k), 0, 10.0e10)[1] * factor # [m]?
    return W
end

# end # end of module