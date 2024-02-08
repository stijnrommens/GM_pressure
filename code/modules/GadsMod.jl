# module Gads_Mod

function U_alpha(t, ah::Float64, W::Float64; kappa::Float64)
    """Calculate the Gibbs adsorption energy of an α-ion

    Args:
        t (float): Distance from G/L-interface [m]
        ah (float): Hydrated radius of the ion [m]
        kappa (float): Inverse Debye-Hückel lenght [1/m]
        W (float): Energy to bring ion out of Gibbs dividing surface [m?]

    Returns:
        U (float): Gibbs adsorption energy [J]
    """
    if t < 1e10ah
        U = 1000
    else
        U = W*1e10ah/t * exp(-2kappa * (1e-10t - ah))
    end
    return U
end

# end # end of module