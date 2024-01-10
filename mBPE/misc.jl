function U_beta1(z, ah, kappa, elc, beta, epsilon_o, epsilon_w)
    if z < 1e10ah # while inside hydrated ion radius
        U = 1/(4pi*epsilon_o) * elc^2*beta / (1e-10z*4epsilon_w) * exp(-2kappa*1e-10z) - 3.05
    else # outside of radius
        U = 1/(4pi*epsilon_o) * elc^2*beta / (1e-10z*4epsilon_w) * exp(-2kappa*1e-10z)
    end
end;
U_H(t) = U_beta1(t, 1.97e-10, kappa, elc, beta, epsilon_o, epsilon_w)
W_H    = W(+1, 1.97e-10);
H    = (+1, U_H);


function U_beta2(z, ah, kappa, W)
    if z < 1e10ah # while inside hydrated ion radius
        U = W * 1e10ah/z * exp(-2kappa * (1e-10z - ah)) - 2.1
    else # outside of radius
        U = W * 1e10ah/z * exp(-2kappa * (1e-10z - ah))
    end
end;
W_ClO4 = W(-1, 2.83e-10);
U_ClO4(t) = U_beta2(t, 2.83e-10, kappa, W_ClO4)
ClO4 = (-1, U_ClO4);