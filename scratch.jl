function Uhard(q, a0)
    return 3 * (q * elc)^2 / (4 * eps_o * eps_w * a0)
end

function Usoft(q, a0)
    return (q * elc) ^ 2 / (2 * log(2) * eps_o * eps_w * a0)
end

function get_delta(a0, q)
    return (Uhard(q, a0) - Usoft(q, a0)) / (kB * T)
end

const elc = 1.6e-19
const eps_o = 8.9e-12
const eps_w = 78
const kB = 1.4e-23

T = 298

q_val = [1, 2]
a_val = [1e-10, 2e-10, 3e-10, 20e-10]

for (q, a) in Base.product(q_val, a_val)
    delta = get_delta(a, q)
    println("a0 = $a, q = $q.\t$delta")
end
