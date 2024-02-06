# C, i = C_range[1], 1
# include("inputs.jl")
# results[i] = sol2[3]

# C, i = C_range[2], 2
# include("inputs.jl")
# results[i] = sol2[3]

# C, i = C_range[3], 3
# include("inputs.jl")
# results[i] = sol2[3]

# C, i = C_range[4], 4
# include("inputs.jl")
# results[i] = sol2[3]

# C, i = C_range[5], 5
# include("inputs.jl")
# results[i] = sol2[3]

# C, i = C_range[6], 6
# include("inputs.jl")
# results[i] = sol2[3]

# C, i = C_range[7], 7
# include("inputs.jl")
# results[i] = sol2[3]

# C, i = C_range[8], 8
# include("inputs.jl")
# results[i] = sol2[3]

# C, i = C_range[9], 9
# include("inputs.jl")
# results[i] = sol2[3]

# C, i = C_range[10], 10
# include("inputs.jl")
# results[i] = sol2[3]



# C_range = range(1e-10, 0.1, 11)
# results = zeros(length(C_range))

# for (i,c) in enumerate(C_range)
#     C, i = C_range[i], i
#     print(C)
#     print("\n\n")
#     print(ionS)
#     print("\n\n")
#     include("inputs.jl")
#     results[i] = sol2[3]
# end

# results = results .* C_range
# print("\n\n")
# print(results)

# scatter(C_range, results,
#     # ylims=(-0.45,0.45),
#     # xlims=(0,1),
#     ylabel="tension [mN/m.M]",
#     xlabel="concentration [M]"
# )

# print("Hello world")

# ion_list = [#[C, +1, 2.50e-10],  # Na, Ionic concentration [M], Charge [-], Hydrated radius [m] from Levin (2009)
            # [C, +1, 2.50e-10],
            #[C, -1, 2.00e-10]]#,  # Cl, Levin (2010)
            # [C, +1, 2.50e-10]]#,  # NH4, Kielland (1937)
            # [C, -2, 3.79e-10], # SO4, Levin (2010)
            # [2C, +1, 3.15e-10]]; # K
            # [C, +1, 3.15e-10]];
            # [C, +2, 3.63e-10]]; # Mg





# function run_file(concentration)
#     C = concentration
#     include("inputs.jl")
# end
# run_file(0.1)
