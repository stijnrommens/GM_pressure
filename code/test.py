import julia
julia_path = "/Applications/Julia-1.9.app/Contents/Resources/julia/bin/julia"
julia.install(julia = julia_path)
jl = julia.Julia(runtime = julia_path, compiled_modules=False)

from julia import Main
script_path = "/Users/stijnrommens/GithubProjects/GM_pressure/code/main.jl"
Main.include(script_path)
Main.main()
