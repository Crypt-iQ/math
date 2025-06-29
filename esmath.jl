module esmath

    export dx, dy, ∇², Laplace, timestep

    include("gradients.jl")
    include("heat.jl")

end
