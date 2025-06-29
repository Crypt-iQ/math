""" These functions are copied from milankl/ShallowWaters.jl """

""" dx of u """
function dx(dudx::Matrix{T},u::Matrix{T}) where {T<:AbstractFloat}
    m,n = size(dudx)
    @boundscheck (m+1,n) == size(u) || throw(BoundsError())

    @inbounds for i ∈ 1:m, j ∈ 1:n
        dudx[i,j] = u[i+1,j] - u[i,j]
    end
end

""" dy of u """
function dy(dudy::Matrix{T},u::Matrix{T}) where {T<:AbstractFloat}
    m,n = size(dudy)
    @boundscheck (m,n+1) == size(u) || throw(BoundsError())

    @inbounds for i ∈ 1:m, j ∈ 1:n
        dudy[i,j] = u[i,j+1] - u[i,j]
    end
end

""" 2nd order centered Laplace-operator d/dx^2 + d/dy^2 
    This is the second partial """
function ∇²(du::Matrix{T},u::Matrix{T}) where {T<:AbstractFloat}
    m,n = size(du)
    @boundscheck (m+2,n+2) == size(u) || throw(BoundsError())

    @inbounds for i ∈ 1:m
        for j ∈ 1:n
            #       1
            #    1 -4 1
            #       1
            ui1j1 = u[i+1,j+1]
            du[i,j] = ((u[i+2,j+1] - ui1j1) - (ui1j1 - u[i,j+1])) + ((u[i+1,j+2] - ui1j1) - (ui1j1 - u[i+1,j]))
        end
    end
end

""" Exposed to be able to use in REPL """
function Laplace(du::Matrix{T},u::Matrix{T}) where {T<:AbstractFloat}
    ∇²(du,u)
end
