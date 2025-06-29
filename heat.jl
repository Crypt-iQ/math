""" The heat equation! 
    dudt = ∇²u
   
    Given a matrix `ic` that describes the temperature at each spatial point for fixed
    time t=0, we can use the heat equation to integrate over time and make a nice
    time-series visual.

    For simplicity, choose u(0,y,t) = u(x,0,t) = u(L,y,t) = u(x,L,t) = 0 as the boundary
    condition for t >= 0.

    This requires the matrix `ic` representing the initial condition to already
    satisfy the boundary condition.

    Currently this is 2d as the Laplacian implementation is 2d.
"""

function timestep(ic::Matrix{T}) where {T<:AbstractFloat}
    m,n = size(ic)
    step = zeros(m-2,n-2)

    # Calculate dudt
    ∇²(step, ic)

    # Add dudt to ic
    @inbounds for i ∈ 1:m, j ∈ 1:n
        if i == 1 || i == m || j == 1 || j == n
            # Do nothing at the boundaries
        else
            ic[i,j] += step[i-1,j-1]
        end
    end
end
