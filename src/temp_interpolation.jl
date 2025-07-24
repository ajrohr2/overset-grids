function eight_node_quadratic(point::Vector{Float64}, x_nodes::Vector{Float64}, y_nodes::Vector{Float64}, function_vals::Vector{Float64}, ξη::Vector{Float64}, N, dN_dξ, dN_dη, dNξ::Vector{Float64}, dNη::Vector{Float64}, R::Vector{Float64}, J::Matrix{Float64}, Δξη::Vector{Float64}, N_vals::Vector{Float64})
    # We assume the nodes are in counterclockwise order in the vectors, with node 1 bottom left, node 2 bottom right, node 3 top right, node 4 top left, node 5 bottom middle, node 6 right middle, node 7 top middle, node 8 left middle

    # Initialize the solved point
    ξη[1] = ξη[2] = 0

    @inbounds for i in 1:8
        dNξ[i] = dN_dξ[i](ξη[1], ξη[2])
        dNη[i] = dN_dη[i](ξη[1], ξη[2])
    end

    x = zero(Float64)
    y = zero(Float64)
    @inbounds for i in 1:8
        x += N[i](ξη[1], ξη[2]) * x_nodes[i]
        y += N[i](ξη[1], ξη[2]) * y_nodes[i]
    end

    R[1] = point[1] - x
    R[2] = point[2] - y

    while sqrt(R[1]^2 + R[2]^2) ≥ 1e-12
        J[1,1] = J[2,1] = J[1,2] = J[2,2] = 0.0
        @inbounds for i in 1:8
            J[1,1] += dNξ[i] * x_nodes[i]
            J[1,2] += dNη[i] * x_nodes[i]
            J[2,1] += dNξ[i] * y_nodes[i]
            J[2,2] += dNη[i] * y_nodes[i]
        end

        detJ = (J[1,1] * J[2,2]) - (J[1,2] * J[2,1])
        Δξη[1] = (J[2,2] * R[1] - J[1,2] * R[2]) / detJ
        Δξη[2] = (-J[2,1] * R[1] + J[1,1] * R[2]) / detJ

        ξη[1] += Δξη[1]
        ξη[2] += Δξη[2]
        ξ, η = ξη[1], ξη[2]

        @inbounds for i in 1:8
            N_vals[i] = N[i](ξ, η)
            dNξ[i] = dN_dξ[i](ξ, η)
            dNη[i] = dN_dη[i](ξ, η)
        end

        x = y = 0.0
        @inbounds for i in 1:8
            x += N_vals[i] * x_nodes[i]
            y += N_vals[i] * y_nodes[i]
        end

        R[1] = point[1] - x
        R[2] = point[2] - y
    end

    f_interp = 0.0
    @inbounds for i in 1:8
        f_interp += N_vals[i] * function_vals[i]
    end
    return f_interp
end
