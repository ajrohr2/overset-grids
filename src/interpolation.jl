using LinearAlgebra, StaticArrays
# DO NOT PASS IN HALO CELLS. ONLY THE VALID DOMAIN

# We first need to map a given interpolation point to computational space of the background mesh
function eight_node_quadratic(point::Vector{Float64}, x_nodes::Vector{Float64}, y_nodes::Vector{Float64}, function_vals::Vector{Float64}, ξη::Vector{Float64}, N, dN_dξ, dN_dη, dNξ::Vector{Float64}, dNη::Vector{Float64}, R::Vector{Float64}, J::Matrix{Float64}, Δξη::Vector{Float64}, N_vals::Vector{Float64})
    # We assume the nodes are in counterclockwise order in the vectors, with node 1 bottom left, node 2 bottom right, node 3 top right, node 4 top left, node 5 bottom middle, node 6 right middle, node 7 top middle, node 8 left middle

    # Initialize the solved point
    ξη[1] = ξη[2] = 0

    @inbounds for i in 1:8
        dNξ[i] = dN_dξ(i, ξη[1], ξη[2])
        dNη[i] = dN_dη(i, ξη[1], ξη[2])
    end

    x = zero(Float64)
    y = zero(Float64)
    @inbounds for i in 1:8
        x += N(i, ξη[1], ξη[2]) * x_nodes[i]
        y += N(i, ξη[1], ξη[2]) * y_nodes[i]
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
            N_vals[i] = N(i, ξ, η)
            dNξ[i] = dN_dξ(i, ξ, η)
            dNη[i] = dN_dη(i, ξ, η)
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
        f_interp += N(i, ξη[1], ξη[2]) * function_vals[i]
    end
    return f_interp
end

function extract_patch!(point, coordinates_x, coordinates_y, kdtree, x_nodes, y_nodes, nodes)
    # First find the closest point in Euclidean space. This uses a KDTree approach
    idxs, dists = knn(kdtree, point, 1)
    i, j = CartesianIndices(coordinates_x)[idxs[1]].I # This is the "center point" of the patch we want

    @inbounds if i == 1 
        if j == 1
            nodes[1] = CartesianIndex(i, j)
            nodes[2] = CartesianIndex(i+2, j)
            nodes[3] = CartesianIndex(i+2, j+2)
            nodes[4] = CartesianIndex(i, j+2)
            nodes[5] = CartesianIndex(i+1, j)
            nodes[6] = CartesianIndex(i+1, j+1)
            nodes[7] = CartesianIndex(i+1, j+2)
            nodes[8] = CartesianIndex(i, j+1)
        elseif j == size(coordinates_x)[2]
            nodes[1] = CartesianIndex(i, j-2)
            nodes[2] = CartesianIndex(i+2, j-2)
            nodes[3] = CartesianIndex(i+2, j)
            nodes[4] = CartesianIndex(i, j)
            nodes[5] = CartesianIndex(i+1, j-2)
            nodes[6] = CartesianIndex(i+2, j-1)
            nodes[7] = CartesianIndex(i+1, j)
            nodes[8] = CartesianIndex(i, j-1)
        else
            nodes[1] = CartesianIndex(i, j-1)
            nodes[2] = CartesianIndex(i+2, j-1)
            nodes[3] = CartesianIndex(i+2, j+1)
            nodes[4] = CartesianIndex(i, j+1)
            nodes[5] = CartesianIndex(i+1, j-1)
            nodes[6] = CartesianIndex(i+1, j)
            nodes[7] = CartesianIndex(i+1, j+1)
            nodes[8] = CartesianIndex(i, j)
        end
    elseif i == size(coordinates_x)[1]
        if j == 1
            nodes[1] = CartesianIndex(i-2, j)
            nodes[2] = CartesianIndex(i, j)
            nodes[3] = CartesianIndex(i, j+2)
            nodes[4] = CartesianIndex(i-2, j+2)
            nodes[5] = CartesianIndex(i-1, j)
            nodes[6] = CartesianIndex(i, j+1)
            nodes[7] = CartesianIndex(i-1, j+2)
            nodes[8] = CartesianIndex(i, j+1)
        elseif j == size(coordinates_x)[2]
            nodes[1] = CartesianIndex(i-2, j-2)
            nodes[2] = CartesianIndex(i, j-2)
            nodes[3] = CartesianIndex(i, j)
            nodes[4] = CartesianIndex(i-2, j)
            nodes[5] = CartesianIndex(i-1, j-2)
            nodes[6] = CartesianIndex(i, j-1)
            nodes[7] = CartesianIndex(i-1, j)
            nodes[8] = CartesianIndex(i-2, j-1)
        else
            nodes[1] = CartesianIndex(i-2, j-1)
            nodes[2] = CartesianIndex(i, j-1)
            nodes[3] = CartesianIndex(i, j+1)
            nodes[4] = CartesianIndex(i-2, j+1)
            nodes[5] = CartesianIndex(i-1, j-1)
            nodes[6] = CartesianIndex(i, j)
            nodes[7] = CartesianIndex(i-1, j+1)
            nodes[8] = CartesianIndex(i-2, j)
        end
    else
        if j == 1
            nodes[1] = CartesianIndex(i-1, j)
            nodes[2] = CartesianIndex(i+1, j)
            nodes[3] = CartesianIndex(i+1, j+2)
            nodes[4] = CartesianIndex(i-1, j+2)
            nodes[5] = CartesianIndex(i, j)
            nodes[6] = CartesianIndex(i+1, j+1)
            nodes[7] = CartesianIndex(i, j+2)
            nodes[8] = CartesianIndex(i-1, j+1)
        elseif j == size(coordinates_x)[2]
            nodes[1] = CartesianIndex(i-1, j-2)
            nodes[2] = CartesianIndex(i+1, j-2)
            nodes[3] = CartesianIndex(i+1, j)
            nodes[4] = CartesianIndex(i-1, j)
            nodes[5] = CartesianIndex(i, j-2)
            nodes[6] = CartesianIndex(i+1, j-1)
            nodes[7] = CartesianIndex(i, j)
            nodes[8] = CartesianIndex(i-1, j-1)
        else
            nodes[1] = CartesianIndex(i-1, j-1)
            nodes[2] = CartesianIndex(i+1, j-1)
            nodes[3] = CartesianIndex(i+1, j+1)
            nodes[4] = CartesianIndex(i-1, j+1)
            nodes[5] = CartesianIndex(i, j-1)
            nodes[6] = CartesianIndex(i+1, j)
            nodes[7] = CartesianIndex(i, j+1)
            nodes[8] = CartesianIndex(i-1, j)
        end
    end

    for k in eachindex(x_nodes)
        x_nodes[k] = coordinates_x[nodes[k]]
        y_nodes[k] = coordinates_y[nodes[k]]
    end
end

"""
    error_estimate(true_val, estimated_val)

Estimate the error between a true value and estimated value.

This operates on indiviual values, hence you should use broadcasting to apply it to matrices.
"""
error_estimate(true_val, estimated_val) = abs(estimated_val - true_val) / abs(true_val)

# DO NOT ALLOW THE FUNCTION MATRICES TO HAVE HALO CELLS
"""
    interpolate_to_grid!(grid1::NamedTuple{(:grid, :component_mesh), Tuple{T, ComponentMesh2D}}, grid2::NamedTuple{(:grid, :component_mesh), Tuple{S, ComponentMesh2D}}, grid1_func_vals::Matrix, grid2_func_vals::Matrix) where {T <: CurvilinearGrids.AbstractCurvilinearGrid2D, S <: CurvilinearGrids.AbstractCurvilinearGrid2D}

Interpolate a solution from `grid1_func_vals` using coordinates from `grid1` to `grid2` using an eight-point-quadratic isoparametric method.

In order to use this function, you must first create your meshes using `create_components`.

Note: This function uses centroid values (cell centers). Due to this, the two matrices `grid1_func_vals` and `grid2_func_vals` should have the same size as the domain of the centroid coordinate arrays for their respective grids. Observe this means the function matrices should NOT have halo cells.
"""
function interpolate_to_grid!(
    grid1::NamedTuple{(:grid, :component_mesh), Tuple{T, ComponentMesh2D}},
    grid2::NamedTuple{(:grid, :component_mesh), Tuple{S, ComponentMesh2D}},
    grid1_func_vals::Matrix,
    grid2_func_vals::Matrix
) where {T <: CurvilinearGrids.AbstractCurvilinearGrid2D, S <: CurvilinearGrids.AbstractCurvilinearGrid2D}

    # We assume to interplate from grid1 to grid2, and that grid2 has interpolation cells marked
    # Pre-allocate some matrices
    point = Vector{Float64}(undef, 2)
    F_val = Vector{Float64}(undef, 8)
    x_nodes = Vector{Float64}(undef, 8)
    y_nodes = Vector{Float64}(undef, 8)
    nodes = Vector{CartesianIndex{2}}(undef, 8)
    coords_x_grid1 = grid1.grid.centroid_coordinates.x[grid1.grid.iterators.cell.domain]
    coords_y_grid1 = grid1.grid.centroid_coordinates.y[grid1.grid.iterators.cell.domain]
    coords_x_grid2 = grid2.grid.centroid_coordinates.x[grid2.grid.iterators.cell.domain]
    coords_y_grid2 = grid2.grid.centroid_coordinates.y[grid2.grid.iterators.cell.domain]
    pts = Vector{SVector{2, Float64}}(undef, length(coords_x_grid1))
    @inbounds for idx in 1:size(pts, 1)
        pts[idx] = SVector(coords_x_grid1[idx], coords_y_grid1[idx])
    end
    kdtree = KDTree(pts)

    # Pre-allocations for the eight_node_quadratic
    ξη = zeros(Float64, 2)
    @inline function eval_N(i::Int, ξ, η)
        if i == 1 
            return 0.25 * (1 - ξ) * (1 - η) * (-ξ - η - 1)
        elseif i == 2 
            return 0.25 * (1 + ξ) * (1 - η) * (ξ - η - 1)
        elseif i == 3 
            return 0.25 * (1 + ξ) * (1 + η) * (ξ + η - 1)
        elseif i == 4 
            return 0.25 * (1 - ξ) * (1 + η) * (-ξ + η - 1)
        elseif i == 5 
            return 0.5 * (1 - ξ^2) * (1 - η)
        elseif i == 6 
            return 0.5 * (1 + ξ) * (1 - η^2)
        elseif i == 7 
            return 0.5 * (1 - ξ^2) * (1 + η)
        elseif i == 8 
            return 0.5 * (1 - ξ) * (1 - η^2)
        else
            return zero(Float64)
        end
    end
    @inline function eval_dN_dξ(i::Int, ξ, η)
        if i == 1 
            return 0.25 * (1 - η) * (2 * ξ + η)
        elseif i == 2 
            return 0.25 * (1 - η) * (2 * ξ - η)
        elseif i == 3 
            return 0.25 * (1 + η) * (2 * ξ + η)
        elseif i == 4 
            return 0.25 * (1 + η) * (2 * ξ - η)
        elseif i == 5 
            return ξ * (η - 1)
        elseif i == 6 
            return 0.5 * (1 - η^2)
        elseif i == 7 
            return -ξ * (η + 1)
        elseif i == 8 
            return 0.5 * (η^2 - 1)
        else
            return zero(Float64)
        end
    end
    @inline function eval_dN_dη(i::Int, ξ, η)
        if i == 1 
            return 0.25 * (1 - ξ) * (ξ + 2 * η)
        elseif i == 2 
            return -0.25 * (1 + ξ) * (ξ - 2 * η)
        elseif i == 3 
            return 0.25 * (1 + ξ) * (ξ + 2 * η)
        elseif i == 4 
            return -0.25 * (1 - ξ) * (ξ - 2 * η)
        elseif i == 5 
            return 0.5 * (ξ^2 - 1)
        elseif i == 6 
            return -(ξ + 1) * η
        elseif i == 7 
            return 0.5 * (1 - ξ^2)
        elseif i == 8 
            return (ξ - 1) * η
        else
            return zero(Float64)
        end
    end
    dNξ = [eval_dN_dξ(i, ξη[1], ξη[2]) for i in 1:8]
    dNη = [eval_dN_dη(i, ξη[1], ξη[2]) for i in 1:8]
    R = zeros(Float64, 2)
    J = zeros(Float64, 2,2)
    Δξη = zeros(Float64, 2)
    N_vals = zeros(Float64, 8)

    @inbounds for c in CartesianIndices(grid2.component_mesh.blank_mask)
        grid2.component_mesh.blank_mask[c] == -1 || continue

        point[1] = coords_x_grid2[c]
        point[2] = coords_y_grid2[c]

        extract_patch!(point, coords_x_grid1, coords_y_grid1, kdtree, x_nodes, y_nodes, nodes)

        @inbounds for i in eachindex(F_val)
            F_val[i] = grid1_func_vals[nodes[i]]
        end

        grid2_func_vals[c] = eight_node_quadratic(point, x_nodes, y_nodes, F_val, ξη, eval_N, eval_dN_dξ, eval_dN_dη, dNξ, dNη, R, J, Δξη, N_vals)
    end
end
# N = (
#     (ξ, η) -> 0.25 * (1 - ξ) * (1 - η) * (-ξ - η - 1),
#     (ξ, η) -> 0.25 * (1 + ξ) * (1 - η) * (ξ - η - 1),
#     (ξ, η) -> 0.25 * (1 + ξ) * (1 + η) * (ξ + η - 1),
#     (ξ, η) -> 0.25 * (1 - ξ) * (1 + η) * (-ξ + η - 1),
#     (ξ, η) -> 0.5 * (1 - ξ^2) * (1 - η),
#     (ξ, η) -> 0.5 * (1 + ξ) * (1 - η^2),
#     (ξ, η) -> 0.5 * (1 - ξ^2) * (1 + η),
#     (ξ, η) -> 0.5 * (1 - ξ) * (1 - η^2)
# )
# dN_dξ = (
#     (ξ, η) -> 0.25 * (1 - η) * (2 * ξ + η),
#     (ξ, η) -> 0.25 * (1 - η) * (2 * ξ - η),
#     (ξ, η) -> 0.25 * (1 + η) * (2 * ξ + η),
#     (ξ, η) -> 0.25 * (1 + η) * (2 * ξ - η),
#     (ξ, η) -> ξ * (η - 1),
#     (ξ, η) -> 0.5 * (1 - η^2),
#     (ξ, η) -> -ξ * (η + 1),
#     (ξ, η) -> 0.5 * (η^2 - 1)
# )
# dN_dη = (
#     (ξ, η) -> 0.25 * (1 - ξ) * (ξ + 2 * η),
#     (ξ, η) -> -0.25 * (1 + ξ) * (ξ - 2 * η),
#     (ξ, η) -> 0.25 * (1 + ξ) * (ξ + 2 * η),
#     (ξ, η) -> -0.25 * (1 - ξ) * (ξ - 2 * η),
#     (ξ, η) -> 0.5 * (ξ^2 - 1),
#     (ξ, η) -> -(ξ + 1) * η,
#     (ξ, η) -> 0.5 * (1 - ξ^2),
#     (ξ, η) -> (ξ - 1) * η
# )
