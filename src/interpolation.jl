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

function ts_node_triquadratic(point::Vector{Float64}, x_nodes::Array{Float64}, y_nodes::Array{Float64}, z_nodes::Array{Float64}, function_vals::Array{Float64}, ξηζ::Vector{Float64}, N, dN_dξ, dN_dη, dN_dζ, dNξ::Array{Float64}, dNη::Array{Float64}, dNζ::Array{Float64}, R::Vector{Float64}, J::Matrix{Float64}, Δξηζ::Vector{Float64}, N_vals::Array{Float64})
    # Initialize the solved point
    ξηζ[1] = ξηζ[2] = ξηζ[3] = 0

    @inbounds for i in CartesianIndices(x_nodes)
        dNξ[i] = dN_dξ(i, ξηζ[1], ξηζ[2], ξηζ[3])
        dNη[i] = dN_dη(i, ξηζ[1], ξηζ[2], ξηζ[3])
        dNζ[i] = dN_dζ(i, ξηζ[1], ξηζ[2], ξηζ[3])
    end

    x = zero(Float64)
    y = zero(Float64)
    z = zero(Float64)
    @inbounds for i in CartesianIndices(x_nodes)
        x += N(i, ξηζ[1], ξηζ[2], ξηζ[3]) * x_nodes[i]
        y += N(i, ξηζ[1], ξηζ[2], ξηζ[3]) * y_nodes[i]
        z += N(i, ξηζ[1], ξηζ[2], ξηζ[3]) * z_nodes[i]
    end

    R[1] = point[1] - x
    R[2] = point[2] - y
    R[3] = point[3] - z

    while sqrt(R[1]^2 + R[2]^2 + R[3]^2) ≥ 1e-12
        @inbounds for i in 1:9
            J[i] = 0.0
        end
        @inbounds for i in CartesianIndices(x_nodes)
            J[1,1] += dNξ[i] * x_nodes[i]
            J[1,2] += dNη[i] * x_nodes[i]
            J[1,3] += dNζ[i] * x_nodes[i]
            J[2,1] += dNξ[i] * y_nodes[i]
            J[2,2] += dNη[i] * y_nodes[i]
            J[2,3] += dNζ[i] * y_nodes[i]
            J[3,1] += dNξ[i] * z_nodes[i]
            J[3,2] += dNη[i] * z_nodes[i]
            J[3,3] += dNζ[i] * z_nodes[i]
        end

        # Access J variables for cleanliness
        a, b, c = J[1,1], J[1,2], J[1,3]
        d, e, f = J[2,1], J[2,2], J[2,3]
        g, h, i = J[3,1], J[3,2], J[3,3]
        #
        detJ = a * (e*i - f*h) - b * (d*i - g*f) + c * (d*h - g*e)
        # Inverse J, but written out
        inv00 = (e*i - f*h); inv01 = -(b*i - c*h); inv02 = (b*f - c*e)
        inv10 = -(d*i - f*g); inv11 = (a*i - c*g); inv12 = -(a*f - c*d)
        inv20 = (d*h - e*g); inv21 = -(a*h - b*g); inv22 = (a*e - b*d)
        # 
        Δξηζ[1] = (inv00*R[1] + inv01*R[2] + inv02*R[3]) / detJ
        Δξηζ[2] = (inv10*R[1] + inv11*R[2] + inv12*R[3]) / detJ
        Δξηζ[3] = (inv20*R[1] + inv21*R[2] + inv22*R[3]) / detJ
        
        ξηζ[1] += Δξηζ[1]
        ξηζ[2] += Δξηζ[2]
        ξηζ[3] += Δξηζ[3]
        ξ, η, ζ = ξηζ[1], ξηζ[2], ξηζ[3]

        @inbounds for i in CartesianIndices(x_nodes)
            N_vals[i] = N(i, ξ, η, ζ)
            dNξ[i] = dN_dξ(i, ξ, η, ζ)
            dNη[i] = dN_dη(i, ξ, η, ζ)
            dNζ[i] = dN_dζ(i, ξ, η, ζ)
        end

        x = y = z = 0.0
        @inbounds for i in CartesianIndices(x_nodes)
            x += N_vals[i] * x_nodes[i]
            y += N_vals[i] * y_nodes[i]
            z += N_vals[i] * z_nodes[i]
        end

        R[1] = point[1] - x
        R[2] = point[2] - y
        R[3] = point[3] - z
    end

    f_interp = 0.0
    @inbounds for i in CartesianIndices(x_nodes)
        f_interp += N(i, ξηζ[1], ξηζ[2], ξηζ[3]) * function_vals[i]
    end
    return f_interp
end

function extract_patch2d!(point, coordinates_x, coordinates_y, kdtree, x_nodes, y_nodes, nodes)
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
function extract_patch3d!(point, offsets, coordinates_x, coordinates_y, coordinates_z, kdtree, x_nodes, y_nodes, z_nodes, nodes)
    # First find the closest point in Euclidean space. This uses a KDTree approach
    idxs, dists = knn(kdtree, point, 1)
    i, j, k = CartesianIndices(coordinates_x)[idxs[1]].I # This is the "center point" of the patch we want
    i0 = min(i, size(coordinates_x)[1] - 2)
    j0 = min(j, size(coordinates_x)[2] - 2)
    k0 = min(k, size(coordinates_x)[3] - 2)

    @inbounds for m in 1:27
        di, dj, dk = offsets[m]
        nodes[m] = CartesianIndex(i0+di, j0+dj, k0+dk)
        x_nodes[m] = coordinates_x[nodes[m]]
        y_nodes[m] = coordinates_y[nodes[m]]
        z_nodes[m] = coordinates_z[nodes[m]]
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
    interpolate_to_grid2d!(grid1::NamedTuple{(:grid, :component_mesh), Tuple{T, ComponentMesh2D}}, grid2::NamedTuple{(:grid, :component_mesh), Tuple{S, ComponentMesh2D}}, grid1_func_vals::Matrix, grid2_func_vals::Matrix) where {T <: CurvilinearGrids.AbstractCurvilinearGrid2D, S <: CurvilinearGrids.AbstractCurvilinearGrid2D}

Interpolate a solution from `grid1_func_vals` using coordinates from `grid1` to `grid2` using an eight-point-quadratic isoparametric method.

In order to use this function, you must first create your meshes using `create_components`.

Note: This function uses centroid values (cell centers). Due to this, the two matrices `grid1_func_vals` and `grid2_func_vals` should have the same size as the domain of the centroid coordinate arrays for their respective grids. Observe this means the function matrices should NOT have halo cells.
"""
function interpolate_to_grid2d!(
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
    kdtree = grid1.component_mesh.kdtree

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

        extract_patch2d!(point, coords_x_grid1, coords_y_grid1, kdtree, x_nodes, y_nodes, nodes)

        @inbounds for i in eachindex(F_val)
            F_val[i] = grid1_func_vals[nodes[i]]
        end

        grid2_func_vals[c] = eight_node_quadratic(point, x_nodes, y_nodes, F_val, ξη, eval_N, eval_dN_dξ, eval_dN_dη, dNξ, dNη, R, J, Δξη, N_vals)
    end
end
"""
    interpolate_to_grid3d!(grid1::NamedTuple{(:grid, :component_mesh), Tuple{T, ComponentMesh2D}}, grid2::NamedTuple{(:grid, :component_mesh), Tuple{S, ComponentMesh2D}}, grid1_func_vals::Matrix, grid2_func_vals::Matrix) where {T <: CurvilinearGrids.AbstractCurvilinearGrid2D, S <: CurvilinearGrids.AbstractCurvilinearGrid2D}

Interpolate a solution from `grid1_func_vals` using coordinates from `grid1` to `grid2` using an eight-point-quadratic isoparametric method.

In order to use this function, you must first create your meshes using `create_components`.

Note: This function uses centroid values (cell centers). Due to this, the two matrices `grid1_func_vals` and `grid2_func_vals` should have the same size as the domain of the centroid coordinate arrays for their respective grids. Observe this means the function matrices should NOT have halo cells.
"""
function interpolate_to_grid3d!(
    grid1::NamedTuple{(:grid, :component_mesh), Tuple{T, ComponentMesh3D}},
    grid2::NamedTuple{(:grid, :component_mesh), Tuple{S, ComponentMesh3D}},
    grid1_func_vals::Array,
    grid2_func_vals::Array
) where {T <: CurvilinearGrids.AbstractCurvilinearGrid3D, S <: CurvilinearGrids.AbstractCurvilinearGrid3D}

    # We assume to interplate from grid1 to grid2, and that grid2 has interpolation cells marked
    # Pre-allocate some matrices
    point = Vector{Float64}(undef, 3)
    F_val = Array{Float64}(undef, 3, 3, 3)
    x_nodes = Array{Float64}(undef, 3, 3, 3)
    y_nodes = Array{Float64}(undef, 3, 3, 3)
    z_nodes = Array{Float64}(undef, 3, 3, 3)
    nodes = Array{CartesianIndex{3}}(undef, 3, 3, 3)
    CI = CartesianIndices(nodes)
    @views begin
        coords_x_grid1 = grid1.grid.centroid_coordinates.x[grid1.grid.iterators.cell.domain]
        coords_y_grid1 = grid1.grid.centroid_coordinates.y[grid1.grid.iterators.cell.domain]
        coords_z_grid1 = grid1.grid.centroid_coordinates.z[grid1.grid.iterators.cell.domain]
        coords_x_grid2 = grid2.grid.centroid_coordinates.x[grid2.grid.iterators.cell.domain]
        coords_y_grid2 = grid2.grid.centroid_coordinates.y[grid2.grid.iterators.cell.domain]
        coords_z_grid2 = grid2.grid.centroid_coordinates.z[grid2.grid.iterators.cell.domain]
    end
    kdtree = grid1.component_mesh.kdtree

    # Pre-allocations for the 27-node-triquadratic
    ξηζ = zeros(Float64, 3)
    # 1D shape functions
    @inline function x_shape(i::Int, x)
        if i == 1
            # minus 1
            return (x * (x - 1)) / 2 
        elseif i == 2 
            # 0
            return 1 - x^2
        else
            # plus 1
            return (x * (x + 1)) / 2
        end
    end
    # 1D shape functions derivatives
    @inline function dx_shape(i::Int, x)
        if i == 1
            # minus 1
            return ((2 * x) - 1) / 2
        elseif i == 2 
            # 0
            return 2 * x
        else
            # plus 1
            return ((2 * x) + 1) / 2
        end
    end
    @inline function eval_N(c::CartesianIndex{3}, ξ, η, ζ)
        return x_shape(c[1], ξ) * x_shape(c[2], η) * x_shape(c[3], ζ)
    end
    @inline function eval_dN_dξ(c::CartesianIndex{3}, ξ, η, ζ)
        return dx_shape(c[1], ξ) * x_shape(c[2], η) * x_shape(c[3], ζ)
    end
    @inline function eval_dN_dη(c::CartesianIndex{3}, ξ, η, ζ)
        return x_shape(c[1], ξ) * dx_shape(c[2], η) * x_shape(c[3], ζ)
    end
    @inline function eval_dN_dζ(c::CartesianIndex{3}, ξ, η, ζ)
        return x_shape(c[1], ξ) * x_shape(c[2], η) * dx_shape(c[3], ζ)
    end
    dNξ = zeros(Float64, 3, 3, 3)
    dNη = zeros(Float64, 3, 3, 3)
    dNζ = zeros(Float64, 3, 3, 3)
    @inbounds for i in CI
        dNξ[i] = eval_dN_dξ(i, ξηζ[1], ξηζ[2], ξηζ[3])
        dNη[i] = eval_dN_dη(i, ξηζ[1], ξηζ[2], ξηζ[3])
        dNζ[i] = eval_dN_dζ(i, ξηζ[1], ξηζ[2], ξηζ[3])
    end
    R = zeros(Float64, 3)
    J = zeros(Float64, 3, 3)
    Δξηζ = zeros(Float64, 3)
    N_vals = zeros(Float64, 3, 3, 3)
    offsets = [(c[1]-1, c[2]-1, c[3]-1) for c in CI]

    @inbounds for c in CartesianIndices(grid2.component_mesh.blank_mask)
        grid2.component_mesh.blank_mask[c] == -1 || continue

        point[1] = coords_x_grid2[c]
        point[2] = coords_y_grid2[c]
        point[3] = coords_z_grid2[c]

        extract_patch3d!(point, offsets, coords_x_grid1, coords_y_grid1, coords_z_grid1, kdtree, x_nodes, y_nodes, z_nodes, nodes)

        @inbounds for i in eachindex(F_val)
            F_val[i] = grid1_func_vals[nodes[i]]
        end

        grid2_func_vals[c] = ts_node_triquadratic(point, x_nodes, y_nodes, z_nodes, F_val, ξηζ, eval_N, eval_dN_dξ, eval_dN_dη, eval_dN_dζ, dNξ, dNη, dNζ, R, J, Δξηζ, N_vals)
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
