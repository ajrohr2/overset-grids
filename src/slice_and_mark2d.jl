function slicer!(meshes::Dict{Int, Vector{ComponentMesh2D}}, grids::Tuple{Vararg{CurvilinearGrids.AbstractCurvilinearGrid2D}}, centroids::Bool, mark_interior::Bool)
    rays = Vector{Line}(undef, 4)
    intersection_list = zeros(Int16, 4)

    maxlen = 0
    for grid in grids
        l = length(CurvilinearGrids.coords(grid)[1])
        if l > maxlen
            maxlen = l
        end
    end
    overlap = Vector{CartesianIndex{2}}(undef, maxlen)
    if mark_interior
        interior = Vector{CartesianIndex{2}}(undef, maxlen)
    end

    @inbounds for i in keys(meshes)
        for mesh_i in meshes[i]
            @inbounds for j in (i+1):length(meshes)
                for mesh_j in meshes[j]
                    # Mesh i should be cutting mesh j, meaning we need to find the CartesianIndices where mesh j will change
                    grid_j = grids[mesh_j.grid_index]

                    boundary_polygon = mesh_i.boundary_polygon

                    if centroids
                        x_mj, y_mj = CurvilinearGrids.centroids(grid_j)
                    else
                        x_mj, y_mj = CurvilinearGrids.coords(grid_j)
                    end

                    # Need to find a bounding box for the combined grids
                    largest_x = max(mesh_i.bounding_box[2], mesh_j.bounding_box[2])
                    smallest_x = min(mesh_i.bounding_box[1], mesh_j.bounding_box[1])
                    
                    largest_y = max(mesh_i.bounding_box[4], mesh_j.bounding_box[4])
                    smallest_y = min(mesh_i.bounding_box[3], mesh_j.bounding_box[3])

                    # Pre-allocate arrays
                    overlap_num = 0
                    interior_num = 0

                    # Here is eventually where we want to implement spiral search
                    @inbounds for c_mj in CartesianIndices((1:size(x_mj)[1], 1:size(x_mj)[2]))
                        if x_mj[c_mj] < mesh_i.bounding_box[1] || x_mj[c_mj] > mesh_i.bounding_box[2] || y_mj[c_mj] < mesh_i.bounding_box[3] || y_mj[c_mj] > mesh_i.bounding_box[4]
                            continue
                        end
                        create_rays!(c_mj, rays, x_mj, y_mj, (smallest_x, smallest_y), (largest_x, largest_y))

                        @inbounds for l in eachindex(rays)
                            for segment in boundary_polygon
                                if determine_intersection(rays[l], segment)
                                    intersection_list[l] += 1
                                end
                            end

                            if intersection_list[l] % 2 == 1 
                                intersection_list[1] = intersection_list[2] = intersection_list[3] = intersection_list[4] = 0
                                overlap_num += 1
                                overlap[overlap_num] = c_mj
                                break
                            end
                        end

                        # Experimentally mark interior cells in the main slicing loop
                        if mark_interior && !any(==(0), intersection_list)
                            intersection_list[1] = intersection_list[2] = intersection_list[3] = intersection_list[4] = 0
                            interior_num += 1
                            interior[interior_num] = c_mj
                        end
                    end

                    # Update the blank matrix for overlap cells
                    @inbounds for c in 1:overlap_num
                        mesh_j.blank_mask[overlap[c]] = 1
                    end
                    if mark_interior
                        # Update the blank matrix for interior cells
                        @inbounds for c in 1:interior_num
                            mesh_j.blank_mask[interior[c]] = 2
                        end
                    end
                end
            end
        end
    end
end

# Function to determine interpolation cells

"""
    mark_interpolation_cells!(meshes::Dict{Int, Vector{ComponentMesh2D}}, num_interp_cells::Int)

Mark the first `num_interp_cells` border overlap cells as interpolation cells.

If your grids have interior cells, this function is unlikely to work as expected unless you first call `slice_interior!`. To learn more about interior cells, see the docstring for `slice_interior!`.

Interpolation cells are marked with a -1 in the `blank_mask` parameter of each mesh.
"""
function mark_interpolation_cells!(meshes::Dict{Int, Vector{ComponentMesh2D}}, num_interp_cells::Int)
    neighbors = ones(Int16, 4)
    for mesh_list in values(meshes)
        for mesh in mesh_list
            interp_points = Vector{CartesianIndex{2}}(undef, length(mesh.boundary_polygon)*num_interp_cells)
            CI = CartesianIndices(mesh.blank_mask)
            for point in CI
                if mesh.blank_mask[point] == 1 
                    if point[1] > 1
                        neighbors[1] = mesh.blank_mask[point[1]-1, point[2]]
                    end
                    if point[1] < size(mesh.blank_mask)[1]
                        neighbors[2] = mesh.blank_mask[point[1]+1, point[2]]
                    end
                    if point[2] > 1
                        neighbors[3] = mesh.blank_mask[point[1], point[2]-1]
                    end
                    if point[2] < size(mesh.blank_mask)[2]
                        neighbors[4] = mesh.blank_mask[point[1], point[2]+1]
                    end

                    num_surrounding_zeros = count(==(0), neighbors)
                    if num_surrounding_zeros >= 1 
                        mesh.blank_mask[point] = -1 
                    end

                    neighbors[1] = neighbors[2] = neighbors[3] = neighbors[4] = 1
                end
            end
            for _ in 1:num_interp_cells-1
                num_interp_points = 0
                for point in CI
                    if mesh.blank_mask[point] == 1 
                        if point[1] > 1
                            neighbors[1] = mesh.blank_mask[point[1]-1, point[2]]
                        end
                        if point[1] < size(mesh.blank_mask)[1]
                            neighbors[2] = mesh.blank_mask[point[1]+1, point[2]]
                        end
                        if point[2] > 1
                            neighbors[3] = mesh.blank_mask[point[1], point[2]-1]
                        end
                        if point[2] < size(mesh.blank_mask)[2]
                            neighbors[4] = mesh.blank_mask[point[1], point[2]+1]
                        end

                        num_surrounding_interps = count(==(-1), neighbors)
                        if num_surrounding_interps >= 1 
                            num_interp_points += 1
                            interp_points[num_interp_points] = point
                        end

                        neighbors[1] = neighbors[2] = neighbors[3] = neighbors[4] = 1
                    end
                end
                for c in 1:num_interp_points
                    mesh.blank_mask[interp_points[c]] = -1
                end
            end
        end
    end
end

function mark_background_interpolation!(meshes::Dict{Int, Vector{ComponentMesh2D}}, grids::Tuple{Vararg{CurvilinearGrids.AbstractCurvilinearGrid2D}}, centroids::Bool, num_interp_points::Int64)
    dict_keys = sort(collect(keys(meshes)))
    nodes = Vector{CartesianIndex{2}}(undef, 8)
    point = Vector{Float64}(undef, 2)
    idxs = Vector{Int}(undef, 1)
    dists = Vector{Float64}(undef, 1)
    @inbounds for i in dict_keys
        for marked_mesh in meshes[i]
            marked_grid = grids[marked_mesh.grid_index]
            if centroids
                x_mj, y_mj = CurvilinearGrids.centroids(marked_grid)
            else
                x_mj, y_mj = CurvilinearGrids.coords(marked_grid)
            end
            @inbounds for j in (i+1):length(meshes)
                for marking_mesh in meshes[j]
                    cis = CartesianIndices(marking_mesh.blank_mask)
                    for interp_point in CartesianIndices(marked_mesh.blank_mask)
                        if (interp_point[1] > num_interp_points && interp_point[2] > num_interp_points && interp_point[1] < size(marked_mesh.blank_mask)[1] - (num_interp_points - 1) && interp_point[2] < size(marked_mesh.blank_mask)[2] - (num_interp_points - 1)) || (marked_mesh.blank_mask[interp_point] == -1)
                            continue
                        end

                        if x_mj[interp_point] < marking_mesh.bounding_box[1] || x_mj[interp_point] > marking_mesh.bounding_box[2] || y_mj[interp_point] < marking_mesh.bounding_box[3] || y_mj[interp_point] > marking_mesh.bounding_box[4]
                            continue
                        end

                        point[1] = x_mj[interp_point]
                        point[2] = y_mj[interp_point]

                        extract_patch_marking!(point, cis, marking_mesh.kdtree, nodes, idxs, dists)

                        num_ones = 0
                        @inbounds for c in nodes
                            if marking_mesh.blank_mask[c] == 1 
                                num_ones += 1
                            end
                        end
                        if num_ones == 0
                            marked_mesh.blank_mask[interp_point] = -1
                        end
                    end
                end
            end
        end
    end
end
function extract_patch_marking!(point, cis, kdtree, nodes, idxs, dists)
    # First find the closest point in Euclidean space. This uses a KDTree approach
    knn!(idxs, dists, kdtree, point, 1)
    i, j = cis[idxs[1]].I # This is the "center point" of the patch we want

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
        elseif j == size(cis)[2]
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
    elseif i == size(cis)[1]
        if j == 1
            nodes[1] = CartesianIndex(i-2, j)
            nodes[2] = CartesianIndex(i, j)
            nodes[3] = CartesianIndex(i, j+2)
            nodes[4] = CartesianIndex(i-2, j+2)
            nodes[5] = CartesianIndex(i-1, j)
            nodes[6] = CartesianIndex(i, j+1)
            nodes[7] = CartesianIndex(i-1, j+2)
            nodes[8] = CartesianIndex(i, j+1)
        elseif j == size(cis)[2]
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
        elseif j == size(cis)[2]
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
end
