function slicer!(meshes::Dict{Int, Vector{<:ComponentMesh2D}}, centroids)
    rays = Vector{Line}(undef, 4)
    intersection_list = zeros(Int16, 4)

    @inbounds for i in keys(meshes)
        for mesh_i in meshes[i]
            @inbounds for j in (i+1):length(meshes)
                for mesh_j in meshes[j]
                    # Mesh i should be cutting mesh j, meaning we need to find the CartesianIndices where mesh j will change
                    grid_j = mesh_j.grid 

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
                    overlap = Vector{CartesianIndex}(undef, length(x_mj))
                    overlap_num = 0

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
                                    # println("Ray $(rays[l]) intersected line $segment ")
                                end
                            end

                            if intersection_list[l] % 2 == 1 
                                intersection_list[1] = intersection_list[2] = intersection_list[3] = intersection_list[4] = 0
                                overlap_num += 1
                                overlap[overlap_num] = c_mj
                                break
                            end
                        end
                    end

                    # Update the iblank matrix
                    @inbounds for c in 1:overlap_num
                        mesh_j.blank_mask[overlap[c]] = 1
                    end
                end
            end
        end
    end
end

# Function to remove interiors of sliced meshes. The idea is if you wrap meshes around an object, the background mesh doesn't need to simulate the area the object takes up. So we mark it for deletion.

function slice_interior!(meshes::Dict{Int, Vector{<:ComponentMesh2D}})
    # Update to use rays rather than array slices
    for mesh in vcat(values(meshes)...)
        for point in CartesianIndices(mesh.blank_mask)
            if mesh.blank_mask[point] == 0 
                left = mesh.blank_mask[1:point[1], point[2]]
                right = mesh.blank_mask[point[1]:end, point[2]]
                top = mesh.blank_mask[point[1], point[2]:end]
                bottom = mesh.blank_mask[point[1], 1:point[2]]

                num_containing_one = count(x -> any(x .== 1), [left, right, top, bottom])
                if num_containing_one >= 3 
                    mesh.blank_mask[point] = 2 
                end
            end
        end
    end
end

# Function to determine interpolation cells

function mark_interpolation_cells!(meshes::Dict{Int, Vector{<:ComponentMesh2D}}, num_interp_cells::Int)
    neighbors = ones(Int16, 4)
    for mesh_list in values(meshes)
        for mesh in mesh_list
            interp_points = Vector{CartesianIndex}(undef, length(mesh.boundary_polygon)*num_interp_cells)
            for point in CartesianIndices(mesh.blank_mask)
                if mesh.blank_mask[point] == 1 
                    try
                        neighbors[1] = mesh.blank_mask[point[1]-1, point[2]]
                    catch e
                        # pass
                    end
                    try
                        neighbors[2] = mesh.blank_mask[point[1]+1, point[2]]
                    catch e
                        # pass
                    end
                    try
                        neighbors[3] = mesh.blank_mask[point[1], point[2]-1]
                    catch e
                        # pass
                    end
                    try
                        neighbors[4] = mesh.blank_mask[point[1], point[2]+1]
                    catch e
                        # pass
                    end

                    num_surrounding_zeros = count(x -> x == 0, neighbors)
                    if num_surrounding_zeros >= 1 
                        mesh.blank_mask[point] = -1 
                    end

                    neighbors[1] = neighbors[2] = neighbors[3] = neighbors[4] = 1
                end
            end
            for _ in 1:num_interp_cells-1
                num_interp_points = 0
                for point in CartesianIndices(mesh.blank_mask)
                    if mesh.blank_mask[point] == 1 
                        try
                            neighbors[1] = mesh.blank_mask[point[1]-1, point[2]]
                        catch e
                            # pass
                        end
                        try
                            neighbors[2] = mesh.blank_mask[point[1]+1, point[2]]
                        catch e
                            # pass
                        end
                        try
                            neighbors[3] = mesh.blank_mask[point[1], point[2]-1]
                        catch e
                            # pass
                        end
                        try
                            neighbors[4] = mesh.blank_mask[point[1], point[2]+1]
                        catch e
                            # pass
                        end

                        num_surrounding_interps = count(x -> x == -1, neighbors)
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
