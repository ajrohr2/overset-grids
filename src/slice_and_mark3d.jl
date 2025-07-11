function slicer!(meshes::Dict{Int, Vector{ComponentMesh3D}}, centroids)
    rays = Vector{Ray}(undef, 6)
    intersection_list = zeros(Int16, 6)

    @inbounds for i in keys(meshes)
        for mesh_i in meshes[i]
            @inbounds for j in (i+1):length(meshes)
                for mesh_j in meshes[j]
                    # Mesh i should be cutting mesh j, meaning we need to find the CartesianIndices where mesh j will change
                    grid_j = mesh_j.grid 

                    boundary_polygon = mesh_i.boundary_polygon

                    if centroids
                        x_mj, y_mj, z_mj = CurvilinearGrids.centroids(grid_j)
                    else
                        x_mj, y_mj, z_mj = CurvilinearGrids.coords(grid_j)
                    end

                    # Need to find a bounding box for the combined grids
                    largest_x = max(mesh_i.bounding_box[2], mesh_j.bounding_box[2])
                    smallest_x = min(mesh_i.bounding_box[1], mesh_j.bounding_box[1])
                    
                    largest_y = max(mesh_i.bounding_box[4], mesh_j.bounding_box[4])
                    smallest_y = min(mesh_i.bounding_box[3], mesh_j.bounding_box[3])

                    largest_z = max(mesh_i.bounding_box[6], mesh_j.bounding_box[6])
                    smallest_z = min(mesh_i.bounding_box[5], mesh_j.bounding_box[5])

                    # Pre-allocate arrays
                    overlap = Vector{CartesianIndex}(undef, length(x_mj))
                    overlap_num = 0

                    # Here is eventually where we want to implement spiral search
                    @inbounds for c_mj in CartesianIndices(size(x_mj))
                        if x_mj[c_mj] < mesh_i.bounding_box[1] || x_mj[c_mj] > mesh_i.bounding_box[2] || y_mj[c_mj] < mesh_i.bounding_box[3] || y_mj[c_mj] > mesh_i.bounding_box[4] || z_mj[c_mj] < mesh_i.bounding_box[5] || z_mj[c_mj] > mesh_i.bounding_box[6]
                            continue
                        end
                        create_rays!(c_mj, rays, x_mj, y_mj, z_mj, (smallest_x, smallest_y, smallest_z), (largest_x, largest_y, largest_z))

                        @inbounds for l in eachindex(rays)
                            for face in boundary_polygon
                                if determine_intersection(face, rays[l])
                                    intersection_list[l] += 1
                                    # println("Ray $(rays[l]) intersected line $segment ")
                                end
                            end

                            if intersection_list[l] % 2 == 1 
                                intersection_list[1] = intersection_list[2] = intersection_list[3] = intersection_list[4] = intersection_list[5] = intersection_list[6] = 0
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

# Function to determine interpolation cells

function mark_interpolation_cells!(meshes::Dict{Int, Vector{ComponentMesh3D}}, num_interp_cells::Int)
    neighbors = ones(Int16, 6)
    for mesh_list in values(meshes)
        for mesh in mesh_list
            interp_points = Vector{CartesianIndex}(undef, length(mesh.boundary_polygon)*num_interp_cells)
            for point in CartesianIndices(mesh.blank_mask)
                if mesh.blank_mask[point] == 1 
                    try
                        neighbors[1] = mesh.blank_mask[point[1]-1, point[2], point[3]]
                    catch e
                        # pass
                    end
                    try
                        neighbors[2] = mesh.blank_mask[point[1]+1, point[2], point[3]]
                    catch e
                        # pass
                    end
                    try
                        neighbors[3] = mesh.blank_mask[point[1], point[2]-1, point[3]]
                    catch e
                        # pass
                    end
                    try
                        neighbors[4] = mesh.blank_mask[point[1], point[2]+1, point[3]]
                    catch e
                        # pass
                    end
                    try
                        neighbors[5] = mesh.blank_mask[point[1], point[2], point[3]-1]
                    catch e
                        # pass
                    end
                    try
                        neighbors[6] = mesh.blank_mask[point[1], point[2], point[3]+1]
                    catch e
                        # pass
                    end

                    num_surrounding_zeros = count(x -> x == 0, neighbors)
                    if num_surrounding_zeros >= 1 
                        mesh.blank_mask[point] = -1 
                    end

                    neighbors[1] = neighbors[2] = neighbors[3] = neighbors[4] = neighbors[5] = neighbors[6] = 1
                end
            end
            for _ in 1:num_interp_cells-1
                num_interp_points = 0
                for point in CartesianIndices(mesh.blank_mask)
                    if mesh.blank_mask[point] == 1 
                        try
                            neighbors[1] = mesh.blank_mask[point[1]-1, point[2], point[3]]
                        catch e
                            # pass
                        end
                        try
                            neighbors[2] = mesh.blank_mask[point[1]+1, point[2], point[3]]
                        catch e
                            # pass
                        end
                        try
                            neighbors[3] = mesh.blank_mask[point[1], point[2]-1, point[3]]
                        catch e
                            # pass
                        end
                        try
                            neighbors[4] = mesh.blank_mask[point[1], point[2]+1, point[3]]
                        catch e
                            # pass
                        end
                        try
                            neighbors[5] = mesh.blank_mask[point[1], point[2], point[3]-1]
                        catch e
                            # pass
                        end
                        try
                            neighbors[6] = mesh.blank_mask[point[1], point[2], point[3]+1]
                        catch e
                            # pass
                        end

                        num_surrounding_interps = count(x -> x == -1, neighbors)
                        if num_surrounding_interps >= 1 
                            num_interp_points += 1
                            interp_points[num_interp_points] = point
                        end

                        neighbors[1] = neighbors[2] = neighbors[3] = neighbors[4] = neighbors[5] = neighbors[6] = 1
                    end
                end
                for c in 1:num_interp_points
                    mesh.blank_mask[interp_points[c]] = -1
                end
            end
        end
    end
end
