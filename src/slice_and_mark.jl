function slicer!(meshes::Dict{Int, Vector{ComponentMesh}}, centroids)
    for i in eachindex(collect(keys(meshes)))
        for mesh_i in meshes[i]
            for mesh_j in vcat([meshes[j] for j in i+1:length(meshes)]...)
                # Mesh i should be cutting mesh j, meaning we need to find the CartesianIndices where mesh j will change
                overlap = []

                grid_i = mesh_i.grid 
                grid_j = mesh_j.grid 

                x_mi, y_mi = CurvilinearGrids.coords(grid_i)
                mi_cartinds = CartesianIndices((1:size(x_mi)[1], 1:size(x_mi)[2]))
                boundary_polygon = create_boundary_polygon(mi_cartinds, x_mi, y_mi)

                if centroids
                    x_mj, y_mj = CurvilinearGrids.centroids(grid_j)
                else
                    x_mj, y_mj = CurvilinearGrids.coords(grid_j)
                end

                # Need to find a bounding box for the combined grids
                largest_x = max(maximum(x_mi), maximum(x_mj))
                smallest_x = min(minimum(x_mi), minimum(x_mj))
                
                largest_y = max(maximum(y_mi), maximum(y_mj))
                smallest_y = min(minimum(y_mi), minimum(y_mj))

                # Here is eventually where we want to implement spiral search
                for c_mj in CartesianIndices((1:size(x_mj)[1], 1:size(x_mj)[2]))
                    rays = create_rays(c_mj, x_mj, y_mj, (smallest_x, smallest_y), (largest_x, largest_y))
                    intersection_list = []

                    for l in eachindex(rays)
                        push!(intersection_list, 0)

                        for k in eachindex(boundary_polygon)
                            for segment in boundary_polygon[k]
                                if determine_intersection(rays[l], segment)
                                    intersection_list[l] += 1
                                    # println("Ray $(rays[l]) intersected line $segment ")
                                end
                            end
                        end

                        if intersection_list[l] % 2 == 1 
                            push!(overlap, c_mj)
                            break
                        end
                    end

                    # if !all(v -> v % 2 == 0, intersection_list)
                    #     push!(overlap, c_mj)
                    # end
                end

                # Update the iblank matrix
                for c in overlap
                    mesh_j.blank_mask[c] = 1
                end
            end
        end
    end
end

# Function to remove interiors of sliced meshes. The idea is if you wrap meshes around an object, the background mesh doesn't need to simulate the area the object takes up. So we mark it for deletion.

function slice_interior!(meshes::Dict{Int, Vector{ComponentMesh}})
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

function mark_interpolation_cells!(meshes::Dict{Int, Vector{ComponentMesh}}, num_interp_cells::Int)
    for mesh in vcat(values(meshes)...)
        for point in CartesianIndices(mesh.blank_mask)[1:end, 1:end]
            if mesh.blank_mask[point] == 1 
                neighbors = []

                try
                    push!(neighbors, mesh.blank_mask[point[1]-1, point[2]])
                catch e
                    # pass
                end
                try
                    push!(neighbors, mesh.blank_mask[point[1]+1, point[2]])
                catch e
                    # pass
                end
                try
                    push!(neighbors, mesh.blank_mask[point[1], point[2]-1])
                catch e
                    # pass
                end
                try
                    push!(neighbors, mesh.blank_mask[point[1], point[2]+1])
                catch e
                    # pass
                end

                num_surrounding_zeros = count(x -> x == 0, neighbors)
                if num_surrounding_zeros >= 1 
                    mesh.blank_mask[point] = -1 
                end
            end
        end
        for _ in 1:num_interp_cells-1
            interp_points = []
            for point in CartesianIndices(mesh.blank_mask)[1:end, 1:end]
                if mesh.blank_mask[point] == 1 
                    neighbors = []
                    try
                        push!(neighbors, mesh.blank_mask[point[1]-1, point[2]])
                    catch e
                        # pass
                    end
                    try
                        push!(neighbors, mesh.blank_mask[point[1]+1, point[2]])
                    catch e
                        # pass
                    end
                    try
                        push!(neighbors, mesh.blank_mask[point[1], point[2]-1])
                    catch e
                        # pass
                    end
                    try
                        push!(neighbors, mesh.blank_mask[point[1], point[2]+1])
                    catch e
                        # pass
                    end

                    num_surrounding_interps = count(x -> x == -1, neighbors)
                    if num_surrounding_interps >= 1 
                        push!(interp_points, point)
                    end
                end
            end
            for c in interp_points
                mesh.blank_mask[c] = -1
            end
        end
    end
end
