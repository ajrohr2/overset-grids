function slicer!(meshes::Dict{Int, Vector{ComponentMesh3D}}, grids::Tuple{Vararg{CurvilinearGrids.AbstractCurvilinearGrid3D}}, centroids::Bool, mark_interior::Bool)
    rays = Vector{Ray}(undef, 6)
    intersection_list = zeros(Int16, 6)

    # Maybe we can define this with the largest size possible.
    maxlen = 0
    for grid in grids
        l = length(CurvilinearGrids.coords(grid)[1])
        if l > maxlen
            maxlen = l
        end
    end
    overlap = Vector{CartesianIndex{3}}(undef, maxlen)
    # interior = Vector{CartesianIndex{3}}(undef, maxlen)
    if mark_interior
        interior = Vector{CartesianIndex{3}}(undef, maxlen)
    end

    @inbounds for i in keys(meshes)
        for mesh_i in meshes[i]
            @inbounds for j in (i+1):length(meshes)
                for mesh_j in meshes[j]
                    # Mesh i should be cutting mesh j, meaning we need to find the CartesianIndices where mesh j will change
                    grid_j = grids[mesh_j.grid_index]

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

                    overlap_num = 0
                    interior_num = 0

                    # Here is eventually where we want to implement spiral search
                    @inbounds for c_mj in CartesianIndices(size(x_mj))
                        if x_mj[c_mj] < mesh_i.bounding_box[1] || x_mj[c_mj] > mesh_i.bounding_box[2] || y_mj[c_mj] < mesh_i.bounding_box[3] || y_mj[c_mj] > mesh_i.bounding_box[4] || z_mj[c_mj] < mesh_i.bounding_box[5] || z_mj[c_mj] > mesh_i.bounding_box[6]
                            continue
                        end
                        create_rays!(c_mj, rays, x_mj, y_mj, z_mj, (smallest_x, smallest_y, smallest_z), (largest_x, largest_y, largest_z))

                        @inbounds for l in eachindex(rays)
                            for face in boundary_polygon
                                if determine_intersection(face, rays[l])
                                    # if c_mj == CartesianIndex(17, 30, 20)
                                    #     println("Face: $(face)")
                                    # end
                                    intersection_list[l] += 1
                                end
                            end

                            # if c_mj == CartesianIndex(17, 30, 20)
                            #     println("Ray: $(rays[l])")
                            #     println(intersection_list)
                            #     println("---")
                            # end

                            if intersection_list[l] % 2 == 1 
                                intersection_list[1] = intersection_list[2] = intersection_list[3] = intersection_list[4] = intersection_list[5] = intersection_list[6] = 0
                                overlap_num += 1
                                overlap[overlap_num] = c_mj
                                break
                            end
                        end

                        # Experimentally mark interior cells in the main slicing loop
                        if mark_interior && !any(==(0), intersection_list)
                            intersection_list[1] = intersection_list[2] = intersection_list[3] = intersection_list[4] = intersection_list[5] = intersection_list[6] = 0
                            interior_num += 1
                            interior[interior_num] = c_mj
                        end
                    end

                    # Update the blank matrix for overlap cells
                    @inbounds for c in 1:overlap_num
                        mesh_j.blank_mask[overlap[c]] = 1
                    end
                    # Update the blank matrix for interior cells
                    # @inbounds for c in 1:interior_num
                    #     mesh_j.blank_mask[interior[c]] = 2
                    # end
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
    mark_interpolation_cells!(meshes::Dict{Int, Vector{ComponentMesh3D}}, num_interp_cells::Int)

Mark the first `num_interp_cells` border overlap cells as interpolation cells.

If your grids have interior cells, this function is unlikely to work as expected unless you first call `slice_interior!`. To learn more about interior cells, see the docstring for `slice_interior!`.

Interpolation cells are marked with a -1 in the `blank_mask` parameter of each mesh.
"""
function mark_interpolation_cells!(meshes::Dict{Int, Vector{ComponentMesh3D}}, num_interp_cells::Int)
    neighbors = ones(Int16, 6)
    for mesh_list in values(meshes)
        for mesh in mesh_list
            CI = CartesianIndices(mesh.blank_mask)
            interp_points = Vector{CartesianIndex{3}}(undef, length(mesh.boundary_polygon)*num_interp_cells)
            for point in CI
                if mesh.blank_mask[point] == 1 
                    if point[1] > 1
                        neighbors[1] = mesh.blank_mask[point[1]-1, point[2], point[3]]
                    end
                    if point[1] < size(mesh.blank_mask)[1]
                        neighbors[2] = mesh.blank_mask[point[1]+1, point[2], point[3]]
                    end
                    if point[2] > 1
                        neighbors[3] = mesh.blank_mask[point[1], point[2]-1, point[3]]
                    end
                    if point[2] < size(mesh.blank_mask)[2]
                        neighbors[4] = mesh.blank_mask[point[1], point[2]+1, point[3]]
                    end
                    if point[3] > 1
                        neighbors[5] = mesh.blank_mask[point[1], point[2], point[3]-1]
                    end
                    if point[3] < size(mesh.blank_mask)[3]
                        neighbors[6] = mesh.blank_mask[point[1], point[2], point[3]+1]
                    end

                    num_surrounding_zeros = count(==(0), neighbors)
                    if num_surrounding_zeros >= 1 
                        mesh.blank_mask[point] = -1 
                    end

                    neighbors[1] = neighbors[2] = neighbors[3] = neighbors[4] = neighbors[5] = neighbors[6] = 1
                end
            end
            for _ in 1:num_interp_cells-1
                num_interp_points = 0
                for point in CI
                    if mesh.blank_mask[point] == 1 
                        if point[1] > 1
                            neighbors[1] = mesh.blank_mask[point[1]-1, point[2], point[3]]
                        end
                        if point[1] < size(mesh.blank_mask)[1]
                            neighbors[2] = mesh.blank_mask[point[1]+1, point[2], point[3]]
                        end
                        if point[2] > 1
                            neighbors[3] = mesh.blank_mask[point[1], point[2]-1, point[3]]
                        end
                        if point[2] < size(mesh.blank_mask)[2]
                            neighbors[4] = mesh.blank_mask[point[1], point[2]+1, point[3]]
                        end
                        if point[3] > 1
                            neighbors[5] = mesh.blank_mask[point[1], point[2], point[3]-1]
                        end
                        if point[3] < size(mesh.blank_mask)[3]
                            neighbors[6] = mesh.blank_mask[point[1], point[2], point[3]+1]
                        end

                        num_surrounding_interps = count(==(-1), neighbors)
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
