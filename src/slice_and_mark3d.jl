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

                    @inbounds for c_mj in CartesianIndices(size(x_mj))
                        if x_mj[c_mj] < mesh_i.bounding_box[1] || x_mj[c_mj] > mesh_i.bounding_box[2] || y_mj[c_mj] < mesh_i.bounding_box[3] || y_mj[c_mj] > mesh_i.bounding_box[4] || z_mj[c_mj] < mesh_i.bounding_box[5] || z_mj[c_mj] > mesh_i.bounding_box[6]
                            continue
                        end
                        create_rays!(c_mj, rays, x_mj, y_mj, z_mj, (smallest_x, smallest_y, smallest_z), (largest_x, largest_y, largest_z))

                        @inbounds for l in eachindex(rays)
                            for face in boundary_polygon
                                if determine_intersection(face, rays[l])
                                    intersection_list[l] += 1
                                end
                            end

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

function mark_background_interpolation!(meshes::Dict{Int, Vector{ComponentMesh3D}}, grids::Tuple{Vararg{CurvilinearGrids.AbstractCurvilinearGrid3D}}, centroids::Bool, num_interp_points::Int64)
    dict_keys = sort(collect(keys(meshes)))
    nodes = Array{CartesianIndex{3}}(undef, 3, 3, 3)
    offsets = [(c[1]-1, c[2]-1, c[3]-1) for c in CartesianIndices(nodes)]
    point = Vector{Float64}(undef, 3)
    idxs = Vector{Int}(undef, 1)
    dists = Vector{Float64}(undef, 1)
    @inbounds for i in dict_keys
        for marked_mesh in meshes[i]
            marked_grid = grids[marked_mesh.grid_index]
            if centroids
                x_mj, y_mj, z_mj = CurvilinearGrids.centroids(marked_grid)
            else
                x_mj, y_mj, z_mj = CurvilinearGrids.coords(marked_grid)
            end
            @inbounds for j in (i+1):length(meshes)
                for marking_mesh in meshes[j]
                    cis = CartesianIndices(marking_mesh.blank_mask)
                    for interp_point in CartesianIndices(marked_mesh.blank_mask)
                        if (interp_point[1] > num_interp_points && interp_point[2] > num_interp_points && interp_point[3] > num_interp_points && interp_point[1] < size(marked_mesh.blank_mask)[1] - (num_interp_points - 1) && interp_point[2] < size(marked_mesh.blank_mask)[2] - (num_interp_points - 1) && interp_point[3] < size(marked_mesh.blank_mask)[3] - (num_interp_points - 1)) || (marked_mesh.blank_mask[interp_point] == -1)
                            continue
                        end

                        if x_mj[interp_point] < marking_mesh.bounding_box[1] || x_mj[interp_point] > marking_mesh.bounding_box[2] || y_mj[interp_point] < marking_mesh.bounding_box[3] || y_mj[interp_point] > marking_mesh.bounding_box[4] || z_mj[interp_point] < marking_mesh.bounding_box[5] || z_mj[interp_point] > marking_mesh.bounding_box[6]
                            continue
                        end

                        point[1] = x_mj[interp_point]
                        point[2] = y_mj[interp_point]
                        point[3] = z_mj[interp_point]

                        extract_patch_marking!(point, offsets, cis, marking_mesh.kdtree, nodes, idxs, dists)

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
function extract_patch_marking!(point, offsets, cis, kdtree, nodes, idxs, dists)
    # First find the closest point in Euclidean space. This uses a KDTree approach
    knn!(idxs, dists, kdtree, point, 1)
    i, j, k = cis[idxs[1]].I # This is the "center point" of the patch we want
    i0 = min(i, size(cis)[1] - 2)
    j0 = min(j, size(cis)[2] - 2)
    k0 = min(k, size(cis)[3] - 2)

    @inbounds for m in 1:27
        di, dj, dk = offsets[m]
        nodes[m] = CartesianIndex(i0+di, j0+dj, k0+dk)
    end
end
# function mark_background_interpolation!(meshes::Dict{Int, Vector{ComponentMesh3D}}, grids::Tuple{Vararg{CurvilinearGrids.AbstractCurvilinearGrid3D}}, centroids::Bool; num_interp_points=5)
#     rays = Vector{Ray}(undef, 6)
#     intersection_list = zeros(Int16, 6)
# 
#     dict_keys = sort(collect(keys(meshes)))
# 
#     maxlen = 0
#     for grid in grids
#         l = length(CurvilinearGrids.coords(grid)[1])
#         if l > maxlen
#             maxlen = l
#         end
#     end
#     overlap = Vector{CartesianIndex{3}}(undef, maxlen * 2)
# 
#     @inbounds for i in dict_keys
#         for marked_mesh in meshes[i]
#             marked_grid = grids[marked_mesh.grid_index]
#             if centroids
#                 x_mj, y_mj, z_mj = CurvilinearGrids.centroids(marked_grid)
#             else
#                 x_mj, y_mj, z_mj = CurvilinearGrids.coords(marked_grid)
#             end
#             @inbounds for j in (i+1):length(meshes)
#                 for marking_mesh in meshes[j]
#                     boundary_polygon = marking_mesh.boundary_polygon
# 
#                     # Need to find a bounding box for the combined grids
#                     largest_x = max(marked_mesh.bounding_box[2], marking_mesh.bounding_box[2])
#                     smallest_x = min(marked_mesh.bounding_box[1], marking_mesh.bounding_box[1])
#                     
#                     largest_y = max(marked_mesh.bounding_box[4], marking_mesh.bounding_box[4])
#                     smallest_y = min(marked_mesh.bounding_box[3], marking_mesh.bounding_box[3])
# 
#                     largest_z = max(marked_mesh.bounding_box[6], marking_mesh.bounding_box[6])
#                     smallest_z = min(marked_mesh.bounding_box[5], marking_mesh.bounding_box[5])
# 
#                     overlap_num = 0
#                     for interp_point in CartesianIndices(marked_mesh.blank_mask)
#                         if (interp_point[1] > num_interp_points && interp_point[2] > num_interp_points && interp_point[3] > num_interp_points && interp_point[1] < size(marked_mesh.blank_mask)[1] - (num_interp_points - 1) && interp_point[2] < size(marked_mesh.blank_mask)[2] - (num_interp_points - 1) && interp_point[3] < size(marked_mesh.blank_mask)[3] - (num_interp_points - 1)) || (marked_mesh.blank_mask[interp_point] == -1)
#                             continue
#                         end
# 
#                         if x_mj[interp_point] < marking_mesh.bounding_box[1] || x_mj[interp_point] > marking_mesh.bounding_box[2] || y_mj[interp_point] < marking_mesh.bounding_box[3] || y_mj[interp_point] > marking_mesh.bounding_box[4] || z_mj[interp_point] < marking_mesh.bounding_box[5] || z_mj[interp_point] > marking_mesh.bounding_box[6]
#                             continue
#                         end
#                         create_rays!(interp_point, rays, x_mj, y_mj, z_mj, (smallest_x, smallest_y, smallest_z), (largest_x, largest_y, largest_z))
# 
#                         @inbounds for l in eachindex(rays)
#                             for face in boundary_polygon
#                                 if determine_intersection(face, rays[l])
#                                     intersection_list[l] += 1
#                                     # println("Ray $(rays[l]) intersected line $segment ")
#                                 end
#                             end
# 
#                             if intersection_list[l] % 2 == 1 
#                                 intersection_list[1] = intersection_list[2] = intersection_list[3] = intersection_list[4] = intersection_list[5] = intersection_list[6] = 0
#                                 overlap_num += 1
#                                 overlap[overlap_num] = interp_point
#                                 break
#                             end
#                         end
#                     end
#                     for k in 1:overlap_num
#                         marked_mesh.blank_mask[overlap[k]] = -1
#                     end
#                 end
#             end
#         end
#     end
# end
