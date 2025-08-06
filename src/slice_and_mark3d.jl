function slicer!(meshes::Dict{Int, Vector{ComponentMesh3D}}, grids::Tuple{Vararg{CurvilinearGrids.AbstractCurvilinearGrid3D}}, centroids::Bool, mark_interior::Bool)
    # First index is beginning point, second is direction
    bvtn_num = Vector{Int64}(undef, 2)
    pert = SVector(1e-4 * ℯ / 100, 1e-4 * ℯ / 100, 1e-4 * ℯ / 100)
    ray_directions = [
        SVector(1.0, 0.0, 0.0) + pert,
        SVector(0.0, 1.0, 0.0) + pert,
        SVector(0.0, 0.0, 1.0) + pert,
        SVector(-1.0, 0.0, 0.0) + pert,
        SVector(0.0, -1.0, 0.0) + pert,
        SVector(0.0, 0.0, -1.0) + pert,
    ]
    intersection_list = zeros(Int16, 6)

    # Maybe we can define this with the largest size possible.
    maxlen = 0
    grid_ind = 0
    for grid in eachindex(grids)
        l = length(CurvilinearGrids.coords(grids[grid])[1])
        if l > maxlen
            maxlen = l
            grid_ind = grid
        end
    end
    maxbvtt = 0
    for mesh_arr in values(meshes)
        for mesh in mesh_arr
            if length(mesh.boundary_polygon) > maxbvtt
                maxbvtt = length(mesh.boundary_polygon)
            end
        end
    end
    bvtt1 = similar(meshes[1][1].bvh.nodes, SVector{2, Int}, 2*maxbvtt)
    bvtt2 = similar(meshes[1][1].bvh.nodes, SVector{2, Int}, 2*maxbvtt)
    overlap = Vector{CartesianIndex{3}}(undef, maxlen)
    if mark_interior
        interior = Vector{CartesianIndex{3}}(undef, maxlen)
    end

    @inbounds for i in keys(meshes)
        for mesh_i in meshes[i]
            @inbounds for j in (i+1):length(meshes)
                for mesh_j in meshes[j]
                    # Mesh i should be cutting mesh j, meaning we need to find the CartesianIndices where mesh j will change
                    boundary_polygon = mesh_i.boundary_polygon

                    x_mj = mesh_j.x
                    y_mj = mesh_j.y
                    z_mj = mesh_j.z

                    overlap_num = 0
                    interior_num = 0

                    @inbounds for c_mj in CartesianIndices(size(x_mj))
                        if x_mj[c_mj] < mesh_i.bounding_box[1] || x_mj[c_mj] > mesh_i.bounding_box[2] || y_mj[c_mj] < mesh_i.bounding_box[3] || y_mj[c_mj] > mesh_i.bounding_box[4] || z_mj[c_mj] < mesh_i.bounding_box[5] || z_mj[c_mj] > mesh_i.bounding_box[6]
                            continue
                        end
                        r = SVector{3}(x_mj[c_mj], y_mj[c_mj], z_mj[c_mj])

                        traverse_rays!(mesh_i.bvh, r, ray_directions, bvtt1, bvtt2, bvtn_num)

                        contacts = bvtn_num[1] == 1 ? bvtt1 : bvtt2

                        @inbounds for idx in 1:bvtn_num[2]
                            contact = contacts[idx]
                            d = ray_directions[contact[2]]
                            if determine_intersection(boundary_polygon[contact[1]], r, d)
                                intersection_list[contact[2]] += 1 
                            end
                        end
                        for intersection in intersection_list
                            if isodd(intersection)
                                fill!(intersection_list, 0)
                                overlap_num += 1
                                overlap[overlap_num] = c_mj
                            end
                        end

                        # Experimentally mark interior cells in the main slicing loop
                        if mark_interior
                            all_nonzero = true
                            for intersection in intersection_list
                                if intersection == 0 
                                    all_nonzero = false
                                    break
                                end
                            end
                            if all_nonzero
                                fill!(intersection_list, 0)
                                interior_num += 1
                                interior[interior_num] = c_mj
                            end
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

# Function to determine interpolation cells | All allocations are from the array "interp_points"

"""
    mark_interpolation_cells!(meshes::Dict{Int, Vector{ComponentMesh3D}}, num_interp_cells::Int)

Mark the first `num_interp_cells` border overlap cells as interpolation cells.

If your grids have interior cells, this function is unlikely to work as expected unless you first call `slice_interior!`. To learn more about interior cells, see the docstring for `slice_interior!`.

Interpolation cells are marked with a -1 in the `blank_mask` parameter of each mesh.
"""
function mark_interpolation_cells!(meshes::Dict{Int, Vector{ComponentMesh3D}}, num_interp_cells::Int)
    neighbors = ones(Int16, 6)
    maxlen = 0
    for mesh_list in values(meshes)
        for mesh in mesh_list
            l = length(mesh.boundary_polygon)
            if l > maxlen
                maxlen = l
            end
        end
    end
    interp_points = Vector{CartesianIndex{3}}(undef, maxlen*num_interp_cells)
    for mesh_list in values(meshes)
        for mesh in mesh_list
            CI = CartesianIndices(mesh.blank_mask)
            for point in CI
                @inbounds if mesh.blank_mask[point] == 1 
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

                    for n in neighbors
                        if n == 0 
                            mesh.blank_mask[point] = -1 
                            break
                        end
                    end

                    fill!(neighbors, 1)
                end
            end
            for _ in 1:num_interp_cells-1
                num_interp_points = 0
                for point in CI
                    @inbounds if mesh.blank_mask[point] == 1 
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

                        for n in neighbors
                            if n == -1 
                                num_interp_points += 1
                                interp_points[num_interp_points] = point
                                break
                            end
                        end

                        fill!(neighbors, 1)
                    end
                end
                @inbounds for c in 1:num_interp_points
                    mesh.blank_mask[interp_points[c]] = -1
                end
            end
        end
    end
end

function mark_background_interpolation!(meshes::Dict{Int, Vector{ComponentMesh3D}}, num_interp_points::Int64)
    dict_keys = sort(collect(keys(meshes)))
    nodes = Array{CartesianIndex{3}}(undef, 3, 3, 3)
    offsets = [(c[1]-1, c[2]-1, c[3]-1) for c in CartesianIndices(nodes)]
    point = Vector{Float64}(undef, 3)
    idxs = Vector{Int}(undef, 1)
    dists = Vector{Float64}(undef, 1)
    @inbounds for i in dict_keys
        for marked_mesh in meshes[i]
            x_mj = marked_mesh.x 
            y_mj = marked_mesh.y 
            z_mj = marked_mesh.z
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
