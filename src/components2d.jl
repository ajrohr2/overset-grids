"""
    create_components(grids::Tuple{Vararg{CurvilinearGrids.AbstractCurvilinearGrid2D}}; centroids=true, check_overlap=true, mark_interpolation=true, num_interp_points=5, mark_interior=true, resolution_type=:min)
    create_components(grids::Tuple{Vararg{CurvilinearGrids.AbstractCurvilinearGrid2D}}, names::Tuple{Vararg{String}}; centroids=true, check_overlap=true, mark_interpolation=true, num_interp_points=5, mark_interior=true, resolution_type=:min)
    create_components(grids::Tuple{Vararg{CurvilinearGrids.AbstractCurvilinearGrid2D}}, z_orders::Tuple{Vararg{Int}}; centroids=true, mark_interpolation=true, num_interp_points=5, mark_interior=true, check_overlap=true)

Create a dictionary of component meshes from a list of 2D `CurvilinearGrids` using either centroid coordinates or node coordinates.

It may be useful to name the component meshes you create. To this end, you may specify a tuple of names when calling this function.

Should you know the z-order of the collection of meshes you have, you may specify a tuples of z-orders when calling this function.

If you don't know the z-order of the collection of meshes you have, the function will attempt to determine it automatically using the resolution of each grid. You may specify the metric it uses to do this by setting `resolution_type` to `:min`, `:avg`, or `:max`. The resultant z-order will prioritize grids that have the smallest minimum, average, or maximum cell area.

You may specify if you wish to mark cells/points as interpolation with the `mark_interpolation` argument. You may also specify the number of cells (from the outside, in) you would like to be interpolation cells with the `num_interp_points` argument.

If your grids have interior points, you should ensure `mark_interior` is true, otherwise the interpolation marking will not work as expected. "Interior" points are defined as regions your current grids may cover that aren't actually a part of the problem you are solving (think the inside of a circular grid that is meant to simulate water flowing around a steel ball). Interior points are represented with a 2 in the `blank_mask` parameter of each mesh.

This function only returns metadata for the tuple of grids you pass in. It is imperative to save the tuple of grids you use to construct the `ComponentMesh`s.

Note: If you have two meshes with the same resolution, marking the `check_overlap` argument `true` will run overlap detection among all of the meshes on the same z-level. This can be computationally taxing, so if you are sure your equivalent-resolution meshes don't overlap, mark this `false`.
"""
function create_components(grids::Tuple{Vararg{CurvilinearGrids.AbstractCurvilinearGrid2D}}; centroids=true, check_overlap=true, mark_interpolation=true, num_interp_points=5, mark_interior=true, resolution_type=:min)
    meshes = Dict{Int, Vector{ComponentMesh2D}}() 
    zs = determine_z_order(resolution_type, grids)
    for (z, grid_list) in zs
        if centroids
            meshes[z] = ComponentMesh2D[ComponentMesh2D(zeros(Int8, (grids[grid_index].nnodes[1]-1, grids[grid_index].nnodes[2]-1)), grid_index, z, nothing, get_boundary(grids[grid_index])..., generate_kdtree(grids[grid_index])) for grid_index in grid_list]
        else
            meshes[z] = ComponentMesh2D[ComponentMesh2D(zeros(Int8, grids[grid_index].nnodes), grid_index, z, nothing, get_boundary(grids[grid_index])..., generate_kdtree(grids[grid_index])) for grid_index in grid_list]
        end
    end

    if check_overlap
        check_illegal_meshes(meshes, centroids, grids)
    end

    slicer!(meshes, grids, centroids, mark_interior)

    if mark_interpolation
        mark_interpolation_cells!(meshes, num_interp_points)
    end

    mark_background_interpolation!(meshes, grids, centroids, num_interp_points)

    return reassociate(meshes, grids)
end
function create_components(grids::Tuple{Vararg{CurvilinearGrids.AbstractCurvilinearGrid2D}}, names::Tuple{Vararg{String}}; centroids=true, check_overlap=true, mark_interpolation=true, num_interp_points=5, mark_interior=true, resolution_type=:min)
    meshes = Dict{Int, Vector{ComponentMesh2D}}() 
    zs = determine_z_order(resolution_type, grids)
    for (z, grid_list) in zs
        if centroids
            meshes[z] = ComponentMesh2D[ComponentMesh2D(zeros(Int8, (grids[grid_index].nnodes[1]-1, grids[grid_index].nnodes[2]-1)), grid_index, z, names[grid_index], get_boundary(grids[grid_index])..., generate_kdtree(grids[grid_index])) for grid_index in grid_list]
        else
            meshes[z] = ComponentMesh2D[ComponentMesh2D(zeros(Int8, grids[grid_index].nnodes), grid_index, z, names[grid_index], get_boundary(grids[grid_index])..., generate_kdtree(grids[grid_index])) for grid_index in grid_list]
        end
    end

    if check_overlap
        check_illegal_meshes(meshes, centroids, grids)
    end

    slicer!(meshes, grids, centroids, mark_interior)

    if mark_interpolation
        mark_interpolation_cells!(meshes, num_interp_points)
    end

    mark_background_interpolation!(meshes, grids, centroids, num_interp_points)

    return reassociate(meshes, grids)
end
function create_components(grids::Tuple{Vararg{CurvilinearGrids.AbstractCurvilinearGrid2D}}, z_orders::Tuple{Vararg{Int}}; centroids=true, mark_interpolation=true, num_interp_points=5, mark_interior=true, check_overlap=true)
    meshes = Dict{Int, Vector{ComponentMesh2D}}() 
    for grid_index in eachindex(grids)
        z = z_orders[grid_index]
        grid = grids[grid_index]
        if z in keys(meshes)
            if centroids
                push!(meshes[z], ComponentMesh2D(zeros(Int8, (grid.nnodes[1]-1, grid.nnodes[2]-1)), grid, z, nothing, get_boundary(grid)..., generate_kdtree(grid)))
            else
                push!(meshes[z], ComponentMesh2D(zeros(Int8, grid.nnodes), grid, z, nothing, get_boundary(grid)..., generate_kdtree(grid)))
            end
        else
            if centroids
                meshes[z] = ComponentMesh2D[ComponentMesh2D(zeros(Int8, (grid.nnodes[1]-1, grid.nnodes[2]-1)), grid_index, z, nothing, get_boundary(grid)..., generate_kdtree(grid))]
            else
                meshes[z] = ComponentMesh2D[ComponentMesh2D(zeros(Int8, grid.nnodes), grid_index, z, nothing, get_boundary(grid)..., generate_kdtree(grid))]
            end
        end
    end

    if check_overlap
        check_illegal_meshes(meshes, centroids, grids)
    end

    slicer!(meshes, grids, centroids, mark_interior)

    if mark_interpolation
        mark_interpolation_cells!(meshes, num_interp_points)
    end

    mark_background_interpolation!(meshes, grids, centroids, num_interp_points)

    return reassociate(meshes, grids)
end

# --- Helper functions --- #
function determine_z_order(resolution_type, grids::Tuple{Vararg{CurvilinearGrids.AbstractCurvilinearGrid2D}})
    order = Dict{Float64, Vector{Int}}() 
    for grid_index in eachindex(grids)
        if resolution_type == :min
            m = trunc(minimum(grids[grid_index].cell_center_metrics.J[grids[grid_index].iterators.cell.domain]), digits=12)
        elseif resolution_type == :avg 
            m = trunc(sum(grids[grid_index].cell_center_metrics.J[grids[grid_index].iterators.cell.domain]) / length(grids[grid_index].cell_center_metrics.J[grids[grid_index].iterators.cell.domain]), digits=12)
        elseif resolution_type == :max
            m = trunc(maximum(grids[grid_index].cell_center_metrics.J[grids[grid_index].iterators.cell.domain]), digits=12)
        else 
            error("Unknown resolution function $(resolution_type).")
        end
        if m in keys(order)
            push!(order[m], grid_index)
        else
            order[m] = Int[grid_index]
        end
    end

    zs = Dict{Int, Vector{Int}}()
    k = sort(collect(keys(order)))
    for idx in eachindex(k)
        zs[idx] = order[k[idx]]
    end

    return zs
end

# Check if two meshes of the same z-order overlap. If two exist, return an error.
function check_illegal_meshes(meshes::Dict{Int, Vector{ComponentMesh2D}}, centroids::Bool, grids::Tuple{Vararg{CurvilinearGrids.AbstractCurvilinearGrid2D}})
    rays = Vector{Line}(undef, 4)
    intersection_list = zeros(Int16, 4)

    for (z, z_level) in meshes
        for mesh_i in z_level
            for mesh_j in setdiff(z_level, [mesh_i])
                
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

                # Here is eventually where we want to implement spiral search
                for c_mj in CartesianIndices((1:size(x_mj)[1], 1:size(x_mj)[2]))
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
                            if mesh_i.name != ""
                                error("Two equivalent z-order meshes overlap! Check the meshes $(mesh_i.name), $(mesh_j.name) in z-order $(z).")
                            else
                                error("Two equivalent z-order meshes overlap!")
                            end
                        end

                        for i in eachindex(intersection_list)
                            intersection_list[i] = 0 
                        end
                    end
                end
            end 
        end 
    end 
end 

function generate_kdtree(grid::CurvilinearGrids.AbstractCurvilinearGrid2D)
    coords_x = grid.centroid_coordinates.x[grid.iterators.cell.domain]
    coords_y = grid.centroid_coordinates.y[grid.iterators.cell.domain]
    pts = Vector{SVector{2, Float64}}(undef, length(coords_x))
    @inbounds for idx in 1:size(pts, 1)
        pts[idx] = SVector(coords_x[idx], coords_y[idx])
    end
    return KDTree(pts)
end
