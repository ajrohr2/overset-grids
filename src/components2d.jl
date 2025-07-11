# Helper struct for intersection detections
struct Line
    point_0::NTuple{2, Float64}
    point_1::NTuple{2, Float64}
end

struct ComponentMesh2D{G}
    blank_mask::Array{Int8, 2}
    grid::G
    z_order::Int
    name::Union{String, Nothing}
    bounding_box::NTuple{4, Float64}
    boundary_polygon::Vector{Line}
end

"""
    create_components(grids::Vararg{CurvilinearGrids.AbstractCurvilinearGrid2D}; centroids=true, check_overlap=true, resolution_type=:min)
    create_components(grids::Vararg{Tuple{CurvilinearGrids.AbstractCurvilinearGrid2D, String}}; centroids=true, check_overlap=true, resolution_type=:min)
    create_components(grids::Vararg{Tuple{CurvilinearGrids.AbstractCurvilinearGrid2D, Int}}; centroids=true, check_overlap=true)

Create a dictionary of component meshes from a list of `CurvilinearGrids` using either centroid coordinates or node coordinates.

It may be useful to name the component meshes you create. To this end, you may specify tuples of `(grid, name)` when calling this function.

Should you know the z-order of the collection of meshes you have, you may specify tuples of `(grid, z-order)` when calling this function.

If you don't know the z-order of the collection of meshes you have, the function will attempt to determine it automatically using the resolution of each grid. You may specify the metric it uses to do this by setting `resolution_type` to `:min`, `:avg`, or `:max`. The resultant z-order will prioritize grids that have the smallest minimum, average, or maximum cell area.

Note: If you have two meshes with the same resolution, marking the `check_overlap` argument `true` will run overlap detection among all of the meshes on the same z-level. This can be computationally taxing, so if you are sure your equivalent-resolution meshes don't overlap, mark this `false`.
"""
function create_components(grids::Vararg{CurvilinearGrids.AbstractCurvilinearGrid2D}; centroids=true, check_overlap=true, resolution_type=:min)
    meshes = Dict{Int, Vector{<:ComponentMesh2D}}() 
    zs = determine_z_order(resolution_type, grids...)
    for (z, grid_list) in zs
        if centroids
            meshes[z] = [ComponentMesh2D{typeof(grid)}(zeros(Int8, (grid.nnodes[1]-1, grid.nnodes[2]-1)), grid, z, nothing, get_boundary!(grid)...) for grid in grid_list]
        else
            meshes[z] = [ComponentMesh2D{typeof(grid)}(zeros(Int8, grid.nnodes), grid, z, nothing, get_boundary!(grid)...) for grid in grid_list]
        end
    end

    if check_overlap
        check_illegal_meshes(meshes, centroids)
    end

    slicer!(meshes, centroids)

    return meshes
end
function create_components(grids::Vararg{Tuple{CurvilinearGrids.AbstractCurvilinearGrid2D, String}}; centroids=true, check_overlap=true, resolution_type=:min)
    meshes = Dict{Int, Vector{<:ComponentMesh2D}}() 
    zs = determine_z_order(resolution_type, grids...)
    for (z, grid_list) in zs
        if centroids
            meshes[z] = [ComponentMesh2D{typeof(grid)}(zeros(Int8, (grid[1].nnodes[1]-1, grid[1].nnodes[2]-1)), grid[1], z, grid[2], get_boundary!(grid[1])...) for grid in grid_list]
        else
            meshes[z] = [ComponentMesh2D{typeof(grid)}(zeros(Int8, grid[1].nnodes), grid[1], z, grid[2], get_boundary!(grid[1])...) for grid in grid_list]
        end
    end

    if check_overlap
        check_illegal_meshes(meshes, centroids)
    end

    slicer!(meshes, centroids)

    return meshes
end
function create_components(grids::Vararg{Tuple{CurvilinearGrids.AbstractCurvilinearGrid2D, Int}}; centroids=true, check_overlap=true)
    meshes = Dict{Int, Vector{<:ComponentMesh2D}}() 
    for (grid, z) in grids
        if z in keys(meshes)
            if centroids
                push!(meshes[z], ComponentMesh2D{typeof(grid)}(zeros(Int8, (grid.nnodes[1]-1, grid.nnodes[2]-1)), grid, z, nothing, get_boundary!(grid)...))
            else
                push!(meshes[z], ComponentMesh2D{typeof(grid)}(zeros(Int8, grid.nnodes), grid, z, nothing, get_boundary!(grid)...))
            end
        else
            if centroids
                meshes[z] = [ComponentMesh2D{typeof(grid)}(zeros(Int8, (grid.nnodes[1]-1, grid.nnodes[2]-1)), grid, z, nothing, get_boundary!(grid)...)]
            else
                meshes[z] = [ComponentMesh2D{typeof(grid)}(zeros(Int8, grid.nnodes), grid, z, nothing, get_boundary!(grid)...)]
            end
        end
    end

    if check_overlap
        check_illegal_meshes(meshes, centroids)
    end

    slicer!(meshes, centroids)

    return meshes
end

# --- Helper functions --- #
function determine_z_order(resolution_type, grids::Vararg{CurvilinearGrids.AbstractCurvilinearGrid2D})
    order = Dict() 
    for grid in grids
        if resolution_type == :min
            m = trunc(minimum(grid.cell_center_metrics.J[grid.iterators.cell.domain]), digits=12)
        elseif resolution_type == :avg 
            m = trunc(sum(grid.cell_center_metrics.J[grid.iterators.cell.domain]) / length(grid.cell_center_metrics.J[grid.iterators.cell.domain]), digits=12)
        elseif resolution_type == :max
            m = trunc(maximum(grid.cell_center_metrics.J[grid.iterators.cell.domain]), digits=12)
        else 
            error("Unknown resolution function $(resolution_type).")
        end
        if m in keys(order)
            push!(order[m], grid)
        else
            order[m] = [grid]
        end
    end

    zs = Dict()
    k = sort(collect(keys(order)))
    for idx in eachindex(k)
        zs[idx] = order[k[idx]]
    end

    return zs
end
function determine_z_order(resolution_type, grids::Vararg{Tuple{CurvilinearGrids.AbstractCurvilinearGrid2D, String}})
    order = Dict() 
    for grid in grids
        if resolution_type == :min
            m = minimum(grid[1].cell_center_metrics.J[grid[1].iterators.cell.domain])
        elseif resolution_type == :avg 
            m = sum(grid[1].cell_center_metrics.J[grid[1].iterators.cell.domain]) / length(grid[1].cell_center_metrics.J[grid[1].iterators.cell.domain])
        elseif resolution_type == :max
            m = maximum(grid[1].cell_center_metrics.J[grid[1].iterators.cell.domain])
        else 
            error("Unknown resolution function $(resolution_type).")
        end
        if m in keys(order)
            push!(order[m], grid)
        else
            order[m] = [grid]
        end
    end

    zs = Dict()
    k = sort(collect(keys(order)))
    for idx in eachindex(k)
        zs[idx] = order[k[idx]]
    end

    return zs
end

# Check if two meshes of the same z-order overlap. If two exist, return an error.
function check_illegal_meshes(meshes::Dict{Int, Vector{<:ComponentMesh2D}}, centroids)
    rays = Vector{Line}(undef, 4)
    intersection_list = zeros(Int16, 4)

    for (z, z_level) in meshes
        for mesh_i in z_level
            for mesh_j in setdiff(z_level, [mesh_i])
                
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
                    end
                end
            end 
        end 
    end 
end 
