struct Plane
    point_0::SVector{3, Float64}
    point_1::SVector{3, Float64}
    point_2::SVector{3, Float64}
    normal::SVector{3, Float64}
end
function Plane(point_0, point_1, point_2)
    normal = normalize(cross((point_1 .- point_0), (point_2 .- point_1)))
    return Plane(point_0, point_1, point_2, normal)
end

struct Ray
    begin_point::SVector{3, Float64}
    end_point::SVector{3, Float64}
    direction::SVector{3, Float64}
end
function Ray(begin_point, end_point)
    direction = normalize(end_point .- begin_point)
    return Ray(begin_point, end_point, direction)
end

struct ComponentMesh3D
    blank_mask::Array{Int8, 3}
    grid::CurvilinearGrids.AbstractCurvilinearGrid3D
    z_order::Int
    name::Union{String, Nothing}
    bounding_box::NTuple{6, Float64}
    boundary_polygon::Vector{Plane}
end

function create_components(grids::Vararg{CurvilinearGrids.AbstractCurvilinearGrid3D}; centroids=true, check_overlap=true, resolution_type=:min)
    meshes = Dict{Int, Vector{ComponentMesh3D}}()
    zs = determine_z_order(resolution_type, grids...)
    for (z, grid_list) in zs
        if centroids
            meshes[z] = [ComponentMesh3D(zeros(Int8, (grid.nnodes[1]-1, grid.nnodes[2]-1, grid.nnodes[3]-1)), grid, z, nothing, get_boundary!(grid)...) for grid in grid_list]
        else
            meshes[z] = [ComponentMesh3D(zeros(Int8, grid.nnodes[1]), grid, z, nothing, get_boundary!(grid)...) for grid in grid_list]
        end
    end

    if check_overlap
        #check_illegal_meshes(meshes, centroids)
    end

    slicer!(meshes, centroids)

    return meshes
end
function create_components(grids::Vararg{Tuple{CurvilinearGrids.AbstractCurvilinearGrid3D, String}}; centroids=true, check_overlap=true, resolution_type=:min)
    meshes = Dict{Int, Vector{ComponentMesh3D}}()
    zs = determine_z_order(resolution_type, grids...)
    for (z, grid_list) in zs
        if centroids
            meshes[z] = [ComponentMesh3D(zeros(Int8, (grid[1].nnodes[1]-1, grid[1].nnodes[2]-1, grid[1].nnodes[3]-1)), grid[1], z, grid[2], get_boundary!(grid[1])...) for grid in grid_list]
        else
            meshes[z] = [ComponentMesh3D(zeros(Int8, grid[1].nnodes), grid[1], z, grid[2], get_boundary!(grid[1])...) for grid in grid_list]
        end
    end

    if check_overlap
        #check_illegal_meshes(meshes, centroids)
    end

    slicer!(meshes, centroids)

    return meshes
end
function create_components(grids::Vararg{Tuple{CurvilinearGrids.AbstractCurvilinearGrid3D, Int}}; centroids=true, check_overlap=true)
    meshes = Dict{Int, Vector{ComponentMesh3D}}()
    for (grid, z) in grids
        if z in keys(meshes)
            if centroids
                push!(meshes[z], ComponentMesh3D(zeros(Int8, (grid.nnodes[1]-1, grid.nnodes[2]-1, grid.nnodes[3]-1)), grid, z, nothing, get_boundary!(grid)...))
            else
                push!(meshes[z], ComponentMesh3D(zeros(Int8, grid.nnodes), grid, z, nothing, get_boundary!(grid)...))
            end
        else
            if centroids
                meshes[z] = [ComponentMesh3D(zeros(Int8, (grid.nnodes[1]-1, grid.nnodes[2]-1, grid.nnodes[3]-1)), grid, z, nothing, get_boundary!(grid)...)]
            else
                meshes[z] = [ComponentMesh3D(zeros(Int8, grid.nnodes), grid, z, nothing, get_boundary!(grid)...)]
            end
        end
    end

    if check_overlap
        #check_illegal_meshes(meshes, centroids)
    end

    slicer!(meshes, centroids)

    return meshes
end

# --- Helper functions --- #
function determine_z_order(resolution_type, grids::Vararg{CurvilinearGrids.AbstractCurvilinearGrid3D})
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
function determine_z_order(resolution_type, grids::Vararg{Tuple{CurvilinearGrids.AbstractCurvilinearGrid3D, String}})
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
