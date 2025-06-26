using CurvilinearGrids

struct ComponentMesh{D}
    blank_mask::Array{Int8, D}
    grid::CurvilinearGrids.AbstractCurvilinearGrid
    z_order::Int
end

function create_componenets(grids::Vararg{CurvilinearGrids.AbstractCurvilinearGrid})
    meshes = Dict{Int, ComponentMesh}() 
    zs = determine_z_order(grids...)
    for (z, grid) in zs
        meshes[z] = ComponentMesh{length(grid.nnodes)}(zeros(Int8, grid.nnodes), grid, z)
    end

    slicer(meshes)

    return meshes
end

# --- Helper functions --- #
function determine_z_order(grids::Vararg{CurvilinearGrids.AbstractCurvilinearGrid})
    order = Dict() 
    for grid in grids
        order[minimum(grid.cell_center_metrics.J[grid.iterators.cell.domain])] = grid
    end
    
    if length(grids) != length(order)
        error("Two grids likely have the same resolution. Unsure how to proceed!")
    end

    zs = Dict()
    k = sort(collect(keys(order)))
    for idx in eachindex(k)
        zs[idx] = order[k[idx]]
    end

    return zs
end

# --- Orientation and ray intersection based overlap detection --- #
struct Line
    point_0::Tuple
    point_1::Tuple
end

function orientation(A::Tuple, B::Tuple, C::Tuple)
    o = (B[1] - A[1])*(C[2] - A[2])-(B[2]-A[2])*(C[1]-A[1])
    if o < 0 
        return -1
    elseif o > 0 
        return 1 
    else
        return 0 
    end
end

function determine_intersection(l1::Line, l2::Line)
    oABC = orientation(l1.point_0, l1.point_1, l2.point_0)
    oABD = orientation(l1.point_0, l1.point_1, l2.point_1)
    oCDA = orientation(l2.point_0, l2.point_1, l1.point_0)
    oCDB = orientation(l2.point_0, l2.point_1, l1.point_1)

    # Maybe if all the points are colinear, then we shouldn't count it as an intersection (rays that run along a boundary don't cross the boundary)

    if oABC == oABD == oCDA == oCDB == 0 
        # Ray runs along boundary
        return false
    end

    # This means some points are colinear
    if oABC == 0 
        if l1.point_0[1] ≤ l2.point_0[1] ≤ l1.point_1[1] && l1.point_0[2] ≤ l2.point_0[2] ≤ l1.point_1[2]
            return true
        end
    elseif oABD == 0 
        if l1.point_0[1] ≤ l2.point_1[1] ≤ l1.point_1[1] && l1.point_0[2] ≤ l2.point_1[2] ≤ l1.point_1[2]
            return true
        end
    elseif oCDA == 0 
        if l2.point_0[1] ≤ l1.point_0[1] ≤ l2.point_1[1] && l2.point_0[2] ≤ l1.point_0[2] ≤ l2.point_1[2]
            return true
        end
    elseif oCDB == 0 
        if l2.point_0[1] ≤ l1.point_1[1] ≤ l2.point_1[1] && l2.point_0[2] ≤ l1.point_1[2] ≤ l2.point_1[2]
            return true
        end
    end

    if oABC != oABD && oCDA != oCDB
        return true
    else
        return false
    end
end

function create_rays(point::CartesianIndex, x_array, y_array, bottom_corner, top_corner)
    # We need the rays to exist in the largest global space
    rays = []
    push!(rays, Line((bottom_corner[1], y_array[1, point[2]]), (x_array[point], y_array[point])))
    push!(rays, Line((x_array[point], y_array[point]), (top_corner[1], y_array[end, point[2]])))
    push!(rays, Line((x_array[point[1], 1], bottom_corner[2]), (x_array[point], y_array[point])))
    push!(rays, Line((x_array[point], y_array[point]), (x_array[point[1], end], top_corner[2])))

    return rays
end

function create_boundary_polygon(cartesian_indices, x_array, y_array)
    boundary1 = cartesian_indices[1,:]
    boundary2 = cartesian_indices[:,1]
    boundary3 = cartesian_indices[end,:]
    boundary4 = cartesian_indices[:,end]

    polygon1 = []
    polygon2 = []
    polygon3 = []
    polygon4 = []

    for i in 1:length(boundary1)-1
        push!(polygon1, Line((x_array[boundary1[i]], y_array[boundary1[i]]), (x_array[boundary1[i+1]], y_array[boundary1[i+1]])))
    end
    for i in 1:length(boundary2)-1
        push!(polygon2, Line((x_array[boundary2[i]], y_array[boundary2[i]]), (x_array[boundary2[i+1]], y_array[boundary2[i+1]])))
    end
    for i in 1:length(boundary3)-1
        push!(polygon3, Line((x_array[boundary3[i]], y_array[boundary3[i]]), (x_array[boundary3[i+1]], y_array[boundary3[i+1]])))
    end
    for i in 1:length(boundary4)-1
        push!(polygon4, Line((x_array[boundary4[i]], y_array[boundary4[i]]), (x_array[boundary4[i+1]], y_array[boundary4[i+1]])))
    end
    
    return [polygon1, polygon2, polygon3, polygon4]
end

function slicer(meshes::Dict{Int, ComponentMesh}; overlap_cell_num=5)
    for i in eachindex(collect(keys(meshes)))
        for j in i+1:length(meshes)
            # Mesh i should be cutting mesh j, meaning we need to find the CartesianIndices where mesh j will change
            overlap = []

            grid_i = meshes[i].grid
            grid_j = meshes[j].grid
            x_mi, y_mi = CurvilinearGrids.coords(grid_i)
            mi_cartinds = CartesianIndices((1:size(x_mi)[1], 1:size(x_mi)[2]))
            boundary_polygon = create_boundary_polygon(mi_cartinds, x_mi, y_mi)
            x_mj, y_mj = CurvilinearGrids.coords(grid_j)

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
                    intersection_list .%= 2
                end

                if !all(v -> v == 0, intersection_list)
                    push!(overlap, c_mj)
                end
            end

            # Update the iblank matrix
            for c in overlap
                meshes[j].blank_mask[c] = 1
            end
        end
    end
end

function save_vtk_with_threshold(mesh, blank, fn)
    @info "Writing to $fn.vti"

    xyz_n = CurvilinearGrids.coords(mesh)
    domain = mesh.iterators.cell.domain

    @views vtk_grid(fn, xyz_n) do vtk
        vtk["J", VTKCellData()] = mesh.cell_center_metrics.J[domain]

        # vtk["volume", VTKCellData()] = CurvilinearGrids.cellvolume.(Ref(mesh), domain)

        vtk["xi", VTKCellData(), component_names=["x1", "x2", "t"]] = (
            mesh.cell_center_metrics.ξ.x₁[domain],
            mesh.cell_center_metrics.ξ.x₂[domain],
            mesh.cell_center_metrics.ξ.t[domain],
        )

        vtk["eta", VTKCellData(), component_names=["x1", "x2", "t"]] = (
            mesh.cell_center_metrics.η.x₁[domain],
            mesh.cell_center_metrics.η.x₂[domain],
            mesh.cell_center_metrics.η.t[domain],
        )

        vtk["dx_di", VTKCellData(), component_names=["xi", "eta"]] = (
            mesh.cell_center_metrics.x₁.ξ[domain], mesh.cell_center_metrics.x₁.η[domain]
        )

        vtk["dy_di", VTKCellData(), component_names=["xi", "eta"]] = (
            mesh.cell_center_metrics.x₂.ξ[domain], mesh.cell_center_metrics.x₂.η[domain]
        )
              
        vtk["blank", VTKPointData()] = blank
    end
end
