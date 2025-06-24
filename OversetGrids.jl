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

function slicer(meshes::Dict{Int, ComponentMesh}; overlap_cell_num=5)
    for i in eachindex(collect(keys(meshes)))
        for j in i+1:length(meshes)
            # Mesh i should be cutting mesh j, meaning we need to find the CartesianIndices where mesh j will change
            overlap = []

            grid_i = meshes[i].grid
            grid_j = meshes[j].grid
            x_mi, y_mi = CurvilinearGrids.coords(grid_i)
            mi_cartinds = CartesianIndices((1:size(x_mi)[1], 1:size(x_mi)[2]))
            mi_boundary = vcat(mi_cartinds[1,:], mi_cartinds[:,1], mi_cartinds[size(x_mi)[1],:], mi_cartinds[:,size(x_mi)[2]])
            x_mj, y_mj = CurvilinearGrids.coords(grid_j)

            # Here is eventually where we want to implement spiral search
            for c_mj in CartesianIndices((1:size(x_mj)[1], 1:size(x_mj)[2]))
                # Here's the logic: The current point is marked for deletion if and only if it's x and y value are between exactly 1 x and y boundary point in the cutting grid
                current_j_point = (x_mj[c_mj], y_mj[c_mj])

                xis_lower = [i for i in mi_boundary if x_mi[i] ≤ current_j_point[1]]
                xis_upper = [i for i in mi_boundary if x_mi[i] > current_j_point[1]]
                yis_lower = [i for i in mi_boundary if y_mi[i] ≤ current_j_point[2]]
                yis_upper = [i for i in mi_boundary if y_mi[i] > current_j_point[2]]

                # Logic here is if the reverse of the CartesianIndex also exists, then the line projected from the selected point crosses two borders, and therefore isn't inside the mesh
                quad3 = _remove_opposite_pairs(intersect(xis_lower, yis_lower), size(x_mi)[1], size(x_mi)[1])
                quad1 = _remove_opposite_pairs(intersect(xis_upper, yis_upper), size(x_mi)[1], size(x_mi)[1])
                quad2 = _remove_opposite_pairs(intersect(xis_lower, yis_upper), size(x_mi)[1], size(x_mi)[1])
                quad4 = _remove_opposite_pairs(intersect(xis_upper, yis_lower), size(x_mi)[1], size(x_mi)[1])

                if quad1 && quad2 && quad3 && quad4
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

function _remove_opposite_pairs(idx_list::Vector, nrows::Int, ncols::Int)
    to_remove = Set{CartesianIndex}()

    for idx in idx_list
        i, j = Tuple(idx)
        # Skip if already marked for removal
        if idx in to_remove
            continue
        end

        # Determine the opposite index
        if i == 1
            opposite = CartesianIndex(nrows, j)
        elseif i == nrows
            opposite = CartesianIndex(1, j)
        elseif j == 1
            opposite = CartesianIndex(i, ncols)
        elseif j == ncols
            opposite = CartesianIndex(i, 1)
        else
            continue  # not on the boundary
        end

        # If the opposite is in the list, mark both for removal
        if opposite in idx_list
            push!(to_remove, idx)
            push!(to_remove, opposite)
        end
    end

    # Remove marked indices from the original list (in-place)
    filter!(x -> x ∉ to_remove, idx_list)

    return length(idx_list) == 0 ? false : true
end
