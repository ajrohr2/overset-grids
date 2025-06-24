# function cut_holes(mesh1, mesh2; overlap_cell_num=5)
#     # For now, we assume that mesh1 is completely contained in mesh2
# 
#     cut_toward_center = true
# 
#     x_m1, y_m1 = CurvilinearGrids.coords(mesh1)
#     x_m2, y_m2 = CurvilinearGrids.coords(mesh2)
# 
#     domain_m1, domain_m2 = mesh1.iterators.cell.domain, mesh2.iterators.cell.domain
# 
#     overlap_points_m2 = []
# 
#     centroid_m1_x = mesh1.centroid_coordinates.x[domain_m1]
#     centroid_m1_y = mesh1.centroid_coordinates.y[domain_m1]
#     boundary_of_object = mesh1.centroid_coordinates.x[domain_m1[1,:]]
# 
#     blanked_cells = shift(domain_m1[1:overlap_cell_num,:], -5)
# 
#     for c_m1 in blanked_cells, c_m2 in CartesianIndices((1:size(x_m2)[1], 1:size(x_m2)[2]))
#         tolerance_x = max(abs(mesh1.cell_center_metrics.x₁.ξ[domain_m1][c_m1]), abs(mesh1.cell_center_metrics.x₁.η[domain_m1][c_m1])) # * 0.6
#         tolerance_y = max(abs(mesh1.cell_center_metrics.x₂.ξ[domain_m1][c_m1]), abs(mesh1.cell_center_metrics.x₂.η[domain_m1][c_m1])) # * 0.6
#         if abs(centroid_m1_x[c_m1]-x_m2[c_m2]) ≤ tolerance_x && abs(centroid_m1_y[c_m1]-y_m2[c_m2]) ≤ tolerance_y
#             push!(overlap_points_m2, c_m2)
#         end
#     end
# 
#     return unique(overlap_points_m2)
# end

function cut_holes(mesh1, mesh2)
    # Mesh i should be cutting mesh j, meaning we need to find the CartesianIndices where mesh j will change
    overlap = []

    grid_i = mesh1
    grid_j = mesh2
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
        quad3 = remove_opposite_pairs!(intersect(xis_lower, yis_lower), size(x_mi)[1], size(x_mi)[1])
        quad1 = remove_opposite_pairs!(intersect(xis_upper, yis_upper), size(x_mi)[1], size(x_mi)[1])
        quad2 = remove_opposite_pairs!(intersect(xis_lower, yis_upper), size(x_mi)[1], size(x_mi)[1])
        quad4 = remove_opposite_pairs!(intersect(xis_upper, yis_lower), size(x_mi)[1], size(x_mi)[1])

        if quad1 && quad2 && quad3 && quad4
            push!(overlap, c_mj)
        end
    end
    return overlap
end

# function remove_opposite_pairs!(idx_list::Vector, nrows::Int, ncols::Int)
#     to_remove = Set{CartesianIndex}()
# 
#     for idx in idx_list
#         i, j = Tuple(idx)
#         # Skip if already marked for removal
#         if idx in to_remove
#             continue
#         end
# 
#         # Determine the opposite index
#         if i == 1
#             opposite = CartesianIndex(nrows, j)
#         elseif i == nrows
#             opposite = CartesianIndex(1, j)
#         elseif j == 1
#             opposite = CartesianIndex(i, ncols)
#         elseif j == ncols
#             opposite = CartesianIndex(i, 1)
#         else
#             continue  # not on the boundary
#         end
# 
#         # If the opposite is in the list, mark both for removal
#         if opposite in idx_list
#             push!(to_remove, idx)
#             push!(to_remove, opposite)
#         end
#     end
# 
#     # Remove marked indices from the original list (in-place)
#     filter!(x -> x ∉ to_remove, idx_list)
# 
#     return idx_list
# end
function remove_opposite_pairs!(idx_list::Vector, nrows::Int, ncols::Int)
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
