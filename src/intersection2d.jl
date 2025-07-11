# --- Orientation and ray intersection based overlap detection --- #
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

function create_rays!(point::CartesianIndex, rays, x_array, y_array, bottom_corner, top_corner)
    # We need the rays to exist in the largest global space
    perturbation_x = abs(x_array[2] - x_array[1]) * ℯ / 100
    perturbation_y = abs(y_array[2] - y_array[1]) * ℯ / 100

    rays[1] = Line((bottom_corner[1], y_array[1, point[2]]+perturbation_y), (x_array[point], y_array[point]))
    rays[2] = Line((x_array[point], y_array[point]), (top_corner[1], y_array[end, point[2]]+perturbation_y))
    rays[3] = Line((x_array[point[1], 1]+perturbation_x, bottom_corner[2]), (x_array[point], y_array[point]))
    rays[4] = Line((x_array[point], y_array[point]), (x_array[point[1], end]+perturbation_x, top_corner[2]))
end

function get_boundary!(grid::CurvilinearGrids.AbstractCurvilinearGrid2D)
    x, y = CurvilinearGrids.coords(grid)
    mi_inds = CartesianIndices((1:size(x,1), 1:size(x,2)))
    boundary_polygon = create_boundary_polygon(mi_inds, x, y)
    return (minimum(x), maximum(x), minimum(y), maximum(y)), boundary_polygon
end
function create_boundary_polygon(cartesian_indices, x_array, y_array)
    boundary1 = cartesian_indices[1,:]
    boundary2 = cartesian_indices[:,1]
    boundary3 = cartesian_indices[end,:]
    boundary4 = cartesian_indices[:,end]

    polygon = Line[]

    for i in 1:length(boundary1)-1
        push!(polygon, Line((x_array[boundary1[i]], y_array[boundary1[i]]), (x_array[boundary1[i+1]], y_array[boundary1[i+1]])))
    end
    for i in 1:length(boundary2)-1
        push!(polygon, Line((x_array[boundary2[i]], y_array[boundary2[i]]), (x_array[boundary2[i+1]], y_array[boundary2[i+1]])))
    end
    for i in 1:length(boundary3)-1
        push!(polygon, Line((x_array[boundary3[i]], y_array[boundary3[i]]), (x_array[boundary3[i+1]], y_array[boundary3[i+1]])))
    end
    for i in 1:length(boundary4)-1
        push!(polygon, Line((x_array[boundary4[i]], y_array[boundary4[i]]), (x_array[boundary4[i+1]], y_array[boundary4[i+1]])))
    end
    
    return polygon
end
