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
