using LinearAlgebra, StaticArrays
# function lagrange_basis(x, x_vals, i)
#     l = 1.0
#     for j in eachindex(x_vals)
#         if j != i
#             l *= (x - x_vals[j]) / (x_vals[i] - x_vals[j])
#         end 
#     end 
#     return l 
# end 
# 
# function lagrange_interpolate_2D(x, y, x_vals, y_vals, function_vals)
#     n = length(x_vals)
#     m = length(y_vals)
#     
#     result = 0.0
# 
#     for i in 1:n
#         for j in 1:m
#             result += function_vals[i,j] * lagrange_basis(x, x_vals, i) * lagrange_basis(y, y_vals, j)
#         end
#     end
#     return result
# end

weight(point, point_i, p) = 1 / norm(point - point_i)^p

function find_k_nearest_2d(point, field_x, field_y, blank; k=5)
    neighbors = Vector{CartesianIndex{2}}()
    for i in 1:k
        neighbor = CartesianIndex(1,1)
        dist = Inf # Should probably be the max distance in the grid.
        for c in CartesianIndices(field_x)
            tmp_dist = norm(point - SVector(field_x[c], field_y[c]))
            if tmp_dist < dist && blank[c] < 1 && c ∉ neighbors
                neighbor = c
                dist = tmp_dist
            end
        end
        push!(neighbors, neighbor)
    end
    if length(neighbors) < k
        error("Not enough neighbors")
    end
    return neighbors
end

function inverse_distance_weighting(point, num_sample_points, x_vals, y_vals, function_vals, blank_mask, p)
    known_cart = find_k_nearest_2d(point, x_vals, y_vals, blank_mask, k=num_sample_points)
    numerator = 0.0
    denominator = 0.0
    for i in 1:num_sample_points
        # Power of 3 is arbitrary
        numerator += weight(point, SVector(x_vals[known_cart[i]], y_vals[known_cart[i]]), p) * function_vals[known_cart[i]]
        denominator += weight(point, SVector(x_vals[known_cart[i]], y_vals[known_cart[i]]), p)
    end
    return numerator / denominator
end
# function find_k_nearest_2d(point, field_x, field_y, blank; k=5)
#     neighbors_x = Vector{CartesianIndex{2}}(undef, k)
#     neighbors_y = Vector{CartesianIndex{2}}(undef, k)
#     for i in 1:k
#         neighbor_x = CartesianIndex(1,1)
#         dist_x = abs(field_x[neighbor_x] - point[1])
#         neighbor_y = CartesianIndex(1,1)
#         dist_y = abs(field_y[neighbor_y] - point[2])
#         for c in CartesianIndices(field_x)
#             tmp_dist_x = abs(field_x[c] - point[1])
#             tmp_dist_y = abs(field_y[c] - point[2])
#             if tmp_dist_x < dist_x && blank[c] != 2 && c ∉ neighbors_x
#                 neighbor_x = c
#                 dist_x = tmp_dist_x
#             end
#             if tmp_dist_y < dist_y && blank[c] != 2 && c ∉ neighbors_y
#                 neighbor_y = c
#                 dist_y = tmp_dist_y
#             end
#         end
#         neighbors_x[i] = neighbor_x
#         neighbors_y[i] = neighbor_y
#     end
#     # While this looks like it returns (y,x), it actually returns (x,y). So when passing the x and y vals to the interpolate function, use ...coordinates.x[c] for c in x
#     return neighbors_y, neighbors_x
# end
