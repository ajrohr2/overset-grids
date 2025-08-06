# Helper struct for intersection detections
abstract type ComponentMesh end

struct Line
    point_0::NTuple{2, Float64}
    point_1::NTuple{2, Float64}
end
struct ComponentMesh2D <: ComponentMesh
    blank_mask::Array{Int8, 2}
    grid_index::Int
    z_order::Int
    name::Union{String, Nothing}
    bounding_box::NTuple{4, Float64}
    boundary_polygon::Vector{Line}
    kdtree::NearestNeighbors.KDTree{SVector{2, Float64}, Euclidean, Float64, SVector{2, Float64}}
end

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
    direction::SVector{3, Float64}
end
struct ComponentMesh3D <: ComponentMesh
    blank_mask::Array{Int8, 3}
    grid_index::Int
    z_order::Int
    name::Union{String, Nothing}
    bounding_box::NTuple{6, Float64}
    boundary_polygon::Vector{Plane}
    bvh::ImplicitBVH.BVH{Int32, Vector{ImplicitBVH.BBox{Float64}}, Vector{ImplicitBVH.BBox{Float64}}, Vector{UInt32}, Vector{Int32}}
    kdtree::NearestNeighbors.KDTree{SVector{3, Float64}, Euclidean, Float64, SVector{3, Float64}}
    x::SubArray{Float64, 3, Array{Float64, 3}, Tuple{UnitRange{Int64}, UnitRange{Int64}, UnitRange{Int64}}, false}
    y::SubArray{Float64, 3, Array{Float64, 3}, Tuple{UnitRange{Int64}, UnitRange{Int64}, UnitRange{Int64}}, false}
    z::SubArray{Float64, 3, Array{Float64, 3}, Tuple{UnitRange{Int64}, UnitRange{Int64}, UnitRange{Int64}}, false}
end
