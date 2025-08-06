@inline function determine_intersection(face::Plane, bp::SVector{3, Float64}, dir::SVector{3, Float64})
    denom = dot(dir, face.normal)

    # Nearly perpendicular, so the ray likely doesn't intersect the face inside a triangle.
    if abs(denom) < 1e-5
        return false
    end

    t = dot((face.point_0 - bp), face.normal) / denom

    if t < 0; return false; end;

    point_of_intersection = bp + (t * dir)

    # Determine if the point of intersection is within the triangular face
    v0 = face.point_2 - face.point_0 
    v1 = face.point_1 - face.point_0 
    v2 = point_of_intersection - face.point_0

    d00 = dot(v0, v0)
    d01 = dot(v0, v1)
    d11 = dot(v1, v1)
    d20 = dot(v2, v0)
    d21 = dot(v2, v1)

    denom_bary = d00 * d11 - d01 * d01
    inv_denom = 1.0 / denom_bary
    
    β = (d11 * d20 - d01 * d21) * inv_denom
    γ = (d00 * d21 - d01 * d20) * inv_denom
    α = 1 - β - γ

    # If one is 1 and others are 0, then the point is a vertex. If one is 0 and the other are in (0, 1), then the point is on an edge. This may cause double intersection issues down the line
    return 0 ≤ α ≤ 1 && 0 ≤ β ≤ 1 && 0 ≤ γ ≤ 1 ? true : false
end

@inline norm_zeros(v::SVector{3, Float64}) = SVector{3, Float64}(ifelse(v[1] == 0.0, 0.0, v[1]), ifelse(v[2] == 0.0, 0.0, v[2]), ifelse(v[3] == 0.0, 0.0, v[3]))
@inline function canonical_points_key(p::Plane)
    a = Tuple(norm_zeros(p.point_0))
    b = Tuple(norm_zeros(p.point_1))
    c = Tuple(norm_zeros(p.point_2))
    return sort(SVector(a,b,c))
end
function get_boundary(grid::CurvilinearGrids.AbstractCurvilinearGrid3D)
    x, y, z = CurvilinearGrids.coords(grid)
    mi_inds = CartesianIndices(size(x))
    boundary_polygon = create_boundary_planes(mi_inds, x, y, z)
    bvh = ImplicitBVH.BVH(map(tri -> ImplicitBVH.BBox(tri.point_0, tri.point_1, tri.point_2), boundary_polygon))
    return (minimum(x), maximum(x), minimum(y), maximum(y), minimum(z), maximum(z)), boundary_polygon, bvh
end
@inbounds function create_boundary_planes(cartesian_indices, x, y, z)
    boundary = [cartesian_indices[1,:,:],cartesian_indices[end,:,:],cartesian_indices[:,1,:],cartesian_indices[:,end,:],cartesian_indices[:,:,1],cartesian_indices[:,:,end]]

    polygon = Plane[]

    d = CartesianIndex(1,0)
    dr = CartesianIndex(1,1)
    r = CartesianIndex(0,1)

    for b in boundary
        for c in CartesianIndices(b)
            if c[2] != size(CartesianIndices(b))[1] && c[1] != size(CartesianIndices(b))[1]
                p0, p1, p2, p3 = SVector{3, Float64}(x[b[c]],y[b[c]],z[b[c]]),SVector{3, Float64}(x[b[c+d]],y[b[c+d]],z[b[c+d]]), SVector{3, Float64}(x[b[c+dr]],y[b[c+dr]],z[b[c+dr]]),SVector{3, Float64}(x[b[c+r]],y[b[c+r]],z[b[c+r]])
                if !is_degenerate(p0, p1, p2)
                    push!(polygon, Plane(p0, p1, p2))
                end
                if !is_degenerate(p0, p3, p2)
                    push!(polygon, Plane(p0, p3, p2))
                end
            end
        end
    end

    return unique(canonical_points_key, polygon)
end
@inline function is_degenerate(p0::SVector{3, Float64}, p1::SVector{3, Float64}, p2::SVector{3, Float64}; atol=1e-10)
    v1 = p1 - p0
    v2 = p2 - p0
    area = norm(cross(v1, v2)) / 2 
    return area < atol
end

