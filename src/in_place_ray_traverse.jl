using ImplicitBVH
using StaticArrays
function traverse_rays!(
    bvh::ImplicitBVH.BVH,
    point::SVector{3, Float64},
    directions::Vector{SVector{3, Float64}},
    bvtt1::Vector{SVector{2, Int}},
    bvtt2::Vector{SVector{2, Int}},
    output::Vector{Int},
    start_level::Int=1,
)
    # Allocate and add all possible BVTT contact pairs to start with
    num_bvtt = initial_bvtt!(bvh, bvtt1, start_level)
    num_checks = num_bvtt

    level = start_level
    while level < bvh.tree.levels
        # We can have maximum 2 new checks per BV-ray-pair; resize destination BVTT accordingly
        num_bvtt = traverse_rays_nodes!(bvh, point, directions,
                                        bvtt1, bvtt2, num_bvtt,
                                        level)
        num_checks += num_bvtt

        # Swap source and destination buffers for next iteration
        tmp = bvtt1
        bvtt1 = bvtt2
        bvtt2 = tmp
        # bvtt1, bvtt2 = bvtt2, bvtt1
        level += 1
    end

    num_bvtt = traverse_rays_leaves!(bvh, point, directions,
                                     bvtt1, bvtt2, num_bvtt,
                                     )
    @inbounds if isodd(bvh.tree.levels - start_level)
        output[1] = 1 
        output[2] = num_bvtt
    else
        output[1] = 2 
        output[2] = num_bvtt 
    end
end


function initial_bvtt!(
    bvh::ImplicitBVH.BVH,
    bvtt1::Vector{SVector{2, Int}},
    start_level::Int,
)
    num_rays = 6

    # Generate all possible contact checks between all nodes at the given start_level and all rays
    level_nodes = pow2(start_level - 1)

    # Number of real nodes at the given start_level and number of checks we'll do
    num_real = level_nodes - bvh.tree.virtual_leaves >> (bvh.tree.levels - start_level)
    level_checks = num_real * num_rays

    # If we're not at leaf-level, allocate enough memory for next BVTT expansion
    if start_level == bvh.tree.levels
        initial_number = level_checks
    else
        initial_number = 2 * level_checks
    end

    # Insert all checks to do at this level
    fill_initial_bvtt_rays!(bvtt1, level_nodes, num_real, num_rays)

    return level_checks
end


@inline function fill_initial_bvtt_rays!(bvtt1::Vector{SVector{2, Int}}, level_nodes::Int, num_real::Int, num_rays::Int)
    num_bvtt = 0
    @inbounds for i in level_nodes:level_nodes + num_real - 1
        # Node-node pair checks
        for j in 1:num_rays
            num_bvtt += 1
            bvtt1[num_bvtt] = SVector(i, j)
        end
    end
end

pow2(n::I) where I <: Integer = one(I) << n

@inline function traverse_rays_nodes!(bvh::ImplicitBVH.BVH, point::SVector{3, Float64}, directions::Vector{SVector{3, Float64}}, src::Vector{SVector{2, Int}}, dst::Vector{SVector{2, Int}}, num_src::Int, level::Int)
    # Traverse nodes when level is above leaves

    # Compute number of virtual elements before this level to skip when computing the memory index
    virtual_nodes_level = bvh.tree.virtual_leaves >> (bvh.tree.levels - (level - 1))
    virtual_nodes_before = 2 * virtual_nodes_level - count_ones(virtual_nodes_level)

    num_dst = traverse_rays_nodes_range!(
        bvh, point, directions,
        src, dst,
        virtual_nodes_before,
        1, num_src,
    )

    return num_dst
end

@inline function traverse_rays_nodes_range!(
    bvh::ImplicitBVH.BVH, point::SVector{3, Float64}, directions::Vector{SVector{3, Float64}}, src::Vector{SVector{2, Int}}, dst::Vector{SVector{2, Int}}, num_skips::Int, ilow::Int, ihi::Int,
)
    # Check src[irange[1]:irange[2]] and write to dst[1:num_dst]; dst should be given as a view
    num_dst = 0

    # For each BVTT node-ray pair, check for intersection
    @inbounds for i in ilow:ihi
        # Extract implicit indices of BVH nodes to test
        implicit, iray = src[i]
        node = bvh.nodes[implicit - num_skips]

        # Extract ray

        # If the node and ray is touching, expand BVTT with new possible contacts - i.e. pair
        if isintersection(node, point[1], point[2], point[3], directions[iray][1], directions[iray][2], directions[iray][3])
            # If a node's right child is virtual, don't add that check. Guaranteed to always have
            # at least one real child

            # BVH node's right child is virtual
            if ImplicitBVH.isvirtual(bvh.tree, 2 * implicit + 1)
                dst[num_dst + 1] = SVector(implicit * 2, iray)
                num_dst += 1
            else
                dst[num_dst + 1] = SVector(implicit * 2, iray)
                dst[num_dst + 2] = SVector(implicit * 2 + 1, iray)
                num_dst += 2
            end
        end
    end

    return num_dst
end

@inline function traverse_rays_leaves!(bvh::ImplicitBVH.BVH, point::SVector{3, Float64}, directions::Vector{SVector{3, Float64}}, src::Vector{SVector{2, Int}}, intersections::Vector{SVector{2, Int}}, num_src::Int)
    # Traverse final level, only doing ray-leaf checks

    num_intersections = traverse_rays_leaves_range!(
        bvh, point, directions,
        src, intersections,
        1, num_src,
    )

    return num_intersections
end

@inline function traverse_rays_leaves_range!(
    bvh::ImplicitBVH.BVH, point::SVector{3, Float64}, directions::Vector{SVector{3, Float64}}, src::Vector{SVector{2, Int}}, intersections::Vector{SVector{2, Int}}, ilow::Int, ihi::Int,
)
    # Check src[irange[1]:irange[2]] and write to dst[1:num_dst]; dst should be given as a view
    num_dst = 0

    # Number of implicit indices above leaf-level
    num_above = pow2(bvh.tree.levels - 1) - 1

    # For each BVTT node-ray pair, check for intersection
    @inbounds for i in ilow:ihi
        # Extract implicit indices of BVH leaves to test
        implicit, iray = src[i]

        iorder = bvh.order[implicit - num_above]
        leaf = bvh.leaves[iorder]

        # If leaf-ray intersection, save in intersections
        if isintersection(leaf, point[1], point[2], point[3], directions[iray][1], directions[iray][2], directions[iray][3])
            intersections[num_dst + 1] = SVector(iorder, iray)
            num_dst += 1
        end
    end

    return num_dst
end

@inline function isintersection(b::ImplicitBVH.BBox, p1::Float64, p2::Float64, p3::Float64, d1::Float64, d2::Float64, d3::Float64)
    T = Float64

    @inbounds begin
        inv_d1 = one(T) / d1
        inv_d2 = one(T) / d2
        inv_d3 = one(T) / d3

        # Set x bounds
        t_bound_x1 = (b.lo[1] - p1) * inv_d1
        t_bound_x2 = (b.up[1] - p1) * inv_d1

        tmin = min(t_bound_x1, t_bound_x2)
        tmax = max(t_bound_x1, t_bound_x2)

        # Set y bounds
        t_bound_y1 = (b.lo[2] - p2) * inv_d2
        t_bound_y2 = (b.up[2] - p2) * inv_d2

        tmin = max(tmin, min(t_bound_y1, t_bound_y2))
        tmax = min(tmax, max(t_bound_y1, t_bound_y2))

        # Set z bounds
        t_bound_z1 = (b.lo[3] - p3) * inv_d3
        t_bound_z2 = (b.up[3] - p3) * inv_d3

        tmin = max(tmin, min(t_bound_z1, t_bound_z2))
        tmax = min(tmax, max(t_bound_z1, t_bound_z2))
    end
        
    # If condition satisfied ray intersects box. tmax >= 0 
    # ensure only forwards intersections are counted
    (tmin <= tmax) && (tmax >= 0)
end
