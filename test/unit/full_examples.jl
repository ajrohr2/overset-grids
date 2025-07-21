@testset "Nested Half-Spheres" begin
    inner_half_sphere = rthetaphi_grid((0, 0, 0), (0.75, π, π), (20, 20, 20), :meg6)
    outer_half_sphere = rthetaphi_grid((0, 0, 0), (1, π, π), (20, 20, 20), :meg6)
    
    grids = (inner_half_sphere, outer_half_sphere)
    zs = (1, 2)

    meshes = create_components(grids, zs)

    @test all(==(0), meshes[1][1].component_mesh.blank_mask)
    
    correct_mask = ones(Int8, (20, 20, 20))
    correct_mask[11:15, :, :] .= -1
    correct_mask[16:end, :, :] .= 0

    @test correct_mask == meshes[2][1].component_mesh.blank_mask
end

@testset "Nested Rectangles" begin
    inner_rectangle = RectilinearGrid2D((0.25, 0.25), (0.75, 0.75), (30, 40), :meg6)
    outer_rectangle = RectilinearGrid2D((0, 0), (1, 1), (20, 30), :meg6)

    grids = (inner_rectangle, outer_rectangle)
    zs = (1, 2)

    meshes = create_components(grids, zs; num_interp_points = 3)

    @test all(==(0), meshes[1][1].component_mesh.blank_mask)

    correct_mask = zeros(Int8, (20, 30))
    correct_mask[6:15, 8:23] .= 1

    for c in CartesianIndices(correct_mask)
        if (c[1] ∈ (6, 7, 8) && c[2] ∈ collect(8:23)) || (c[1] ∈ (13, 14, 15) && c[2] ∈ collect(8:23)) || (c[1] ∈ collect(6:15) && c[2] ∈ (8, 9, 10)) || (c[1] ∈ collect(6:15) && c[2] ∈ (21, 22, 23))
            correct_mask[c] = -1
        end
    end

    @test meshes[2][1].component_mesh.blank_mask == correct_mask
end
