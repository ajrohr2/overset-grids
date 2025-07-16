@testset "Intersection" begin
    # Intersecting plane, inside triangle
    in_face = OversetGrids.Plane(SVector(0, -1, 0), SVector(1, 0, 0), SVector(0, 1, 0))
    in_ray = OversetGrids.Ray(SVector(0.5, 0, -1), SVector(0.5, 0, 1))

    @test OversetGrids.determine_intersection(in_face, in_ray) == true

    # Intersecting plane, outside triangle
    out_face = OversetGrids.Plane(SVector(0, -1, 0), SVector(1, 0, 0), SVector(0, 1, 0))
    out_ray = OversetGrids.Ray(SVector(2, -1, 0), SVector(2, 1, 0))
    
    @test OversetGrids.determine_intersection(out_face, out_ray) == false

    # Ray and plane colinear
    col_face = OversetGrids.Plane(SVector(0, -1, 0), SVector(1, 0, 0), SVector(0, 1, 0))
    col_ray = OversetGrids.Ray(SVector(1, 0, 0), SVector(2, 0, 0))

    @test OversetGrids.determine_intersection(col_face, col_ray) == false
end

@testset "Boundary Creation" begin
    grid = UniformGrid3D((0, 0, 0), (1, 1, 1), 0.1, :meg6)

    bounds, polygon = OversetGrids.get_boundary(grid)
    @test bounds == (0, 1, 0, 1, 0, 1)

    correct_bottom = []
    
    for i in 0:0.1:0.9, j in 0:0.1:0.9
        push!(correct_bottom, OversetGrids.Plane(SVector(i, j, 0), SVector(i+0.1, j, 0), SVector(i+0.1, j+0.1, 0)))
        push!(correct_bottom, OversetGrids.Plane(SVector(i, j, 0), SVector(i, j+0.1, 0), SVector(i+0.1, j+0.1, 0)))
    end

    for face in correct_bottom
        @test any(f -> isapprox(f.point_1[1], face.point_1[1]), polygon) && any(f -> isapprox(f.point_1[2], face.point_1[2]), polygon) && any(f -> isapprox(f.point_1[3], face.point_1[3]), polygon) && any(f -> isapprox(f.point_0[1], face.point_0[1]), polygon) && any(f -> isapprox(f.point_0[2], face.point_0[2]), polygon) && any(f -> isapprox(f.point_0[3], face.point_0[3]), polygon) && any(f -> isapprox(f.point_2[1], face.point_2[1]), polygon) && any(f -> isapprox(f.point_2[2], face.point_2[2]), polygon) && any(f -> isapprox(f.point_2[3], face.point_2[3]), polygon)

    end
end
